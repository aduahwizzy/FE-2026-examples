import os
from glob import glob

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from idmtools_calibra.utilities.ll_calculators import beta_binomial

import manifest


user = os.getlogin()
expt_name = f"{user}_FE_example_pickup_proj"
output_dir = os.path.join(manifest.job_directory, "my_outputs")
experiment_output_dir = os.path.join(output_dir, expt_name)
prevalence_file = os.path.join(manifest.input_dir, "prevalence.csv")

PARAMETER_VARIABLES = [
    "x_Temporary_Larval_Habitat",
    "scale_constant",
    "scale_temp_rain",
    "scale_water_veg",
]
RUN_VARIABLE = "Run_Number"
YEAR_COLUMN = "Year"


def load_data():
    preferred_file = os.path.join(
        experiment_output_dir, "U5_PfPR_ClinicalIncidence.csv"
    )
    simulation_files = (
        [preferred_file]
        if os.path.isfile(preferred_file)
        else sorted(
            glob(os.path.join(experiment_output_dir, "*PfPR_ClinicalIncidence*.csv"))
        )
    )
    if not simulation_files:
        raise FileNotFoundError(
            "No PfPR_ClinicalIncidence CSV outputs found in "
            f"{experiment_output_dir}"
        )
    sim_df = pd.concat(
        [pd.read_csv(filename) for filename in simulation_files],
        ignore_index=True,
    )
    if "PfPR U5" in sim_df.columns:
        sim_df = sim_df.rename(columns={"PfPR U5": "PfPR"})
    elif "PfPR_U5" in sim_df.columns:
        sim_df = sim_df.rename(columns={"PfPR_U5": "PfPR"})

    if "agebin" in sim_df.columns:
        # With age_bins=[0.25, 5, 115], the bin labeled 5 is age 0.25-5.
        u5_mask = pd.to_numeric(sim_df["agebin"], errors="coerce").eq(5)
        if not u5_mask.any():
            raise ValueError(
                "Simulation CSV has age-bin rows but no bin labeled 5, expected "
                "for the 0.25-5 year age bin."
            )
        sim_df = sim_df.loc[u5_mask].copy()
    elif not any("U5" in os.path.basename(filename) for filename in simulation_files):
        raise ValueError(
            "The simulation CSV does not identify its PfPR as under-five. "
            "Use U5_PfPR_ClinicalIncidence.csv, a PfPR U5 column, or an "
            "agebin column with the 0.25-5 bin labeled 5."
        )

    target_df = pd.read_csv(prevalence_file)
    if "year" in sim_df.columns and YEAR_COLUMN not in sim_df.columns:
        sim_df = sim_df.rename(columns={"year": YEAR_COLUMN})

    required_sim_columns = set(
        PARAMETER_VARIABLES
        + [RUN_VARIABLE, YEAR_COLUMN, "month", "PfPR", "Pop"]
    )
    missing_sim_columns = required_sim_columns.difference(sim_df.columns)
    if missing_sim_columns:
        raise ValueError(
            "Simulation CSV is missing required columns: "
            f"{', '.join(sorted(missing_sim_columns))}"
        )

    required_target_columns = {YEAR_COLUMN, "DHS_n", "DHS_pos"}
    missing_target_columns = required_target_columns.difference(target_df.columns)
    if missing_target_columns:
        raise ValueError(
            "Prevalence target CSV is missing required columns: "
            f"{', '.join(sorted(missing_target_columns))}"
        )

    sim_df[YEAR_COLUMN] = pd.to_numeric(sim_df[YEAR_COLUMN], errors="raise").astype(
        int
    )
    for column in ["month", "PfPR", "Pop"] + PARAMETER_VARIABLES + [RUN_VARIABLE]:
        sim_df[column] = pd.to_numeric(sim_df[column], errors="raise")
    target_df[YEAR_COLUMN] = pd.to_numeric(
        target_df[YEAR_COLUMN], errors="raise"
    ).astype(int)
    for column in ["DHS_n", "DHS_pos"]:
        target_df[column] = pd.to_numeric(target_df[column], errors="raise")

    if not sim_df["PfPR"].between(0, 1).all():
        raise ValueError("Simulated PfPR values must be proportions between 0 and 1.")
    if (sim_df["Pop"] <= 0).any():
        raise ValueError("Simulated population values must be greater than zero.")
    if target_df[YEAR_COLUMN].duplicated().any():
        raise ValueError("Prevalence target CSV must contain at most one row per year.")
    if (target_df["DHS_n"] <= 0).any():
        raise ValueError("DHS_n sample sizes must be greater than zero.")
    if not np.equal(target_df[["DHS_n", "DHS_pos"]] % 1, 0).all().all():
        raise ValueError("DHS_n and DHS_pos must contain integer counts.")
    if ((target_df["DHS_pos"] < 0) | (target_df["DHS_pos"] > target_df["DHS_n"])).any():
        raise ValueError("DHS_pos values must be between zero and DHS_n.")

    return sim_df, target_df


def make_annual_u5_prevalence(sim_df):
    """Average monthly U5 PfPR and population values within each simulation year."""
    return (
        sim_df.groupby(
            PARAMETER_VARIABLES + [RUN_VARIABLE, YEAR_COLUMN],
            as_index=False,
        )[["PfPR", "Pop"]]
        .mean()
        .rename(columns={"PfPR": "Simulated_Prevalence"})
    )


def score_parameter_sets(run_annual, target_df):
    target_years = set(target_df[YEAR_COLUMN])
    available_years = set(run_annual[YEAR_COLUMN])
    missing_years = target_years.difference(available_years)
    if missing_years:
        raise ValueError(
            "Simulation output has no annual prevalence for target year(s): "
            + ", ".join(str(year) for year in sorted(missing_years))
        )

    joined = run_annual.merge(
        target_df[[YEAR_COLUMN, "DHS_n", "DHS_pos"]],
        on=YEAR_COLUMN,
        how="inner",
        validate="many_to_one",
    )
    expected_year_count = len(target_years)
    scored_year_counts = joined.groupby(
        PARAMETER_VARIABLES + [RUN_VARIABLE]
    )[YEAR_COLUMN].nunique()
    if (scored_year_counts != expected_year_count).any():
        raise ValueError(
            "Some parameter-set/run combinations do not have simulation output "
            "for every target year."
        )
    joined["npos"] = (joined["Simulated_Prevalence"] * joined["Pop"]).round()

    run_scores = []
    group_variables = PARAMETER_VARIABLES + [RUN_VARIABLE]
    for group_key, candidate_df in joined.groupby(group_variables, sort=False):
        likelihood = beta_binomial(
            candidate_df["DHS_n"].to_numpy(),
            candidate_df["Pop"].to_numpy(),
            candidate_df["DHS_pos"].to_numpy(),
            candidate_df["npos"].to_numpy(),
        )
        run_scores.append(
            {
                **dict(zip(group_variables, group_key)),
                "Log_Likelihood": float(likelihood),
            }
        )
    run_scores = pd.DataFrame(run_scores)

    scores = (
        run_scores.groupby(PARAMETER_VARIABLES, as_index=False)
        .agg(
            Mean_Log_Likelihood=("Log_Likelihood", "mean"),
            Runs_Scored=(RUN_VARIABLE, "nunique"),
        )
        .sort_values("Mean_Log_Likelihood", ascending=False)
        .reset_index(drop=True)
    )
    if scores.empty:
        raise ValueError("No parameter sets could be scored against the target years.")
    scores.insert(0, "Rank", np.arange(1, len(scores) + 1))
    scores.insert(1, "Candidate", [f"candidate_{i:03d}" for i in scores["Rank"]])
    return run_scores, scores


def plot_selection(annual, target_df, scores):
    top_candidates = scores.head(min(10, len(scores)))
    top_annual = annual.merge(
        top_candidates[["Candidate"] + PARAMETER_VARIABLES],
        on=PARAMETER_VARIABLES,
        how="inner",
        validate="many_to_many",
    ).groupby(["Candidate", YEAR_COLUMN], as_index=False)[
        "Simulated_Prevalence"
    ].mean()
    observed = target_df.copy()
    observed["Observed_Prevalence"] = observed["DHS_pos"] / observed["DHS_n"]

    fig, axes = plt.subplots(1, 2, figsize=(15, 6))
    for candidate, candidate_df in top_annual.groupby("Candidate"):
        axes[0].plot(
            candidate_df[YEAR_COLUMN],
            candidate_df["Simulated_Prevalence"],
            marker="o",
            linewidth=1,
            alpha=0.85,
            label=candidate,
        )
    axes[0].scatter(
        observed[YEAR_COLUMN],
        observed["Observed_Prevalence"],
        color="black",
        marker="x",
        s=70,
        label="Observed",
        zorder=5,
    )
    axes[0].set_title("Observed vs simulated annual under-five prevalence")
    axes[0].set_xlabel("Year")
    axes[0].set_ylabel("Under-five prevalence")
    axes[0].set_ylim(bottom=0)
    axes[0].legend(fontsize="small")

    best = scores.iloc[0]
    axes[1].plot(
        scores["Rank"], scores["Mean_Log_Likelihood"], marker="o", linewidth=1
    )
    axes[1].scatter(
        [best["Rank"]],
        [best["Mean_Log_Likelihood"]],
        color="red",
        s=60,
        label="Best fit",
        zorder=5,
    )
    axes[1].set_title("Parameter combinations ranked by beta-binomial likelihood")
    axes[1].set_xlabel("Candidate rank (see ranked CSV for parameters)")
    axes[1].set_ylabel("Mean log-likelihood (higher is better)")
    axes[1].legend()

    fig.tight_layout()
    fig.savefig(
        os.path.join(experiment_output_dir, "selection_project_likelihood.png"),
        dpi=150,
    )
    plt.close(fig)


if __name__ == "__main__":
    sim_df, target_df = load_data()
    annual = make_annual_u5_prevalence(sim_df)
    run_scores, scores = score_parameter_sets(annual, target_df)

    scores.to_csv(
        os.path.join(experiment_output_dir, "calibration_likelihood_scores.csv"),
        index=False,
    )
    run_scores.to_csv(
        os.path.join(experiment_output_dir, "calibration_likelihood_by_run.csv"),
        index=False,
    )
    annual.to_csv(
        os.path.join(
            experiment_output_dir, "annual_simulated_prevalence_likelihood.csv"
        ),
        index=False,
    )
    plot_selection(annual, target_df, scores)

    print("Best-fitting habitat parameter combination (highest mean log-likelihood):")
    print(scores.iloc[0].to_string())
    print(
        "\nSaved ranked likelihood scores, per-run scores, annual simulated "
        f"prevalence, and plot to {experiment_output_dir}"
    )
