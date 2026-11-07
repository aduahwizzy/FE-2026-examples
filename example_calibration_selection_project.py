import os
from glob import glob

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

import manifest


user = os.getlogin()
expt_name = f"{user}_FE_example_pickup_proj"
output_dir = os.path.join(manifest.job_directory, "my_outputs")
experiment_output_dir = os.path.join(output_dir, expt_name)

prevalence_file = os.path.join(manifest.input_dir, "malaria_prevalence.csv")

PARAMETER_VARIABLES = [
    "x_Temporary_Larval_Habitat",
    "scale_constant",
    "scale_temp_rain",
    "scale_water_veg",
]
RUN_VARIABLE = "Run_Number"
YEAR_COLUMN = "Year"
TARGET_PREVALENCE_COLUMN = "Prevalence"


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
        PARAMETER_VARIABLES + [RUN_VARIABLE, YEAR_COLUMN, "month", "PfPR"]
    )
    missing_sim_columns = required_sim_columns.difference(sim_df.columns)
    if missing_sim_columns:
        raise ValueError(
            "Simulation CSV is missing required sweep metadata columns: "
            f"{', '.join(sorted(missing_sim_columns))}. "
            "Rerun analyzer_calibration.py with all habitat sweep variables "
            "included in sweep_variables, then rerun this selector. Existing "
            "CSV files cannot recover sweep values that were not exported."
        )

    required_target_columns = {YEAR_COLUMN, TARGET_PREVALENCE_COLUMN}
    missing_target_columns = required_target_columns.difference(target_df.columns)
    if missing_target_columns:
        raise ValueError(
            f"Prevalence target CSV is missing required columns: "
            f"{', '.join(sorted(missing_target_columns))}"
        )

    sim_df[YEAR_COLUMN] = pd.to_numeric(sim_df[YEAR_COLUMN], errors="raise")
    sim_df[YEAR_COLUMN] = sim_df[YEAR_COLUMN].astype(int)
    sim_df["month"] = pd.to_numeric(sim_df["month"], errors="raise")
    sim_df["PfPR"] = pd.to_numeric(sim_df["PfPR"], errors="raise")
    target_df[YEAR_COLUMN] = pd.to_numeric(
        target_df[YEAR_COLUMN], errors="raise"
    ).astype(int)
    target_df[TARGET_PREVALENCE_COLUMN] = (
        pd.to_numeric(target_df[TARGET_PREVALENCE_COLUMN], errors="raise") / 100.0
    )

    for variable in PARAMETER_VARIABLES + [RUN_VARIABLE]:
        sim_df[variable] = pd.to_numeric(sim_df[variable], errors="raise")

    if not sim_df["PfPR"].between(0, 1).all():
        raise ValueError("Simulated PfPR values must be proportions between 0 and 1.")
    if not target_df[TARGET_PREVALENCE_COLUMN].between(0, 1).all():
        raise ValueError(
            "Observed prevalence percentages must be between 0 and 100."
        )

    return sim_df, target_df


def make_annual_u5_prevalence(sim_df):
    """Average the monthly under-five PfPR values within each simulation year."""
    return (
        sim_df.groupby(
            PARAMETER_VARIABLES + [RUN_VARIABLE, YEAR_COLUMN],
            as_index=False,
        )["PfPR"]
        .mean()
        .rename(columns={"PfPR": "Simulated_Prevalence"})
    )


def score_parameter_sets(run_annual, target_df):
    # Average stochastic replicates for each complete habitat parameter set.
    annual = (
        run_annual.groupby(
            PARAMETER_VARIABLES + [YEAR_COLUMN], as_index=False
        )["Simulated_Prevalence"]
        .mean()
    )
    target_years = set(target_df[YEAR_COLUMN])
    available_years = set(annual[YEAR_COLUMN])
    missing_years = target_years.difference(available_years)
    if missing_years:
        raise ValueError(
            "Simulation output has no annual prevalence for target year(s): "
            + ", ".join(str(year) for year in sorted(missing_years))
        )

    joined = annual.merge(
        target_df[[YEAR_COLUMN, TARGET_PREVALENCE_COLUMN]],
        on=YEAR_COLUMN,
        how="inner",
        validate="many_to_one",
    )
    joined["Squared_Error"] = (
        joined["Simulated_Prevalence"] - joined[TARGET_PREVALENCE_COLUMN]
    ) ** 2
    scores = (
        joined.groupby(PARAMETER_VARIABLES, as_index=False)
        .agg(
            RMSE=("Squared_Error", lambda errors: np.sqrt(errors.mean())),
            Mean_Squared_Error=("Squared_Error", "mean"),
            Years_Scored=(YEAR_COLUMN, "nunique"),
        )
        .sort_values("RMSE", ascending=True)
        .reset_index(drop=True)
    )
    if scores.empty:
        raise ValueError("No parameter sets could be scored against the target years.")
    scores.insert(0, "Rank", np.arange(1, len(scores) + 1))
    scores.insert(1, "Candidate", [f"candidate_{i:03d}" for i in scores["Rank"]])
    return annual, scores


def plot_selection(annual, target_df, scores):
    best = scores.iloc[0]
    top_candidates = scores.head(min(10, len(scores)))
    top_annual = annual.merge(
        top_candidates[["Candidate"] + PARAMETER_VARIABLES],
        on=PARAMETER_VARIABLES,
        how="inner",
        validate="many_to_one",
    )

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
        target_df[YEAR_COLUMN],
        target_df[TARGET_PREVALENCE_COLUMN],
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

    axes[1].plot(scores["Rank"], scores["RMSE"], marker="o", linewidth=1)
    axes[1].scatter(
        [best["Rank"]],
        [best["RMSE"]],
        color="red",
        s=60,
        label="Best fit",
        zorder=5,
    )
    axes[1].set_title("Parameter combinations ranked by RMSE")
    axes[1].set_xlabel("Candidate rank (see ranked CSV for parameters)")
    axes[1].set_ylabel("RMSE (under-five prevalence proportion)")
    axes[1].legend()

    fig.tight_layout()
    fig.savefig(os.path.join(experiment_output_dir, "selection_project.png"), dpi=150)
    plt.close(fig)


if __name__ == "__main__":
    sim_df, target_df = load_data()
    run_annual = make_annual_u5_prevalence(sim_df)
    annual, scores = score_parameter_sets(run_annual, target_df)

    scores.to_csv(
        os.path.join(experiment_output_dir, "calibration_scores.csv"),
        index=False,
    )
    annual.to_csv(
        os.path.join(experiment_output_dir, "annual_simulated_prevalence.csv"),
        index=False,
    )
    plot_selection(annual, target_df, scores)

    print("Best-fitting habitat parameter combination (lowest RMSE):")
    print(scores.iloc[0].to_string())
    print(
        "\nSaved ranked scores, annual simulated prevalence, and plot to "
        f"{experiment_output_dir}"
    )
