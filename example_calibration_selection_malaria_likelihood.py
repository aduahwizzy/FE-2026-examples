from math import lgamma
import os
from glob import glob

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

import manifest


user = os.getlogin()
expt_name = f"{user}_FE_example_pickup_proj_test"
output_dir = os.path.join(manifest.job_directory, "my_outputs")
experiment_output_dir = os.path.join(output_dir, expt_name, 'smc_simulations')
prevalence_file = os.path.join(manifest.input_dir, "malaria_prevalence.csv")

PARAMETER_VARIABLES = [
    "x_Temporary_Larval_Habitat",
    "scale_constant",
    "scale_water_vegetation",
    "scale_temp_rain",
]
WATER_VEGETATION_ALIAS = "scale_water_veg"
RUN_VARIABLE = "Run_Number"
DATE_COLUMNS = ["year", "month"]


def beta_binomial(raw_nobs, sim_nobs, raw_data, sim_data):
    """Return idmtools-calibra's mean beta-binomial log-likelihood."""
    likelihood = 0.0
    observation_count = 0
    for raw_n, sim_n, raw_pos, sim_pos in zip(
        raw_nobs, sim_nobs, raw_data, sim_data
    ):
        likelihood += lgamma(raw_n + 1)
        likelihood += lgamma(sim_n + 2)
        likelihood -= lgamma(raw_n + sim_n + 2)
        likelihood += lgamma(raw_pos + sim_pos + 1)
        likelihood += lgamma(raw_n - raw_pos + sim_n - sim_pos + 1)
        likelihood -= lgamma(raw_pos + 1)
        likelihood -= lgamma(raw_n - raw_pos + 1)
        likelihood -= lgamma(sim_pos + 1)
        likelihood -= lgamma(sim_n - sim_pos + 1)
        observation_count += 1

    return likelihood / observation_count if observation_count else likelihood


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
    if "year" not in sim_df.columns and "Year" in sim_df.columns:
        sim_df = sim_df.rename(columns={"Year": "year"})
    column_aliases = {
        "PfPR U5": "PfPR",
        "PfPR_U5": "PfPR",
        "Pop U5": "Pop",
        "Pop_U5": "Pop",
    }
    sim_df = sim_df.rename(
        columns={
            alias: canonical
            for alias, canonical in column_aliases.items()
            if canonical not in sim_df.columns and alias in sim_df.columns
        }
    )

    if "agebin" in sim_df.columns:
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

    if (
        "scale_water_vegetation" not in sim_df.columns
        and WATER_VEGETATION_ALIAS in sim_df.columns
    ):
        sim_df = sim_df.rename(
            columns={WATER_VEGETATION_ALIAS: "scale_water_vegetation"}
        )
    elif (
        "scale_water_vegetation" in sim_df.columns
        and WATER_VEGETATION_ALIAS in sim_df.columns
    ):
        raise ValueError(
            "Simulation CSV contains both 'scale_water_vegetation' and "
            f"'{WATER_VEGETATION_ALIAS}'; remove the ambiguous duplicate column."
        )

    target_df = pd.read_csv(prevalence_file)
    if "year" not in target_df.columns and "Year" in target_df.columns:
        target_df = target_df.rename(columns={"Year": "year"})

    required_sim_columns = set(
        PARAMETER_VARIABLES
        + [RUN_VARIABLE, "year", "month", "PfPR", "Pop"]
    )
    missing_sim_columns = required_sim_columns.difference(sim_df.columns)
    if missing_sim_columns:
        raise ValueError(
            "Simulation CSV is missing required columns: "
            f"{', '.join(sorted(missing_sim_columns))}. Ensure the analyzer "
            "exports all four parameter values and Run_Number."
        )

    required_target_columns = {"year", "month", "DHS_n", "DHS_pos"}
    missing_target_columns = required_target_columns.difference(target_df.columns)
    if missing_target_columns:
        raise ValueError(
            "Prevalence target CSV is missing required columns: "
            f"{', '.join(sorted(missing_target_columns))}"
        )

    for column in ["year", "month"]:
        sim_df[column] = pd.to_numeric(sim_df[column], errors="raise").astype(int)
        target_df[column] = pd.to_numeric(target_df[column], errors="raise").astype(
            int
        )
    for column in ["PfPR", "Pop"] + PARAMETER_VARIABLES + [RUN_VARIABLE]:
        sim_df[column] = pd.to_numeric(sim_df[column], errors="raise")
    for column in ["DHS_n", "DHS_pos"]:
        target_df[column] = pd.to_numeric(target_df[column], errors="raise")

    if not sim_df["PfPR"].between(0, 1).all():
        raise ValueError("Simulated PfPR values must be proportions between 0 and 1.")
    if (sim_df["Pop"] <= 0).any():
        raise ValueError("Simulated population values must be greater than zero.")
    if target_df.duplicated(DATE_COLUMNS).any():
        raise ValueError(
            "Prevalence target CSV must contain at most one row per year and month."
        )
    if not target_df["month"].between(1, 12).all():
        raise ValueError("Prevalence target months must be between 1 and 12.")
    if (target_df["DHS_n"] <= 0).any():
        raise ValueError("DHS_n sample sizes must be greater than zero.")
    if not np.equal(target_df[["DHS_n", "DHS_pos"]] % 1, 0).all().all():
        raise ValueError("DHS_n and DHS_pos must contain integer counts.")
    if ((target_df["DHS_pos"] < 0) | (target_df["DHS_pos"] > target_df["DHS_n"])).any():
        raise ValueError("DHS_pos values must be between zero and DHS_n.")

    return sim_df, target_df


def score_parameter_sets(sim_df, target_df):
    parameter_set_variables = PARAMETER_VARIABLES + [RUN_VARIABLE]
    target_dates = target_df[DATE_COLUMNS].drop_duplicates()
    joined = sim_df.merge(
        target_df[DATE_COLUMNS + ["DHS_n", "DHS_pos"]],
        on=DATE_COLUMNS,
        how="inner",
        validate="many_to_one",
    )

    matched_dates = joined[DATE_COLUMNS].drop_duplicates()
    missing_dates = target_dates.merge(
        matched_dates, on=DATE_COLUMNS, how="left", indicator=True
    )
    if (missing_dates["_merge"] == "left_only").any():
        dates = missing_dates.loc[missing_dates["_merge"] == "left_only", DATE_COLUMNS]
        raise ValueError(
            "Simulation output has no prevalence for target date(s): "
            + ", ".join(
                f"{row.year}-{row.month:02d}" for row in dates.itertuples(index=False)
            )
        )

    expected_date_count = len(target_dates)
    scored_date_counts = joined.groupby(parameter_set_variables)[DATE_COLUMNS].apply(
        lambda dates: len(dates.drop_duplicates())
    )
    if (scored_date_counts != expected_date_count).any():
        raise ValueError(
            "Some parameter-set/run combinations do not have output for every "
            "target survey month."
        )

    joined["npos"] = (joined["PfPR"] * joined["Pop"]).round()
    run_scores = []
    for group_key, candidate_df in joined.groupby(
        parameter_set_variables, sort=False
    ):
        likelihood = beta_binomial(
            candidate_df["DHS_n"].to_numpy(),
            candidate_df["Pop"].to_numpy(),
            candidate_df["DHS_pos"].to_numpy(),
            candidate_df["npos"].to_numpy(),
        )
        run_scores.append(
            {
                **dict(zip(parameter_set_variables, group_key)),
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
        raise ValueError("No parameter sets could be scored against the target months.")
    scores.insert(0, "Rank", np.arange(1, len(scores) + 1))
    scores.insert(1, "Candidate", [f"candidate_{rank:03d}" for rank in scores["Rank"]])
    return run_scores, scores


def plot_selection(sim_df, target_df, run_scores, scores):
    best = scores.iloc[0]
    best_parameters = best[PARAMETER_VARIABLES]
    best_parameter_mask = np.logical_and.reduce(
        [sim_df[column].eq(best_parameters[column]) for column in PARAMETER_VARIABLES]
    )
    selected_simulations = sim_df.loc[best_parameter_mask].copy()
    parameter_mask = np.logical_and.reduce(
        [
            run_scores[column].eq(best_parameters[column])
            for column in PARAMETER_VARIABLES
        ]
    )
    selected_runs = run_scores.loc[parameter_mask]
    if selected_runs.empty:
        raise ValueError("No per-run likelihoods found for the best-fit parameters.")
    best_run = selected_runs.sort_values("Log_Likelihood", ascending=False).iloc[0][
        RUN_VARIABLE
    ]

    best_fit = selected_simulations.loc[
        selected_simulations[RUN_VARIABLE].eq(best_run)
    ].groupby(DATE_COLUMNS, as_index=False)["PfPR"].mean()
    best_fit["date"] = pd.to_datetime(
        dict(year=best_fit["year"], month=best_fit["month"], day=1)
    )

    run_prevalence = (
        selected_simulations.groupby([RUN_VARIABLE] + DATE_COLUMNS, as_index=False)[
            "PfPR"
        ]
        .mean()
    )
    mean_ci = (
        run_prevalence.groupby(DATE_COLUMNS)["PfPR"]
        .agg(["mean", "std", "count"])
        .reset_index()
    )
    mean_ci["ci"] = 1.96 * mean_ci["std"].fillna(0) / np.sqrt(mean_ci["count"])
    mean_ci["date"] = pd.to_datetime(
        dict(year=mean_ci["year"], month=mean_ci["month"], day=1)
    )

    observed = target_df.copy()
    observed["Observed_Prevalence"] = observed["DHS_pos"] / observed["DHS_n"]
    observed["date"] = pd.to_datetime(
        dict(year=observed["year"], month=observed["month"], day=1)
    )

    fig, axis = plt.subplots(figsize=(11, 6))
    axis.plot(
        best_fit["date"],
        best_fit["PfPR"],
        marker="o",
        linewidth=2,
        color="#C0392B",
        label=f"Best-fit simulation (Run_Number={best_run})",
    )
    axis.scatter(
        observed["date"],
        observed["Observed_Prevalence"],
        color="black",
        marker="x",
        s=70,
        label="Observed DHS",
        zorder=5,
    )
    axis.set_title("Best-fit simulation vs observed under-five PfPR")
    axis.set_ylabel("Under 5 PfPR")
    axis.set_ylim(bottom=0)
    axis.legend()
    fig.tight_layout()
    fig.savefig(
        os.path.join(experiment_output_dir, "selection_malaria_best_fit.png"),
        dpi=150,
    )
    plt.close(fig)

    simulation_profiles = (
        sim_df.groupby(
            PARAMETER_VARIABLES + [RUN_VARIABLE] + DATE_COLUMNS,
            as_index=False,
        )["PfPR"]
        .mean()
    )
    simulation_profiles["date"] = pd.to_datetime(
        dict(
            year=simulation_profiles["year"],
            month=simulation_profiles["month"],
            day=1,
        )
    )
    fig, axis = plt.subplots(figsize=(11, 6))
    background_labeled = False
    for group_key, profile in simulation_profiles.groupby(
        PARAMETER_VARIABLES + [RUN_VARIABLE], sort=False
    ):
        parameter_values = group_key[: len(PARAMETER_VARIABLES)]
        run_number = group_key[-1]
        is_best_simulation = (
            all(
                np.isclose(value, best_parameters[column])
                for column, value in zip(PARAMETER_VARIABLES, parameter_values)
            )
            and np.isclose(run_number, best_run)
        )
        if is_best_simulation:
            continue
        axis.plot(
            profile["date"],
            profile[" Under 5 PfPR"],
            color="0.7",
            alpha=0.35,
            linewidth=0.8,
            label="Other simulations" if not background_labeled else None,
        )
        background_labeled = True

    axis.plot(
        best_fit["date"],
        best_fit["PfPR"],
        color="#C0392B",
        linewidth=2.2,
        label=f"Best-fit simulation (Run_Number={best_run})",
        zorder=3,
    )
    axis.scatter(
        observed["date"],
        observed["Observed_Prevalence"],
        color="black",
        marker="x",
        s=70,
        label="Observed DHS",
        zorder=5,
    )
    axis.set_title("Under-five PfPR across all simulations")
    axis.set_ylabel("PfPR")
    axis.set_ylim(bottom=0)
    axis.legend()
    fig.tight_layout()
    fig.savefig(
        os.path.join(experiment_output_dir, "selection_malaria_all_simulations.png"),
        dpi=150,
    )
    plt.close(fig)

    fig, axis = plt.subplots(figsize=(11, 6))
    axis.plot(
        mean_ci["date"],
        mean_ci["mean"],
        marker="o",
        linewidth=2,
        color="#2471A3",
        label="Mean simulated prevalence",
    )
    axis.fill_between(
        mean_ci["date"],
        (mean_ci["mean"] - mean_ci["ci"]).clip(lower=0),
        (mean_ci["mean"] + mean_ci["ci"]).clip(upper=1),
        color="#2471A3",
        alpha=0.2,
        label="95% CI of mean across Run_Number replicates",
    )
    axis.scatter(
        observed["date"],
        observed["Observed_Prevalence"],
        color="black",
        marker="x",
        s=70,
        label="Observed DHS",
        zorder=5,
    )
    axis.set_title("Mean simulated prevalence and 95% CI for best-fit parameters")
    axis.set_ylabel("Prevalence")
    axis.set_ylim(bottom=0)
    axis.legend()
    fig.tight_layout()
    fig.savefig(
        os.path.join(experiment_output_dir, "selection_malaria_mean_ci.png"),
        dpi=150,
    )
    plt.close(fig)


if __name__ == "__main__":
    sim_df, target_df = load_data()
    run_scores, scores = score_parameter_sets(sim_df, target_df)

    scores.to_csv(
        os.path.join(experiment_output_dir, "malaria_likelihood_scores.csv"),
        index=False,
    )
    run_scores.to_csv(
        os.path.join(experiment_output_dir, "malaria_likelihood_by_run.csv"),
        index=False,
    )
    plot_selection(sim_df, target_df, run_scores, scores)

    print("Best-fitting parameter combination (highest mean log-likelihood):")
    print(scores.iloc[0].to_string())
    print(
        "\nSaved ranked likelihood scores, per-run scores, and plot to "
        f"{experiment_output_dir}"
    )
