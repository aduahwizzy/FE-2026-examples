#!/usr/bin/env python3
"""Extract age-specific malaria outcomes and compare the scenario experiments."""

import argparse
import os
import re
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

import manifest
from idmtools.entities import IAnalyzer
from idmtools.entities.simulation import Simulation


ROOT = Path(__file__).resolve().parent
SCENARIO_FILE = ROOT / "inputs" / "smc_vaccine_scenarios.csv"
REPORT_YEARS = range(2026, 2028)
REPORT_FILES = [
    f"output/MalariaSummaryReport_Scenario_{year}.json"
    for year in REPORT_YEARS
]
AGE_GROUP_BINS = {
    "under5": {0.25, 0.5, 5},
    "age3_59m": {0.5, 5},
    "age6_59m": {5},
    "all_ages": {0.25, 0.5, 5, 115},
}
AGE_LABELS = {
    "under5": "Under 5",
    "age3_59m": "3-59 months",
    "age6_59m": "6-59 months",
    "all_ages": "All ages",
}
RATE_CHANNELS = ["PfPR", "Cases", "Severe_cases"]
INCIDENCE_DISPLAY_FACTOR = 0.1
CASE_SCENARIO_COLORS = {
    "S0": "#34495E",
    "S1": "#2471A3",
    "S2": "#E67E22",
    "S3": "#C0392B",
}
CASE_SCENARIO_LABELS = {
    "S0": "S0: None",
    "S1": "S1: SMC only",
    "S2": "S2: Vaccine only",
    "S3": "S3: Combined SMC + vaccine",
}


class MonthlyScenarioAnalyzer(IAnalyzer):
    def __init__(self, output_dir):
        super().__init__(working_dir=str(output_dir), filenames=REPORT_FILES)
        self.output_dir = Path(output_dir)

    def map(self, data, simulation: Simulation):
        frames = []
        for filename in self.filenames:
            if filename not in data:
                raise FileNotFoundError(
                    f"Simulation {simulation.id} is missing {filename}."
                )
            report = data[filename]
            year_match = re.search(r"Scenario_(\d{4})", filename)
            if not year_match:
                raise ValueError(f"Unable to read reporting year from {filename}.")
            year = int(year_match.group(1))
            age_bins = report["Metadata"]["Age Bins"]
            report_data = report["DataByTimeAndAgeBins"]
            channels = {
                "PfPR": "PfPR by Age Bin",
                "Cases": "Annual Clinical Incidence by Age Bin",
                "Severe_cases": "Annual Severe Incidence by Age Bin",
                "New_infections": "New Infections by Age Bin",
                "Pop": "Average Population by Age Bin",
            }
            values = {
                column: report_data[channel][:-1]
                for column, channel in channels.items()
            }
            month_count = len(values["PfPR"])
            if month_count != 12 or any(
                len(channel_values) != month_count
                for channel_values in values.values()
            ):
                raise ValueError(
                    f"{filename} for simulation {simulation.id} should contain "
                    "12 monthly records plus the final annual summary."
                )
            for age_index, age_bin in enumerate(age_bins):
                frame = pd.DataFrame(
                    {
                        "Year": year,
                        "month": np.arange(1, month_count + 1),
                        "agebin": float(age_bin),
                        **{
                            column: [row[age_index] for row in values[column]]
                            for column in channels
                        },
                    }
                )
                frames.append(frame)

        result = pd.concat(frames, ignore_index=True)
        for tag in (
            "Run_Number",
            "Scenario",
            "Stage",
            "SMC_Future",
            "Vaccine_Future",
            "SMC_Coverage",
            "Vaccine_Completion_Coverage",
            "x_Temporary_Larval_Habitat",
            "scale_constant",
            "scale_temp_rain",
            "scale_water_veg",
        ):
            if tag in simulation.tags:
                result[tag] = simulation.tags[tag]
        if "Scenario" not in result.columns or "Run_Number" not in result.columns:
            raise ValueError(
                f"Simulation {simulation.id} is missing scenario or Run_Number tags."
            )
        return result

    def reduce(self, all_data):
        selected = [data for data in all_data.values() if data is not None]
        if not selected:
            raise ValueError("No malaria summary reports were analyzed.")
        self.output_dir.mkdir(parents=True, exist_ok=True)
        output_path = self.output_dir / "Monthly_Scenario_Outcomes.csv"
        pd.concat(selected, ignore_index=True).to_csv(output_path, index=False)
        print(f"Saved monthly scenario outcomes to {output_path}")


def parse_args():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--experiment-id", required=True)
    parser.add_argument(
        "--output-dir",
        type=Path,
        help="Output directory (defaults to experiments/my_outputs/forward_scenarios/<id>).",
    )
    return parser.parse_args()


def run_analyzer(experiment_id, output_dir):
    from idmtools.analysis.analyze_manager import AnalyzeManager
    from idmtools.core import ItemType
    from idmtools.core.platform_factory import Platform

    output_dir.mkdir(parents=True, exist_ok=True)
    with Platform(
        manifest.default_platform,
        job_directory=manifest.job_directory,
        docker_image=manifest.plat_image,
    ):
        analyzer = MonthlyScenarioAnalyzer(output_dir)
        manager = AnalyzeManager(
            configuration={},
            ids=[(experiment_id, ItemType.EXPERIMENT)],
            analyzers=[analyzer],
            partial_analyze_ok=False,
        )
        manager.analyze()


def make_age_specific_monthly(data):
    numeric_columns = [
        "Year",
        "month",
        "Run_Number",
        "agebin",
        "PfPR",
        "Cases",
        "Severe_cases",
        "New_infections",
        "Pop",
    ]
    for column in numeric_columns:
        data[column] = pd.to_numeric(data[column], errors="raise")
    data["Year"] = data["Year"].astype(int)
    data["month"] = data["month"].astype(int)
    data["Run_Number"] = data["Run_Number"].astype(int)

    key_columns = ["Scenario", "Run_Number", "Year", "month"]
    output = []
    for age_group, bins in AGE_GROUP_BINS.items():
        selected = data.loc[
            data["agebin"].apply(lambda value: any(np.isclose(value, b) for b in bins))
        ].copy()
        observed_bins = set(selected["agebin"].unique())
        if any(not any(np.isclose(value, bin_value) for value in observed_bins)
               for bin_value in bins):
            raise ValueError(
                f"Report data is missing age bins {sorted(bins)} required for "
                f"{AGE_LABELS[age_group]}."
            )

        selected["weighted_PfPR"] = selected["PfPR"] * selected["Pop"]
        selected["weighted_Cases"] = selected["Cases"] * selected["Pop"]
        selected["weighted_Severe_cases"] = (
            selected["Severe_cases"] * selected["Pop"]
        )
        summary = selected.groupby(key_columns, as_index=False, sort=False).agg(
            Pop=("Pop", "sum"),
            weighted_PfPR=("weighted_PfPR", "sum"),
            weighted_Cases=("weighted_Cases", "sum"),
            weighted_Severe_cases=("weighted_Severe_cases", "sum"),
            New_infections=("New_infections", "sum"),
        )
        summary["PfPR"] = summary["weighted_PfPR"] / summary["Pop"]
        summary["Cases"] = summary["weighted_Cases"] / summary["Pop"]
        summary["Severe_cases"] = summary["weighted_Severe_cases"] / summary["Pop"]
        summary["Age_Group"] = age_group
        output.append(
            summary.drop(
                columns=[
                    "weighted_PfPR",
                    "weighted_Cases",
                    "weighted_Severe_cases",
                ]
            )
        )
    return pd.concat(output, ignore_index=True)


def run_level_summary(monthly):
    run_year = (
        monthly.groupby(["Scenario", "Run_Number", "Age_Group", "Year"], as_index=False)[
            RATE_CHANNELS
        ]
        .mean()
    )
    run_average = (
        run_year.groupby(["Scenario", "Run_Number", "Age_Group"], as_index=False)[
            RATE_CHANNELS
        ]
        .mean()
    )
    run_average["Clinical_Incidence_per_1000_per_year"] = (
        run_average["Cases"] * 1000
    )
    run_average["Severe_Incidence_per_1000_per_year"] = (
        run_average["Severe_cases"] * 1000
    )
    run_average["PfPR_percent"] = run_average["PfPR"] * 100
    return run_average


def burden_category(row):
    incidence = row["Clinical_Incidence_per_1000_per_year_mean"]
    prevalence = row["PfPR_percent_mean"]
    incidence_category = (
        "Very low" if incidence < 50 else
        "Low" if incidence < 200 else
        "Moderate" if incidence < 500 else
        "High"
    )
    prevalence_category = (
        "Very low" if prevalence < 1 else
        "Low" if prevalence <= 5 else
        "Moderate" if prevalence < 15 else
        "High"
    )
    return pd.Series(
        {
            "all_age_incidence_reference_category": incidence_category,
            "all_age_prevalence_reference_category": prevalence_category,
        }
    )


def summarize_scenarios(run_average, output_dir):
    metrics = [
        "Clinical_Incidence_per_1000_per_year",
        "Severe_Incidence_per_1000_per_year",
        "PfPR_percent",
    ]
    grouped = run_average.groupby(["Scenario", "Age_Group"])[metrics].agg(
        ["mean", "std", "count"]
    )
    summary = grouped.reset_index()
    summary.columns = [
        "_".join(str(part) for part in column if part)
        if isinstance(column, tuple)
        else column
        for column in summary.columns
    ]
    summary_path = output_dir / "Scenario_Summary_5year.csv"
    summary.to_csv(summary_path, index=False)

    all_age = summary.loc[summary["Age_Group"] == "all_ages"].copy()
    if not all_age.empty:
        categories = all_age.apply(burden_category, axis=1)
        all_age = pd.concat([all_age.reset_index(drop=True), categories], axis=1)
        all_age.to_csv(output_dir / "All_Age_Burden_Reference.csv", index=False)
    return summary


def paired_stage1_effects(run_average, output_dir):
    metrics = [
        "Clinical_Incidence_per_1000_per_year",
        "Severe_Incidence_per_1000_per_year",
        "PfPR_percent",
    ]
    effects = []
    for age_group in AGE_GROUP_BINS:
        subset = run_average.loc[run_average["Age_Group"].eq(age_group)]
        for metric in metrics:
            wide = subset.pivot(
                index="Run_Number", columns="Scenario", values=metric
            )
            required = {"S0", "S1", "S2", "S3"}
            if not required.issubset(wide.columns):
                raise ValueError(
                    f"Missing stage-1 scenarios for {age_group}: "
                    f"{sorted(required.difference(wide.columns))}"
                )
            paired = pd.DataFrame(
                {
                    "SMC_only_vs_reference_percent": (
                        np.where(
                            wide["S0"].ne(0),
                            (wide["S0"] - wide["S1"]) / wide["S0"] * 100,
                            np.nan,
                        )
                    ),
                    "Vaccine_only_vs_reference_percent": (
                        np.where(
                            wide["S0"].ne(0),
                            (wide["S0"] - wide["S2"]) / wide["S0"] * 100,
                            np.nan,
                        )
                    ),
                    "Combined_vs_reference_percent": (
                        np.where(
                            wide["S0"].ne(0),
                            (wide["S0"] - wide["S3"]) / wide["S0"] * 100,
                            np.nan,
                        )
                    ),
                    "Vaccine_added_to_SMC_percent": (
                        np.where(
                            wide["S1"].ne(0),
                            (wide["S1"] - wide["S3"]) / wide["S1"] * 100,
                            np.nan,
                        )
                    ),
                }
            )
            paired["Age_Group"] = age_group
            paired["Outcome"] = metric
            paired["Run_Number"] = paired.index
            effects.append(paired.reset_index(drop=True))
    per_run = pd.concat(effects, ignore_index=True)
    per_run.to_csv(output_dir / "Stage1_Paired_Effects_by_Run.csv", index=False)
    per_run.groupby(["Age_Group", "Outcome"], as_index=False).agg(
        **{
            f"{column}_{stat}": (column, stat)
            for column in (
                "SMC_only_vs_reference_percent",
                "Vaccine_only_vs_reference_percent",
                "Combined_vs_reference_percent",
                "Vaccine_added_to_SMC_percent",
            )
            for stat in ("mean", "std", "count")
        }
    ).to_csv(output_dir / "Stage1_Paired_Effects_Summary.csv", index=False)


def plot_stage1_monthly(monthly, output_dir):
    stage1 = monthly.loc[
        monthly["Scenario"].isin(["S0", "S1", "S2", "S3"])
        & monthly["Age_Group"].eq("under5")
    ]
    grouped = (
        stage1.groupby(["Scenario", "Year", "month"], as_index=False)[
            ["Cases", "PfPR"]
        ]
        .mean()
    )
    fig, axes = plt.subplots(2, 1, figsize=(12, 8), sharex=True)
    for scenario in ("S0", "S1", "S2", "S3"):
        frame = grouped.loc[grouped["Scenario"].eq(scenario)].copy()
        frame["date"] = pd.to_datetime(
            dict(year=frame["Year"], month=frame["month"], day=15)
        )
        axes[0].plot(
            frame["date"],
            frame["Cases"] * 1000 * INCIDENCE_DISPLAY_FACTOR,
            linewidth=1.5,
            color=CASE_SCENARIO_COLORS[scenario],
            label=CASE_SCENARIO_LABELS[scenario],
        )
        axes[1].plot(
            frame["date"],
            frame["PfPR"] * 100,
            linewidth=1.5,
            color=CASE_SCENARIO_COLORS[scenario],
            label=CASE_SCENARIO_LABELS[scenario],
        )
    axes[0].set_ylabel("Clinical incidence / 1,000 / year")
    axes[0].set_title("Stage 1: under-five monthly outcomes")
    axes[1].set_ylabel("PfPR (%)")
    axes[1].set_xlabel("Month")
    for axis in axes:
        axis.grid(alpha=0.2)
        axis.legend(title="Scenario")
    fig.tight_layout()
    fig.savefig(output_dir / "Stage1_Under5_Monthly.png", dpi=150)
    plt.close(fig)


def plot_coverage_heatmaps(run_average, scenarios, output_dir):
    stage2 = run_average.loc[
        run_average["Scenario"].isin(
            scenarios.loc[scenarios["stage"].eq("stage2"), "scenario_id"]
        )
    ]
    coverage = scenarios.loc[
        scenarios["stage"].eq("stage2"),
        [
            "scenario_id",
            "smc_coverage",
            "vaccine_completion_coverage",
        ],
    ]
    mean_values = (
        stage2.groupby(["Scenario", "Age_Group"], as_index=False)[
            ["Clinical_Incidence_per_1000_per_year", "PfPR_percent"]
        ]
        .mean()
        .merge(coverage, left_on="Scenario", right_on="scenario_id", validate="many_to_one")
    )
    for age_group in ("under5", "age6_59m"):
        age_data = mean_values.loc[mean_values["Age_Group"].eq(age_group)]
        for metric, short_name, label in (
            (
                "Clinical_Incidence_per_1000_per_year",
                "incidence",
                "Clinical incidence / 1,000 / year",
            ),
            ("PfPR_percent", "pfpr", "PfPR (%)"),
        ):
            matrix = age_data.pivot(
                index="smc_coverage",
                columns="vaccine_completion_coverage",
                values=metric,
            ).sort_index(ascending=True)
            matrix = matrix.reindex(sorted(matrix.columns), axis=1)
            display_matrix = (
                matrix * INCIDENCE_DISPLAY_FACTOR
                if metric == "Clinical_Incidence_per_1000_per_year"
                else matrix
            )
            fig, axis = plt.subplots(figsize=(7, 5))
            image = axis.imshow(
                display_matrix.to_numpy(), cmap="YlOrRd", aspect="auto"
            )
            axis.set_xticks(np.arange(len(matrix.columns)))
            axis.set_xticklabels([f"{value:.0%}" for value in matrix.columns])
            axis.set_yticks(np.arange(len(matrix.index)))
            axis.set_yticklabels([f"{value:.0%}" for value in matrix.index])
            axis.set_xlabel("Four-dose vaccine completion")
            axis.set_ylabel("SMC coverage per cycle")
            axis.set_title(f"Stage 2: {label} - {AGE_LABELS[age_group]}")
            for row_index in range(matrix.shape[0]):
                for column_index in range(matrix.shape[1]):
                    axis.text(
                        column_index,
                        row_index,
                        f"{display_matrix.iloc[row_index, column_index]:.1f}",
                        ha="center",
                        va="center",
                        color="black",
                    )
            fig.colorbar(image, ax=axis, label=label)
            fig.tight_layout()
            fig.savefig(
                output_dir / f"Stage2_{age_group}_{short_name}_heatmap.png",
                dpi=150,
            )
            plt.close(fig)
    mean_values.to_csv(output_dir / "Stage2_Coverage_Grid_Summary.csv", index=False)


def process_results(monthly_path, output_dir):
    output_dir.mkdir(parents=True, exist_ok=True)
    data = pd.read_csv(monthly_path)
    if data.empty:
        raise ValueError(f"No analyzed data in {monthly_path}.")
    age_monthly = make_age_specific_monthly(data)
    age_monthly.to_csv(output_dir / "Monthly_Age_Group_Outcomes.csv", index=False)
    run_average = run_level_summary(age_monthly)
    run_average.to_csv(output_dir / "Run_Level_5year_Outcomes.csv", index=False)
    scenarios = pd.read_csv(SCENARIO_FILE)
    summary = summarize_scenarios(run_average, output_dir)
    paired_stage1_effects(run_average, output_dir)
    plot_stage1_monthly(age_monthly, output_dir)
    plot_coverage_heatmaps(run_average, scenarios, output_dir)
    print(f"Saved scenario summaries and figures to {output_dir}")
    print(
        "For Ghana's regional low-burden reference, use All_Age_Burden_Reference.csv; "
        "the national incidence bands are not under-five thresholds."
    )


def main():
    args = parse_args()
    output_dir = args.output_dir or (
        Path(manifest.job_directory)
        / "my_outputs"
        / "forward_scenarios"
        / args.experiment_id
    )
    run_analyzer(args.experiment_id, output_dir)
    process_results(output_dir / "Monthly_Scenario_Outcomes.csv", output_dir)


if __name__ == "__main__":
    main()
