#!/usr/bin/env python3
"""Plot monthly simulated PfPR and clinical-incidence seasonality."""

import argparse
import calendar
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


PARAMETER_COLUMNS = [
    "x_Temporary_Larval_Habitat",
    "scale_constant",
    "scale_temp_rain",
    "scale_water_veg",
]
# Edit these values to choose a parameter combination from the analyzer CSVs.
PARAMETER_FILTERS = {
    "x_Temporary_Larval_Habitat": 6.812921,
    "scale_constant": 1.0,
    "scale_temp_rain": 1.737801,
    "scale_water_veg": 1.023293,
}
OBSERVED_INCIDENCE_FILE = (
    Path(__file__).resolve().parent
    / "inputs"
    / "northern_region_malaria_incidence.csv"
)
# Diagnostic visualization only; does not adjust simulation data or imply
# that the model's clinical incidence should be calibrated with this factor.
SIMULATED_INCIDENCE_DISPLAY_FACTOR = 0.1
METRICS = {
    "PfPR": "simulated_pfpr_seasonality.png",
    "Cases": "simulated_incidence_seasonality.png",
}
MONTH_NAMES = [
    "Jan", "Feb", "Mar", "Apr", "May", "Jun",
    "Jul", "Aug", "Sep", "Oct", "Nov", "Dec",
]


def parse_args():
    parser = argparse.ArgumentParser(
        description=(
            "Plot monthly PfPR and clinical-incidence seasonality from "
            "PfPR_ClinicalIncidence_Monthly_*.csv analyzer outputs."
        )
    )
    parser.add_argument(
        "--input-dir",
        required=True,
        type=Path,
        help="Directory containing the monthly analyzer CSV files.",
    )
    parser.add_argument(
        "--output-dir",
        type=Path,
        help="Directory for plots (defaults to --input-dir).",
    )
    parser.add_argument(
        "--compare-dir",
        type=Path,
        help=(
            "Optional second scenario directory. Adds PfPR and incidence "
            "seasonality comparison plots using matched Run_Number/year pairs."
        ),
    )
    parser.add_argument("--input-label", help="Label for --input-dir in comparison plots.")
    parser.add_argument("--compare-label", help="Label for --compare-dir in comparison plots.")
    parser.add_argument(
        "--age-group",
        choices=["under5", "age5plus"],
        default="under5",
        help="Age group to plot; default combines report bins 0.25 and 5.",
    )
    parser.add_argument(
        "--filter",
        action="append",
        default=[],
        metavar="COLUMN=VALUE",
        help=(
            "Select one parameter combination, e.g. "
            "--filter scale_temp_rain=1.0. Repeat for multiple parameters."
        ),
    )
    parser.add_argument(
        "--start-year", type=int, help="First calendar year to include."
    )
    parser.add_argument(
        "--end-year", type=int, help="Last calendar year to include."
    )
    parser.add_argument(
        "--observed-incidence-file",
        type=Path,
        default=OBSERVED_INCIDENCE_FILE,
        help=(
            "District-level Northern Region incidence CSV (defaults to "
            "inputs/northern_region_malaria_incidence.csv)."
        ),
    )
    return parser.parse_args()


def load_data(input_dir):
    files = sorted(input_dir.glob("PfPR_ClinicalIncidence_Monthly_*.csv"))
    if not files:
        raise FileNotFoundError(
            f"No PfPR_ClinicalIncidence_Monthly_*.csv files found in {input_dir}"
        )

    frames = [pd.read_csv(filename) for filename in files]
    data = pd.concat(frames, ignore_index=True)
    if "Year" not in data.columns and "year" in data.columns:
        data = data.rename(columns={"year": "Year"})

    required = {"Year", "month", "Run_Number", "PfPR", "Cases", "Pop"}
    missing = required.difference(data.columns)
    if missing:
        raise ValueError(
            f"Monthly analyzer data is missing required columns: {', '.join(sorted(missing))}"
        )

    numeric_columns = ["Year", "month", "Run_Number", "PfPR", "Cases", "Pop"]
    if "agebin" in data.columns:
        numeric_columns.append("agebin")
    numeric_columns.extend(
        column for column in PARAMETER_COLUMNS if column in data.columns
    )
    for column in numeric_columns:
        data[column] = pd.to_numeric(data[column], errors="raise")

    if not data["month"].between(1, 12).all():
        raise ValueError("Month values must be between 1 and 12.")
    if not data["PfPR"].between(0, 1).all():
        raise ValueError("PfPR values must be proportions between zero and one.")
    if (data["Pop"] <= 0).any():
        raise ValueError("Population values must be greater than zero.")

    data["Year"] = data["Year"].astype(int)
    data["month"] = data["month"].astype(int)
    return data


def apply_filters(data, filters):
    for specification in filters:
        column, separator, value = specification.partition("=")
        if not separator or not column or not value:
            raise ValueError(
                f"Invalid filter {specification!r}; use COLUMN=VALUE."
            )
        if column not in PARAMETER_COLUMNS or column not in data.columns:
            raise ValueError(
                f"Filter column {column!r} is not a parameter column in the input."
            )
        try:
            numeric_value = float(value)
            numeric_column = pd.to_numeric(data[column], errors="raise")
        except ValueError:
            mask = data[column].astype(str).eq(value)
        else:
            mask = np.isclose(numeric_column, numeric_value)
        data = data.loc[mask].copy()
        if data.empty:
            raise ValueError(f"Filter {specification!r} selected no rows.")
    return data


def require_one_parameter_set(data):
    parameter_columns = [
        column for column in PARAMETER_COLUMNS if column in data.columns
    ]
    if not parameter_columns:
        return
    combinations = data[parameter_columns].drop_duplicates()
    if len(combinations) > 1:
        available = combinations.to_dict(orient="records")
        raise ValueError(
            "Input contains multiple parameter combinations. Edit "
            "PARAMETER_FILTERS near the top of the script or select one with "
            "--filter COLUMN=VALUE arguments. "
            f"Available combinations: {available}"
        )


def summarize_age_group(data, age_group):
    group_columns = ["Year", "month", "Run_Number"] + [
        column for column in PARAMETER_COLUMNS if column in data.columns
    ]

    if "agebin" not in data.columns:
        return data

    if age_group == "under5":
        chosen = data.loc[
            np.isclose(data["agebin"], 0.25)
            | np.isclose(data["agebin"], 5)
        ].copy()
        available_bins = set(chosen["agebin"].unique())
        if not {0.25, 5}.issubset(available_bins):
            raise ValueError(
                "Under-five aggregation requires age bins 0.25 and 5 in the "
                "monthly analyzer CSV."
            )
    else:
        chosen = data.loc[np.isclose(data["agebin"], 115)].copy()
        if chosen.empty:
            raise ValueError("Age 5+ aggregation requires an age bin labeled 115.")

    chosen["weighted_PfPR"] = chosen["PfPR"] * chosen["Pop"]
    chosen["weighted_Cases"] = chosen["Cases"] * chosen["Pop"]
    summary = chosen.groupby(group_columns, as_index=False, sort=False).agg(
        Pop=("Pop", "sum"),
        weighted_PfPR=("weighted_PfPR", "sum"),
        weighted_Cases=("weighted_Cases", "sum"),
    )
    summary["PfPR"] = summary["weighted_PfPR"] / summary["Pop"]
    summary["Cases"] = summary["weighted_Cases"] / summary["Pop"]
    return summary.drop(columns=["weighted_PfPR", "weighted_Cases"])


def select_complete_years(data):
    months_per_run_year = data.groupby(["Run_Number", "Year"])["month"].nunique()
    complete = months_per_run_year.eq(12)
    incomplete = months_per_run_year.loc[~complete]
    if not incomplete.empty:
        print(
            "Excluding incomplete Run_Number/year records: "
            + ", ".join(f"{run}/{year}" for run, year in incomplete.index)
        )
    keep = complete[complete].index
    complete_keys = pd.MultiIndex.from_frame(data[["Run_Number", "Year"]])
    data = data.loc[complete_keys.isin(keep)].copy()
    if data.empty:
        raise ValueError("No complete 12-month Run_Number/year records remain.")
    return data


def plot_metric(data, metric, ylabel, title, output_path):
    run_month = (
        data.groupby(["Run_Number", "Year", "month"], as_index=False)[metric]
        .mean()
    )
    monthly_stats = run_month.groupby(["Year", "month"])[metric].agg(
        mean="mean",
        low=lambda values: values.quantile(0.10),
        high=lambda values: values.quantile(0.90),
    ).reset_index()
    monthly_stats["date"] = pd.to_datetime(
        dict(
            year=monthly_stats["Year"],
            month=monthly_stats["month"],
            day=15,
        )
    )

    yearly_profiles = (
        run_month.groupby(["Year", "month"], as_index=False)[metric].mean()
    )
    seasonal_stats = yearly_profiles.groupby("month")[metric].agg(
        mean="mean",
        low=lambda values: values.quantile(0.10),
        high=lambda values: values.quantile(0.90),
    ).reindex(range(1, 13))

    scale = (
        100
        if metric == "PfPR"
        else SIMULATED_INCIDENCE_DISPLAY_FACTOR
    )
    incidence_display_note = (
        " (displayed at 0.1×)" if metric == "Cases" else ""
    )
    fig, (time_axis, seasonal_axis) = plt.subplots(
        2, 1, figsize=(11, 8), gridspec_kw={"height_ratios": [1.2, 1]}
    )

    time_axis.plot(
        monthly_stats["date"],
        monthly_stats["mean"] * scale,
        color="#2471A3",
        linewidth=1.6,
        label=f"Mean across Run_Number{incidence_display_note}",
    )
    time_axis.set_title("Monthly values over simulation years")
    time_axis.set_ylabel(ylabel)
    time_axis.legend()
    time_axis.grid(alpha=0.2)

    for _, yearly_data in yearly_profiles.groupby("Year"):
        seasonal_axis.plot(
            yearly_data["month"],
            yearly_data[metric] * scale,
            color="0.55",
            alpha=0.35,
            linewidth=0.8,
        )
    months = np.arange(1, 13)
    seasonal_axis.plot(
        months,
        seasonal_stats["mean"].to_numpy() * scale,
        color="#C0392B",
        marker="o",
        linewidth=2.2,
        label=f"Mean seasonal profile across years{incidence_display_note}",
    )
    seasonal_axis.fill_between(
        months,
        seasonal_stats["low"].to_numpy() * scale,
        seasonal_stats["high"].to_numpy() * scale,
        color="#C0392B",
        alpha=0.16,
        label=f"10th-90th percentile across years{incidence_display_note}",
    )
    seasonal_axis.set_title("Month-of-year profile (thin lines: individual years)")
    seasonal_axis.set_ylabel(ylabel)
    seasonal_axis.set_xlabel("Month")
    seasonal_axis.set_xticks(months)
    seasonal_axis.set_xticklabels(MONTH_NAMES)
    seasonal_axis.legend()
    seasonal_axis.grid(alpha=0.2)

    first_year, last_year = int(data["Year"].min()), int(data["Year"].max())
    fig.suptitle(f"{title} ({first_year}-{last_year})")
    fig.tight_layout()
    fig.savefig(output_path, dpi=150)
    plt.close(fig)


def load_observed_incidence(filename):
    observed = pd.read_csv(filename)
    required = {"Region", "District", "year", "month", "Pop_5", "conf_u5"}
    missing = required.difference(observed.columns)
    if missing:
        raise ValueError(
            f"Observed incidence CSV {filename} is missing columns: "
            f"{', '.join(sorted(missing))}"
        )

    observed = observed.loc[
        observed["Region"].astype(str).str.strip().str.casefold().eq("northern")
    ].copy()
    if observed.empty:
        raise ValueError(
            f"No Northern Region rows were found in observed incidence CSV {filename}."
        )
    if observed["District"].isna().any():
        raise ValueError("Observed incidence data contains rows without a district.")
    if observed.duplicated(["District", "year", "month"]).any():
        raise ValueError(
            "Observed incidence data must have at most one row per district, "
            "year, and month."
        )

    observed["year"] = pd.to_numeric(observed["year"], errors="raise").astype(int)
    observed["Pop_5"] = pd.to_numeric(observed["Pop_5"], errors="raise")
    observed["conf_u5"] = pd.to_numeric(observed["conf_u5"], errors="raise")
    if (observed["Pop_5"] <= 0).any():
        raise ValueError("Observed Pop_5 values must be greater than zero.")
    if (observed["conf_u5"] < 0).any():
        raise ValueError("Observed conf_u5 values cannot be negative.")

    month_lookup = {
        name.casefold(): index
        for index, name in enumerate(calendar.month_name)
        if name
    }
    month_lookup.update(
        {
            abbreviation.casefold(): index
            for index, abbreviation in enumerate(MONTH_NAMES, start=1)
        }
    )
    normalized_months = observed["month"].astype(str).str.strip().str.casefold()
    invalid_months = ~normalized_months.isin(month_lookup)
    if invalid_months.any():
        raise ValueError(
            "Observed incidence CSV has unrecognized month values: "
            f"{sorted(normalized_months.loc[invalid_months].unique())}"
        )
    observed["month"] = normalized_months.map(month_lookup).astype(int)

    monthly = (
        observed.groupby(["year", "month"], as_index=False)
        .agg(
            population=("Pop_5", "sum"),
            confirmed_cases=("conf_u5", "sum"),
        )
        .sort_values(["year", "month"])
    )
    month_counts = monthly.groupby("year")["month"].nunique()
    if not month_counts.eq(12).all():
        incomplete_years = month_counts.loc[~month_counts.eq(12)].index.tolist()
        raise ValueError(
            "Observed incidence data must include all 12 months for each year; "
            f"incomplete years: {incomplete_years}"
        )

    month_start = pd.to_datetime(
        dict(year=monthly["year"], month=monthly["month"], day=1)
    )
    days_in_month = month_start.dt.days_in_month
    monthly["Incidence_per_1000_person_years"] = (
        monthly["confirmed_cases"]
        / monthly["population"]
        * (365 / days_in_month)
        * 1000
    )
    monthly["date"] = month_start + pd.to_timedelta(14, unit="D")
    return monthly


def plot_observed_vs_simulated_incidence(data, observed_file, output_path):
    observed = load_observed_incidence(observed_file)
    simulated_by_run = (
        data.groupby(["Run_Number", "Year", "month"], as_index=False)["Cases"]
        .mean()
    )
    simulated = (
        simulated_by_run.groupby(["Year", "month"], as_index=False)["Cases"]
        .mean()
    )
    simulated["Cases"] *= 1000 * SIMULATED_INCIDENCE_DISPLAY_FACTOR
    simulated["date"] = pd.to_datetime(
        dict(year=simulated["Year"], month=simulated["month"], day=15)
    )

    comparison = simulated.merge(
        observed[["year", "month", "Incidence_per_1000_person_years"]],
        left_on=["Year", "month"],
        right_on=["year", "month"],
        how="inner",
        validate="one_to_one",
    )
    if comparison.empty:
        raise ValueError(
            "Observed incidence and simulated data have no overlapping "
            "year-months. Check --start-year, --end-year, and the input directory."
        )

    fig, axis = plt.subplots(figsize=(12, 6))
    axis.plot(
        comparison["date"],
        comparison["Incidence_per_1000_person_years"],
        color="#D68910",
        marker="o",
        linewidth=1.6,
        markersize=3.5,
        label="Observed confirmed cases",
        zorder=3,
    )
    axis.plot(
        comparison["date"],
        comparison["Cases"],
        color="#2471A3",
        linewidth=1.8,
        label=(
            "Simulated clinical incidence (mean across runs, "
            "displayed at 0.1×)"
        ),
    )
    axis.set_title("Observed and simulated under-five malaria incidence")
    axis.set_ylabel(
        "Annualized incidence per 1,000 under-five person-years"
    )
    axis.set_xlabel("Month")
    axis.grid(alpha=0.2)
    axis.legend()
    fig.tight_layout()
    fig.savefig(output_path, dpi=150)
    plt.close(fig)
    print(
        "Observed/simulated incidence comparison uses overlapping months "
        f"{comparison['Year'].min()}-{comparison['Year'].max()}."
    )


def plot_observed_vs_simulated_seasonality(data, observed_file, output_path):
    observed = load_observed_incidence(observed_file)
    observed_profile = (
        observed.groupby("month", as_index=False)[
            "Incidence_per_1000_person_years"
        ]
        .mean()
        .sort_values("month")
        .set_index("month")["Incidence_per_1000_person_years"]
        .reindex(range(1, 13))
    )

    simulated = (
        data.groupby(["Run_Number", "Year", "month"], as_index=False)["Cases"]
        .mean()
    )
    observed_year_month = observed[["year", "month"]].drop_duplicates()
    simulated = simulated.merge(
        observed_year_month,
        left_on=["Year", "month"],
        right_on=["year", "month"],
        how="inner",
        validate="many_to_one",
    )
    if simulated.empty:
        raise ValueError(
            "Observed and simulated incidence have no overlapping calendar months "
            "for the seasonal profile."
        )

    simulated_profile = (
        simulated.groupby(["Year", "month"], as_index=False)["Cases"]
        .mean()
        .groupby("month")["Cases"]
        .mean()
        .reindex(range(1, 13))
        * 1000
        * SIMULATED_INCIDENCE_DISPLAY_FACTOR
    )

    months = np.arange(1, 13)
    fig, axis = plt.subplots(figsize=(10, 6))
    axis.plot(
        months,
        observed_profile.to_numpy(),
        color="#D68910",
        marker="o",
        linewidth=2,
        label="Observed confirmed cases",
    )
    axis.plot(
        months,
        simulated_profile.to_numpy(),
        color="#2471A3",
        marker="o",
        linewidth=2,
        label="Simulated clinical incidence (mean, displayed at 0.1×)",
    )
    axis.set_title("Observed and simulated under-five incidence seasonality")
    axis.set_ylabel(
        "Annualized incidence per 1,000 under-five person-years "
    )
    axis.set_xlabel("Month")
    axis.set_xticks(months)
    axis.set_xticklabels(MONTH_NAMES)
    axis.grid(alpha=0.2)
    axis.legend()
    fig.tight_layout()
    fig.savefig(output_path, dpi=150)
    plt.close(fig)


def prepare_data(input_dir, args):
    data = load_data(input_dir)
    filters = args.filter or [
        f"{column}={value}" for column, value in PARAMETER_FILTERS.items()
    ]
    data = apply_filters(data, filters)
    require_one_parameter_set(data)

    if args.start_year is not None:
        data = data.loc[data["Year"] >= args.start_year].copy()
    if args.end_year is not None:
        data = data.loc[data["Year"] <= args.end_year].copy()
    if data.empty:
        raise ValueError("The selected year range contains no data.")

    data = summarize_age_group(data, args.age_group)
    return select_complete_years(data)


def align_scenario_pairs(data, compare_data):
    data_pairs = set(
        data[["Run_Number", "Year"]].drop_duplicates().itertuples(index=False, name=None)
    )
    compare_pairs = set(
        compare_data[["Run_Number", "Year"]]
        .drop_duplicates()
        .itertuples(index=False, name=None)
    )
    common_pairs = data_pairs & compare_pairs
    if not common_pairs:
        raise ValueError(
            "The two scenarios have no matching complete Run_Number/year pairs."
        )
    if common_pairs != data_pairs or common_pairs != compare_pairs:
        print(
            "Scenario comparison will use only matched Run_Number/year pairs "
            f"({len(common_pairs)} pairs)."
        )

    pair_frame = pd.DataFrame(
        sorted(common_pairs), columns=["Run_Number", "Year"]
    )
    data = data.merge(pair_frame, on=["Run_Number", "Year"], how="inner")
    compare_data = compare_data.merge(
        pair_frame, on=["Run_Number", "Year"], how="inner"
    )
    return data, compare_data


def mean_seasonal_profile(data, metric):
    run_month = (
        data.groupby(["Run_Number", "Year", "month"], as_index=False)[metric]
        .mean()
    )
    year_month = run_month.groupby(["Year", "month"], as_index=False)[metric].mean()
    return year_month.groupby("month")[metric].mean().reindex(range(1, 13))


def scenario_label(input_dir, explicit_label):
    if explicit_label:
        return explicit_label
    directory_name = input_dir.name.lower()
    if "non_smc" in directory_name or "no_smc" in directory_name:
        return "Without SMC"
    if "smc" in directory_name:
        return "With SMC"
    return input_dir.name


def plot_scenario_comparison(
    data,
    compare_data,
    metric,
    ylabel,
    title,
    input_label,
    compare_label,
    output_path,
):
    first_profile = mean_seasonal_profile(data, metric)
    second_profile = mean_seasonal_profile(compare_data, metric)
    months = np.arange(1, 13)
    scale = (
        100
        if metric == "PfPR"
        else SIMULATED_INCIDENCE_DISPLAY_FACTOR
    )
    incidence_display_note = (
        " (simulated values displayed at 0.1×)" if metric == "Cases" else ""
    )

    fig, axis = plt.subplots(figsize=(10, 6))
    axis.plot(
        months,
        first_profile.to_numpy() * scale,
        color="#2471A3",
        marker="o",
        linewidth=2,
        label=f"{input_label}{incidence_display_note}",
    )
    axis.plot(
        months,
        second_profile.to_numpy() * scale,
        color="#C0392B",
        marker="o",
        linewidth=2,
        label=f"{compare_label}{incidence_display_note}",
    )
    axis.set_title(title)
    axis.set_ylabel(ylabel)
    axis.set_xlabel("Month")
    axis.set_xticks(months)
    axis.set_xticklabels(MONTH_NAMES)
    axis.legend()
    axis.grid(alpha=0.2)
    fig.tight_layout()
    fig.savefig(output_path, dpi=150)
    plt.close(fig)


def main():
    args = parse_args()
    data = prepare_data(args.input_dir, args)
    output_dir = args.output_dir or args.input_dir
    output_dir.mkdir(parents=True, exist_ok=True)
    age_label = "Under 5" if args.age_group == "under5" else "Age 5+"

    plot_metric(
        data,
        "PfPR",
        "PfPR (%)",
        f"Simulated PfPR seasonality - {age_label}",
        output_dir / METRICS["PfPR"],
    )
    plot_metric(
        data,
        "Cases",
        "Annualized clinical incidence (episodes per person-year; displayed at 0.1×)",
        f"Simulated clinical-incidence seasonality - {age_label}",
        output_dir / METRICS["Cases"],
    )

    if args.age_group == "under5":
        plot_observed_vs_simulated_incidence(
            data,
            args.observed_incidence_file,
            output_dir / "observed_vs_simulated_incidence.png",
        )
        plot_observed_vs_simulated_seasonality(
            data,
            args.observed_incidence_file,
            output_dir / "observed_vs_simulated_seasonality.png",
        )
    else:
        print(
            "Skipping observed incidence comparison: the observed file is "
            "under-five, but --age-group is age5plus."
        )

    if args.compare_dir:
        compare_data = prepare_data(args.compare_dir, args)
        data, compare_data = align_scenario_pairs(data, compare_data)
        input_label = scenario_label(args.input_dir, args.input_label)
        compare_label = scenario_label(args.compare_dir, args.compare_label)
        plot_scenario_comparison(
            data,
            compare_data,
            "PfPR",
            "PfPR (%)",
            f"Mean PfPR seasonal profile - {age_label}",
            input_label,
            compare_label,
            output_dir / "simulated_pfpr_scenario_comparison.png",
        )
        plot_scenario_comparison(
            data,
            compare_data,
            "Cases",
            "Annualized clinical incidence (episodes per person-year; displayed at 0.1×)",
            f"Mean clinical-incidence seasonal profile - {age_label}",
            input_label,
            compare_label,
            output_dir / "simulated_incidence_scenario_comparison.png",
        )
        print(f"Saved scenario comparison plots to {output_dir}")

    print(f"Saved PfPR and clinical-incidence plots to {output_dir}")
    print(
        "Clinical incidence is plotted from the analyzer's 'Cases' column "
        "(Annual Clinical Incidence by Age Bin). The observed monthly case "
        "counts are population-denominated and annualized to the same per-year scale."
    )
    print(
        "All simulated incidence plots display simulated values at "
        f"{SIMULATED_INCIDENCE_DISPLAY_FACTOR:g}× for visual inspection only; "
        "the analyzer CSV and observed incidence values are unchanged."
    )


if __name__ == "__main__":
    main()
