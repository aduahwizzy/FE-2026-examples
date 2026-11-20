#!/usr/bin/env python3
"""Extend Northern Region climate through 2030 using the 2016-2025 daily mean."""

import csv
from datetime import date, timedelta
from pathlib import Path

import pandas as pd

from build_weather import build_region


ROOT = Path(__file__).resolve().parent
CLIMATE_DIR = ROOT / "climate"
HISTORICAL_PATH = CLIMATE_DIR / "northern_region_19630101_20251231_daily_raw.csv"
OBSERVED_PATH = CLIMATE_DIR / "northern_region_19810101_20251231_daily_raw.csv"
PROJECTED_PATH = CLIMATE_DIR / "northern_region_19630101_20301231_daily_raw.csv"
WEATHER_DIR = CLIMATE_DIR / "emod_weather_northern_region_2016_2025_climatology"
CLIMATOLOGY_START_YEAR = 2016
CLIMATOLOGY_END_YEAR = 2025
FORECAST_START_YEAR = 2026
FORECAST_END_YEAR = 2030
REGION = "northern_region"
LATITUDE = 9.462781
LONGITUDE = -0.285867
VALUE_COLUMNS = ("airtemp_C", "rh_pct", "rain_mm")


def build_projection():
    history = pd.read_csv(HISTORICAL_PATH, parse_dates=["date"])
    observed = (
        pd.read_csv(OBSERVED_PATH, index_col=0, parse_dates=True)
        .rename_axis("date")
        .reset_index()
    )
    if history["date"].duplicated().any() or observed["date"].duplicated().any():
        raise ValueError("Climate source CSVs contain duplicate dates.")
    if history[list(VALUE_COLUMNS)].isna().any().any():
        raise ValueError("Historical climate contains missing values.")
    if observed[list(VALUE_COLUMNS)].isna().any().any():
        raise ValueError("Observed climate contains missing values.")

    climatology_source = observed.loc[
        observed["date"].dt.year.between(
            CLIMATOLOGY_START_YEAR, CLIMATOLOGY_END_YEAR
        )
        & ~((observed["date"].dt.month == 2) & (observed["date"].dt.day == 29))
    ].copy()
    climatology_source["month_day"] = climatology_source["date"].dt.strftime("%m-%d")
    counts = climatology_source.groupby("month_day")["date"].nunique()
    if len(counts) != 365 or not counts.eq(10).all():
        raise ValueError(
            "The 2016-2025 climate source must provide 10 observations for "
            "each of the 365 non-leap month/day values."
        )
    climatology = (
        climatology_source.groupby("month_day")[list(VALUE_COLUMNS)]
        .mean()
        .to_dict(orient="index")
    )

    history["date"] = history["date"].dt.date
    expected_history = []
    current = date(1963, 1, 1)
    history_end = date(2025, 12, 31)
    while current <= history_end:
        if (current.month, current.day) != (2, 29):
            expected_history.append(current)
        current += timedelta(days=1)
    actual_history = set(history["date"])
    missing = set(expected_history).difference(actual_history)
    if missing:
        raise ValueError(
            f"Historical climate is missing {len(missing)} dates; "
            f"first missing date is {min(missing)}."
        )

    projected_rows = history.loc[
        history["date"].between(date(1963, 1, 1), history_end),
        ["date", *VALUE_COLUMNS],
    ].to_dict(orient="records")
    current = date(FORECAST_START_YEAR, 1, 1)
    forecast_end = date(FORECAST_END_YEAR, 12, 31)
    while current <= forecast_end:
        if (current.month, current.day) != (2, 29):
            values = climatology[current.strftime("%m-%d")]
            projected_rows.append({"date": current, **values})
        current += timedelta(days=1)

    PROJECTED_PATH.parent.mkdir(parents=True, exist_ok=True)
    with PROJECTED_PATH.open("w", newline="") as climate_file:
        writer = csv.DictWriter(climate_file, fieldnames=["date", *VALUE_COLUMNS])
        writer.writeheader()
        for row in projected_rows:
            writer.writerow(
                {
                    "date": row["date"].isoformat(),
                    **{column: f"{row[column]:.10g}" for column in VALUE_COLUMNS},
                }
            )

    expected_rows = (2030 - 1963 + 1) * 365
    if len(projected_rows) != expected_rows:
        raise ValueError(
            f"Expected {expected_rows} projected climate rows, "
            f"wrote {len(projected_rows)}."
        )

    weather_dir, _ = build_region(
        REGION,
        LATITUDE,
        LONGITUDE,
        date(1963, 1, 1),
        date(2030, 12, 31),
        output_dir=str(WEATHER_DIR),
        provenance=(
            "NASA POWER historical daily climate through 2025; "
            "2026-2030 projected from the 2016-2025 daily climatology, "
            "Northern Region Ghana"
        ),
    )
    print(f"Wrote projected daily climate CSV: {PROJECTED_PATH}")
    print(
        f"Built EMOD weather assets in {weather_dir}; 2026-2030 use the "
        "2016-2025 day-of-year mean."
    )


if __name__ == "__main__":
    build_projection()
