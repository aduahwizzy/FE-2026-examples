"""Prefix 1981-2025 daily climate data with a 1981-1990 daily climatology.

The generated file covers 1963-2025 and is suitable for build_weather.py.
February 29 is omitted so every model year contains 365 days.
"""
import argparse
import csv
from collections import defaultdict
from datetime import date, timedelta
from pathlib import Path


CLIMATE_DIR = Path(__file__).resolve().parent / "climate"
SOURCE_START = date(1981, 1, 1)
SOURCE_END = date(2025, 12, 31)
EXTENDED_START = date(1963, 1, 1)
CLIMATOLOGY_END_YEAR = 1990
VALUE_COLUMNS = ("airtemp_C", "rh_pct", "rain_mm")


def read_daily_csv(path):
    with path.open(newline="") as climate_file:
        reader = csv.DictReader(climate_file)
        if not reader.fieldnames or len(reader.fieldnames) < 4:
            raise ValueError(f"{path} does not have a date column and 3 values.")

        date_column = reader.fieldnames[0]
        missing_columns = set(VALUE_COLUMNS).difference(reader.fieldnames)
        if missing_columns:
            raise ValueError(
                f"{path} is missing columns: {', '.join(sorted(missing_columns))}"
            )

        records = {}
        for line_number, row in enumerate(reader, start=2):
            try:
                day = date.fromisoformat(row[date_column])
                values = {
                    column: float(row[column]) for column in VALUE_COLUMNS
                }
            except (TypeError, ValueError) as error:
                raise ValueError(
                    f"Invalid climate data at {path}:{line_number}"
                ) from error
            if day in records:
                raise ValueError(f"Duplicate climate date {day} in {path}.")
            records[day] = values
    return records


def build_climatology(source_records):
    values_by_month_day = defaultdict(lambda: {column: [] for column in VALUE_COLUMNS})
    for day, values in source_records.items():
        if SOURCE_START.year <= day.year <= CLIMATOLOGY_END_YEAR:
            if (day.month, day.day) != (2, 29):
                for column in VALUE_COLUMNS:
                    values_by_month_day[(day.month, day.day)][column].append(
                        values[column]
                    )

    climatology = {}
    for month_day, values in values_by_month_day.items():
        counts = {len(series) for series in values.values()}
        if counts != {10}:
            raise ValueError(
                f"Expected 10 values per year for {month_day}, found {counts}."
            )
        climatology[month_day] = {
            column: sum(series) / len(series)
            for column, series in values.items()
        }

    if len(climatology) != 365:
        raise ValueError(
            f"Expected 365 climatology days from 1981-1990, "
            f"found {len(climatology)}."
        )
    return climatology


def extend_climate(source_path, output_path):
    source_records = read_daily_csv(source_path)
    climatology = build_climatology(source_records)

    output_path.parent.mkdir(parents=True, exist_ok=True)
    expected_source_days = set()
    day = SOURCE_START
    while day <= SOURCE_END:
        if (day.month, day.day) != (2, 29):
            expected_source_days.add(day)
        day += timedelta(days=1)

    missing_days = expected_source_days.difference(source_records)
    if missing_days:
        first_missing = min(missing_days)
        raise ValueError(
            f"Source climate data is missing {len(missing_days)} dates, "
            f"starting with {first_missing}."
        )

    with output_path.open("w", newline="") as climate_file:
        writer = csv.writer(climate_file)
        writer.writerow(("date", *VALUE_COLUMNS))
        day = EXTENDED_START
        while day <= SOURCE_END:
            if (day.month, day.day) == (2, 29):
                day += timedelta(days=1)
                continue

            if day < SOURCE_START:
                values = climatology[(day.month, day.day)]
            else:
                values = source_records[day]
            writer.writerow(
                (day.isoformat(), *(f"{values[column]:.10g}" for column in VALUE_COLUMNS))
            )
            day += timedelta(days=1)

    return len(expected_source_days), output_path


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--input",
        type=Path,
        default=CLIMATE_DIR / "northern_region_19810101_20251231_daily_raw.csv",
        help="Fetched 1981-2025 daily climate CSV.",
    )
    parser.add_argument(
        "--output",
        type=Path,
        default=CLIMATE_DIR / "northern_region_19630101_20251231_daily_raw.csv",
        help="Output CSV for the 1963-2025 daily climate series.",
    )
    args = parser.parse_args()

    source_days, output_path = extend_climate(args.input, args.output)
    total_days = (SOURCE_END.year - EXTENDED_START.year + 1) * 365
    print(
        f"Wrote {total_days} daily rows to {output_path}; "
        f"preserved {source_days} observed days from 1981-2025 and "
        f"filled 1963-1980 using the 1981-1990 day-of-year climatology."
    )


if __name__ == "__main__":
    main()
