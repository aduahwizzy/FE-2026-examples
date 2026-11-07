"""
Build EMOD climate/weather files (.bin + .bin.json) from a requested NASA POWER
daily date range.

Each region is an independent single-node model (node id 1), so each gets its own
weather set in climate/emod_weather_<region>/.

Run climate_fetch.py first with the same date range to produce the daily CSV.
Dates must cover one or more complete 365-day windows after leap days are
removed, as expected by EMOD's daily climate indexing.

Example:
    python build_weather.py --start-date 2019-01-01 --end-date 2023-12-31
"""
import argparse
from datetime import date
import os
import json
import pandas as pd

from emodpy_malaria.weather import csv_to_weather, WeatherVariable, WeatherAttributes

REGIONS = {
    "northern_region": dict(lat=9.462781, lon=-0.285867),
}

NODE_ID = 1
DEFAULT_START_DATE = date(1965, 1, 1)
DEFAULT_END_DATE = date(2025, 12, 31)
DEMOG_ID_REF = "Gridded world grump2.5arcmin"    # must match demographics IdReference
CLIM_DIR = os.path.join(os.path.dirname(__file__), "climate")

FILE_NAMES = {
    WeatherVariable.AIR_TEMPERATURE: "air_temperature.bin",
    WeatherVariable.RELATIVE_HUMIDITY: "relative_humidity.bin",
    WeatherVariable.RAINFALL: "rainfall.bin",
    WeatherVariable.LAND_TEMPERATURE: "land_temperature.bin",
}


def build_region(region, lat, lon, start_date, end_date):
    date_key = f"{start_date:%Y%m%d}_{end_date:%Y%m%d}"
    raw_path = os.path.join(CLIM_DIR, f"{region}_{date_key}_daily_raw.csv")
    daily = pd.read_csv(raw_path, index_col=0, parse_dates=True).sort_index()

    expected_dates = pd.date_range(start=start_date, end=end_date, freq="D")
    expected_dates = expected_dates[~((expected_dates.month == 2) &
                                      (expected_dates.day == 29))]
    if len(expected_dates) < 365 or len(expected_dates) % 365 != 0:
        raise ValueError(
            f"{region}: requested range must contain a positive whole number "
            "of 365-day years after Feb 29 is removed; "
            f"found {len(expected_dates)} days")

    daily = daily.reindex(expected_dates)
    if daily[["airtemp_C", "rh_pct", "rain_mm"]].isna().any().any():
        missing_dates = daily.index[
            daily[["airtemp_C", "rh_pct", "rain_mm"]].isna().any(axis=1)
        ]
        raise ValueError(
            f"{region}: climate data is missing values for "
            f"{', '.join(str(day.date()) for day in missing_dates[:5])}")

    airtemp = daily["airtemp_C"].to_numpy()
    humidity = (daily["rh_pct"] / 100.0).clip(0, 1).to_numpy()
    rainfall = daily["rain_mm"].clip(lower=0).to_numpy()
    n = len(daily)
    df = pd.DataFrame({
        "nodes": NODE_ID, "steps": range(n),
        "airtemp": airtemp, "humidity": humidity,
        "rainfall": rainfall, "landtemp": airtemp,
    })

    attrs = WeatherAttributes(
        reference=DEMOG_ID_REF,
        start_year=start_date.year, end_year=end_date.year,
        start_doy=start_date.timetuple().tm_yday,
        lat_min=lat, lat_max=lat, lon_min=lon, lon_max=lon,
        provenance=f"NASA POWER daily {start_date} through {end_date}, {region}",
        resolution="daily", update_freq="CLIMATE_UPDATE_DAY",
    )

    wdir = os.path.join(CLIM_DIR, f"emod_weather_{region}")
    os.makedirs(wdir, exist_ok=True)
    csv_to_weather(
        csv_data=df, node_column="nodes", step_column="steps",
        weather_columns={WeatherVariable.AIR_TEMPERATURE: "airtemp",
                         WeatherVariable.RELATIVE_HUMIDITY: "humidity",
                         WeatherVariable.RAINFALL: "rainfall",
                         WeatherVariable.LAND_TEMPERATURE: "landtemp"},
        attributes=attrs, weather_dir=wdir, weather_file_names=FILE_NAMES,
    )
    # strip WeatherSchemaVersion so this Eradication uses the legacy parser
    for f in os.listdir(wdir):
        if f.endswith(".bin.json"):
            p = os.path.join(wdir, f)
            meta = json.load(open(p))
            meta.get("Metadata", {}).pop("WeatherSchemaVersion", None)
            json.dump(meta, open(p, "w"), indent=4)
    return wdir, n


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--start-date", type=date.fromisoformat,
                        default=DEFAULT_START_DATE, metavar="YYYY-MM-DD")
    parser.add_argument("--end-date", type=date.fromisoformat,
                        default=DEFAULT_END_DATE, metavar="YYYY-MM-DD")
    args = parser.parse_args()
    if args.start_date > args.end_date:
        parser.error("--start-date must be on or before --end-date")

    for region, coords in REGIONS.items():
        wdir, n = build_region(region, coords["lat"], coords["lon"],
                               args.start_date, args.end_date)
        print(f"{region:12s} -> {wdir}  ({n} days, node {NODE_ID}, "
              f"{args.start_date} through {args.end_date})")


if __name__ == "__main__":
    main()
