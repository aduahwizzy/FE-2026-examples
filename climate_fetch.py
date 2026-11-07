"""
Fetch daily climate for the northern-Ghana study regions from NASA POWER and
save the requested daily record plus a day-of-year climatology per region.

Regions (representative site = regional capital / research site):
  upper_east -> Navrongo (10.90 N, 1.09 W)

NASA POWER daily point API (no key): T2M (air temp C), RH2M (rel humidity %),
PRECTOTCORR (bias-corrected precip mm/day).

Outputs per region in climate/:
    <region>_<start>_<end>_daily_raw.csv   - requested daily record
    <region>_<start>_<end>_climatology.csv - day-of-year means

Dates use ISO format (YYYY-MM-DD). For example:
    python climate_fetch.py --start-date 1965-01-01 --end-date 2025-12-31
"""
import argparse
import os
from datetime import date

import numpy as np
import pandas as pd
import requests

SITES = {
    "northern_region": dict(name="tamale", lat=9.462781, lon=-0.285867),
}

DEFAULT_START_DATE = date(1965, 1, 1)
DEFAULT_END_DATE = date(2025, 12, 31)
OUTDIR = os.path.join(os.path.dirname(__file__), "climate")
PARAMS = ["T2M", "RH2M", "PRECTOTCORR"]


def fetch(lat, lon, start_date, end_date):
    url = ("https://power.larc.nasa.gov/api/temporal/daily/point"
           f"?parameters={','.join(PARAMS)}&community=AG"
           f"&longitude={lon}&latitude={lat}"
           f"&start={start_date:%Y%m%d}&end={end_date:%Y%m%d}&format=JSON")
    r = requests.get(url, timeout=120)
    r.raise_for_status()
    p = r.json()["properties"]["parameter"]
    df = pd.DataFrame({k: p[k] for k in PARAMS})
    df.index = pd.to_datetime(df.index, format="%Y%m%d")
    df = df.replace(-999, np.nan).dropna()
    return df.rename(columns={"T2M": "airtemp_C", "RH2M": "rh_pct", "PRECTOTCORR": "rain_mm"})


def climatology(df):
    d = df.copy()
    d["doy"] = d.index.dayofyear
    d = d[d["doy"] <= 365]
    clim = d.groupby("doy").mean(numeric_only=True).reset_index()
    clim["airtemp"] = clim["airtemp_C"]
    clim["humidity"] = (clim["rh_pct"] / 100.0).clip(0, 1)
    clim["rainfall"] = clim["rain_mm"].clip(lower=0)
    return clim[["doy", "airtemp", "humidity", "rainfall"]]


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--start-date", type=date.fromisoformat,
                        default=DEFAULT_START_DATE, metavar="YYYY-MM-DD")
    parser.add_argument("--end-date", type=date.fromisoformat,
                        default=DEFAULT_END_DATE, metavar="YYYY-MM-DD")
    args = parser.parse_args()
    if args.start_date > args.end_date:
        parser.error("--start-date must be on or before --end-date")

    os.makedirs(OUTDIR, exist_ok=True)
    date_key = f"{args.start_date:%Y%m%d}_{args.end_date:%Y%m%d}"
    print(f"Fetching {args.start_date} through {args.end_date}")
    print(f"{'region':12s} {'site':10s} {'ann.rain':>9s} {'T mean':>7s} {'RH min':>7s}")
    for region, s in SITES.items():
        df = fetch(s["lat"], s["lon"], args.start_date, args.end_date)
        df.to_csv(os.path.join(OUTDIR, f"{region}_{date_key}_daily_raw.csv"))
        clim = climatology(df)
        clim.to_csv(os.path.join(OUTDIR, f"{region}_{date_key}_climatology.csv"), index=False)
        print(f"{region:12s} {s['name']:10s} {clim.rainfall.sum():8.0f}mm "
              f"{clim.airtemp.mean():6.1f}C {clim.humidity.min():6.2f}")
    print("\nSaved date-ranged CSVs to", OUTDIR)


if __name__ == "__main__":
    main()
