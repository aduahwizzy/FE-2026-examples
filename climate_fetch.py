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

Set START_DATE and END_DATE in this script before running it. NASA POWER daily
data starts on 1981-01-01.
"""
import os
from datetime import date

import numpy as np
import pandas as pd
import requests

SITES = {
    "Osun": dict(name="Osun", lat=7.534717, lon=-4.536133),
}

START_DATE = date(1981, 1, 1)
END_DATE = date(2025, 12, 31)
MIN_START_DATE = date(1981, 1, 1)
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
    if START_DATE > END_DATE:
        raise ValueError("START_DATE must be on or before END_DATE")
    if START_DATE < MIN_START_DATE:
        raise ValueError(
            f"NASA POWER daily data starts on {MIN_START_DATE}; "
            "choose a START_DATE on or after that date"
        )

    os.makedirs(OUTDIR, exist_ok=True)
    date_key = f"{START_DATE:%Y%m%d}_{END_DATE:%Y%m%d}"
    print(f"Fetching {START_DATE} through {END_DATE}")
    print(f"{'region':12s} {'site':10s} {'ann.rain':>9s} {'T mean':>7s} {'RH min':>7s}")
    for region, s in SITES.items():
        df = fetch(s["lat"], s["lon"], START_DATE, END_DATE)
        df.to_csv(os.path.join(OUTDIR, f"{region}_{date_key}_daily_raw.csv"))
        clim = climatology(df)
        clim.to_csv(os.path.join(OUTDIR, f"{region}_{date_key}_climatology.csv"), index=False)
        print(f"{region:12s} {s['name']:10s} {clim.rainfall.sum():8.0f}mm "
              f"{clim.airtemp.mean():6.1f}C {clim.humidity.min():6.2f}")
    print("\nSaved date-ranged CSVs to", OUTDIR)


if __name__ == "__main__":
    main()
