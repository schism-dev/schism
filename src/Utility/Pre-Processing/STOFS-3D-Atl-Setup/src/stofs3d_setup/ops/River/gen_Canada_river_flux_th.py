import os
import io
import requests
import pandas as pd
from datetime import datetime, timezone
import pytz  # pip install pytz


def gen_Canada_river_flux_th(
    start_datetime=datetime(2018, 12, 1, tzinfo=pytz.UTC),
    end_datetime=datetime(2020, 1, 2, tzinfo=pytz.UTC)
):
    """
    Generate Canada river flux_th file for SCHISM from Environment Canada daily discharge data.

    The output file is named "02OA016_daily_discharge_YYYY-MM-DD_YYYY-MM-DD.csv" where the dates
    correspond to the start and end dates of the data retrieved.

    The discharge values are negated to convert from inflow to outflow as per SCHISM convention.

    Parameters
    ----------
    start_datetime : datetime, optional
        Start datetime for data retrieval (inclusive). Default is 2022-12-01 UTC.
    end_datetime : datetime, optional
        End datetime for data retrieval (inclusive). Default is 2024-01-02 UTC
    """

    STATION = "02OA016"
    start_str = start_datetime.strftime("%Y-%m-%d")
    end_str = end_datetime.strftime("%Y-%m-%d")
    EPOCH_UTC = datetime(start_datetime.year, start_datetime.month, start_datetime.day, tzinfo=timezone.utc)

    date_col = "Date"
    val_col = "Value/Valeur"

    url = (
        "https://wateroffice.ec.gc.ca/services/daily_data/csv/inline"
        f"?stations[]={STATION}&parameters[]=flow&start_date={start_str}&end_date={end_str}"
    )

    r = requests.get(url, timeout=60)
    r.raise_for_status()

    # Parse CSV. The file may include a header block; pandas handles it well.
    df = pd.read_csv(io.StringIO(r.text))

    # Parse dates as naive midnights, then treat as UTC midnights (OK for daily means)
    ts_utc = pd.to_datetime(df[date_col], format="%Y-%m-%d", errors="coerce").dt.tz_localize("UTC")

    # Convert to seconds since epoch
    secs = (ts_utc - EPOCH_UTC).dt.total_seconds().astype("int64")

    # Build output (drop NaNs)
    q = pd.to_numeric(df[val_col], errors="coerce")
    out = pd.DataFrame({
        f"time_seconds_since_{start_str}T00:00:00Z": secs,
        "discharge_m3s": -q  # negate to follow SCHISM convention (negative for inflow)
    }).dropna()

    # Save
    out_name = f"{STATION}_daily_discharge_{start_str}_{end_str}.csv"
    out.to_csv(out_name, index=False, sep=' ', header=False)
    os.system(f"ln -s {out_name} flux.th")
    print(out.head())
    print("Done")


if __name__ == "__main__":
    year = 2017
    gen_Canada_river_flux_th(
        start_datetime=datetime(year-1, 12, 1, tzinfo=pytz.UTC),
        end_datetime=datetime(year+1, 1, 2, tzinfo=pytz.UTC)
    )
