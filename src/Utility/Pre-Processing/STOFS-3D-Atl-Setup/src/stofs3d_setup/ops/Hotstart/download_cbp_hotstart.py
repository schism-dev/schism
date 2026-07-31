"""
Small Chesapeake Bay Program downloader for STOFS-3D-Atl hotstart files.

This uses the current CBP DataHub WaterQuality API and writes the legacy
hotstart observation format expected by gen_subregion_ic_stofs3d:

    lon lat value station_name
"""

from __future__ import annotations

import datetime
from pathlib import Path
from typing import Any

import pandas as pd
import requests


CBP_WATER_QUALITY_URL = "https://datahub.chesapeakebay.net/api/WaterQuality"

CBP_VARIABLE_IDS = {
    "SALINITY": 83,
    "WTEMP": 123,
}

CBP_VAR_NAMES = {
    "sal": "SALINITY",
    "tem": "WTEMP",
}

CBP_STATION_IDS = {
    "CB1.1": "1149",
    "CB2.1": "1150",
    "CB2.2": "1151",
    "CB3.1": "1152",
    "CB3.2": "1153",
    "CB3.3C": "1154",
    "CB4.1C": "1160",
    "CB4.2C": "1163",
    "CB4.3C": "1166",
    "CB4.4": "1169",
    "CB5.1": "1170",
    "CB5.2": "1172",
    "CB5.3": "1173",
    "CB5.4": "1174",
    "CB5.5": "1176",
    "CB7.1S": "1183",
    "CB7.2": "1184",
    "CB7.3": "1186",
    "CB7.4": "1188",
    "LE5.1": "1331",
    "LE5.2": "1332",
    "LE5.3": "1335",
    "LE5.4": "1336",
    "LE5.5-W": "1338",
    "LE5.6": "1341",
}

DEFAULT_CBP_STATIONS = list(CBP_STATION_IDS)


class GenericObsData:
    __slots__ = ["station_info", "df"]

    def __init__(self, station_info: dict[str, Any] | None = None, df: pd.DataFrame | None = None):
        self.station_info = station_info
        self.df = df


def _request_cbp_json(url: str) -> Any:
    headers = {
        "User-Agent": "stofs3d_setup hotstart CBP downloader",
        "Accept": "application/json,text/plain,*/*",
    }
    response = requests.get(url, headers=headers, timeout=300, allow_redirects=True)
    print(f"[INFO] CBP response status: {response.status_code}", flush=True)
    print(f"[INFO] CBP final URL: {response.url}", flush=True)
    print(f"[INFO] CBP response length: {len(response.content)} bytes", flush=True)
    response.raise_for_status()
    try:
        return response.json()
    except ValueError as exc:
        print("[ERROR] CBP response was not valid JSON.", flush=True)
        print(response.text[:1000], flush=True)
        raise RuntimeError("CBP DataHub returned a non-JSON response.") from exc


def get_cbp(
    stations: list[str] | None = None,
    sample_time: str | pd.Timestamp | None = None,
    varname: str = "SALINITY",
) -> dict[str, GenericObsData]:
    if stations is None:
        stations = DEFAULT_CBP_STATIONS
    if sample_time is None:
        sample_time = "2011-08-24"
    if varname not in CBP_VARIABLE_IDS:
        raise ValueError(f"Unsupported CBP variable: {varname}. Available variables are {list(CBP_VARIABLE_IDS)}.")

    sample_time = pd.to_datetime(sample_time, format="%Y-%m-%d", errors="raise")
    search_window = [
        (sample_time - datetime.timedelta(days=360)).strftime("%m-%d-%Y"),
        (sample_time + datetime.timedelta(days=50)).strftime("%m-%d-%Y"),
    ]

    missing_stations = [station for station in stations if station not in CBP_STATION_IDS]
    if missing_stations:
        raise KeyError("CBP station IDs are not defined for: " + ", ".join(missing_stations))

    station_ids = [CBP_STATION_IDS[station] for station in stations]
    var_id = CBP_VARIABLE_IDS[varname]
    url = (
        f"{CBP_WATER_QUALITY_URL}/{search_window[0]}/{search_window[1]}"
        f"/0,1/2,6,4/35,13,12,34,7,23,3,2/Station/{','.join(station_ids)}/{var_id}"
    )

    print("[INFO] Requesting CBP WaterQuality data", flush=True)
    print(f"[INFO] Variable: {varname}", flush=True)
    print(f"[INFO] Target time: {sample_time}", flush=True)
    print(f"[INFO] Search window: {search_window[0]} to {search_window[1]}", flush=True)
    print(f"[INFO] Number of requested stations: {len(stations)}", flush=True)

    data_json = _request_cbp_json(url)
    if isinstance(data_json, dict):
        for key in ["data", "results", "value"]:
            if key in data_json:
                data_json = data_json[key]
                break

    df_all = pd.DataFrame(data_json)
    if df_all.empty:
        raise RuntimeError(f"CBP returned no {varname} data for {search_window[0]} to {search_window[1]}.")

    df_all = df_all.rename(columns={
        "station": "Station",
        "latitude": "Latitude",
        "longitude": "Longitude",
        "sampleDepthM": "Depth",
        "measureValue": "MeasureValue",
        "parameter": "Parameter",
    })
    required_columns = ["Station", "Longitude", "Latitude", "SampleDate", "Depth", "MeasureValue", "Parameter"]
    missing_columns = [column for column in required_columns if column not in df_all.columns]
    if missing_columns:
        raise KeyError(
            "Required columns are missing from the CBP response: "
            + ", ".join(missing_columns)
            + f". Available columns: {df_all.columns.tolist()}"
        )

    df_all["SampleDate"] = pd.to_datetime(df_all["SampleDate"], errors="coerce")
    df_all["Depth"] = pd.to_numeric(df_all["Depth"], errors="coerce")
    df_all["MeasureValue"] = pd.to_numeric(df_all["MeasureValue"], errors="coerce")
    df_all = df_all.dropna(subset=["Station", "Longitude", "Latitude", "SampleDate", "Depth", "MeasureValue"]).copy()

    cbp_dict: dict[str, GenericObsData] = {}
    for station in stations:
        df = df_all[df_all["Station"].astype(str).str.strip() == station].copy()
        if df.empty:
            print(f"[WARNING] No CBP {varname} data found for {station}; skipping.", flush=True)
            continue

        lon = float(df["Longitude"].iloc[0])
        lat = float(df["Latitude"].iloc[0])
        time_difference = (df["SampleDate"] - sample_time).abs()
        df = df[time_difference == time_difference.min()].copy()
        selected_date = df["SampleDate"].iloc[0]
        df = df.sort_values("Depth")[["Station", "Depth", "MeasureValue", "SampleDate", "Parameter"]].copy()
        df = df.drop_duplicates(subset=["Depth"], keep="last")
        print(f"[INFO] {station}: selected {selected_date.strftime('%Y-%m-%d')}, {len(df)} vertical levels", flush=True)

        cbp_dict[station] = GenericObsData(
            station_info={"station_name": station, "lon": lon, "lat": lat},
            df=df,
        )

    if not cbp_dict:
        raise RuntimeError(f"No valid CBP {varname} station profiles were produced.")

    print(f"[INFO] Number of valid CBP {varname} stations: {len(cbp_dict)}", flush=True)
    return cbp_dict


def get_cbp_obs_for_stofs3d(
    outdir: str | Path | None = None,
    sample_time: str = "2015-09-18",
    varname: list[str] | None = None,
) -> None:
    if varname is None:
        varname = ["sal"]
    if outdir is None:
        outdir = "."
    outdir = Path(outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    for var in varname:
        cbp_var = CBP_VAR_NAMES[var]
        observations = get_cbp(sample_time=sample_time, varname=cbp_var)
        with (outdir / f"mean_{var}_xyz_{sample_time}").open("a", encoding="utf-8") as fout:
            for station, obs in observations.items():
                x = obs.station_info["lon"]
                y = obs.station_info["lat"]
                z = obs.df["MeasureValue"].values[-1]
                st = obs.station_info["station_name"]
                fout.write(f"{x} {y} {z} {st}\n")
