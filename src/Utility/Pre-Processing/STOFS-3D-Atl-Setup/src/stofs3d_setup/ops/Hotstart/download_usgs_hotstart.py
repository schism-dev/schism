"""
Small USGS Water Data OGC downloader for STOFS-3D-Atl hotstart files.

This module intentionally implements only the hotstart observation products
needed by tweak_stofs3d_hotstart.py:

    mean_temperature_xyz_YYYY-MM-DD
    mean_salinity_xyz_YYYY-MM-DD
    mean_conductance_xyz_YYYY-MM-DD
    mean_salinity_cond_xyz_YYYY-MM-DD
    mean_tem_xyz_YYYY-MM-DD -> mean_temperature_xyz_YYYY-MM-DD
    mean_sal_xyz_YYYY-MM-DD -> mean_salinity_cond_xyz_YYYY-MM-DD

The older implementation used climata/NWIS. This uses the modern USGS Water
Data OGC API and keeps the output format identical: lon lat mean_value site_id.
"""

from __future__ import annotations

import os
import time
import hashlib
from math import cos, radians
from datetime import datetime, timedelta
from pathlib import Path
from typing import Any
from urllib.parse import urljoin

import numpy as np
import pandas as pd
import requests


USGS_COLLECTION_URL = "https://api.waterdata.usgs.gov/ogcapi/v0/collections/{collection}/items"

USGS_PARAMETER_CODES = {
    "temperature": "00010",
    "salinity": "00480",
    "conductance": "00095",
}

USGS_ECGC_STATE_BBOXES = {
    "ME": (-71.1, 42.9, -66.8, 47.6),
    "NH": (-72.6, 42.6, -70.6, 45.4),
    "MA": (-73.6, 41.2, -69.8, 42.9),
    "RI": (-71.9, 41.1, -71.1, 42.1),
    "CT": (-73.8, 40.9, -71.7, 42.1),
    "NY": (-79.8, 40.4, -71.8, 45.1),
    "NJ": (-75.6, 38.8, -73.8, 41.4),
    "DE": (-75.9, 38.4, -75.0, 39.9),
    "PA": (-80.6, 39.6, -74.6, 42.3),
    "MD": (-79.6, 37.8, -75.0, 39.8),
    "DC": (-77.2, 38.8, -76.9, 39.0),
    "VA": (-83.7, 36.5, -75.2, 39.5),
    "NC": (-84.4, 33.8, -75.4, 36.7),
    "SC": (-83.4, 32.0, -78.5, 35.3),
    "GA": (-85.7, 30.3, -80.8, 35.1),
    "FL": (-87.8, 24.4, -80.0, 31.1),
    "AL": (-88.6, 30.1, -84.9, 35.1),
    "MI": (-90.5, 41.7, -82.1, 48.4),
    "LA": (-94.1, 28.8, -88.8, 33.1),
    "TX": (-106.7, 25.8, -93.5, 36.6),
}

DEFAULT_STATES = list(USGS_ECGC_STATE_BBOXES)

USGS_STATE_NAMES = {
    "ME": "Maine",
    "NH": "New Hampshire",
    "MA": "Massachusetts",
    "RI": "Rhode Island",
    "CT": "Connecticut",
    "NY": "New York",
    "NJ": "New Jersey",
    "DE": "Delaware",
    "PA": "Pennsylvania",
    "MD": "Maryland",
    "DC": "District of Columbia",
    "VA": "Virginia",
    "NC": "North Carolina",
    "SC": "South Carolina",
    "GA": "Georgia",
    "FL": "Florida",
    "AL": "Alabama",
    "MI": "Michigan",
    "LA": "Louisiana",
    "TX": "Texas",
}


def make_session() -> requests.Session:
    session = requests.Session()
    session.headers.update({
        "User-Agent": "stofs3d_setup hotstart USGS downloader",
        "Accept": "application/json, */*",
    })
    api_key = os.environ.get("USGS_API_KEY")
    if api_key:
        session.headers["X-Api-Key"] = api_key
    return session


def parse_time(value: Any) -> pd.Timestamp:
    return pd.to_datetime(value, utc=True)


def iso_interval(start: Any, end: Any) -> str:
    start_ts = parse_time(start).isoformat().replace("+00:00", "Z")
    end_ts = parse_time(end).isoformat().replace("+00:00", "Z")
    return f"{start_ts}/{end_ts}"


def request_json(
    session: requests.Session,
    url: str,
    params: dict[str, Any] | None = None,
    max_429_retries: int = 8,
) -> dict[str, Any]:
    for attempt in range(max_429_retries + 1):
        response = session.get(url, params=params, timeout=60)
        if response.status_code == 429 and "X-Api-Key" in session.headers:
            print("USGS API key is over rate limit; retrying anonymously", flush=True)
            session.headers.pop("X-Api-Key", None)
            continue
        if response.status_code != 429:
            response.raise_for_status()
            return response.json()
        if attempt == max_429_retries:
            response.raise_for_status()

        wait_seconds = min(120.0, 10.0 * 2**attempt)
        print(f"USGS rate limit hit; sleeping {wait_seconds:g} seconds", flush=True)
        time.sleep(wait_seconds)

    raise RuntimeError("unreachable")


def iter_ogc_features(
    session: requests.Session,
    collection: str,
    params: dict[str, Any],
    limit: int = 500,
    max_429_retries: int = 8,
) -> list[dict[str, Any]]:
    url = USGS_COLLECTION_URL.format(collection=collection)
    params = {"f": "json", "limit": limit, **params}
    features: list[dict[str, Any]] = []

    while url:
        payload = request_json(session, url, params=params, max_429_retries=max_429_retries)
        features.extend(payload.get("features", []))
        next_url = None
        for link in payload.get("links", []):
            if link.get("rel") == "next" and link.get("href"):
                next_url = urljoin(url, link["href"])
                break
        url = next_url
        params = None

    return features


def feature_lon_lat(feature: dict[str, Any]) -> tuple[float | None, float | None]:
    coords = feature.get("geometry", {}).get("coordinates")
    if isinstance(coords, list) and len(coords) >= 2:
        return float(coords[0]), float(coords[1])
    return None, None


def normalize_station_id(station_id: Any) -> str:
    station_id = str(station_id).strip()
    return station_id[5:] if station_id.upper().startswith("USGS-") else station_id


def series_overlaps(properties: dict[str, Any], start: Any, end: Any) -> bool:
    series_start = pd.to_datetime(
        properties.get("begin_utc") or properties.get("begin"), utc=True, errors="coerce"
    )
    series_end = pd.to_datetime(
        properties.get("end_utc") or properties.get("end"), utc=True, errors="coerce"
    )
    request_start = parse_time(start)
    request_end = parse_time(end)
    return not (
        pd.notna(series_start) and series_start > request_end
        or pd.notna(series_end) and series_end < request_start
    )


def _bbox_string(bbox: tuple[float, float, float, float]) -> str:
    return ",".join(f"{value:.12g}" for value in bbox)


def _read_hgrid(hgrid_path: str | os.PathLike[str] | None) -> Any | None:
    if hgrid_path is None:
        return None
    hgrid_path = Path(hgrid_path)
    if not hgrid_path.exists():
        print(f"USGS domain filter: {hgrid_path} not found; using unfiltered state candidates", flush=True)
        return None

    from pylib import schism_grid

    print(f"USGS domain filter: reading {hgrid_path}", flush=True)
    hgrid = schism_grid(str(hgrid_path))
    print(
        "USGS domain filter: "
        f"bbox=({hgrid.x.min():.4f}, {hgrid.y.min():.4f}, {hgrid.x.max():.4f}, {hgrid.y.max():.4f})",
        flush=True,
    )
    return hgrid


def _expanded_hgrid_bbox(hgrid: Any, buffer_km: float) -> tuple[float, float, float, float]:
    lat_mid = 0.5 * (float(hgrid.y.min()) + float(hgrid.y.max()))
    lat_pad = buffer_km / 111.0
    lon_pad = buffer_km / max(111.0 * cos(radians(lat_mid)), 1.0)
    return (
        float(hgrid.x.min()) - lon_pad,
        float(hgrid.x.max()) + lon_pad,
        float(hgrid.y.min()) - lat_pad,
        float(hgrid.y.max()) + lat_pad,
    )


def _boundary_points(hgrid: Any) -> np.ndarray:
    if not hasattr(hgrid, "bndinfo") or not hasattr(hgrid.bndinfo, "nb"):
        hgrid.compute_bnd()

    boundary_node_ids = np.unique(
        np.concatenate([
            np.asarray(hgrid.bndinfo.ibn[i], dtype=int)
            for i in range(int(hgrid.bndinfo.nb))
        ])
    )
    return np.c_[hgrid.x[boundary_node_ids], hgrid.y[boundary_node_ids]]


def _xy_km(points: np.ndarray, reference_lat: float) -> np.ndarray:
    km_per_lon_degree = max(111.0 * cos(radians(reference_lat)), 1.0)
    return np.c_[points[:, 0] * km_per_lon_degree, points[:, 1] * 111.0]


def _within_boundary_buffer(points: np.ndarray, hgrid: Any, buffer_km: float) -> np.ndarray:
    if buffer_km <= 0.0 or len(points) == 0:
        return np.zeros(len(points), dtype=bool)

    try:
        from scipy.spatial import cKDTree
    except ImportError:
        print(
            "USGS domain filter: scipy is unavailable; using inside-grid stations only",
            flush=True,
        )
        return np.zeros(len(points), dtype=bool)

    reference_lat = float(np.nanmean(points[:, 1]))
    boundary_tree = cKDTree(_xy_km(_boundary_points(hgrid), reference_lat))
    distance_km, _ = boundary_tree.query(_xy_km(points, reference_lat), k=1)
    return distance_km <= buffer_km


def filter_series_by_hgrid(
    series: pd.DataFrame,
    hgrid: Any | None,
    domain_buffer_km: float = 100.0,
) -> pd.DataFrame:
    if hgrid is None or series.empty:
        return series

    xmin, xmax, ymin, ymax = _expanded_hgrid_bbox(hgrid, domain_buffer_km)
    in_bbox = (
        series["lon"].between(xmin, xmax)
        & series["lat"].between(ymin, ymax)
    )
    bbox_filtered = series.loc[in_bbox].copy()
    if bbox_filtered.empty:
        return bbox_filtered

    points = bbox_filtered[["lon", "lat"]].to_numpy(dtype=float)
    inside = hgrid.inside_grid(points).astype(bool)
    buffered = _within_boundary_buffer(points, hgrid, domain_buffer_km)
    filtered = bbox_filtered.loc[inside | buffered].reset_index(drop=True)
    print(
        "USGS domain filter: "
        f"{len(filtered)} of {len(series)} candidate stations are inside hgrid "
        f"or within {domain_buffer_km:g} km of its boundary",
        flush=True,
    )
    return filtered


def discover_parameter_series(
    parameter: str,
    start: Any,
    end: Any,
    states: list[str] | None = None,
    session: requests.Session | None = None,
    max_429_retries: int = 8,
) -> pd.DataFrame:
    session = session or make_session()
    states = DEFAULT_STATES if states is None else states
    rows = []

    for state in states:
        state_code = state.upper()
        state_name = USGS_STATE_NAMES.get(state_code, state)
        print(f"USGS parameter {parameter}: discovering {state_code} stations", flush=True)
        features = iter_ogc_features(
            session,
            "time-series-metadata",
            {
                "state_name": state_name,
                "parameter_code": parameter,
            },
            limit=10000,
            max_429_retries=max_429_retries,
        )
        n_before = len(rows)
        for feature in features:
            properties = feature.get("properties", {})
            if str(properties.get("computation_identifier", "")).lower() != "instantaneous":
                continue
            if not series_overlaps(properties, start, end):
                continue
            lon, lat = feature_lon_lat(feature)
            station_id = properties.get("monitoring_location_id")
            if not station_id or lon is None or lat is None:
                continue
            rows.append({
                "station_id": station_id,
                "time_series_id": feature.get("id"),
                "state": state_code,
                "state_name": state_name,
                "lon": lon,
                "lat": lat,
                "parameter": properties.get("parameter_code") or parameter,
                "unit": properties.get("unit_of_measure"),
                "primary": properties.get("primary"),
            })
        print(
            f"USGS parameter {parameter}: {state.upper()} added {len(rows) - n_before} candidate series",
            flush=True,
        )

    series = pd.DataFrame(rows)
    if series.empty:
        return pd.DataFrame(
            columns=[
                "station_id", "time_series_id", "state", "state_name", "lon",
                "lat", "parameter", "unit", "primary",
            ]
        )

    primary = series["primary"].fillna("").astype(str).str.lower() == "primary"
    stations_with_primary = set(series.loc[primary, "station_id"])
    series = series[primary | ~series["station_id"].isin(stations_with_primary)]
    return series.drop_duplicates(["station_id", "parameter"]).sort_values("station_id").reset_index(drop=True)


def download_station_parameter(
    station_id: str,
    parameter: str,
    start: Any,
    end: Any,
    session: requests.Session,
) -> pd.DataFrame:
    features = iter_ogc_features(
        session,
        "continuous",
        {
            "monitoring_location_id": station_id,
            "parameter_code": parameter,
            "datetime": iso_interval(start, end),
        },
    )
    rows = []
    for feature in features:
        properties = feature.get("properties", {})
        lon, lat = feature_lon_lat(feature)
        rows.append({
            "station_id": properties.get("monitoring_location_id", station_id),
            "station_name": properties.get("monitoring_location_name"),
            "lon": lon,
            "lat": lat,
            "time": properties.get("time"),
            "value": properties.get("value"),
            "unit": properties.get("unit_of_measure"),
        })
    observations = pd.DataFrame(rows)
    if observations.empty:
        return observations
    observations["time"] = pd.to_datetime(observations["time"], utc=True, errors="coerce")
    observations["value"] = pd.to_numeric(observations["value"], errors="coerce")
    observations["lon"] = pd.to_numeric(observations["lon"], errors="coerce")
    observations["lat"] = pd.to_numeric(observations["lat"], errors="coerce")
    return observations.dropna(subset=["time", "value"])


def download_state_parameter(
    state_name: str,
    parameter: str,
    start: Any,
    end: Any,
    session: requests.Session,
) -> pd.DataFrame:
    features = iter_ogc_features(
        session,
        "continuous",
        {
            "state_name": state_name,
            "parameter_code": parameter,
            "datetime": iso_interval(start, end),
        },
    )
    rows = []
    for feature in features:
        properties = feature.get("properties", {})
        lon, lat = feature_lon_lat(feature)
        rows.append({
            "station_id": properties.get("monitoring_location_id"),
            "station_name": properties.get("monitoring_location_name"),
            "lon": lon,
            "lat": lat,
            "time": properties.get("time"),
            "value": properties.get("value"),
            "unit": properties.get("unit_of_measure"),
        })
    observations = pd.DataFrame(rows)
    if observations.empty:
        return observations
    observations["station_id"] = observations["station_id"].astype(str)
    observations["time"] = pd.to_datetime(observations["time"], utc=True, errors="coerce")
    observations["value"] = pd.to_numeric(observations["value"], errors="coerce")
    observations["lon"] = pd.to_numeric(observations["lon"], errors="coerce")
    observations["lat"] = pd.to_numeric(observations["lat"], errors="coerce")
    return observations.dropna(subset=["station_id", "time", "value"])


def download_station_group_parameter(
    station_ids: list[str],
    parameter: str,
    start: Any,
    end: Any,
    session: requests.Session,
    max_429_retries: int = 2,
) -> pd.DataFrame:
    features = iter_ogc_features(
        session,
        "continuous",
        {
            "monitoring_location_id": ",".join(station_ids),
            "parameter_code": parameter,
            "datetime": iso_interval(start, end),
        },
        limit=10000,
        max_429_retries=max_429_retries,
    )
    rows = []
    for feature in features:
        properties = feature.get("properties", {})
        lon, lat = feature_lon_lat(feature)
        rows.append({
            "station_id": properties.get("monitoring_location_id"),
            "station_name": properties.get("monitoring_location_name"),
            "lon": lon,
            "lat": lat,
            "time": properties.get("time"),
            "value": properties.get("value"),
            "unit": properties.get("unit_of_measure"),
        })
    observations = pd.DataFrame(rows)
    if observations.empty:
        return observations
    observations["station_id"] = observations["station_id"].astype(str)
    observations["time"] = pd.to_datetime(observations["time"], utc=True, errors="coerce")
    observations["value"] = pd.to_numeric(observations["value"], errors="coerce")
    observations["lon"] = pd.to_numeric(observations["lon"], errors="coerce")
    observations["lat"] = pd.to_numeric(observations["lat"], errors="coerce")
    return observations.dropna(subset=["station_id", "time", "value"])


def practical_salinity_from_conductance_25c(conductance_ms_cm: float) -> float:
    """Convert conductivity in mS/cm to practical salinity at 25 C and 0 dbar."""
    if np.isnan(conductance_ms_cm):
        return np.nan

    t = 25.0
    c3515 = 42.9140
    rt35 = (
        0.6766097
        + 0.0200564 * t
        + 0.0001104259 * t**2
        - 0.00000069698 * t**3
        + 0.0000000010031 * t**4
    )
    rt = max((conductance_ms_cm / c3515) / rt35, 0.0)
    rtx = rt**0.5

    a0, a1, a2, a3, a4, a5 = 0.0080, -0.1692, 25.3851, 14.0941, -7.0261, 2.7081
    b0, b1, b2, b3, b4, b5 = 0.0005, -0.0056, -0.0066, -0.0375, 0.0636, -0.0144
    k = 0.0162
    return (
        a0
        + (a1 + (a2 + (a3 + (a4 + a5 * rtx) * rtx) * rtx) * rtx) * rtx
        + ((t - 15.0) / (1.0 + k * (t - 15.0)))
        * (b0 + (b1 + (b2 + (b3 + (b4 + b5 * rtx) * rtx) * rtx) * rtx) * rtx)
    )


def _mean_value(observations: pd.DataFrame, parameter: str) -> float:
    value = pd.to_numeric(observations["value"], errors="coerce").mean()
    if parameter == USGS_PARAMETER_CODES["conductance"]:
        value = practical_salinity_from_conductance_25c(value * 0.001)
    return value


def _write_station_mean(
    fout: Any,
    station_id: str,
    observations: pd.DataFrame,
    station_info: pd.Series,
    parameter: str,
) -> bool:
    value = _mean_value(observations, parameter)
    lon = observations["lon"].dropna().iloc[0] if observations["lon"].notna().any() else station_info["lon"]
    lat = observations["lat"].dropna().iloc[0] if observations["lat"].notna().any() else station_info["lat"]
    if np.isnan(value) or pd.isna(lon) or pd.isna(lat):
        return False
    fout.write(f"{lon} {lat} {value} {normalize_station_id(station_id)}\n")
    return True


def _csv_has_columns(path: Path, columns: set[str]) -> bool:
    if not path.exists() or path.stat().st_size == 0:
        return False
    try:
        return columns.issubset(set(pd.read_csv(path, nrows=0).columns))
    except pd.errors.EmptyDataError:
        return False


def _chunks(values: list[str], size: int) -> list[list[str]]:
    return [values[i:i + size] for i in range(0, len(values), size)]


def _load_or_discover_state_series(
    outdir: Path,
    parameter: str,
    state: str,
    start: Any,
    end: Any,
    session: requests.Session,
    expected_columns: set[str],
    max_429_retries: int = 8,
) -> pd.DataFrame:
    state_code = state.upper()
    state_cache = outdir / f"usgs_series_{parameter}_{state_code}_{start}_{end}.csv"
    if _csv_has_columns(state_cache, expected_columns):
        return pd.read_csv(
            state_cache,
            dtype={"station_id": str, "time_series_id": str, "parameter": str},
        )

    state_series = discover_parameter_series(
        parameter,
        start,
        end,
        states=[state_code],
        session=session,
        max_429_retries=max_429_retries,
    )
    state_series.to_csv(state_cache, index=False)
    return state_series


def write_parameter_mean_file(
    parameter: str,
    outfilename: Path,
    start: Any,
    end: Any,
    states: list[str] | None = None,
    hgrid: Any | None = None,
    domain_buffer_km: float = 100.0,
    request_interval_s: float = 1.0,
    stations_per_request: int = 25,
    obs_max_429_retries: int = 2,
    metadata_max_429_retries: int = 8,
    session: requests.Session | None = None,
) -> None:
    if outfilename.exists():
        print(f"{outfilename} exists, skipping ...")
        return

    session = session or make_session()
    outfilename.parent.mkdir(parents=True, exist_ok=True)
    states_label = "_".join(states or DEFAULT_STATES)
    series_cache = outfilename.parent / f"usgs_series_{parameter}_{states_label}_{start}_{end}.csv"
    filtered_cache = outfilename.parent / f"usgs_series_{parameter}_{states_label}_{start}_{end}_filtered.csv"
    expected_series_columns = {"station_id", "time_series_id", "state", "state_name", "lon", "lat", "parameter"}
    if _csv_has_columns(filtered_cache, expected_series_columns):
        series = pd.read_csv(filtered_cache, dtype={"station_id": str, "time_series_id": str, "parameter": str})
        print(f"USGS parameter {parameter}: using cached filtered station list {filtered_cache}", flush=True)
    else:
        if _csv_has_columns(series_cache, expected_series_columns):
            series = pd.read_csv(series_cache, dtype={"station_id": str, "time_series_id": str, "parameter": str})
            print(f"USGS parameter {parameter}: using cached station list {series_cache}", flush=True)
        else:
            discovered = []
            for state in (states or DEFAULT_STATES):
                try:
                    discovered.append(
                        _load_or_discover_state_series(
                            outfilename.parent,
                            parameter,
                            state,
                            start,
                            end,
                            session,
                            expected_series_columns,
                            max_429_retries=metadata_max_429_retries,
                        )
                    )
                except Exception as exc:
                    print(
                        f"Warning: failed discovering {state} {parameter} stations: {exc}",
                        flush=True,
                    )
            discovered = [df for df in discovered if not df.empty]
            series = pd.concat(discovered, ignore_index=True) if discovered else pd.DataFrame()
            series.to_csv(series_cache, index=False)
        series = filter_series_by_hgrid(series, hgrid, domain_buffer_km=domain_buffer_km)
        series.to_csv(filtered_cache, index=False)
    print(f"USGS parameter {parameter}: {len(series)} candidate instantaneous stations", flush=True)

    tmp_path = outfilename.with_suffix(outfilename.suffix + ".tmp")
    if series.empty:
        tmp_path.write_text("", encoding="utf-8")
        os.replace(tmp_path, outfilename)
        print(f"USGS parameter {parameter}: wrote 0 station means to {outfilename}", flush=True)
        return

    completed_station_ids: set[str] = set()
    if tmp_path.exists():
        previous = pd.read_csv(
            tmp_path,
            sep=r"\s+",
            names=["lon", "lat", "value", "station_id"],
            dtype={"station_id": str},
        )
        completed_station_ids = set(previous["station_id"].astype(str))
        print(
            f"USGS parameter {parameter}: resuming from {tmp_path} "
            f"with {len(completed_station_ids)} completed station means",
            flush=True,
        )

    n_written = 0
    with tmp_path.open("a", encoding="utf-8") as fout:
        for state, state_series in series.groupby("state", sort=True):
            station_map = state_series.set_index("station_id")
            state_station_ids = state_series["station_id"].astype(str).tolist()
            for chunk_index, station_chunk in enumerate(
                _chunks(state_station_ids, stations_per_request),
                start=1,
            ):
                pending_station_ids = [
                    station_id for station_id in station_chunk
                    if normalize_station_id(station_id) not in completed_station_ids
                ]
                if not pending_station_ids:
                    continue
                chunk_key = hashlib.sha1(",".join(station_chunk).encode("utf-8")).hexdigest()[:12]
                obs_cache = outfilename.parent / (
                    f"usgs_obs_{parameter}_{state}_{chunk_index:03d}_{chunk_key}_{start}_{end}.csv"
                )
                try:
                    if obs_cache.exists():
                        observations = pd.read_csv(obs_cache, dtype={"station_id": str})
                        print(
                            f"USGS parameter {parameter}: using cached {state} "
                            f"chunk {chunk_index} observations",
                            flush=True,
                        )
                    else:
                        print(
                            f"USGS parameter {parameter}: bulk downloading {state} "
                            f"chunk {chunk_index} ({len(station_chunk)} stations)",
                            flush=True,
                        )
                        observations = download_station_group_parameter(
                            station_chunk,
                            parameter,
                            start,
                            end,
                            session=session,
                            max_429_retries=obs_max_429_retries,
                        )
                        observations.to_csv(obs_cache, index=False)
                except Exception as exc:
                    print(
                        f"Warning: failed bulk downloading {state} {parameter} "
                        f"chunk {chunk_index}: {exc}",
                        flush=True,
                    )
                    continue
                if observations.empty:
                    continue
                selected = set(pending_station_ids)
                observations = observations[observations["station_id"].astype(str).isin(selected)]
                for station_id, station_observations in observations.groupby("station_id", sort=True):
                    if _write_station_mean(
                        fout,
                        station_id,
                        station_observations,
                        station_map.loc[station_id],
                        parameter,
                    ):
                        n_written += 1
                fout.flush()
                print(
                    f"USGS parameter {parameter}: wrote {n_written} new station means after "
                    f"{state} chunk {chunk_index}",
                    flush=True,
                )
                if request_interval_s > 0.0:
                    time.sleep(request_interval_s)

    os.replace(tmp_path, outfilename)
    print(f"USGS parameter {parameter}: wrote {n_written} station means to {outfilename}", flush=True)


def _replace_symlink(link_path: Path, target_name: str) -> None:
    if link_path.exists() or link_path.is_symlink():
        link_path.unlink()
    link_path.symlink_to(target_name)


def get_usgs_obs_for_stofs3d(
    outdir: str | os.PathLike[str],
    start_date_str: str = "2015-09-18",
    end_date_str: str | None = None,
    vars: list[str] | None = None,
    states: list[str] | None = None,
    hgrid_path: str | os.PathLike[str] | None = None,
    domain_buffer_km: float = 100.0,
    request_interval_s: float = 1.0,
    stations_per_request: int = 25,
    obs_max_429_retries: int = 2,
    metadata_max_429_retries: int = 8,
) -> None:
    if vars is None:
        vars = ["temperature", "salinity", "conductance"]
    if end_date_str is None:
        end_date_str = (
            datetime.strptime(start_date_str, "%Y-%m-%d") + timedelta(days=1)
        ).strftime("%Y-%m-%d")

    outdir = Path(outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    session = make_session()
    hgrid = _read_hgrid(hgrid_path)
    outfilenames = {}

    for var in vars:
        parameter = USGS_PARAMETER_CODES[var]
        outfilenames[var] = outdir / f"mean_{var}_xyz_{start_date_str}"
        write_parameter_mean_file(
            parameter=parameter,
            outfilename=outfilenames[var],
            start=start_date_str,
            end=end_date_str,
            states=states,
            hgrid=hgrid,
            domain_buffer_km=domain_buffer_km,
            request_interval_s=request_interval_s,
            stations_per_request=stations_per_request,
            obs_max_429_retries=obs_max_429_retries,
            metadata_max_429_retries=metadata_max_429_retries,
            session=session,
        )

    salinity_cond = outdir / f"mean_salinity_cond_xyz_{start_date_str}"
    with salinity_cond.open("w", encoding="utf-8") as fout:
        for var in ["salinity", "conductance"]:
            path = outfilenames.get(var)
            if path is not None and path.exists():
                fout.write(path.read_text(encoding="utf-8"))

    if "temperature" in outfilenames:
        _replace_symlink(outdir / f"mean_tem_xyz_{start_date_str}", outfilenames["temperature"].name)
    _replace_symlink(outdir / f"mean_sal_xyz_{start_date_str}", salinity_cond.name)
