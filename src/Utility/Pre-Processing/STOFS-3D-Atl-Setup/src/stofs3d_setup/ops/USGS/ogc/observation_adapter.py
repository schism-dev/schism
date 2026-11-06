"""Adapt USGS OGC station observations for source/sink workflows."""

from __future__ import annotations

import json
import os
import time
from pathlib import Path
from typing import Any, Callable, Iterable, List, Mapping, Optional

import numpy as np
import pandas as pd

from stofs3d_setup.ops.USGS.iv.download_usgs import (
    GenericObsData,
    RequiredStationDownloadError,
    convert_to_ObsData,
    detect_data_gap,
    df_to_generic_obs,
    find_data_gaps,
    load_cache_parquet,
    save_cache_csv,
    save_cache_parquet,
)
from stofs3d_setup.ops.USGS.ogc.download_policy import (
    USGS_DOWNLOAD_POLICY,
    USGS_MAX_MISSING_GAP_DAYS,
)
from stofs3d_setup.ops.USGS.ogc.download_usgs import (
    MIN_WATER_TEMPERATURE_C,
    USGS_STATE_NAMES,
    discover_monitoring_locations,
    discover_station_ids_with_parameter_in_bbox,
    download_station_group_parameter,
    log_bulk_diagnostic,
    make_session,
    normalize_station_id,
)


usgs_var_dict = {
    "streamflow": {"id": "00060", "unit_conv": 0.028316846592},
    "salinity": {"id": "00480", "unit_conv": 1.0},
    "gauge height": {"id": "00065", "unit_conv": 1.0},
    "temperature": {"id": "00010", "unit_conv": 1.0},
    "conductance": {"id": "00095", "unit_conv": 0.001},
}

ecgc_states = list(USGS_STATE_NAMES)
CACHE_QUALITY_VERSION = 1
WindowResultCallback = Callable[
    [list[str], str, Any, Any, pd.DataFrame | None, Exception | None], None
]
WindowStartCallback = Callable[[list[str], str, Any, Any], None]


class ObservationPersistenceError(RuntimeError):
    """Raised when a successful OGC response cannot be safely persisted."""


def _normalize_station_ids(station_ids: Iterable[str] | None) -> list[str]:
    normalized = {
        normalize_station_id(station_id)
        for station_id in (station_ids or ())
        if station_id is not None and str(station_id).strip()
    }
    return sorted(normalized)


def get_usgs_stations_from_state(
    states: List[str] | None = None,
    data_type: str = "iv",
    parameter: Optional[str] = None,
    start: Any = "1800-01-01T00:00:00Z",
    end: Any = "2200-01-01T00:00:00Z",
    locations: pd.DataFrame | None = None,
) -> pd.DataFrame:
    """List locations with the requested instantaneous series in the period."""
    del data_type
    states = ecgc_states if states is None else states
    parameter = parameter or usgs_var_dict["streamflow"]["id"]
    locations = (
        discover_monitoring_locations(states=states)
        if locations is None
        else locations
    )
    locations = filter_station_inventory_by_parameter(
        locations,
        parameter=parameter,
        start=start,
        end=end,
    )
    return pd.DataFrame(
        {
            "site_no": locations["station_id"].map(normalize_station_id),
            "station_nm": locations["station_name"],
            "state": locations["state"],
            "dec_long_va": pd.to_numeric(locations["lon"], errors="coerce"),
            "dec_lat_va": pd.to_numeric(locations["lat"], errors="coerce"),
            "begin_date": pd.NaT,
            "end_date": pd.NaT,
        }
    ).reset_index(drop=True)


def filter_station_inventory_by_parameter(
    locations: pd.DataFrame,
    parameter: str,
    start: Any,
    end: Any,
) -> pd.DataFrame:
    """Keep inventory locations with a matching series for the run period."""
    if locations.empty:
        return locations.copy()

    matching_ids = set()
    for state, state_locations in locations.groupby("state"):
        bbox = (
            float(state_locations["lon"].min()),
            float(state_locations["lat"].min()),
            float(state_locations["lon"].max()),
            float(state_locations["lat"].max()),
        )
        print(f"USGS parameter {parameter}: screening {state} by bounding box", flush=True)
        matching_ids.update(discover_station_ids_with_parameter_in_bbox(
            bbox=bbox,
            parameter=parameter,
            start=start,
            end=end,
        ))
    normalized_ids = locations["station_id"].map(normalize_station_id)
    return locations.loc[normalized_ids.isin(matching_ids)].reset_index(drop=True)


def _as_utc(value) -> pd.Timestamp:
    timestamp = pd.Timestamp(value)
    if timestamp.tzinfo is None:
        return timestamp.tz_localize("UTC")
    return timestamp.tz_convert("UTC")


def _iter_download_windows(start, end, chunk_days: int):
    """Yield adjacent OGC windows, with inclusive endpoints and no gaps."""
    window_start = _as_utc(start)
    end = _as_utc(end)
    if end < window_start:
        raise ValueError("USGS download end time precedes start time")
    if chunk_days <= 0:
        raise ValueError("USGS chunk size must be positive")

    while window_start <= end:
        window_end = min(window_start + pd.Timedelta(days=chunk_days), end)
        yield window_start, window_end
        if window_end == end:
            break
        window_start = window_end + pd.Timedelta(nanoseconds=1)


def _manifest_path(cache_path: str | Path) -> Path:
    return Path(f"{cache_path}.manifest.json")


def _cached_station_ids_for_period(
    cache_path: str | Path,
    parameter: str,
    start,
    end,
) -> set[str] | None:
    manifest_path = _manifest_path(cache_path)
    if not manifest_path.is_file():
        return None
    try:
        manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
    except (OSError, json.JSONDecodeError):
        return None
    compatible_format = (
        manifest.get("api") == "usgs-waterdata-ogc-v1"
        or manifest.get("version") == 1
    )
    stations = manifest.get("stations")
    try:
        same_period = (
            _as_utc(manifest["start"]) == _as_utc(start)
            and _as_utc(manifest["end"]) == _as_utc(end)
        )
    except (KeyError, TypeError, ValueError):
        return None
    if not (
        compatible_format
        and manifest.get("quality_version") == CACHE_QUALITY_VERSION
        and manifest.get("parameter") == str(parameter)
        and same_period
        and isinstance(stations, list)
    ):
        return None
    return set(_normalize_station_ids(stations))


def _write_cache_manifest(
    cache_path: str | Path,
    station_ids: list[str],
    parameter: str,
    start,
    end,
) -> None:
    manifest = {
        "api": "usgs-waterdata-ogc-v1",
        "quality_version": CACHE_QUALITY_VERSION,
        "parameter": str(parameter),
        "stations": station_ids,
        "start": _as_utc(start).isoformat(),
        "end": _as_utc(end).isoformat(),
    }
    _manifest_path(cache_path).write_text(
        json.dumps(manifest, indent=2, sort_keys=True),
        encoding="utf-8",
    )


def _read_cache(cache_path: str | Path) -> list[GenericObsData]:
    cache_path = str(cache_path)
    if cache_path.endswith((".parquet", ".pq")):
        return df_to_generic_obs(load_cache_parquet(cache_path))
    if cache_path.endswith(".csv"):
        frame = pd.read_csv(cache_path, parse_dates=["date"])
        frame["date"] = pd.to_datetime(frame["date"], utc=True)
        return df_to_generic_obs(frame)
    raise ValueError(f"Unsupported USGS cache format: {cache_path}")


def _filter_temperature_records(
    records: list[GenericObsData], parameter: str
) -> list[GenericObsData]:
    """Keep the physical temperature floor at the cache boundary too."""
    if str(parameter) != "00010":
        return records
    filtered = []
    for record in records:
        values = pd.to_numeric(record.df["value"], errors="coerce")
        filtered.append(GenericObsData(
            station_info=record.station_info,
            df=record.df.loc[values >= MIN_WATER_TEMPERATURE_C].copy(),
        ))
    return filtered


def _write_cache(records: list[GenericObsData], cache_path: str | Path) -> None:
    cache_path = str(cache_path)
    if cache_path.endswith((".parquet", ".pq")):
        save_cache_parquet(records, cache_path)
    elif cache_path.endswith(".csv"):
        save_cache_csv(records, cache_path)
    else:
        raise ValueError(f"Unsupported USGS cache format: {cache_path}")


def _record_failures(
    records: list[GenericObsData],
    station_ids: Iterable[str],
    start: pd.Timestamp,
    end: pd.Timestamp,
    max_gap: pd.Timedelta | None,
) -> dict[str, str]:
    records_by_id = {
        normalize_station_id(record.station_info.get("id")): record
        for record in records
    }
    failures = {}
    for station_id in sorted(set(station_ids)):
        record = records_by_id.get(station_id)
        if record is None or record.df.empty:
            failures[station_id] = "no observations returned"
            continue
        if max_gap is None:
            continue
        data = record.df.copy()
        data.index = pd.DatetimeIndex(pd.to_datetime(data["date"], utc=True))
        gap = detect_data_gap(
            data,
            start=start,
            end=end,
            freq=pd.Timedelta(hours=1),
        )
        if gap > max_gap:
            failures[station_id] = (
                f"largest missing duration {gap} exceeds {max_gap}"
            )
    return failures


def _records_from_observations(
    observations: pd.DataFrame,
    parameter: str,
) -> list[GenericObsData]:
    records = []
    if observations.empty:
        return records

    for station_id, station_data in observations.groupby(
        "station_id", sort=True
    ):
        station_data = station_data.sort_values("time").drop_duplicates(
            "time", keep="last"
        )
        info = {
            "id": normalize_station_id(station_id),
            "name": station_data["station_name"].dropna().iloc[0]
            if station_data["station_name"].notna().any()
            else None,
            "lon": station_data["lon"].dropna().iloc[0]
            if station_data["lon"].notna().any()
            else np.nan,
            "lat": station_data["lat"].dropna().iloc[0]
            if station_data["lat"].notna().any()
            else np.nan,
            "var_name": parameter,
            "var_code": parameter,
            "unit": station_data["unit"].dropna().iloc[0]
            if station_data["unit"].notna().any()
            else None,
        }
        frame = pd.DataFrame(
            {
                "date": pd.to_datetime(station_data["time"], utc=True),
                "value": pd.to_numeric(station_data["value"], errors="coerce"),
            }
        ).dropna(subset=["date", "value"])
        records.append(GenericObsData(station_info=info, df=frame))
    return records


def download_stations(
    param_id: str,
    station_ids: List[str],
    cache_fname: Optional[str] = None,
    datelist: pd.DatetimeIndex = pd.date_range(
        start="2000-01-01", end="2000-01-02"
    ),
    stations_chunk_size: int = 10,
    days_per_chunk: Optional[int] = None,
    required_station_ids: Optional[Iterable[str]] = None,
    optional_station_ids: Optional[Iterable[str]] = None,
    download_policy: Optional[Mapping[str, Any]] = None,
    max_missing_gap_days: Optional[float] = None,
    validation_start=None,
    validation_end=None,
    warning_log_path: str | Path | None = None,
    warning_input_files: Iterable[str] = (),
    warning_forcing_name: str | None = None,
    allow_adaptive_retry: bool = True,
    raise_on_fetch_error: bool = False,
    validate_coverage: bool = True,
    observation_window_callback: WindowResultCallback | None = None,
    observation_window_start_callback: WindowStartCallback | None = None,
) -> List[GenericObsData]:
    """Download station observations from the OGC continuous collection.

    The start callback marks a station/window pending immediately before its
    request; the result callback receives each completed or failed window,
    including successful empty responses.
    """
    station_ids = _normalize_station_ids(station_ids)
    required = set(_normalize_station_ids(required_station_ids))
    optional = set(_normalize_station_ids(optional_station_ids))
    forced = required | optional
    if forced - set(station_ids):
        raise ValueError(
            "Forced USGS stations must be included in station_ids: "
            f"{sorted(forced - set(station_ids))}"
        )
    if required & optional:
        raise ValueError("USGS stations cannot be both required and optional")

    policy = dict(download_policy or USGS_DOWNLOAD_POLICY)
    chunk_days = days_per_chunk or policy["primary_chunk_days"]
    max_missing_gap_days = (
        USGS_MAX_MISSING_GAP_DAYS
        if max_missing_gap_days is None and validate_coverage
        else max_missing_gap_days
    )
    max_gap = (
        None
        if max_missing_gap_days is None
        else pd.Timedelta(days=max_missing_gap_days)
    )
    request_start = pd.to_datetime(
        datelist[0] if validation_start is None else validation_start,
        utc=True,
    )
    request_end = pd.to_datetime(
        datelist[-1] if validation_end is None else validation_end,
        utc=True,
    )
    fetch_start = _as_utc(datelist[0])
    fetch_end = _as_utc(datelist[-1])
    if request_end < request_start:
        raise ValueError("USGS validation end precedes start")

    reusable_records: list[GenericObsData] = []
    stations_to_fetch = station_ids
    if cache_fname is not None and os.path.isfile(cache_fname):
        cached_station_ids = _cached_station_ids_for_period(
            cache_fname, param_id, fetch_start, fetch_end
        )
        if cached_station_ids is not None:
            requested = set(station_ids)
            cached = [
                record for record in _filter_temperature_records(
                    _read_cache(cache_fname), param_id
                )
                if normalize_station_id(record.station_info["id"]) in requested
            ]
            cached_record_ids = {
                normalize_station_id(record.station_info["id"])
                for record in cached
            }
            failures = _record_failures(
                cached,
                forced | cached_record_ids,
                request_start,
                request_end,
                max_gap,
            )
            if cached_station_ids == requested and not failures:
                return cached
            reusable_records = [
                record for record in cached
                if normalize_station_id(record.station_info["id"])
                not in failures
            ]
            reusable_ids = {
                normalize_station_id(record.station_info["id"])
                for record in reusable_records
            }
            stations_to_fetch = [
                station_id for station_id in station_ids
                if station_id not in reusable_ids
            ]
            print(
                f"OGC cache: reusing {len(reusable_records)} station(s); "
                f"fetching {len(stations_to_fetch)} station(s)",
                flush=True,
            )

    session = make_session()
    observation_frames: dict[str, list[pd.DataFrame]] = {}
    failed_request_station_ids: set[str] = set()
    for station_offset in range(0, len(stations_to_fetch), stations_chunk_size):
        station_group = stations_to_fetch[
            station_offset:station_offset + stations_chunk_size
        ]
        for window_start, window_end in _iter_download_windows(
            fetch_start,
            fetch_end,
            chunk_days,
        ):
            started = time.monotonic()
            log_bulk_diagnostic(
                "window_start",
                group_start=station_offset + 1,
                group_size=len(station_group),
                stations=",".join(station_group),
                parameter=param_id,
                start=window_start.isoformat(),
                end=window_end.isoformat(),
            )
            if observation_window_start_callback is not None:
                observation_window_start_callback(
                    station_group,
                    param_id,
                    window_start,
                    window_end,
                )
            try:
                observations = download_station_group_parameter(
                    station_ids=station_group,
                    parameter=param_id,
                    start=window_start,
                    end=window_end,
                    session=session,
                    primary_only=True,
                    page_size=policy.get("page_size", 1000),
                )
            except Exception as exc:
                if observation_window_callback is not None:
                    observation_window_callback(
                        station_group,
                        param_id,
                        window_start,
                        window_end,
                        None,
                        exc,
                    )
                log_bulk_diagnostic(
                    "window_error",
                    group_start=station_offset + 1,
                    error=type(exc).__name__,
                    elapsed_s=f"{time.monotonic() - started:.1f}",
                )
                if raise_on_fetch_error:
                    raise
                print(
                    f"OGC download failed for {station_group} / {param_id} "
                    f"from {window_start} to {window_end}: {exc}"
                )
                failed_request_station_ids.update(station_group)
                continue
            if observation_window_callback is not None:
                observation_window_callback(
                    station_group,
                    param_id,
                    window_start,
                    window_end,
                    observations,
                    None,
                )
            log_bulk_diagnostic(
                "window_complete",
                group_start=station_offset + 1,
                records=len(observations),
                elapsed_s=f"{time.monotonic() - started:.1f}",
            )
            if observations.empty:
                continue
            for station_id, station_data in observations.groupby(
                "station_id", sort=False
            ):
                normalized_id = normalize_station_id(station_id)
                observation_frames.setdefault(normalized_id, []).append(
                    station_data
                )

    records = reusable_records + _records_from_observations(
        pd.concat(
            [frame for pieces in observation_frames.values() for frame in pieces],
            ignore_index=True,
        )
        if observation_frames
        else pd.DataFrame(),
        param_id,
    )
    records = _filter_temperature_records(records, param_id)

    if allow_adaptive_retry and station_ids:
        from stofs3d_setup.ops.USGS.ogc.series import download_usgs_series

        records_by_id = {
            normalize_station_id(record.station_info["id"]): record
            for record in records
        }
        retry_windows = {
            station_id: [(request_start, request_end)]
            for station_id in failed_request_station_ids
        }
        records_by_id = {
            normalize_station_id(record.station_info["id"]): record
            for record in records
        }
        for station_id in forced - records_by_id.keys():
            retry_windows[station_id] = [(request_start, request_end)]
        for station_id, record in records_by_id.items():
            record_times = pd.DatetimeIndex(
                pd.to_datetime(record.df["date"], utc=True)
            )
            gaps = (
                []
                if max_gap is None
                else find_data_gaps(
                    record_times,
                    request_start,
                    request_end,
                    max_gap,
                )
            )
            if gaps:
                retry_windows[station_id] = gaps

        for station_id, windows in sorted(retry_windows.items()):
            if cache_fname is None:
                raise ValueError(
                    "cache_fname is required for adaptive OGC retries"
                )
            retry_pieces = []
            for window_start, window_end in windows:
                try:
                    series = download_usgs_series(
                        station_id=station_id,
                        parameter_id=str(param_id),
                        start_time=window_start,
                        end_time=window_end,
                        usgs_cache_folder=Path(cache_fname).parent,
                        policy=policy,
                        observation_window_callback=(
                            observation_window_callback
                        ),
                        observation_window_start_callback=(
                            observation_window_start_callback
                        ),
                    )
                except Exception as exc:
                    if isinstance(exc, ObservationPersistenceError):
                        raise
                    print(f"Adaptive OGC retry failed for {station_id}: {exc}")
                    continue
                if series is not None and not series.empty:
                    retry_pieces.append(series)
            if not retry_pieces:
                continue

            old_record = records_by_id.get(station_id)
            if old_record is not None:
                old_series = pd.Series(
                    old_record.df["value"].to_numpy(dtype=float),
                    index=pd.DatetimeIndex(
                        pd.to_datetime(old_record.df["date"], utc=True)
                    ),
                )
                retry_pieces.insert(0, old_series)
            merged = pd.concat(retry_pieces).sort_index()
            merged = merged[~merged.index.duplicated(keep="last")]
            station_info = (
                old_record.station_info.copy()
                if old_record is not None
                else {
                    "id": station_id,
                    "name": None,
                    "lon": np.nan,
                    "lat": np.nan,
                    "var_name": param_id,
                    "var_code": param_id,
                    "unit": None,
                }
            )
            records = [
                record
                for record in records
                if normalize_station_id(record.station_info["id"])
                != station_id
            ]
            records.append(
                GenericObsData(
                    station_info=station_info,
                    df=pd.DataFrame(
                        {"date": merged.index, "value": merged.to_numpy()}
                    ),
                )
            )
            records_by_id[station_id] = records[-1]

    records = _filter_temperature_records(records, param_id)
    checked_station_ids = forced | {
        normalize_station_id(record.station_info["id"])
        for record in records
    }
    failures = _record_failures(
        records,
        checked_station_ids,
        request_start,
        request_end,
        max_gap,
    )
    required_failures = {
        station_id: reason
        for station_id, reason in failures.items()
        if station_id in required
    }
    if required_failures:
        raise RequiredStationDownloadError(
            f"Required USGS station(s) failed OGC validation for parameter "
            f"{param_id} after adaptive retries: {required_failures}"
        )

    optional_failures = {
        station_id: reason
        for station_id, reason in failures.items()
        if station_id in optional
    }
    if optional_failures and warning_log_path is not None:
        from stofs3d_setup.utils.warning_log import append_input_warning

        for station_id, reason in optional_failures.items():
            append_input_warning(
                warning_log_path,
                message=(
                    f"Optional Forced USGS "
                    f"{warning_forcing_name or param_id} replacement "
                    f"station {station_id} failed OGC validation after "
                    f"adaptive retries: {reason}"
                ),
                input_files=warning_input_files,
            )

    unusable = set(failures) - required
    records = [
        record
        for record in records
        if normalize_station_id(record.station_info["id"]) not in unusable
    ]

    if cache_fname is not None:
        Path(cache_fname).parent.mkdir(parents=True, exist_ok=True)
        _write_cache(records, cache_fname)
        _write_cache_manifest(
            cache_fname,
            station_ids,
            param_id,
            fetch_start,
            fetch_end,
        )
    return records
