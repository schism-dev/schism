"""High-level USGS time-series retrieval and interpolation helpers."""

from pathlib import Path

import numpy as np
import pandas as pd

from stofs3d_setup.ops.Source_sink.Replace_with_USGS.download_usgs import (
    download_stations,
)
from stofs3d_setup.ops.Source_sink.Replace_with_USGS.station_mappings import (
    CFS_TO_CMS,
    USGS_FLOW_PARAMETER_ID,
    USGS_FLOW_STATION_BY_NAME,
    USGS_TEMPERATURE_PARAMETER_ID,
    USGS_TEMPERATURE_STATION_BY_NAME,
)

MAX_USGS_INTERPOLATION_GAP = pd.Timedelta("6 hours")
HUDSON_USGS_DOWNLOAD_CHUNK_DAYS = 100
HUDSON_USGS_RETRY_CHUNK_DAYS = 20
HUDSON_USGS_MIN_CHUNK_DAYS = 5


def _normalize_usgs_station_id(value) -> str:
    """Normalize a USGS station ID while preserving leading zeros."""
    if value is None:
        return ""

    value = str(value).strip()

    if value.endswith(".0"):
        value = value[:-2]

    if value.isdigit():
        value = value.zfill(8)

    return value


def _as_utc_timestamp(value) -> pd.Timestamp:
    """Convert a datetime-like value to a timezone-aware UTC Timestamp."""
    ts = pd.Timestamp(value)
    if ts.tzinfo is None:
        return ts.tz_localize("UTC")
    return ts.tz_convert("UTC")


def _model_datetimes(start_time, model_time: np.ndarray) -> pd.DatetimeIndex:
    """Convert SCHISM elapsed seconds to UTC datetimes."""
    start = _as_utc_timestamp(start_time)
    seconds = np.asarray(model_time, dtype=float).reshape(-1)
    return pd.DatetimeIndex(start + pd.to_timedelta(seconds, unit="s"))


def _station_id_from_record(record) -> str:
    """Return a normalized USGS station ID from a download record."""
    info = getattr(record, "station_info", {}) or {}
    return _normalize_usgs_station_id(info.get("id", ""))


def _extract_usgs_values(record) -> pd.Series:
    """
    Return a cleaned, time-indexed USGS value series.

    download_stations() normally returns a dataframe with ``date`` and
    ``value``. A numeric fallback is retained for minor API-format changes.
    """
    df = record.df.copy()
    if "date" not in df.columns:
        raise ValueError("USGS record does not contain a 'date' column")

    index = pd.DatetimeIndex(pd.to_datetime(df["date"], utc=True))

    if "value" in df.columns:
        values = pd.to_numeric(df["value"], errors="coerce")
    else:
        candidate_columns = [
            col for col in df.columns
            if col != "date" and pd.api.types.is_numeric_dtype(
                pd.to_numeric(df[col], errors="coerce")
            )
        ]
        if not candidate_columns:
            raise ValueError("USGS record does not contain a numeric value column")
        values = pd.to_numeric(df[candidate_columns[0]], errors="coerce")

    series = pd.Series(values.to_numpy(dtype=float), index=index)
    series = series.replace([np.inf, -np.inf], np.nan).dropna()
    series = series[~series.index.duplicated(keep="first")].sort_index()
    return series




def _download_usgs_series_standard(
    station_id: str,
    parameter_id: str,
    start_time,
    end_time,
    usgs_cache_folder,
) -> pd.Series | None:
    """Standard full-period download used for all locations except Hudson."""
    cache_dir = Path(usgs_cache_folder)
    cache_dir.mkdir(parents=True, exist_ok=True)

    station_id = _normalize_usgs_station_id(station_id)
    start = _as_utc_timestamp(start_time)
    end = _as_utc_timestamp(end_time)
    padded_start = start - pd.Timedelta(days=1)
    padded_end = end + pd.Timedelta(days=1)

    cache_file = cache_dir / (
        f"artificial_island_usgs_{station_id}_{parameter_id}_"
        f"{padded_start.strftime('%Y%m%d')}_"
        f"{padded_end.strftime('%Y%m%d')}.pq"
    )

    try:
        records = download_stations(
            param_id=str(parameter_id),
            station_ids=[station_id],
            cache_fname=str(cache_file),
            datelist=pd.date_range(
                start=padded_start.tz_localize(None),
                end=padded_end.tz_localize(None),
                freq="D",
            ),
        )
    except Exception as exc:
        print(
            f"[ARTIFICIAL ISLAND PATCH] warning: USGS download failed "
            f"for station {station_id}, parameter {parameter_id}: {exc}"
        )
        return None

    matching = [
        record for record in records
        if _station_id_from_record(record) == station_id
    ]
    if not matching:
        print(
            f"[ARTIFICIAL ISLAND PATCH] warning: USGS station "
            f"{station_id}, parameter {parameter_id} returned no data."
        )
        return None

    try:
        series = _extract_usgs_values(matching[0])
    except Exception as exc:
        print(
            f"[ARTIFICIAL ISLAND PATCH] warning: failed to parse USGS "
            f"station {station_id}, parameter {parameter_id}: {exc}"
        )
        return None

    return None if series.empty else series


def _iter_hudson_chunks(start_time, end_time, chunk_days=HUDSON_USGS_DOWNLOAD_CHUNK_DAYS):
    """Yield Hudson-only non-overlapping date windows."""
    start = _as_utc_timestamp(start_time)
    end = _as_utc_timestamp(end_time)

    chunk_start = start
    n = 1
    while chunk_start <= end:
        chunk_end = min(
            chunk_start + pd.Timedelta(days=chunk_days - 1),
            end,
        )
        yield n, chunk_start, chunk_end
        chunk_start = chunk_end + pd.Timedelta(days=1)
        n += 1


def _merge_usgs_series(series_list):
    """Merge successful USGS pieces preserving original observation times."""
    usable = [s for s in series_list if s is not None and not s.empty]
    if not usable:
        return None
    merged = pd.concat(usable).sort_index()
    merged = merged[~merged.index.duplicated(keep="first")]
    return merged.replace([np.inf, -np.inf], np.nan).dropna()


def _download_hudson_window(
    station_id,
    parameter_id,
    window_start,
    window_end,
    cache_dir,
):
    """Download one Hudson window."""
    station_id = _normalize_usgs_station_id(station_id)

    cache_file = cache_dir / (
        f"hudson_usgs_{station_id}_{parameter_id}_"
        f"{window_start.strftime('%Y%m%d')}_"
        f"{window_end.strftime('%Y%m%d')}.pq"
    )

    try:
        records = download_stations(
            param_id=str(parameter_id),
            station_ids=[station_id],
            cache_fname=str(cache_file),
            datelist=pd.date_range(
                start=window_start.tz_localize(None),
                end=window_end.tz_localize(None),
                freq="D",
            ),
        )
    except Exception as exc:
        print(
            f"[HUDSON USGS] warning: request failed "
            f"{window_start:%Y-%m-%d} to {window_end:%Y-%m-%d}: {exc}"
        )
        return None

    matching = [
        record for record in records
        if _station_id_from_record(record) == station_id
    ]
    if not matching:
        return None

    try:
        series = _extract_usgs_values(matching[0])
    except Exception:
        return None

    return None if series.empty else series


def _download_hudson_window_adaptive(
    station_id,
    parameter_id,
    window_start,
    window_end,
    cache_dir,
    retry_chunk_days=HUDSON_USGS_RETRY_CHUNK_DAYS,
    min_chunk_days=HUDSON_USGS_MIN_CHUNK_DAYS,
):
    """Retry failed Hudson windows as 20-day, then 10-day, then 5-day pieces."""
    series = _download_hudson_window(
        station_id,
        parameter_id,
        window_start,
        window_end,
        cache_dir,
    )
    if series is not None:
        return [series]

    span_days = (window_end.normalize() - window_start.normalize()).days + 1
    if span_days <= min_chunk_days:
        print(
            f"[HUDSON USGS] giving up: "
            f"{window_start:%Y-%m-%d} to {window_end:%Y-%m-%d}"
        )
        return []

    if span_days > retry_chunk_days:
        pieces = []
        for _, s, e in _iter_hudson_chunks(
            window_start,
            window_end,
            chunk_days=retry_chunk_days,
        ):
            pieces.extend(
                _download_hudson_window_adaptive(
                    station_id,
                    parameter_id,
                    s,
                    e,
                    cache_dir,
                    retry_chunk_days,
                    min_chunk_days,
                )
            )
        return pieces

    left_days = span_days // 2
    left_end = window_start + pd.Timedelta(days=left_days - 1)
    right_start = left_end + pd.Timedelta(days=1)

    pieces = _download_hudson_window_adaptive(
        station_id,
        parameter_id,
        window_start,
        left_end,
        cache_dir,
        retry_chunk_days,
        min_chunk_days,
    )
    pieces.extend(
        _download_hudson_window_adaptive(
            station_id,
            parameter_id,
            right_start,
            window_end,
            cache_dir,
            retry_chunk_days,
            min_chunk_days,
        )
    )
    return pieces


def _download_hudson_usgs_series(
    station_id,
    parameter_id,
    start_time,
    end_time,
    usgs_cache_folder,
):
    """Hudson-only 100d -> 20d -> 10d -> 5d adaptive downloader."""
    cache_dir = Path(usgs_cache_folder)
    cache_dir.mkdir(parents=True, exist_ok=True)

    start = _as_utc_timestamp(start_time) - pd.Timedelta(days=1)
    end = _as_utc_timestamp(end_time) + pd.Timedelta(days=1)

    pieces = []
    chunks = list(_iter_hudson_chunks(start, end))

    print(
        f"[HUDSON USGS] adaptive download for {station_id}/{parameter_id}: "
        f"{len(chunks)} primary chunk(s)"
    )

    for _, chunk_start, chunk_end in chunks:
        pieces.extend(
            _download_hudson_window_adaptive(
                station_id,
                parameter_id,
                chunk_start,
                chunk_end,
                cache_dir,
            )
        )

    return _merge_usgs_series(pieces)


def _interpolate_usgs_with_original_fallback(
    series: pd.Series | None,
    target_time: pd.DatetimeIndex,
    original_values: np.ndarray,
    scale: float = 1.0,
    max_interp_gap: pd.Timedelta = MAX_USGS_INTERPOLATION_GAP,
) -> tuple[np.ndarray, np.ndarray]:
    """
    Blend USGS observations with original forcing.

    Exact USGS observations and interpolation across short observation gaps
    are used. Model times outside the observation period, inside long gaps,
    or otherwise unsupported retain original forcing.

    Returns
    -------
    blended
        Final values containing USGS where supported and original forcing
        elsewhere.
    use_usgs
        Boolean mask identifying records supplied by USGS.
    """
    original_values = np.asarray(original_values, dtype=float).reshape(-1)
    target_time = pd.DatetimeIndex(target_time)

    if len(original_values) != len(target_time):
        raise ValueError(
            "Original forcing length does not match model datetime length"
        )

    blended = original_values.copy()
    use_usgs = np.zeros(len(target_time), dtype=bool)

    if series is None:
        return blended, use_usgs

    series = series.copy()
    series = series.replace([np.inf, -np.inf], np.nan).dropna()
    series = series[~series.index.duplicated(keep="first")].sort_index()

    if series.empty:
        return blended, use_usgs

    if series.index.tz is None:
        series.index = series.index.tz_localize("UTC")
    else:
        series.index = series.index.tz_convert("UTC")

    if target_time.tz is None:
        target_time = target_time.tz_localize("UTC")
    else:
        target_time = target_time.tz_convert("UTC")

    obs_time_ns = series.index.asi8
    target_time_ns = target_time.asi8
    obs_values = series.to_numpy(dtype=float) * float(scale)

    right_idx = np.searchsorted(
        obs_time_ns,
        target_time_ns,
        side="left",
    )
    left_idx = right_idx - 1

    safe_right = np.minimum(right_idx, len(obs_time_ns) - 1)
    exact = (
        (right_idx < len(obs_time_ns))
        & (obs_time_ns[safe_right] == target_time_ns)
    )

    if np.any(exact):
        exact_obs_idx = right_idx[exact]
        blended[exact] = obs_values[exact_obs_idx]
        use_usgs[exact] = True

    between = (
        (left_idx >= 0)
        & (right_idx < len(obs_time_ns))
        & (~exact)
    )

    candidate_idx = np.where(between)[0]
    if len(candidate_idx) > 0:
        left_obs_idx = left_idx[candidate_idx]
        right_obs_idx = right_idx[candidate_idx]

        gap_ns = (
            obs_time_ns[right_obs_idx]
            - obs_time_ns[left_obs_idx]
        )
        short_gap = gap_ns <= int(max_interp_gap.value)
        valid_target_idx = candidate_idx[short_gap]

        if len(valid_target_idx) > 0:
            left_obs_idx = left_idx[valid_target_idx]
            right_obs_idx = right_idx[valid_target_idx]

            denominator = (
                obs_time_ns[right_obs_idx]
                - obs_time_ns[left_obs_idx]
            )
            fraction = (
                target_time_ns[valid_target_idx]
                - obs_time_ns[left_obs_idx]
            ) / denominator

            interpolated = (
                obs_values[left_obs_idx]
                + fraction
                * (
                    obs_values[right_obs_idx]
                    - obs_values[left_obs_idx]
                )
            )

            blended[valid_target_idx] = interpolated
            use_usgs[valid_target_idx] = True

    return blended, use_usgs



def _get_usgs_forcing(
    name: str,
    start_time,
    model_time: np.ndarray,
    usgs_cache_folder,
) -> tuple[
    pd.Series | None,
    pd.Series | None,
    str | None,
    str | None,
]:
    """Return explicit USGS flow/temperature; chunk only Hudson River."""
    flow_station_id = USGS_FLOW_STATION_BY_NAME.get(name)
    temperature_station_id = USGS_TEMPERATURE_STATION_BY_NAME.get(name)

    target_time = _model_datetimes(start_time, model_time)
    period_start = target_time[0]
    period_end = target_time[-1]

    flow_series = None
    temperature_series = None

    downloader = (
        _download_hudson_usgs_series
        if name == "Hudson River"
        else _download_usgs_series_standard
    )

    if flow_station_id:
        flow_station_id = _normalize_usgs_station_id(flow_station_id)
        flow_series = downloader(
            flow_station_id,
            USGS_FLOW_PARAMETER_ID,
            period_start,
            period_end,
            usgs_cache_folder,
        )
    else:
        print(
            f"[ARTIFICIAL ISLAND PATCH] {name}: "
            "no USGS flow station configured."
        )

    if temperature_station_id:
        temperature_station_id = _normalize_usgs_station_id(
            temperature_station_id
        )
        temperature_series = downloader(
            temperature_station_id,
            USGS_TEMPERATURE_PARAMETER_ID,
            period_start,
            period_end,
            usgs_cache_folder,
        )
    else:
        print(
            f"[ARTIFICIAL ISLAND PATCH] {name}: "
            "no USGS temperature station configured."
        )

    return (
        flow_series,
        temperature_series,
        flow_station_id,
        temperature_station_id,
    )
