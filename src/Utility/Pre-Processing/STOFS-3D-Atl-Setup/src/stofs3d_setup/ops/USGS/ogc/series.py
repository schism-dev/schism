"""High-level USGS time-series retrieval and interpolation helpers."""

from dataclasses import dataclass
from pathlib import Path
from typing import Any, Mapping

import numpy as np
import pandas as pd

from stofs3d_setup.ops.USGS.ogc.observation_adapter import (
    ObservationPersistenceError,
    download_stations,
)
from stofs3d_setup.ops.USGS.iv.download_usgs import find_data_gaps
from stofs3d_setup.ops.USGS.ogc.download_policy import (
    USGS_DOWNLOAD_POLICY,
    USGS_MAX_MISSING_GAP_DAYS,
)
from stofs3d_setup.ops.Source_sink.Replace_with_USGS.station_mappings import (
    CFS_TO_CMS,
    USGS_FLOW_PARAMETER_ID,
    USGS_TEMPERATURE_PARAMETER_ID,
)

MAX_USGS_INTERPOLATION_GAP = pd.Timedelta("6 hours")


@dataclass(frozen=True)
class UsgsDownloadPolicy:
    """Window sizes and cache naming for one-station USGS retrieval."""

    primary_chunk_days: int | None = None
    retry_chunk_days: tuple[int, ...] = ()
    cache_prefix: str = "usgs_adaptive"
    page_size: int = 1000

    def __post_init__(self) -> None:
        chunk_days = (
            ([] if self.primary_chunk_days is None else [self.primary_chunk_days])
            + list(self.retry_chunk_days)
        )
        if any(days <= 0 for days in chunk_days):
            raise ValueError("USGS download chunk sizes must be positive")
        if any(
            later >= earlier
            for earlier, later in zip(chunk_days, chunk_days[1:])
        ):
            raise ValueError(
                "USGS download chunk sizes must be strictly decreasing"
            )
        if not isinstance(self.cache_prefix, str):
            raise ValueError("USGS download cache_prefix must be a string")
        if not self.cache_prefix.strip():
            raise ValueError("USGS download cache_prefix must not be empty")
        if Path(self.cache_prefix).name != self.cache_prefix:
            raise ValueError(
                "USGS download cache_prefix must not contain a path"
            )
        if type(self.page_size) is not int or self.page_size <= 0:
            raise ValueError("USGS OGC page_size must be a positive integer")


DEFAULT_USGS_DOWNLOAD_POLICY = UsgsDownloadPolicy(
    primary_chunk_days=USGS_DOWNLOAD_POLICY["primary_chunk_days"],
    retry_chunk_days=tuple(USGS_DOWNLOAD_POLICY["retry_chunk_days"]),
    cache_prefix=USGS_DOWNLOAD_POLICY["cache_prefix"],
    page_size=USGS_DOWNLOAD_POLICY.get("page_size", 1000),
)


def _as_usgs_download_policy(
    value: UsgsDownloadPolicy | Mapping[str, Any] | None,
) -> UsgsDownloadPolicy:
    """Normalize an optional configured USGS download policy."""
    if value is None:
        return DEFAULT_USGS_DOWNLOAD_POLICY
    if isinstance(value, UsgsDownloadPolicy):
        return value
    policy = dict(value)
    supported_keys = {
        "primary_chunk_days",
        "retry_chunk_days",
        "cache_prefix",
        "page_size",
    }
    unknown_keys = set(policy) - supported_keys
    if unknown_keys:
        raise ValueError(
            "Unsupported USGS download policy setting(s): "
            f"{sorted(unknown_keys)}"
        )
    retry_chunk_days = policy.get("retry_chunk_days", ())
    return UsgsDownloadPolicy(
        primary_chunk_days=policy.get("primary_chunk_days"),
        retry_chunk_days=tuple(retry_chunk_days or ()),
        cache_prefix=policy.get(
            "cache_prefix",
            DEFAULT_USGS_DOWNLOAD_POLICY.cache_prefix,
        ),
        page_size=policy.get(
            "page_size", DEFAULT_USGS_DOWNLOAD_POLICY.page_size
        ),
    )


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




def _iter_date_windows(start_time, end_time, chunk_days: int | None):
    """Yield numbered, non-overlapping inclusive date windows."""
    start = _as_utc_timestamp(start_time)
    end = _as_utc_timestamp(end_time)
    if end < start:
        raise ValueError("USGS download end time precedes start time")
    if chunk_days is None:
        yield 1, start, end
        return
    if chunk_days <= 0:
        raise ValueError("USGS download chunk_days must be positive")

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


def _download_usgs_window(
    station_id: str,
    parameter_id: str,
    window_start: pd.Timestamp,
    window_end: pd.Timestamp,
    cache_dir: Path,
    cache_prefix: str,
    observation_window_callback=None,
    observation_window_start_callback=None,
) -> "UsgsWindowResult":
    """Download one window and distinguish empty data from request failure."""
    station_id = _normalize_usgs_station_id(station_id)

    cache_file = cache_dir / (
        f"{cache_prefix}_{station_id}_{parameter_id}_"
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
            allow_adaptive_retry=False,
            raise_on_fetch_error=True,
            validate_coverage=False,
            observation_window_callback=observation_window_callback,
            observation_window_start_callback=(
                observation_window_start_callback
            ),
        )
    except Exception as exc:
        if isinstance(exc, ObservationPersistenceError):
            raise
        print(
            f"[USGS SERIES] warning: request failed "
            f"for station {station_id}, parameter {parameter_id}, "
            f"{window_start:%Y-%m-%d} to {window_end:%Y-%m-%d}: {exc}"
        )
        return UsgsWindowResult(series=None, retryable=True)

    matching = [
        record for record in records
        if _station_id_from_record(record) == station_id
    ]
    if not matching:
        return UsgsWindowResult(series=None, retryable=False)

    try:
        series = _extract_usgs_values(matching[0])
    except Exception as exc:
        print(
            f"[USGS SERIES] warning: failed to parse station {station_id}, "
            f"parameter {parameter_id}: {exc}"
        )
        return UsgsWindowResult(series=None, retryable=True)

    if series.empty:
        return UsgsWindowResult(series=None, retryable=False)
    return UsgsWindowResult(series=series, retryable=False)


@dataclass(frozen=True)
class UsgsWindowResult:
    """Observation outcome for one window; empty and failed are distinct."""

    series: pd.Series | None
    retryable: bool


def _download_usgs_window_adaptive(
    station_id: str,
    parameter_id: str,
    window_start: pd.Timestamp,
    window_end: pd.Timestamp,
    cache_dir: Path,
    cache_prefix: str,
    retry_chunk_days: tuple[int, ...],
    observation_window_callback=None,
    observation_window_start_callback=None,
) -> list[pd.Series]:
    """Retry request failures using successively smaller configured windows."""
    if (
        observation_window_callback is None
        and observation_window_start_callback is None
    ):
        result = _download_usgs_window(
            station_id,
            parameter_id,
            window_start,
            window_end,
            cache_dir,
            cache_prefix,
        )
    else:
        result = _download_usgs_window(
            station_id,
            parameter_id,
            window_start,
            window_end,
            cache_dir,
            cache_prefix,
            observation_window_callback=observation_window_callback,
            observation_window_start_callback=(
                observation_window_start_callback
            ),
        )
    if result.series is not None:
        return [result.series]
    if not result.retryable:
        return []

    span_days = (
        window_end.normalize() - window_start.normalize()
    ).days + 1
    smaller_retry_windows = tuple(
        days for days in retry_chunk_days if days < span_days
    )
    if not smaller_retry_windows:
        print(
            f"[USGS SERIES] giving up: station {station_id}, parameter "
            f"{parameter_id}, "
            f"{window_start:%Y-%m-%d} to {window_end:%Y-%m-%d}"
        )
        return []

    next_chunk_days, *remaining_chunk_days = smaller_retry_windows
    pieces = []
    for _, retry_start, retry_end in _iter_date_windows(
        window_start,
        window_end,
        next_chunk_days,
    ):
        pieces.extend(
            _download_usgs_window_adaptive(
                station_id=station_id,
                parameter_id=parameter_id,
                window_start=retry_start,
                window_end=retry_end,
                cache_dir=cache_dir,
                cache_prefix=cache_prefix,
                retry_chunk_days=tuple(remaining_chunk_days),
                observation_window_callback=observation_window_callback,
                observation_window_start_callback=(
                    observation_window_start_callback
                ),
            )
        )
    return pieces


def download_usgs_series(
    station_id: str,
    parameter_id: str,
    start_time,
    end_time,
    usgs_cache_folder,
    policy: UsgsDownloadPolicy | Mapping[str, Any] | None = None,
    observation_window_callback=None,
    observation_window_start_callback=None,
) -> pd.Series | None:
    """Download one station using an optional adaptive chunk policy."""
    policy = _as_usgs_download_policy(policy)
    cache_dir = Path(usgs_cache_folder)
    cache_dir.mkdir(parents=True, exist_ok=True)

    start = _as_utc_timestamp(start_time) - pd.Timedelta(days=1)
    end = _as_utc_timestamp(end_time) + pd.Timedelta(days=1)

    pieces = []
    windows = list(_iter_date_windows(start, end, policy.primary_chunk_days))

    print(
        f"[USGS SERIES] download for {station_id}/{parameter_id}: "
        f"{len(windows)} primary window(s), retry windows "
        f"{policy.retry_chunk_days or 'disabled'}"
    )

    for _, window_start, window_end in windows:
        window_pieces = _download_usgs_window_adaptive(
            station_id=station_id,
            parameter_id=parameter_id,
            window_start=window_start,
            window_end=window_end,
            cache_dir=cache_dir,
            cache_prefix=policy.cache_prefix,
            retry_chunk_days=policy.retry_chunk_days,
            observation_window_callback=observation_window_callback,
            observation_window_start_callback=(
                observation_window_start_callback
            ),
        )
        pieces.extend(window_pieces)

    merged = _merge_usgs_series(pieces)
    if merged is None or merged.empty:
        return None

    gap_windows = find_data_gaps(
        merged.index,
        start=start,
        end=end,
        max_gap=pd.Timedelta(days=USGS_MAX_MISSING_GAP_DAYS),
        freq=pd.Timedelta(hours=1),
    )
    for gap_start, gap_end in gap_windows:
        for _, window_start, window_end in _iter_date_windows(
            gap_start,
            gap_end,
            policy.primary_chunk_days,
        ):
            pieces.extend(
                _download_usgs_window_adaptive(
                    station_id=station_id,
                    parameter_id=parameter_id,
                    window_start=window_start,
                    window_end=window_end,
                    cache_dir=cache_dir,
                    cache_prefix=policy.cache_prefix,
                    retry_chunk_days=policy.retry_chunk_days,
                    observation_window_callback=observation_window_callback,
                    observation_window_start_callback=(
                        observation_window_start_callback
                    ),
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



def _download_configured_usgs_series(
    station_id: str | None,
    parameter_id: str,
    start_time,
    model_time: np.ndarray,
    usgs_cache_folder,
    download_policy: UsgsDownloadPolicy | Mapping[str, Any] | None = None,
) -> tuple[pd.Series | None, str | None]:
    """Read or fetch one configured station series through the shared store."""
    if not station_id:
        return None, None

    station_id = _normalize_usgs_station_id(station_id)
    target_time = _model_datetimes(start_time, model_time)
    from stofs3d_setup.ops.USGS.ogc.store import (
        DEFAULT_USGS_CACHE_FOLDER,
        get_usgs_observations,
    )

    policy = _as_usgs_download_policy(download_policy)
    store_policy = {
        "primary_chunk_days": policy.primary_chunk_days,
        "retry_chunk_days": list(policy.retry_chunk_days),
        "cache_prefix": policy.cache_prefix,
        "page_size": policy.page_size,
    }
    data_start = target_time[0].normalize() - pd.Timedelta(days=1)
    data_end = target_time[-1].normalize() + pd.Timedelta(days=1)
    records = get_usgs_observations(
        param_id=str(parameter_id),
        station_ids=[station_id],
        database_dir=usgs_cache_folder or DEFAULT_USGS_CACHE_FOLDER,
        datelist=pd.date_range(data_start, data_end, freq="D"),
        download_policy=store_policy,
        validation_start=target_time[0],
        validation_end=target_time[-1],
    )
    matching = [
        record for record in records
        if _station_id_from_record(record) == station_id
    ]
    if not matching:
        return None, station_id
    series = _extract_usgs_values(matching[0])
    return series, station_id


def _get_usgs_flow(
    name: str,
    station_id: str | None,
    start_time,
    model_time: np.ndarray,
    usgs_cache_folder,
    download_policy: UsgsDownloadPolicy | Mapping[str, Any] | None = None,
) -> tuple[pd.Series | None, str | None]:
    """Return streamflow observations for the station configured by caller."""
    if not station_id:
        print(
            f"[USGS SERIES] {name}: "
            "no USGS flow station configured."
        )
    return _download_configured_usgs_series(
        station_id=station_id,
        parameter_id=USGS_FLOW_PARAMETER_ID,
        start_time=start_time,
        model_time=model_time,
        usgs_cache_folder=usgs_cache_folder,
        download_policy=download_policy,
    )


def _get_usgs_temperature(
    name: str,
    station_id: str | None,
    start_time,
    model_time: np.ndarray,
    usgs_cache_folder,
    download_policy: UsgsDownloadPolicy | Mapping[str, Any] | None = None,
) -> tuple[pd.Series | None, str | None]:
    """Return temperature observations for the station configured by caller."""
    if not station_id:
        print(
            f"[USGS SERIES] {name}: "
            "no USGS temperature station configured."
        )
    return _download_configured_usgs_series(
        station_id=station_id,
        parameter_id=USGS_TEMPERATURE_PARAMETER_ID,
        start_time=start_time,
        model_time=model_time,
        usgs_cache_folder=usgs_cache_folder,
        download_policy=download_policy,
    )


def _get_usgs_forcing(
    name: str,
    flow_station_id: str | None,
    temperature_station_id: str | None,
    start_time,
    model_time: np.ndarray,
    usgs_cache_folder,
    download_policy: UsgsDownloadPolicy | Mapping[str, Any] | None = None,
) -> tuple[
    pd.Series | None,
    pd.Series | None,
    str | None,
    str | None,
]:
    """Return explicitly configured USGS flow and temperature series."""
    flow_series, flow_station_id = _get_usgs_flow(
        name=name,
        station_id=flow_station_id,
        start_time=start_time,
        model_time=model_time,
        usgs_cache_folder=usgs_cache_folder,
        download_policy=download_policy,
    )
    temperature_series, temperature_station_id = _get_usgs_temperature(
        name=name,
        station_id=temperature_station_id,
        start_time=start_time,
        model_time=model_time,
        usgs_cache_folder=usgs_cache_folder,
        download_policy=download_policy,
    )

    return (
        flow_series,
        temperature_series,
        flow_station_id,
        temperature_station_id,
    )
