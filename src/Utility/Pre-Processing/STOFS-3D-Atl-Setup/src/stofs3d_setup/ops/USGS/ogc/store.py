"""Persistent, resumable storage for USGS observation downloads.

The store is an immutable-part Parquet dataset with small JSON coverage
records. It is intentionally file-based so readers never need to hold a
database-wide lock on the shared NFS cache.
"""

from __future__ import annotations

from collections import defaultdict
from datetime import timedelta
import hashlib
import json
import os
from pathlib import Path
import tempfile
from typing import Any, Iterable, Mapping
import uuid

import numpy as np
import pandas as pd

from stofs3d_setup.ops.USGS.iv.download_usgs import (
    GenericObsData,
    RequiredStationDownloadError,
)
from stofs3d_setup.ops.USGS.ogc.download_policy import (
    USGS_DOWNLOAD_POLICY,
    USGS_MAX_MISSING_GAP_DAYS,
)
from stofs3d_setup.ops.USGS.ogc.download_usgs import normalize_station_id


DEFAULT_USGS_CACHE_FOLDER = Path(
    "/sciclone/schism10/Hgrid_projects/USGS_Cache"
)
STORE_FORMAT_VERSION = 1
LEGACY_IMPORT_VERSION = 1
PENDING_REQUEST_TIMEOUT = timedelta(hours=3)
_OBSERVATION_COLUMNS = (
    "station_id",
    "parameter",
    "time",
    "value",
    "unit",
    "qualifier",
    "time_series_id",
    "station_name",
    "lon",
    "lat",
    "retrieved_at",
    "ingest_source",
)


def _as_utc(value: Any) -> pd.Timestamp:
    timestamp = pd.Timestamp(value)
    if timestamp.tzinfo is None:
        return timestamp.tz_localize("UTC")
    return timestamp.tz_convert("UTC")


def _day(value: Any) -> pd.Timestamp:
    return _as_utc(value).normalize()


def _time_key(value: Any) -> int:
    return int(_as_utc(value).value)


def _station_ids(values: Iterable[str] | None) -> list[str]:
    return sorted({
        normalize_station_id(value)
        for value in (values or ())
        if value is not None and str(value).strip()
    })


def _atomic_json(path: Path, payload: Mapping[str, Any]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_name(f".{path.name}.{uuid.uuid4().hex}.tmp")
    temporary.write_text(
        json.dumps(payload, indent=2, sort_keys=True), encoding="utf-8"
    )
    os.replace(temporary, path)


def _atomic_parquet(frame: pd.DataFrame, path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_name(f".{path.name}.{uuid.uuid4().hex}.tmp")
    try:
        frame.to_parquet(
            temporary, engine="pyarrow", compression="zstd", index=False
        )
        os.replace(temporary, path)
    finally:
        temporary.unlink(missing_ok=True)


def _empty_observations() -> pd.DataFrame:
    return pd.DataFrame(columns=_OBSERVATION_COLUMNS)


class USGSObservationStore:
    """Append-only observations and resumable per-station coverage on disk."""

    def __init__(self, root: str | Path = DEFAULT_USGS_CACHE_FOLDER):
        self.root = Path(root).expanduser()
        self.observation_root = self.root / "observations"
        self.coverage_root = self.root / "coverage"
        self.import_root = self.root / "imports"

    def _station_observation_root(
        self, parameter: str, station_id: str
    ) -> Path:
        return (
            self.observation_root
            / f"parameter={parameter}"
            / f"station={station_id}"
        )

    def _station_coverage_root(
        self, parameter: str, station_id: str
    ) -> Path:
        return (
            self.coverage_root
            / f"parameter={parameter}"
            / f"station={station_id}"
        )

    def _coverage_path(
        self,
        parameter: str,
        station_id: str,
        start: Any,
        end: Any,
    ) -> Path:
        start_day = _day(start)
        end_day = _day(end)
        return self._station_coverage_root(parameter, station_id) / (
            f"{start_day:%Y%m%d}_{end_day:%Y%m%d}.json"
        )

    def record_window(
        self,
        station_ids: Iterable[str],
        parameter: str,
        start: Any,
        end: Any,
        observations: pd.DataFrame | None,
        error: BaseException | None = None,
        *,
        source: str = "ogc-api",
    ) -> None:
        """Persist one successful or failed station/window request atomically."""
        station_ids = _station_ids(station_ids)
        parameter = str(parameter)
        start_utc = _as_utc(start)
        end_utc = _as_utc(end)
        if end_utc < start_utc:
            raise ValueError("USGS store window end precedes its start")

        frame = self._normalize_observations(
            observations, parameter=parameter, source=source
        )
        for station_id in station_ids:
            station_frame = frame.loc[frame["station_id"] == station_id]
            if not station_frame.empty:
                self._write_station_window(
                    station_frame,
                    parameter=parameter,
                    station_id=station_id,
                    start=start_utc,
                    end=end_utc,
                )
            state = "failed" if error is not None else "complete"
            payload = {
                "format_version": STORE_FORMAT_VERSION,
                "parameter": parameter,
                "station_id": station_id,
                "start": _day(start_utc).isoformat(),
                "end": _day(end_utc).isoformat(),
                "status": state,
                "attempted_at": pd.Timestamp.now(tz="UTC").isoformat(),
                "source": source,
            }
            if error is not None:
                # Avoid persisting request URLs, headers, or API credentials.
                payload["error_type"] = type(error).__name__
            _atomic_json(
                self._coverage_path(
                    parameter, station_id, start_utc, end_utc
                ),
                payload,
            )

    def mark_pending(
        self,
        station_ids: Iterable[str],
        parameter: str,
        start: Any,
        end: Any,
    ) -> None:
        """Record an in-flight request before waiting for the API response."""
        for station_id in _station_ids(station_ids):
            self._write_coverage(
                parameter=str(parameter),
                station_id=station_id,
                start=start,
                end=end,
                status="pending",
                source="ogc-api",
            )

    @staticmethod
    def _normalize_observations(
        observations: pd.DataFrame | None,
        parameter: str,
        source: str,
    ) -> pd.DataFrame:
        if observations is None or observations.empty:
            return _empty_observations()

        frame = observations.copy()
        if "station_id" not in frame and "site_no" in frame:
            frame = frame.rename(columns={"site_no": "station_id"})
        if "time" not in frame and "date" in frame:
            frame = frame.rename(columns={"date": "time"})
        if "station_id" not in frame or "time" not in frame:
            raise ValueError(
                "USGS observations require station_id and time columns"
            )
        frame["station_id"] = frame["station_id"].map(normalize_station_id)
        frame["parameter"] = str(parameter)
        frame["time"] = pd.to_datetime(frame["time"], utc=True, errors="coerce")
        frame["value"] = pd.to_numeric(frame.get("value"), errors="coerce")
        for column in _OBSERVATION_COLUMNS:
            if column not in frame:
                frame[column] = None
        frame["retrieved_at"] = pd.Timestamp.now(tz="UTC")
        frame["ingest_source"] = source
        frame["qualifier"] = frame["qualifier"].map(
            lambda value: " | ".join(map(str, value))
            if isinstance(value, (list, tuple, set))
            else value
        )
        frame["lon"] = pd.to_numeric(frame["lon"], errors="coerce")
        frame["lat"] = pd.to_numeric(frame["lat"], errors="coerce")
        return frame.loc[
            frame["time"].notna() & frame["value"].notna(),
            list(_OBSERVATION_COLUMNS),
        ].reset_index(drop=True)

    def _write_station_window(
        self,
        frame: pd.DataFrame,
        parameter: str,
        station_id: str,
        start: pd.Timestamp,
        end: pd.Timestamp,
    ) -> None:
        start_key = _time_key(start)
        end_key = _time_key(end)
        part_id = uuid.uuid4().hex
        by_month = frame.groupby(
            [frame["time"].dt.year, frame["time"].dt.month], sort=True
        )
        for (year, month), month_frame in by_month:
            path = (
                self._station_observation_root(parameter, station_id)
                / f"year={int(year):04d}"
                / f"month={int(month):02d}"
                / f"part_{start_key}_{end_key}_{part_id}.parquet"
            )
            _atomic_parquet(month_frame, path)

    def _coverage_records(
        self, parameter: str, station_id: str
    ) -> list[dict[str, Any]]:
        folder = self._station_coverage_root(parameter, station_id)
        records = []
        if not folder.is_dir():
            return records
        for path in folder.glob("*.json"):
            try:
                record = json.loads(path.read_text(encoding="utf-8"))
            except (OSError, json.JSONDecodeError):
                continue
            if (
                record.get("format_version") == STORE_FORMAT_VERSION
                and record.get("parameter") == str(parameter)
                and record.get("station_id") == station_id
            ):
                records.append(record)
        return records

    def covered_intervals(
        self,
        parameter: str,
        station_id: str,
        *,
        retry_failed_after: timedelta | None = None,
        retry_pending_after: timedelta | None = None,
        now: pd.Timestamp | None = None,
    ) -> list[tuple[pd.Timestamp, pd.Timestamp]]:
        """Return completed ranges and, optionally, recently failed ranges."""
        now = _as_utc(now or pd.Timestamp.now(tz="UTC"))
        intervals = []
        for record in self._coverage_records(
            str(parameter), normalize_station_id(station_id)
        ):
            include = record.get("status") == "complete"
            if (
                not include
                and record.get("status") == "failed"
                and retry_failed_after is not None
            ):
                try:
                    attempted_at = _as_utc(record["attempted_at"])
                except (KeyError, TypeError, ValueError):
                    attempted_at = pd.Timestamp.min.tz_localize("UTC")
                include = now - attempted_at < retry_failed_after
            if (
                not include
                and record.get("status") == "pending"
                and retry_pending_after is not None
            ):
                try:
                    attempted_at = _as_utc(record["attempted_at"])
                except (KeyError, TypeError, ValueError):
                    attempted_at = pd.Timestamp.min.tz_localize("UTC")
                include = now - attempted_at < retry_pending_after
            if not include:
                continue
            try:
                intervals.append((_day(record["start"]), _day(record["end"])))
            except (KeyError, TypeError, ValueError):
                continue
        return intervals

    def missing_intervals(
        self,
        parameter: str,
        station_id: str,
        start: Any,
        end: Any,
        *,
        retry_failed_after: timedelta | None = None,
        retry_pending_after: timedelta | None = None,
        now: pd.Timestamp | None = None,
    ) -> list[tuple[pd.Timestamp, pd.Timestamp]]:
        """Find calendar-day gaps not already completed or recently attempted."""
        start_day = _day(start)
        end_day = _day(end)
        if end_day < start_day:
            raise ValueError("USGS requested coverage end precedes start")
        intervals = sorted(
            self.covered_intervals(
                parameter,
                station_id,
                retry_failed_after=retry_failed_after,
                retry_pending_after=retry_pending_after,
                now=now,
            ),
            key=lambda interval: interval[0],
        )
        merged: list[list[pd.Timestamp]] = []
        for interval_start, interval_end in intervals:
            interval_start = max(start_day, interval_start)
            interval_end = min(end_day, interval_end)
            if interval_end < interval_start:
                continue
            if not merged or interval_start > merged[-1][1] + pd.Timedelta(days=1):
                merged.append([interval_start, interval_end])
            else:
                merged[-1][1] = max(merged[-1][1], interval_end)

        missing = []
        cursor = start_day
        for interval_start, interval_end in merged:
            if interval_start > cursor:
                missing.append((cursor, interval_start - pd.Timedelta(days=1)))
            cursor = max(cursor, interval_end + pd.Timedelta(days=1))
        if cursor <= end_day:
            missing.append((cursor, end_day))
        return missing

    def has_coverage(
        self,
        parameter: str,
        station_id: str,
        start: Any,
        end: Any,
    ) -> bool:
        return not self.missing_intervals(parameter, station_id, start, end)

    def read_observations(
        self,
        parameter: str,
        station_ids: Iterable[str],
        start: Any,
        end: Any,
    ) -> pd.DataFrame:
        """Read a selected parameter/station/date slice from Parquet parts."""
        start_utc = _as_utc(start)
        end_utc = _as_utc(end)
        frames = []
        for station_id in _station_ids(station_ids):
            station_root = self._station_observation_root(
                str(parameter), station_id
            )
            months = pd.period_range(
                start_utc.tz_localize(None).to_period("M"),
                end_utc.tz_localize(None).to_period("M"),
                freq="M",
            )
            for month in months:
                month_root = (
                    station_root
                    / f"year={month.year:04d}"
                    / f"month={month.month:02d}"
                )
                for path in month_root.glob("*.parquet"):
                    try:
                        frame = pd.read_parquet(path, engine="pyarrow")
                    except (OSError, ValueError, ImportError):
                        continue
                    frame["time"] = pd.to_datetime(
                        frame["time"], utc=True, errors="coerce"
                    )
                    frames.append(frame)
        if not frames:
            return _empty_observations()

        observations = pd.concat(frames, ignore_index=True)
        observations = observations.loc[
            observations["time"].between(start_utc, end_utc, inclusive="both")
        ]
        observations["_source_rank"] = observations["ingest_source"].map(
            {"legacy-cache": 0, "ogc-api": 1}
        ).fillna(0)
        observations = observations.sort_values(
            ["_source_rank", "retrieved_at", "time"], kind="stable"
        ).drop_duplicates(
            ["station_id", "parameter", "time_series_id", "time"],
            keep="last",
        )
        return observations.drop(columns="_source_rank").reset_index(drop=True)

    def import_legacy_caches(
        self,
        parameter: str,
        station_ids: Iterable[str],
        start: Any,
        end: Any,
    ) -> int:
        """Import compatible flat cache files without removing or changing them.

        A recognized manifest can establish coverage only for station series
        that pass the current value and gap checks. Older or incomplete cache
        records may contribute observations without establishing coverage.
        """
        requested_ids = set(_station_ids(station_ids))
        if not requested_ids or not self.root.is_dir():
            return 0

        from stofs3d_setup.ops.USGS.ogc.observation_adapter import (
            _filter_temperature_records,
            _read_cache,
            _record_failures,
        )

        imported = 0
        for cache_path in self.root.glob("*.pq"):
            manifest_path = Path(f"{cache_path}.manifest.json")
            if not manifest_path.is_file():
                continue
            try:
                manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
                cache_parameter = str(manifest.get("parameter"))
                cache_start = _as_utc(manifest["start"])
                cache_end = _as_utc(manifest["end"])
                manifest_ids = set(_station_ids(manifest.get("stations")))
            except (OSError, json.JSONDecodeError, KeyError, TypeError, ValueError):
                continue
            if cache_parameter != str(parameter):
                continue
            if cache_end < _as_utc(start) or cache_start > _as_utc(end):
                continue
            selected_ids = requested_ids & manifest_ids
            if not selected_ids:
                continue

            try:
                stat = cache_path.stat()
                receipt_id = hashlib.sha256(
                    str(cache_path.resolve()).encode("utf-8")
                ).hexdigest()
                receipt_path = self.import_root / f"{receipt_id}.json"
                receipt = json.loads(receipt_path.read_text(encoding="utf-8"))
            except (OSError, json.JSONDecodeError):
                stat = cache_path.stat()
                receipt = {}
                receipt_path = self.import_root / (
                    hashlib.sha256(
                        str(cache_path.resolve()).encode("utf-8")
                    ).hexdigest()
                    + ".json"
                )
            already_imported = set(
                _station_ids(receipt.get("station_ids", ()))
            ) if (
                receipt.get("importer_version") == LEGACY_IMPORT_VERSION
                and receipt.get("size") == stat.st_size
                and receipt.get("mtime_ns") == stat.st_mtime_ns
            ) else set()
            ids_to_import = selected_ids - already_imported
            if not ids_to_import:
                continue

            try:
                records = _read_cache(cache_path)
            except (OSError, ValueError, ImportError):
                continue
            recognized_manifest = (
                manifest.get("api") == "usgs-waterdata-ogc-v1"
                or manifest.get("version") == 1
            )
            imported_ids = set(already_imported)
            for record in records:
                station_id = normalize_station_id(
                    record.station_info.get("id", "")
                )
                if station_id not in ids_to_import:
                    continue
                frame = record.df.copy()
                if "date" not in frame and frame.index.name == "date":
                    frame = frame.reset_index()
                frame = frame.rename(columns={"date": "time"})
                frame["station_id"] = station_id
                frame["station_name"] = record.station_info.get("name")
                frame["lon"] = record.station_info.get("lon")
                frame["lat"] = record.station_info.get("lat")
                frame["unit"] = record.station_info.get("unit")
                frame["time_series_id"] = None
                frame["qualifier"] = None
                self._write_station_window(
                    self._normalize_observations(
                        frame, parameter=str(parameter), source="legacy-cache"
                    ),
                    parameter=str(parameter),
                    station_id=station_id,
                    start=cache_start,
                    end=cache_end,
                )
                valid_records = _filter_temperature_records(
                    [record], str(parameter)
                )
                if (
                    recognized_manifest
                    and valid_records
                    and not _record_failures(
                        valid_records,
                        [station_id],
                        cache_start,
                        cache_end,
                        pd.Timedelta(days=USGS_MAX_MISSING_GAP_DAYS),
                    )
                ):
                    self._write_coverage(
                        parameter=str(parameter),
                        station_id=station_id,
                        start=cache_start,
                        end=cache_end,
                        status="complete",
                        source="legacy-cache",
                    )
                imported_ids.add(station_id)
                imported += 1

            _atomic_json(
                receipt_path,
                {
                    "importer_version": LEGACY_IMPORT_VERSION,
                    "path": str(cache_path),
                    "size": stat.st_size,
                    "mtime_ns": stat.st_mtime_ns,
                    "station_ids": sorted(imported_ids),
                },
            )
        return imported

    def _write_coverage(
        self,
        parameter: str,
        station_id: str,
        start: Any,
        end: Any,
        status: str,
        source: str,
    ) -> None:
        payload = {
            "format_version": STORE_FORMAT_VERSION,
            "parameter": str(parameter),
            "station_id": normalize_station_id(station_id),
            "start": _day(start).isoformat(),
            "end": _day(end).isoformat(),
            "status": status,
            "attempted_at": pd.Timestamp.now(tz="UTC").isoformat(),
            "source": source,
        }
        _atomic_json(
            self._coverage_path(parameter, station_id, start, end), payload
        )

    def status(self) -> dict[str, Any]:
        """Return compact storage counts for the standalone downloader CLI."""
        observation_parts = list(self.observation_root.rglob("*.parquet"))
        legacy_cache_files = list(self.root.glob("*.pq"))
        coverage_paths = list(self.coverage_root.rglob("*.json"))
        status_counts: dict[str, int] = defaultdict(int)
        for path in coverage_paths:
            try:
                record = json.loads(path.read_text(encoding="utf-8"))
            except (OSError, json.JSONDecodeError):
                continue
            status_counts[str(record.get("status", "unknown"))] += 1
        return {
            "root": str(self.root),
            "observation_parts": len(observation_parts),
            "observation_bytes": sum(
                path.stat().st_size for path in observation_parts
            ),
            "legacy_cache_files": len(legacy_cache_files),
            "legacy_cache_bytes": sum(
                path.stat().st_size for path in legacy_cache_files
            ),
            "coverage_windows": dict(sorted(status_counts.items())),
        }


def _warning_for_incomplete_station(
    station_id: str,
    parameter: str,
    warning_log_path: str | Path | None,
    warning_input_files: Iterable[str],
    warning_forcing_name: str | None,
    *,
    required: bool,
    missing: list[tuple[pd.Timestamp, pd.Timestamp]],
) -> None:
    description = ", ".join(
        f"{start:%Y-%m-%d} to {end:%Y-%m-%d}" for start, end in missing
    )
    message = (
        f"USGS {warning_forcing_name or parameter} station {station_id} has "
        f"incomplete cached coverage ({description})."
    )
    if required:
        raise RequiredStationDownloadError(message)
    print(f"[USGS STORE] warning: {message}", flush=True)
    if warning_log_path is not None:
        from stofs3d_setup.utils.warning_log import append_input_warning

        append_input_warning(
            warning_log_path,
            message=message,
            input_files=warning_input_files,
        )


def get_usgs_observations(
    param_id: str,
    station_ids: Iterable[str],
    database_dir: str | Path = DEFAULT_USGS_CACHE_FOLDER,
    datelist: pd.DatetimeIndex | None = None,
    *,
    required_station_ids: Iterable[str] | None = None,
    optional_station_ids: Iterable[str] | None = None,
    mode: str = "fetch_missing",
    retry_failed_after_hours: float = 24.0,
    stations_chunk_size: int = 10,
    days_per_chunk: int | None = None,
    download_policy: Mapping[str, Any] | None = None,
    max_missing_gap_days: float | None = None,
    validation_start: Any | None = None,
    validation_end: Any | None = None,
    warning_log_path: str | Path | None = None,
    warning_input_files: Iterable[str] = (),
    warning_forcing_name: str | None = None,
) -> list[GenericObsData]:
    """Return cached USGS data, fetching only uncovered station/date windows.

    ``offline`` mode only reads the persistent store. Successful empty queries
    count as covered; failed or rate-limited queries remain retryable after the
    configured cooldown. Coverage describes API query completion, not the
    consumer's later gap/interpolation policy.
    """
    from stofs3d_setup.ops.USGS.ogc import observation_adapter
    from stofs3d_setup.ops.USGS.ogc.observation_adapter import (
        ObservationPersistenceError,
    )

    if mode not in {"fetch_missing", "offline"}:
        raise ValueError("mode must be 'fetch_missing' or 'offline'")
    station_ids = _station_ids(station_ids)
    required = set(_station_ids(required_station_ids))
    optional = set(_station_ids(optional_station_ids))
    forced = required | optional
    if forced - set(station_ids):
        raise ValueError(
            "Forced USGS stations must be included in station_ids: "
            f"{sorted(forced - set(station_ids))}"
        )
    if required & optional:
        raise ValueError("USGS stations cannot be both required and optional")
    if not station_ids:
        return []
    if datelist is None or len(datelist) == 0:
        raise ValueError("datelist must contain at least one timestamp")
    if not pd.DatetimeIndex(datelist).is_monotonic_increasing:
        raise ValueError("datelist must be sorted in increasing order")

    database = USGSObservationStore(database_dir)
    parameter = str(param_id)
    fetch_start = _as_utc(datelist[0])
    fetch_end = _as_utc(datelist[-1])
    validation_start = fetch_start if validation_start is None else _as_utc(
        validation_start
    )
    validation_end = fetch_end if validation_end is None else _as_utc(
        validation_end
    )
    if fetch_end < fetch_start or validation_end < validation_start:
        raise ValueError("USGS request end time precedes start time")

    if mode == "fetch_missing":
        database.import_legacy_caches(
            parameter, station_ids, fetch_start, fetch_end
        )
    policy = dict(download_policy or USGS_DOWNLOAD_POLICY)
    chunk_days = (
        policy["primary_chunk_days"]
        if days_per_chunk is None
        else days_per_chunk
    )
    if type(chunk_days) is not int or chunk_days <= 0:
        raise ValueError("days_per_chunk must be a positive integer")
    if type(stations_chunk_size) is not int or stations_chunk_size <= 0:
        raise ValueError("stations_chunk_size must be a positive integer")
    if retry_failed_after_hours < 0:
        raise ValueError("retry_failed_after_hours must not be negative")
    if (
        max_missing_gap_days is not None
        and max_missing_gap_days <= 0
    ):
        raise ValueError("max_missing_gap_days must be positive")
    retry_delay = timedelta(hours=max(0.0, retry_failed_after_hours))
    now = pd.Timestamp.now(tz="UTC")
    missing_by_window: dict[
        tuple[pd.Timestamp, pd.Timestamp], list[str]
    ] = defaultdict(list)
    for station_id in station_ids:
        missing = database.missing_intervals(
            parameter,
            station_id,
            fetch_start,
            fetch_end,
            retry_failed_after=retry_delay,
            retry_pending_after=PENDING_REQUEST_TIMEOUT,
            now=now,
        )
        for window in missing:
            missing_by_window[window].append(station_id)

    def persist_window(
        requested_ids,
        requested_parameter,
        window_start,
        window_end,
        observations,
        error,
    ):
        try:
            database.record_window(
                requested_ids,
                requested_parameter,
                window_start,
                window_end,
                observations,
                error,
            )
        except Exception as exc:
            raise ObservationPersistenceError(
                "Unable to persist USGS observations or coverage"
            ) from exc

    def mark_window_pending(
        requested_ids,
        requested_parameter,
        window_start,
        window_end,
    ):
        try:
            database.mark_pending(
                requested_ids,
                requested_parameter,
                window_start,
                window_end,
            )
        except Exception as exc:
            raise ObservationPersistenceError(
                "Unable to persist pending USGS request state"
            ) from exc

    if mode == "fetch_missing":
        for (window_start, window_end), window_station_ids in sorted(
            missing_by_window.items(), key=lambda item: item[0]
        ):
            print(
                f"[USGS STORE] fetching {len(window_station_ids)} station(s), "
                f"parameter {parameter}, {window_start:%Y-%m-%d} to "
                f"{window_end:%Y-%m-%d}",
                flush=True,
            )
            with tempfile.TemporaryDirectory(
                prefix="stofs-usgs-fill-"
            ) as temp_dir:
                observation_adapter.download_stations(
                    param_id=parameter,
                    station_ids=window_station_ids,
                    cache_fname=str(Path(temp_dir) / "fill.pq"),
                    datelist=pd.date_range(
                        window_start.tz_localize(None),
                        window_end.tz_localize(None),
                        freq="D",
                    ),
                    stations_chunk_size=stations_chunk_size,
                    days_per_chunk=chunk_days,
                    download_policy=policy,
                    max_missing_gap_days=None,
                    validation_start=window_start,
                    validation_end=window_end,
                    allow_adaptive_retry=True,
                    validate_coverage=False,
                    observation_window_callback=persist_window,
                    observation_window_start_callback=mark_window_pending,
                )

    full_coverage_ids = {
        station_id
        for station_id in station_ids
        if database.has_coverage(
            parameter, station_id, fetch_start, fetch_end
        )
    }
    for station_id in station_ids:
        if station_id in full_coverage_ids:
            continue
        _warning_for_incomplete_station(
            station_id,
            parameter,
            warning_log_path,
            warning_input_files,
            warning_forcing_name,
            required=station_id in required,
            missing=database.missing_intervals(
                parameter, station_id, fetch_start, fetch_end
            ),
        )

    validation_ids = sorted(full_coverage_ids)
    if not validation_ids:
        return []

    records = _records_from_store(
        database.read_observations(
            parameter, validation_ids, fetch_start, fetch_end
        ),
        parameter,
    )
    if mode == "fetch_missing":
        _retry_quality_gaps(
            records=records,
            forced_station_ids=forced & full_coverage_ids,
            parameter=parameter,
            start=validation_start,
            end=validation_end,
            max_missing_gap_days=max_missing_gap_days,
            download_policy=policy,
            observation_window_callback=persist_window,
            observation_window_start_callback=mark_window_pending,
        )
        records = _records_from_store(
            database.read_observations(
                parameter, validation_ids, fetch_start, fetch_end
            ),
            parameter,
        )

    return _validate_cached_records(
        records=records,
        parameter=parameter,
        required_station_ids=required & full_coverage_ids,
        optional_station_ids=optional & full_coverage_ids,
        validation_start=validation_start,
        validation_end=validation_end,
        max_missing_gap_days=max_missing_gap_days,
        warning_log_path=warning_log_path,
        warning_input_files=warning_input_files,
        warning_forcing_name=warning_forcing_name,
    )


def _retry_quality_gaps(
    records: list[GenericObsData],
    forced_station_ids: Iterable[str],
    parameter: str,
    start: Any,
    end: Any,
    max_missing_gap_days: float | None,
    download_policy: Mapping[str, Any],
    observation_window_callback,
    observation_window_start_callback,
) -> None:
    """Retry only long gaps (or empty forced series), preserving cached data."""
    from stofs3d_setup.ops.USGS.iv.download_usgs import find_data_gaps
    from stofs3d_setup.ops.USGS.ogc.series import download_usgs_series

    gap_days = (
        USGS_MAX_MISSING_GAP_DAYS
        if max_missing_gap_days is None
        else max_missing_gap_days
    )
    records_by_id = {
        normalize_station_id(record.station_info.get("id")): record
        for record in records
    }
    retry_windows: dict[
        str, list[tuple[pd.Timestamp, pd.Timestamp]]
    ] = {}
    for station_id, record in records_by_id.items():
        record_times = pd.DatetimeIndex(
            pd.to_datetime(record.df["date"], utc=True)
        )
        gaps = find_data_gaps(
            record_times,
            _as_utc(start),
            _as_utc(end),
            pd.Timedelta(days=gap_days),
        )
        if gaps:
            retry_windows[station_id] = gaps

    for station_id in set(_station_ids(forced_station_ids)) - set(records_by_id):
        retry_windows[station_id] = [(_as_utc(start), _as_utc(end))]

    if not retry_windows:
        return

    with tempfile.TemporaryDirectory(prefix="stofs-usgs-gap-retry-") as temp_dir:
        for station_id, windows in sorted(retry_windows.items()):
            for window_start, window_end in windows:
                try:
                    download_usgs_series(
                        station_id=station_id,
                        parameter_id=parameter,
                        start_time=window_start,
                        end_time=window_end,
                        usgs_cache_folder=temp_dir,
                        policy=download_policy,
                        observation_window_callback=observation_window_callback,
                        observation_window_start_callback=(
                            observation_window_start_callback
                        ),
                    )
                except Exception as exc:
                    from stofs3d_setup.ops.USGS.ogc.observation_adapter import (
                        ObservationPersistenceError,
                    )

                    if isinstance(exc, ObservationPersistenceError):
                        raise
                    print(
                        f"[USGS STORE] targeted retry failed for "
                        f"{station_id}/{parameter}: {type(exc).__name__}",
                        flush=True,
                    )


def _validate_cached_records(
    records: list[GenericObsData],
    parameter: str,
    required_station_ids: Iterable[str],
    optional_station_ids: Iterable[str],
    validation_start: Any,
    validation_end: Any,
    max_missing_gap_days: float | None,
    warning_log_path: str | Path | None,
    warning_input_files: Iterable[str],
    warning_forcing_name: str | None,
) -> list[GenericObsData]:
    """Apply current model-quality checks without making any network calls."""
    from stofs3d_setup.ops.USGS.ogc.observation_adapter import (
        _filter_temperature_records,
        _record_failures,
    )

    required = set(_station_ids(required_station_ids))
    optional = set(_station_ids(optional_station_ids))
    checked = required | optional | {
        normalize_station_id(record.station_info.get("id"))
        for record in records
    }
    gap_days = (
        USGS_MAX_MISSING_GAP_DAYS
        if max_missing_gap_days is None
        else max_missing_gap_days
    )
    records = _filter_temperature_records(records, parameter)
    failures = _record_failures(
        records,
        checked,
        _as_utc(validation_start),
        _as_utc(validation_end),
        pd.Timedelta(days=gap_days),
    )
    required_failures = {
        station_id: reason
        for station_id, reason in failures.items()
        if station_id in required
    }
    if required_failures:
        raise RequiredStationDownloadError(
            f"Required USGS station(s) failed stored-data validation for "
            f"parameter {parameter}: {required_failures}"
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
                    f"{warning_forcing_name or parameter} replacement "
                    f"station {station_id} failed stored-data validation: "
                    f"{reason}"
                ),
                input_files=warning_input_files,
            )
    unusable = set(failures) - required
    return [
        record for record in records
        if normalize_station_id(record.station_info.get("id"))
        not in unusable
    ]


def _records_from_store(
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
        def first_value(column: str, default=None):
            if column not in station_data:
                return default
            values = station_data[column].dropna()
            return values.iloc[0] if not values.empty else default

        info = {
            "id": normalize_station_id(station_id),
            "name": first_value("station_name"),
            "lon": first_value("lon", np.nan),
            "lat": first_value("lat", np.nan),
            "var_name": str(parameter),
            "var_code": str(parameter),
            "unit": first_value("unit"),
        }
        frame = pd.DataFrame({
            "date": pd.to_datetime(station_data["time"], utc=True),
            "value": pd.to_numeric(station_data["value"], errors="coerce"),
        }).dropna(subset=["date", "value"])
        records.append(GenericObsData(station_info=info, df=frame))
    return records
