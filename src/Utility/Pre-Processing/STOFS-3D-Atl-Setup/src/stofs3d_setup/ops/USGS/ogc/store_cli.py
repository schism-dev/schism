"""Command-line access to the resumable USGS observation store."""

from __future__ import annotations

import argparse
import json
from pathlib import Path
import time

import pandas as pd

from stofs3d_setup.ops.USGS.iv.download_usgs import RequiredStationDownloadError
from stofs3d_setup.ops.USGS.ogc.store import (
    DEFAULT_USGS_CACHE_FOLDER,
    USGSObservationStore,
    get_usgs_observations,
)


def _parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description="Fetch or inspect the persistent USGS observation store."
    )
    subparsers = parser.add_subparsers(dest="command", required=True)

    fetch = subparsers.add_parser(
        "fetch", help="fetch missing station/parameter/date coverage"
    )
    fetch.add_argument(
        "--database-dir",
        type=Path,
        default=DEFAULT_USGS_CACHE_FOLDER,
    )
    fetch.add_argument("--parameter", required=True, help="USGS parameter code")
    fetch.add_argument(
        "--station", action="append", default=[], help="station ID (repeatable)"
    )
    fetch.add_argument(
        "--stations-file",
        type=Path,
        help="optional text file with one station ID per line",
    )
    fetch.add_argument(
        "--required-station",
        action="append",
        default=[],
        help="station that must pass the existing quality checks",
    )
    fetch.add_argument("--start", required=True, help="inclusive UTC date/time")
    fetch.add_argument("--end", required=True, help="inclusive UTC date/time")
    fetch.add_argument("--days-per-chunk", type=int)
    fetch.add_argument(
        "--retry-failed-after-hours", type=float, default=24.0
    )
    fetch.add_argument(
        "--offline",
        action="store_true",
        help="read only; do not issue OGC requests",
    )
    fetch.add_argument(
        "--watch-seconds",
        type=float,
        default=0.0,
        help="repeat until interrupted, retrying failed coverage after cooldown",
    )

    status = subparsers.add_parser("status", help="summarize store contents")
    status.add_argument(
        "--database-dir",
        type=Path,
        default=DEFAULT_USGS_CACHE_FOLDER,
    )
    return parser


def _station_list(args: argparse.Namespace) -> list[str]:
    station_ids = list(args.station)
    if args.stations_file is not None:
        station_ids.extend(
            line.strip()
            for line in args.stations_file.read_text(
                encoding="utf-8"
            ).splitlines()
            if line.strip() and not line.lstrip().startswith("#")
        )
    return sorted(set(station_ids))


def _run_fetch(args: argparse.Namespace) -> int:
    station_ids = _station_list(args)
    if not station_ids:
        raise ValueError("provide --station or --stations-file")
    start = pd.Timestamp(args.start)
    end = pd.Timestamp(args.end)
    if start.tzinfo is None:
        start = start.tz_localize("UTC")
    if end.tzinfo is None:
        end = end.tz_localize("UTC")
    if end < start:
        raise ValueError("--end must not precede --start")

    mode = "offline" if args.offline else "fetch_missing"
    while True:
        try:
            records = get_usgs_observations(
                param_id=args.parameter,
                station_ids=station_ids,
                database_dir=args.database_dir,
                datelist=pd.date_range(start, end, freq="D"),
                required_station_ids=args.required_station,
                mode=mode,
                retry_failed_after_hours=args.retry_failed_after_hours,
                days_per_chunk=args.days_per_chunk,
            )
            observation_count = sum(len(record.df) for record in records)
            print(
                f"USGS store query returned {len(records)} station series and "
                f"{observation_count} observation(s).",
                flush=True,
            )
        except RequiredStationDownloadError as exc:
            print(f"USGS required station is not ready: {exc}", flush=True)
            if args.watch_seconds <= 0:
                return 2
        if args.watch_seconds <= 0:
            return 0
        time.sleep(max(1.0, args.watch_seconds))


def main(argv: list[str] | None = None) -> int:
    args = _parser().parse_args(argv)
    if args.command == "status":
        print(json.dumps(USGSObservationStore(args.database_dir).status(), indent=2))
        return 0
    return _run_fetch(args)


if __name__ == "__main__":
    raise SystemExit(main())
