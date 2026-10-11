#!/usr/bin/env python3
"""Live-check USGS station inventory and flow/temperature metadata screening."""

from __future__ import annotations

import sys
from pathlib import Path


REPOSITORY_ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPOSITORY_ROOT / "src"))

from stofs3d_setup.ops.USGS.ogc.download_usgs import (
    discover_monitoring_locations,
    discover_station_ids_with_parameter,
)


STATE = "NY"
START = "2017-01-01"
END = "2017-12-31"
STATION_CHECKS = (
    ("01358000", "00060", "Hudson flow"),
    ("01359139", "00010", "Hudson temperature"),
)


def main() -> int:
    print(f"Loading {STATE} USGS monitoring-location inventory...")
    locations = discover_monitoring_locations(states=[STATE])
    if locations.empty:
        raise RuntimeError(
            f"USGS returned no parsed locations for {STATE}; "
            "check the inventory state_code filter and response fields"
        )
    location_ids = set(locations["station_id"].astype(str))
    print(f"Loaded {len(locations)} USGS locations.")

    for station_id, parameter, label in STATION_CHECKS:
        if station_id not in location_ids:
            raise RuntimeError(
                f"{label}: station {station_id} is missing from the {STATE} inventory"
            )

        matches = discover_station_ids_with_parameter(
            [station_id],
            parameter=parameter,
            start=START,
            end=END,
        )
        if station_id not in matches:
            raise RuntimeError(
                f"{label}: no parameter {parameter} series overlaps {START} to {END}"
            )

        print(
            f"PASS {label}: station {station_id} has parameter {parameter} "
            f"metadata overlapping {START} to {END}"
        )

    print("USGS inventory and station-screening checks passed.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
