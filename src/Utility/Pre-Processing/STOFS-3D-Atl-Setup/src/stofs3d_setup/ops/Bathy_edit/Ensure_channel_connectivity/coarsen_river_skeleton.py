#!/usr/bin/env python3
"""Remove manual line densification while preserving river meanders."""

from __future__ import annotations

import argparse
import json
from pathlib import Path
from typing import Sequence

import geopandas as gpd
import numpy as np
from shapely import count_coordinates, force_2d, hausdorff_distance
from shapely.geometry import box

try:
    from stofs3d_setup.ops.Bathy_edit.Ensure_channel_connectivity.hydrofabric_match import (
        ANALYSIS_GEOGRAPHIC_CRS,
        MATCH_CRS,
    )
except ModuleNotFoundError:  # Permit direct execution from this directory.
    from hydrofabric_match import (  # type: ignore[no-redef]
        ANALYSIS_GEOGRAPHIC_CRS,
        MATCH_CRS,
    )


def coarsen_skeleton(
    skeleton: gpd.GeoDataFrame,
    tolerance_m: float = 1.0,
    footprint_lonlat: tuple[float, float, float, float] | None = None,
) -> tuple[gpd.GeoDataFrame, dict[str, object]]:
    """Simplify densified lines in a metric CRS without smoothing real bends.

    Douglas--Peucker simplification retains every displacement larger than the
    requested tolerance and preserves each LineString's endpoints. The default
    1 m tolerance is deliberately small relative to both the approximately
    50 m manual vertex spacing and the channel-mapping search distances.
    """
    if skeleton.crs is None:
        raise ValueError("River skeleton has no CRS")
    if tolerance_m <= 0:
        raise ValueError("tolerance_m must be positive")
    selected = skeleton.copy()
    input_feature_count = len(selected)
    if footprint_lonlat is not None:
        footprint = gpd.GeoSeries(
            [box(*footprint_lonlat)], crs=ANALYSIS_GEOGRAPHIC_CRS
        ).to_crs(selected.crs).iloc[0]
        selected = selected.loc[selected.geometry.intersects(footprint)].copy()
    selected = selected.loc[~selected.geometry.is_empty & selected.geometry.notna()]
    metric = selected.to_crs(MATCH_CRS).reset_index(drop=True)
    metric.geometry = force_2d(metric.geometry.array)
    original_geometry = metric.geometry.array.copy()
    original_lengths = metric.length.to_numpy(dtype=float)
    original_coordinates = int(count_coordinates(original_geometry))

    metric.geometry = metric.geometry.simplify(
        tolerance_m, preserve_topology=False
    )
    simplified_lengths = metric.length.to_numpy(dtype=float)
    simplified_coordinates = int(count_coordinates(metric.geometry.array))
    deviations = np.asarray(
        hausdorff_distance(original_geometry, metric.geometry.array), dtype=float
    )
    start_preserved = []
    end_preserved = []
    for original, simplified in zip(original_geometry, metric.geometry.array):
        original_coordinates_xy = np.asarray(original.coords, dtype=float)[:, :2]
        simplified_coordinates_xy = np.asarray(simplified.coords, dtype=float)[:, :2]
        start_preserved.append(
            np.array_equal(original_coordinates_xy[0], simplified_coordinates_xy[0])
        )
        end_preserved.append(
            np.array_equal(original_coordinates_xy[-1], simplified_coordinates_xy[-1])
        )
    metric.insert(0, "river_idx", np.arange(len(metric), dtype=np.int64))
    metric["n_stations"] = np.asarray(
        [len(geometry.coords) for geometry in metric.geometry], dtype=np.int64
    )

    length_loss = original_lengths - simplified_lengths
    summary: dict[str, object] = {
        "tolerance_m": tolerance_m,
        "input_features": input_feature_count,
        "selected_features": len(metric),
        "original_coordinates": original_coordinates,
        "simplified_coordinates": simplified_coordinates,
        "coordinate_reduction_fraction": (
            1.0 - simplified_coordinates / original_coordinates
            if original_coordinates
            else 0.0
        ),
        "original_total_length_m": float(original_lengths.sum()),
        "simplified_total_length_m": float(simplified_lengths.sum()),
        "total_length_loss_m": float(length_loss.sum()),
        "relative_total_length_loss": (
            float(length_loss.sum() / original_lengths.sum())
            if original_lengths.sum()
            else 0.0
        ),
        "max_feature_length_loss_m": float(length_loss.max(initial=0.0)),
        "max_hausdorff_distance_m": float(deviations.max(initial=0.0)),
        "p99_hausdorff_distance_m": (
            float(np.percentile(deviations, 99)) if len(deviations) else 0.0
        ),
        "all_startpoints_preserved": bool(np.all(start_preserved)),
        "all_endpoints_preserved": bool(np.all(end_preserved)),
        "output_crs": MATCH_CRS,
    }
    return metric, summary


def get_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--skeleton", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--tolerance-m", type=float, default=1.0)
    parser.add_argument(
        "--footprint",
        nargs=4,
        type=float,
        default=None,
        metavar=("MIN_LON", "MIN_LAT", "MAX_LON", "MAX_LAT"),
        help="Select features intersecting a longitude/latitude mesh footprint.",
    )
    parser.add_argument("--summary", type=Path, default=None)
    return parser


def main(argv: Sequence[str] | None = None) -> int:
    args = get_parser().parse_args(argv)
    skeleton = gpd.read_file(args.skeleton)
    footprint = tuple(args.footprint) if args.footprint is not None else None
    simplified, summary = coarsen_skeleton(
        skeleton, tolerance_m=args.tolerance_m, footprint_lonlat=footprint
    )
    args.output.parent.mkdir(parents=True, exist_ok=True)
    simplified.to_parquet(args.output, index=False, compression="zstd")
    summary_path = args.summary or args.output.with_suffix(".summary.json")
    summary_path.write_text(json.dumps(summary, indent=2, sort_keys=True) + "\n")
    print(json.dumps(summary, indent=2, sort_keys=True), flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
