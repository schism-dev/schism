#!/usr/bin/env python3
"""Compare one hydrofabric dredging run with a legacy connectivity hgrid."""

from __future__ import annotations

import argparse
import json
from pathlib import Path
from typing import Sequence

import geopandas as gpd
import numpy as np
import pandas as pd

try:
    from stofs3d_setup.ops.Bathy_edit.Ensure_channel_connectivity.hydrofabric_dredge import (
        ANALYSIS_GEOGRAPHIC_CRS,
        read_hgrid_nodes,
    )
except ModuleNotFoundError:  # Permit direct execution from this directory.
    from hydrofabric_dredge import (  # type: ignore[no-redef]
        ANALYSIS_GEOGRAPHIC_CRS,
        read_hgrid_nodes,
    )


def _point_layer(frame: pd.DataFrame) -> gpd.GeoDataFrame:
    return gpd.GeoDataFrame(
        frame,
        geometry=gpd.points_from_xy(frame.x, frame.y),
        crs=ANALYSIS_GEOGRAPHIC_CRS,
    )


def compare(
    run_dir: Path,
    cache_dir: Path,
    legacy_hgrid: Path,
    difference_threshold_m: float = 0.5,
) -> tuple[dict[str, object], dict[str, pd.DataFrame]]:
    requested = pd.read_parquet(run_dir / "dredge_requested_mesh_nodes.parquet")
    changed = pd.read_parquet(run_dir / "dredged_mesh_nodes.parquet")
    if not requested.mesh_node_id.is_unique or not changed.mesh_node_id.is_unique:
        raise RuntimeError("New-run mesh-node IDs must be unique")

    node_id = np.load(cache_dir / "node_id.npy", mmap_mode="r")
    x = np.load(cache_dir / "x.npy", mmap_mode="r")
    y = np.load(cache_dir / "y.npy", mmap_mode="r")
    base_dp = np.load(cache_dir / "dp.npy", mmap_mode="r")
    legacy = read_hgrid_nodes(legacy_hgrid)
    if not np.array_equal(legacy.node_id, node_id):
        raise RuntimeError("Legacy and cached base hgrid node IDs differ")
    if not np.allclose(legacy.x, x, rtol=0.0, atol=1.0e-12) or not np.allclose(
        legacy.y, y, rtol=0.0, atol=1.0e-12
    ):
        raise RuntimeError("Legacy and cached base hgrid coordinates differ")

    legacy_changed_mask = legacy.dp > base_dp + 1.0e-12
    legacy_ids = set(legacy.node_id[legacy_changed_mask].astype(int))
    new_ids = set(changed.mesh_node_id.astype(int))
    position_by_id = pd.Series(
        np.arange(len(node_id), dtype=np.int64), index=np.asarray(node_id)
    )

    def legacy_frame(ids: set[int]) -> pd.DataFrame:
        ordered = np.asarray(sorted(ids), dtype=np.int64)
        positions = position_by_id.loc[ordered].to_numpy(dtype=np.int64)
        return pd.DataFrame(
            {
                "mesh_node_id": ordered,
                "x": np.asarray(x[positions]),
                "y": np.asarray(y[positions]),
                "original_dp": np.asarray(base_dp[positions]),
                "legacy_final_dp": legacy.dp[positions],
                "legacy_delta_m": legacy.dp[positions]
                - np.asarray(base_dp[positions]),
            }
        )

    new_only = changed.loc[~changed.mesh_node_id.isin(legacy_ids)].copy()
    legacy_only = legacy_frame(legacy_ids - new_ids)
    shared_ids = np.asarray(sorted(new_ids & legacy_ids), dtype=np.int64)
    shared_positions = position_by_id.loc[shared_ids].to_numpy(dtype=np.int64)
    new_by_id = changed.set_index("mesh_node_id")
    shared = pd.DataFrame(
        {
            "mesh_node_id": shared_ids,
            "x": np.asarray(x[shared_positions]),
            "y": np.asarray(y[shared_positions]),
            "original_dp": np.asarray(base_dp[shared_positions]),
            "legacy_final_dp": legacy.dp[shared_positions],
            "new_final_dp": new_by_id.loc[shared_ids].final_dp.to_numpy(),
        }
    )
    shared["new_minus_legacy_m"] = shared.new_final_dp - shared.legacy_final_dp
    shared["abs_final_difference_m"] = shared.new_minus_legacy_m.abs()
    shared_large_difference = shared.loc[
        shared.abs_final_difference_m > difference_threshold_m
    ].copy()

    layers = {
        "new_only_changed": new_only,
        "legacy_only_changed": legacy_only,
        "shared_diff_gt0p5m": shared_large_difference,
        "new_capped_at6m": requested.loc[requested.dredging_delta_capped].copy(),
        "new_raw_delta_gt10m": requested.loc[
            requested.requested_dredging_delta_m > 10.0
        ].copy(),
        "new_intersection_only": changed.loc[
            changed.intersection_target & ~changed.forward_requested
        ].copy(),
    }
    abs_difference = shared.abs_final_difference_m.to_numpy(dtype=float)
    summary: dict[str, object] = {
        "run_dir": str(run_dir),
        "legacy_hgrid": str(legacy_hgrid),
        "base_hgrid_cache": str(cache_dir),
        "difference_threshold_m": difference_threshold_m,
        "counts": {
            "new_requested": len(requested),
            "new_changed": len(changed),
            "legacy_changed": int(legacy_changed_mask.sum()),
            "shared_changed": len(shared_ids),
            "new_only_changed": len(new_only),
            "legacy_only_changed": len(legacy_only),
            "shared_difference_gt_threshold": len(shared_large_difference),
            "new_capped_at6m": len(layers["new_capped_at6m"]),
            "new_raw_delta_gt10m": len(layers["new_raw_delta_gt10m"]),
            "new_intersection_only": len(layers["new_intersection_only"]),
        },
        "shared_abs_depth_difference_m": {
            "median": float(np.median(abs_difference)) if len(abs_difference) else None,
            "p95": float(np.percentile(abs_difference, 95))
            if len(abs_difference)
            else None,
            "p99": float(np.percentile(abs_difference, 99))
            if len(abs_difference)
            else None,
            "max": float(np.max(abs_difference)) if len(abs_difference) else None,
        },
        "gpkg_layers": {name: len(frame) for name, frame in layers.items()},
    }
    return summary, layers


def get_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--run-dir", type=Path, required=True)
    parser.add_argument("--cache-dir", type=Path, required=True)
    parser.add_argument("--legacy-hgrid", type=Path, required=True)
    parser.add_argument("--difference-threshold-m", type=float, default=0.5)
    return parser


def main(argv: Sequence[str] | None = None) -> int:
    args = get_parser().parse_args(argv)
    summary, layers = compare(
        args.run_dir,
        args.cache_dir,
        args.legacy_hgrid,
        args.difference_threshold_m,
    )
    gpkg_path = args.run_dir / "legacy_comparison.gpkg"
    gpkg_path.unlink(missing_ok=True)
    for name, frame in layers.items():
        if not frame.empty:
            _point_layer(frame).to_file(
                gpkg_path, layer=name, driver="GPKG", mode="w"
            )
    summary["comparison_gpkg"] = str(gpkg_path)
    summary_path = args.run_dir / "legacy_comparison_summary.json"
    summary_path.write_text(json.dumps(summary, indent=2) + "\n")
    print(json.dumps(summary, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
