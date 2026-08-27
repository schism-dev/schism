#!/usr/bin/env python3
"""Compare the three v32e hydrofabric dredging diagnostics with legacy output."""

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


DEFAULT_RUN_ROOT = Path(
    "/sciclone/schism10/feiye/CIROH/Channel_Geometry/"
    "Connectivity_test/dredge_full_runs_v32e"
)
DEFAULT_CACHE = Path(
    "/sciclone/schism10/feiye/CIROH/Channel_Geometry/"
    "Connectivity_test/dredge_full_v32e/input_cache/hgrid_nodes"
)
DEFAULT_LEGACY_HGRID = Path(
    "/sciclone/schism10/Hgrid_projects/STOFS3D-v7.4/v32e/"
    "Bathy_edit_channel_variations/Channel_variations_verified_pre_channel/"
    "lower_bank_2m/Ensure_channel_connectivity/hgrid_channel_lower_bank_2m.gr3"
)


def _read_product(run_dir: Path, name: str) -> pd.DataFrame:
    return pd.read_parquet(run_dir / name)


def _point_layer(frame: pd.DataFrame) -> gpd.GeoDataFrame:
    return gpd.GeoDataFrame(
        frame,
        geometry=gpd.points_from_xy(frame.x, frame.y),
        crs=ANALYSIS_GEOGRAPHIC_CRS,
    )


def compare(
    run_root: Path,
    cache_dir: Path,
    legacy_hgrid: Path,
) -> tuple[dict[str, object], dict[str, pd.DataFrame]]:
    run_dirs = {
        "test1": run_root / "01_ml_lowerbank",
        "test2": run_root / "02_constant2_lowerbank",
        "test3": run_root / "03_constant2_lowerbank_forward",
    }
    requested = {
        key: _read_product(path, "dredge_requested_mesh_nodes.parquet")
        for key, path in run_dirs.items()
    }
    changed = {
        key: _read_product(path, "dredged_mesh_nodes.parquet")
        for key, path in run_dirs.items()
    }
    for key in run_dirs:
        if not requested[key].mesh_node_id.is_unique:
            raise RuntimeError(f"{key} requested mesh-node IDs are not unique")
        if not changed[key].mesh_node_id.is_unique:
            raise RuntimeError(f"{key} changed mesh-node IDs are not unique")

    test1 = requested["test1"].set_index("mesh_node_id")
    test2 = requested["test2"].set_index("mesh_node_id")
    if not test1.index.equals(test2.index):
        raise RuntimeError("Test 1 and Test 2 requested-node indices differ")
    comparison = pd.DataFrame(
        {
            "mesh_node_id": test1.index,
            "x": test1.x,
            "y": test1.y,
            "original_dp": test1.original_dp,
            "test1_requested_dp": test1.requested_dp,
            "test1_final_dp": test1.final_dp,
            "test1_requested_delta": test1.requested_dredging_delta_m,
            "test1_applied_delta": test1.dredging_delta_m,
            "test1_status": test1.depth_screen_status,
            "test2_requested_dp": test2.requested_dp,
            "test2_final_dp": test2.final_dp,
            "test2_requested_delta": test2.requested_dredging_delta_m,
            "test2_applied_delta": test2.dredging_delta_m,
            "test2_status": test2.depth_screen_status,
        }
    ).reset_index(drop=True)
    comparison["test1_minus_test2_final_m"] = (
        comparison.test1_final_dp - comparison.test2_final_dp
    )
    comparison["abs_final_difference_m"] = (
        comparison.test1_minus_test2_final_m.abs()
    )

    changed_ids = {
        key: set(frame.mesh_node_id.astype(int)) for key, frame in changed.items()
    }
    test1_only = comparison.loc[
        comparison.mesh_node_id.isin(changed_ids["test1"] - changed_ids["test2"])
    ].copy()
    test2_only = comparison.loc[
        comparison.mesh_node_id.isin(changed_ids["test2"] - changed_ids["test1"])
    ].copy()
    test2_intersection_added = changed["test2"].loc[
        ~changed["test2"].mesh_node_id.isin(changed_ids["test3"])
    ].copy()

    base_node_id = np.load(cache_dir / "node_id.npy", mmap_mode="r")
    base_x = np.load(cache_dir / "x.npy", mmap_mode="r")
    base_y = np.load(cache_dir / "y.npy", mmap_mode="r")
    base_dp = np.load(cache_dir / "dp.npy", mmap_mode="r")
    legacy = read_hgrid_nodes(legacy_hgrid)
    if not np.array_equal(legacy.node_id, base_node_id):
        raise RuntimeError("Legacy and cached base hgrid node IDs differ")
    legacy_changed_mask = legacy.dp > base_dp + 1.0e-12
    legacy_changed_ids = set(legacy.node_id[legacy_changed_mask].astype(int))
    id_to_position = pd.Series(
        np.arange(len(base_node_id), dtype=np.int64), index=np.asarray(base_node_id)
    )

    def base_frame(node_ids: set[int]) -> pd.DataFrame:
        ids = np.asarray(sorted(node_ids), dtype=np.int64)
        positions = id_to_position.loc[ids].to_numpy(dtype=np.int64)
        return pd.DataFrame(
            {
                "mesh_node_id": ids,
                "x": np.asarray(base_x[positions]),
                "y": np.asarray(base_y[positions]),
                "original_dp": np.asarray(base_dp[positions]),
                "legacy_final_dp": legacy.dp[positions],
                "legacy_delta_m": legacy.dp[positions]
                - np.asarray(base_dp[positions]),
            }
        )

    legacy_only = base_frame(legacy_changed_ids - changed_ids["test3"])
    test3_only = changed["test3"].loc[
        ~changed["test3"].mesh_node_id.isin(legacy_changed_ids)
    ].copy()
    shared_ids = np.asarray(
        sorted(legacy_changed_ids & changed_ids["test3"]), dtype=np.int64
    )
    shared_positions = id_to_position.loc[shared_ids].to_numpy(dtype=np.int64)
    test3_by_id = changed["test3"].set_index("mesh_node_id")
    legacy_test3_difference = pd.DataFrame(
        {
            "mesh_node_id": shared_ids,
            "x": np.asarray(base_x[shared_positions]),
            "y": np.asarray(base_y[shared_positions]),
            "original_dp": np.asarray(base_dp[shared_positions]),
            "legacy_final_dp": legacy.dp[shared_positions],
            "test3_final_dp": test3_by_id.loc[shared_ids].final_dp.to_numpy(),
        }
    )
    legacy_test3_difference["test3_minus_legacy_m"] = (
        legacy_test3_difference.test3_final_dp
        - legacy_test3_difference.legacy_final_dp
    )
    legacy_test3_difference["abs_final_difference_m"] = (
        legacy_test3_difference.test3_minus_legacy_m.abs()
    )
    legacy_test3_difference = legacy_test3_difference.loc[
        legacy_test3_difference.abs_final_difference_m > 1.0e-9
    ].copy()

    layers = {
        "t1_t2_diff_gt0p5m": comparison.loc[
            comparison.abs_final_difference_m > 0.5
        ].copy(),
        "t1_only_changed": test1_only,
        "t2_only_changed": test2_only,
        "t2_intersection_added": test2_intersection_added,
        "t1_capped_at6m": requested["test1"].loc[
            requested["test1"].dredging_delta_capped
        ].copy(),
        "t2_capped_at6m": requested["test2"].loc[
            requested["test2"].dredging_delta_capped
        ].copy(),
        "t3_capped_at6m": requested["test3"].loc[
            requested["test3"].dredging_delta_capped
        ].copy(),
        "t1_raw_delta_gt10m": requested["test1"].loc[
            requested["test1"].requested_dredging_delta_m > 10.0
        ].copy(),
        "t2_raw_delta_gt10m": requested["test2"].loc[
            requested["test2"].requested_dredging_delta_m > 10.0
        ].copy(),
        "t3_raw_delta_gt10m": requested["test3"].loc[
            requested["test3"].requested_dredging_delta_m > 10.0
        ].copy(),
        "legacy_only_changed": legacy_only,
        "t3_only_changed": test3_only,
        "legacy_t3_depth_diff": legacy_test3_difference,
    }

    test2_changed_mask = test2.would_deepen.astype(bool)
    test1_minus_test2 = test1.final_dp - test2.final_dp
    additional_depth_counts = {
        str(threshold): int(
            (test2_changed_mask & (test1_minus_test2 > threshold)).sum()
        )
        for threshold in (1.0e-12, 0.5, 1.0, 2.0, 5.0)
    }
    summary: dict[str, object] = {
        "runs": {key: str(value) for key, value in run_dirs.items()},
        "legacy_comparison_hgrid": str(legacy_hgrid),
        "counts": {
            **{
                f"{key}_requested": len(requested[key]) for key in run_dirs
            },
            **{f"{key}_changed": len(changed[key]) for key in run_dirs},
            **{
                f"{key}_capped_at6m": int(
                    requested[key].dredging_delta_capped.sum()
                )
                for key in run_dirs
            },
            "test1_test2_same_requested_node_set": True,
            "test1_test2_final_difference_gt0p5m": len(
                layers["t1_t2_diff_gt0p5m"]
            ),
            "test1_only_changed": len(test1_only),
            "test2_only_changed": len(test2_only),
            "test2_intersection_added_changed": len(test2_intersection_added),
            "test2_changed_receiving_additional_test1_depth": (
                additional_depth_counts
            ),
            "legacy_changed": int(legacy_changed_mask.sum()),
            "legacy_test3_changed_overlap": len(
                legacy_changed_ids & changed_ids["test3"]
            ),
            "legacy_only_changed": len(legacy_only),
            "test3_only_changed": len(test3_only),
            "legacy_test3_shared_depth_difference": len(
                legacy_test3_difference
            ),
        },
        "gpkg_layers": {key: len(value) for key, value in layers.items()},
    }
    return summary, layers


def get_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--run-root", type=Path, default=DEFAULT_RUN_ROOT)
    parser.add_argument("--cache-dir", type=Path, default=DEFAULT_CACHE)
    parser.add_argument("--legacy-hgrid", type=Path, default=DEFAULT_LEGACY_HGRID)
    return parser


def main(argv: Sequence[str] | None = None) -> int:
    args = get_parser().parse_args(argv)
    summary, layers = compare(args.run_root, args.cache_dir, args.legacy_hgrid)
    gpkg_path = args.run_root / "comparison_review.gpkg"
    gpkg_path.unlink(missing_ok=True)
    for name, frame in layers.items():
        if frame.empty:
            continue
        _point_layer(frame).to_file(
            gpkg_path, layer=name, driver="GPKG", mode="w"
        )
    summary["comparison_gpkg"] = str(gpkg_path)
    summary_path = args.run_root / "comparison_summary.json"
    summary_path.write_text(json.dumps(summary, indent=2) + "\n")
    print(json.dumps(summary, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
