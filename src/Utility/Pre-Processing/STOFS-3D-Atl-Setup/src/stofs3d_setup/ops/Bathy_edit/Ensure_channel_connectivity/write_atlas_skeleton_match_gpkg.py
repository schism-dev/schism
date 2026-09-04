#!/usr/bin/env python3
"""Write a compact GIS review package for Atlas skeleton-vertex matching."""

from __future__ import annotations

import argparse
import json
from pathlib import Path
from typing import Sequence

import geopandas as gpd
import numpy as np
import pandas as pd
from scipy.spatial import cKDTree
from shapely import STRtree, hausdorff_distance, intersects_xy, points, prepare
from shapely.geometry import LineString

try:
    from stofs3d_setup.ops.Bathy_edit.Ensure_channel_connectivity.hydrofabric_match import (
        MATCH_CRS,
    )
except ModuleNotFoundError:  # Permit direct execution from this directory.
    from hydrofabric_match import MATCH_CRS  # type: ignore[no-redef]


def _inner_arc_vertex_arrays(
    arcs: gpd.GeoDataFrame,
) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    """Explode inner arcs without constructing millions of Python records."""
    selected = arcs.loc[
        (arcs.local_arc_idx > 0) & (arcs.local_arc_idx < arcs.n_arcs - 1)
    ]
    coordinates: list[np.ndarray] = []
    river_ids: list[np.ndarray] = []
    arc_ids: list[np.ndarray] = []
    station_ids: list[np.ndarray] = []
    for arc in selected.itertuples(index=False):
        xy = np.asarray(arc.geometry.coords, dtype=float)[:, :2]
        coordinates.append(xy)
        river_ids.append(np.full(len(xy), int(arc.river_idx), dtype=np.int64))
        arc_ids.append(np.full(len(xy), int(arc.local_arc_idx), dtype=np.int16))
        station_ids.append(np.arange(len(xy), dtype=np.int32))
    if not coordinates:
        empty_xy = np.empty((0, 2), dtype=float)
        return (
            empty_xy,
            np.empty(0, dtype=np.int64),
            np.empty(0, dtype=np.int16),
            np.empty(0, dtype=np.int32),
        )
    return (
        np.concatenate(coordinates),
        np.concatenate(river_ids),
        np.concatenate(arc_ids),
        np.concatenate(station_ids),
    )


def build_strict_inner_arc_catalog(
    arcs: gpd.GeoDataFrame,
    mesh_xy: np.ndarray,
    river_ids: set[int] | None = None,
    selection_geometry=None,
    selection_buffer_m: float = 500.0,
    max_arc_node_distance_m: float = 50.0,
) -> tuple[pd.DataFrame, dict[str, int]]:
    """Map only RiverMapper inner-arc vertices to their nearest mesh nodes."""
    metric_arcs = arcs.to_crs(MATCH_CRS)
    if river_ids is not None:
        metric_arcs = metric_arcs.loc[metric_arcs.river_idx.isin(river_ids)].copy()
    if selection_geometry is not None:
        xmin, ymin, xmax, ymax = selection_geometry.bounds
        bounds = metric_arcs.geometry.bounds
        metric_arcs = metric_arcs.loc[
            (bounds.maxx >= xmin - selection_buffer_m)
            & (bounds.minx <= xmax + selection_buffer_m)
            & (bounds.maxy >= ymin - selection_buffer_m)
            & (bounds.miny <= ymax + selection_buffer_m)
        ].copy()
    inner_xy, inner_river, inner_arc, inner_station = _inner_arc_vertex_arrays(
        metric_arcs
    )
    if len(inner_xy) == 0:
        raise ValueError("No inner RiverMapper arc vertices are available")
    mesh_tree = cKDTree(mesh_xy)
    inner_distance, inner_position = mesh_tree.query(inner_xy, k=1, workers=-1)
    # Node exclusivity is intentionally not imposed. At a confluence an inner
    # tributary arc can legitimately share a node reached by another bank arc.
    valid = inner_distance <= max_arc_node_distance_m
    catalog = pd.DataFrame(
        {
            "river_idx": inner_river[valid],
            "local_arc_idx": inner_arc[valid],
            "station_idx": inner_station[valid],
            "arc_x_m": inner_xy[valid, 0],
            "arc_y_m": inner_xy[valid, 1],
            "mesh_position": inner_position[valid].astype(np.int64),
            "arc_to_mesh_m": inner_distance[valid],
        }
    )
    catalog.sort_values(
        ["river_idx", "mesh_position", "arc_to_mesh_m"], inplace=True
    )
    catalog.drop_duplicates(["river_idx", "mesh_position"], inplace=True)
    stats = {
        "river_arcs_in_selection": len(metric_arcs),
        "inner_arc_vertices": len(inner_xy),
        "inner_vertices_beyond_node_guard": int(
            (inner_distance > max_arc_node_distance_m).sum()
        ),
        "strict_inner_river_node_pairs": len(catalog),
        "strict_inner_mesh_nodes": int(catalog.mesh_position.nunique()),
    }
    return catalog.reset_index(drop=True), stats


def match_skeleton_centerlines_to_rivermapper(
    skeleton: gpd.GeoDataFrame,
    rivermapper: gpd.GeoDataFrame,
    candidate_radius_m: float = 500.0,
    max_hausdorff_m: float = 1_000.0,
) -> pd.DataFrame:
    """Choose one whole-line RiverMapper counterpart for each skeleton line."""
    left = skeleton.to_crs(MATCH_CRS).reset_index(drop=True)
    right = rivermapper.to_crs(MATCH_CRS).reset_index(drop=True)
    left_geometry = left.geometry.to_numpy()
    right_geometry = right.geometry.to_numpy()
    skeleton_idx, rivermapper_idx = STRtree(right_geometry).query(
        left_geometry, predicate="dwithin", distance=candidate_radius_m
    )
    candidates = pd.DataFrame(
        {
            "_skeleton_position": skeleton_idx.astype(np.int64),
            "_rivermapper_position": rivermapper_idx.astype(np.int64),
        }
    )
    if len(candidates):
        candidates["hausdorff_m"] = np.asarray(
            hausdorff_distance(
                left_geometry[candidates._skeleton_position],
                right_geometry[candidates._rivermapper_position],
            ),
            dtype=float,
        )
        candidates["length_difference_m"] = np.abs(
            left.length.to_numpy()[candidates._skeleton_position]
            - right.length.to_numpy()[candidates._rivermapper_position]
        )
        candidates.sort_values(
            ["_skeleton_position", "hausdorff_m", "length_difference_m"],
            inplace=True,
        )
        best = candidates.drop_duplicates("_skeleton_position").set_index(
            "_skeleton_position"
        )
    else:
        best = pd.DataFrame()
    records = pd.DataFrame(
        {
            "skeleton_river_idx": left.river_idx.to_numpy(dtype=np.int64),
            "source_arc_river_idx": np.full(len(left), -1, dtype=np.int64),
            "river_match_hausdorff_m": np.full(len(left), np.nan),
            "river_match_length_difference_m": np.full(len(left), np.nan),
            "river_match_status": np.full(len(left), "unmatched", dtype=object),
        }
    )
    if len(best):
        positions = best.index.to_numpy(dtype=np.int64)
        scores = best.hausdorff_m.to_numpy(dtype=float)
        accepted = np.isfinite(scores) & (scores <= max_hausdorff_m)
        positions = positions[accepted]
        chosen = best.loc[positions]
        records.loc[positions, "source_arc_river_idx"] = right.river_idx.to_numpy(
            dtype=np.int64
        )[chosen._rivermapper_position.to_numpy(dtype=np.int64)]
        records.loc[positions, "river_match_hausdorff_m"] = chosen.hausdorff_m.to_numpy()
        records.loc[positions, "river_match_length_difference_m"] = (
            chosen.length_difference_m.to_numpy()
        )
        records.loc[positions, "river_match_status"] = "matched"
    return records


def explode_skeleton_vertices(
    centerlines: gpd.GeoDataFrame, correspondence: pd.DataFrame
) -> pd.DataFrame:
    """Return every retained skeleton coordinate with stable feature indices."""
    frames = []
    for row in centerlines.itertuples(index=False):
        xy = np.asarray(row.geometry.coords, dtype=float)[:, :2]
        frames.append(
            pd.DataFrame(
                {
                    "skeleton_river_idx": int(row.river_idx),
                    "vertex_idx": np.arange(len(xy), dtype=np.int32),
                    "skeleton_x_m": xy[:, 0],
                    "skeleton_y_m": xy[:, 1],
                }
            )
        )
    vertices = pd.concat(frames, ignore_index=True)
    return vertices.merge(
        correspondence, on="skeleton_river_idx", how="left", validate="many_to_one"
    )


def map_skeleton_vertices_to_strict_inner_nodes(
    vertices: pd.DataFrame,
    coverage_geometry,
    catalog: pd.DataFrame,
    mesh_xy: np.ndarray,
    node_id: np.ndarray,
    hgrid_x: np.ndarray,
    hgrid_y: np.ndarray,
    hgrid_dp: np.ndarray,
    association_radius_m: float = 500.0,
) -> pd.DataFrame:
    """Use strict inner nodes for associated skeleton vertices, else nearest."""
    output = vertices.copy()
    source_xy = output[["skeleton_x_m", "skeleton_y_m"]].to_numpy(dtype=float)
    prepare(coverage_geometry)
    inside = np.asarray(
        intersects_xy(coverage_geometry, source_xy[:, 0], source_xy[:, 1]),
        dtype=bool,
    )
    arc_distance = np.full(len(output), np.inf)
    catalog_idx = np.full(len(output), -1, dtype=np.int64)
    catalog_groups = {
        int(river_idx): group
        for river_idx, group in catalog.groupby("river_idx", sort=False)
    }
    for source_river_idx, positions in output.groupby(
        "source_arc_river_idx", sort=False
    ).groups.items():
        if int(source_river_idx) < 0:
            continue
        river_catalog = catalog_groups.get(int(source_river_idx))
        if river_catalog is None:
            continue
        row_positions = np.asarray(positions, dtype=np.int64)
        distances, local_idx = cKDTree(
            river_catalog[["arc_x_m", "arc_y_m"]].to_numpy(dtype=float)
        ).query(source_xy[row_positions], k=1, workers=1)
        arc_distance[row_positions] = distances
        catalog_idx[row_positions] = river_catalog.index.to_numpy(dtype=np.int64)[
            local_idx
        ]
    _, nearest_mesh_position = cKDTree(mesh_xy).query(
        source_xy, k=1, workers=-1
    )
    associated = (
        inside & (catalog_idx >= 0) & (arc_distance <= association_radius_m)
    )
    chosen_position = nearest_mesh_position.astype(np.int64)
    chosen_position[associated] = catalog.loc[
        catalog_idx[associated], "mesh_position"
    ].to_numpy(dtype=np.int64)
    chosen_position[~inside] = -1
    valid = chosen_position >= 0
    output["inside_hgrid"] = inside
    output["river_associated"] = associated
    output["mapping_method"] = np.select(
        [~inside, associated],
        ["outside_hgrid", "strict_inner_river_node"],
        default="nearest_mesh_node",
    )
    output["nearest_inner_arc_vertex_m"] = arc_distance
    output["mesh_position"] = chosen_position
    output["mesh_node_id"] = -1
    output["mesh_x"] = np.nan
    output["mesh_y"] = np.nan
    output["mesh_dp"] = np.nan
    output["distance_m"] = np.nan
    output.loc[valid, "mesh_node_id"] = node_id[chosen_position[valid]]
    output.loc[valid, "mesh_x"] = hgrid_x[chosen_position[valid]]
    output.loc[valid, "mesh_y"] = hgrid_y[chosen_position[valid]]
    output.loc[valid, "mesh_dp"] = hgrid_dp[chosen_position[valid]]
    output.loc[valid, "distance_m"] = np.linalg.norm(
        source_xy[valid] - mesh_xy[chosen_position[valid]], axis=1
    )
    output["source_local_arc_idx"] = -1
    output["source_station_idx"] = -1
    output.loc[associated, "source_local_arc_idx"] = catalog.loc[
        catalog_idx[associated], "local_arc_idx"
    ].to_numpy()
    output.loc[associated, "source_station_idx"] = catalog.loc[
        catalog_idx[associated], "station_idx"
    ].to_numpy()
    return output


def get_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--matches-gpkg", type=Path, required=True)
    parser.add_argument("--river-arcs", type=Path, required=True)
    parser.add_argument("--river-centerlines", type=Path, required=True)
    parser.add_argument("--hgrid-cache", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--association-radius-m", type=float, default=500.0)
    parser.add_argument("--max-arc-node-distance-m", type=float, default=50.0)
    parser.add_argument("--river-match-radius-m", type=float, default=500.0)
    parser.add_argument("--max-river-hausdorff-m", type=float, default=1_000.0)
    return parser


def main(argv: Sequence[str] | None = None) -> int:
    args = get_parser().parse_args(argv)
    node_id = np.load(args.hgrid_cache / "node_id.npy", mmap_mode="r")
    hgrid_x = np.load(args.hgrid_cache / "x.npy", mmap_mode="r")
    hgrid_y = np.load(args.hgrid_cache / "y.npy", mmap_mode="r")
    hgrid_dp = np.load(args.hgrid_cache / "dp.npy", mmap_mode="r")
    mesh_xy = np.load(args.hgrid_cache / "metric_xy.npy", mmap_mode="r")
    coverage = gpd.read_parquet(args.hgrid_cache / "hgrid_coverage.parquet")
    coverage_geometry = coverage.geometry.iloc[0]
    selected = gpd.read_file(args.matches_gpkg, layer="selected_intervals").to_crs(
        MATCH_CRS
    )
    centerlines = gpd.read_file(args.matches_gpkg, layer="river_centerlines").to_crs(
        MATCH_CRS
    )
    arcs = gpd.read_parquet(args.river_arcs)
    rivermapper_centerlines = gpd.read_parquet(args.river_centerlines)
    correspondence = match_skeleton_centerlines_to_rivermapper(
        centerlines,
        rivermapper_centerlines,
        candidate_radius_m=args.river_match_radius_m,
        max_hausdorff_m=args.max_river_hausdorff_m,
    )
    matched_river_ids = set(
        correspondence.loc[
            correspondence.river_match_status == "matched", "source_arc_river_idx"
        ].astype(int)
    )
    catalog, catalog_stats = build_strict_inner_arc_catalog(
        arcs,
        mesh_xy,
        river_ids=matched_river_ids,
        selection_geometry=coverage_geometry,
        selection_buffer_m=args.association_radius_m,
        max_arc_node_distance_m=args.max_arc_node_distance_m,
    )
    vertex_input = explode_skeleton_vertices(centerlines, correspondence)
    mapped = map_skeleton_vertices_to_strict_inner_nodes(
        vertex_input,
        coverage_geometry,
        catalog,
        mesh_xy,
        node_id,
        hgrid_x,
        hgrid_y,
        hgrid_dp,
        association_radius_m=args.association_radius_m,
    )
    mapped_points = gpd.GeoDataFrame(
        mapped,
        geometry=points(mapped.skeleton_x_m, mapped.skeleton_y_m),
        crs=MATCH_CRS,
    )

    valid = mapped.mesh_position.to_numpy(dtype=np.int64) >= 0
    valid_rows = mapped.loc[valid].copy()
    source_xy = valid_rows[["skeleton_x_m", "skeleton_y_m"]].to_numpy(dtype=float)
    target_xy = mesh_xy[valid_rows.mesh_position.to_numpy(dtype=np.int64)]
    links = gpd.GeoDataFrame(
        valid_rows.drop(columns=["skeleton_x_m", "skeleton_y_m"]),
        geometry=[LineString((source, target)) for source, target in zip(source_xy, target_xy)],
        crs=MATCH_CRS,
    )

    node_groups = catalog.groupby("mesh_position", sort=True)
    unique_positions = np.asarray(list(node_groups.groups), dtype=np.int64)
    in_river = gpd.GeoDataFrame(
        {
            "mesh_position": unique_positions,
            "mesh_node_id": node_id[unique_positions],
            "mesh_dp": hgrid_dp[unique_positions],
            "river_count": [group.river_idx.nunique() for _, group in node_groups],
            "river_idx_list": [
                ",".join(group.river_idx.astype(str).drop_duplicates())
                for _, group in node_groups
            ],
        },
        geometry=points(mesh_xy[unique_positions, 0], mesh_xy[unique_positions, 1]),
        crs=MATCH_CRS,
    )
    if args.output.exists():
        args.output.unlink()
    layers = {
        "hgrid_coverage": coverage,
        "river_centerlines": centerlines,
        "selected_intervals": selected,
        "skeleton_rivermapper_correspondence": centerlines.merge(
            correspondence,
            left_on="river_idx",
            right_on="skeleton_river_idx",
            how="left",
            validate="one_to_one",
        ),
        "in_river_mesh_nodes": in_river,
        "skeleton_vertex_mapping": mapped_points,
        "mapping_links": links,
        "outside_hgrid_vertices": mapped_points.loc[~mapped_points.inside_hgrid],
        "long_links_gt_500m": links.loc[links.distance_m > 500.0],
    }
    for layer_name, frame in layers.items():
        if not frame.empty:
            frame.to_file(args.output, layer=layer_name, driver="GPKG", mode="w")
    summary = {
        "skeleton_centerlines": len(centerlines),
        "skeleton_vertices": len(mapped_points),
        "in_river_mesh_nodes": len(in_river),
        "association_radius_m": args.association_radius_m,
        "max_arc_node_distance_m": args.max_arc_node_distance_m,
        "river_match_radius_m": args.river_match_radius_m,
        "max_river_hausdorff_m": args.max_river_hausdorff_m,
        "matched_skeleton_centerlines": int(
            (correspondence.river_match_status == "matched").sum()
        ),
        "unmatched_skeleton_centerlines": int(
            (correspondence.river_match_status != "matched").sum()
        ),
        "watershed_filter_applied": False,
        "strict_inner_catalog": catalog_stats,
        "mapping_methods": mapped.mapping_method.value_counts().astype(int).to_dict(),
        "outside_hgrid_vertices": int((~mapped.inside_hgrid).sum()),
        "long_links_gt_500m": int((links.distance_m > 500.0).sum()),
        "layers": {name: len(frame) for name, frame in layers.items()},
    }
    args.output.with_suffix(".summary.json").write_text(
        json.dumps(summary, indent=2, sort_keys=True) + "\n"
    )
    print(json.dumps(summary, indent=2, sort_keys=True), flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
