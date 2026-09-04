#!/usr/bin/env python3
"""Map hydrofabric reaches back to nearby in-river SCHISM mesh nodes.

The default mapping unit is one complete ``COMID`` geometry.  It returns every
in-river mesh node within the requested point-to-line distance.  The optional
vertex mode instead returns one nearest in-river mesh node for every coordinate
of every LineString part.  When either query has no viable in-river candidate,
the nearest node in the complete mesh is returned and explicitly marked as a
fallback.

All spatial decisions are made in the workflow's metric matching CRS.  The
COMID implementation searches line segments rather than only line vertices,
so nodes beside a long, sparsely digitized reach are not missed.
"""

from __future__ import annotations

import math
from typing import Iterator, Literal

import geopandas as gpd
import numpy as np
import pandas as pd
from scipy.spatial import cKDTree
from shapely import distance as shapely_distance
from shapely import points as shapely_points
from shapely.geometry import LineString, MultiLineString
from shapely.ops import unary_union

try:
    from stofs3d_setup.ops.Bathy_edit.Ensure_channel_connectivity.hydrofabric_dredge import (
        HgridNodes,
        project_lonlat,
    )
    from stofs3d_setup.ops.Bathy_edit.Ensure_channel_connectivity.hydrofabric_match import (
        MATCH_CRS,
    )
except ModuleNotFoundError:  # Permit direct execution from this directory.
    from hydrofabric_dredge import HgridNodes, project_lonlat  # type: ignore[no-redef]
    from hydrofabric_match import MATCH_CRS  # type: ignore[no-redef]


MappingGranularity = Literal["comid", "vertex"]

INVERSE_MAPPING_COLUMNS = [
    "COMID",
    "granularity",
    "part_idx",
    "vertex_idx",
    "mesh_position",
    "mesh_node_id",
    "mesh_x",
    "mesh_y",
    "mesh_dp",
    "distance_m",
    "match_method",
    "is_fallback",
    "candidate_count",
    "target_x_m",
    "target_y_m",
]


def _iter_lines(geometry) -> Iterator[LineString]:
    """Yield non-empty LineStrings from a linear or collection geometry."""
    if geometry is None or geometry.is_empty:
        return
    if isinstance(geometry, LineString):
        yield geometry
        return
    if isinstance(geometry, MultiLineString) or hasattr(geometry, "geoms"):
        for part in geometry.geoms:
            yield from _iter_lines(part)


def _group_comid_geometries(hydrofabric: gpd.GeoDataFrame) -> list[tuple[object, object]]:
    records: list[tuple[object, object]] = []
    for comid, group in hydrofabric.groupby("COMID", sort=True, dropna=False):
        if pd.isna(comid):
            continue
        geometry = unary_union(group.geometry.to_numpy())
        if geometry is not None and not geometry.is_empty:
            records.append((comid, geometry))
    return records


def _segment_radius_candidates(
    geometry,
    mesh_xy: np.ndarray,
    tree: cKDTree,
    radius_m: float,
) -> np.ndarray:
    """Find all mesh positions within an exact distance of a line geometry."""
    candidates: set[int] = set()
    for line in _iter_lines(geometry):
        coordinates = np.asarray(line.coords, dtype=float)[:, :2]
        for start, end in zip(coordinates[:-1], coordinates[1:]):
            length = float(np.linalg.norm(end - start))
            pieces = max(1, int(math.ceil(length / radius_m)))
            edges = np.linspace(start, end, pieces + 1)
            for sub_start, sub_end in zip(edges[:-1], edges[1:]):
                midpoint = (sub_start + sub_end) / 2.0
                half_length = float(np.linalg.norm(sub_end - sub_start)) / 2.0
                candidates.update(
                    tree.query_ball_point(midpoint, radius_m + half_length)
                )
    if not candidates:
        return np.empty(0, dtype=np.int64)
    positions = np.fromiter(candidates, dtype=np.int64)
    points = shapely_points(mesh_xy[positions, 0], mesh_xy[positions, 1])
    distances = np.asarray(shapely_distance(geometry, points), dtype=float)
    return positions[distances <= radius_m + 1.0e-9]


def _nearest_mesh_position_to_line(
    geometry,
    mesh_xy: np.ndarray,
    tree: cKDTree,
) -> tuple[int, float]:
    """Return the exact nearest mesh point to a linear geometry.

    For each line segment, the nearest node to its midpoint establishes an
    upper distance bound.  Every point that could improve that bound lies in a
    ball of radius ``upper_bound + half_segment_length`` around the midpoint.
    Exact point-to-segment distances within those balls therefore recover the
    global nearest mesh node without constructing millions of Shapely points.
    """
    best_position = -1
    best_distance = math.inf
    for line in _iter_lines(geometry):
        coordinates = np.asarray(line.coords, dtype=float)[:, :2]
        if len(coordinates) == 1:
            distance, position = tree.query(coordinates[0], k=1)
            if float(distance) < best_distance:
                best_position, best_distance = int(position), float(distance)
            continue
        for start, end in zip(coordinates[:-1], coordinates[1:]):
            midpoint = (start + end) / 2.0
            _, seed_position = tree.query(midpoint, k=1)
            seed_position = int(seed_position)
            segment = LineString((start, end))
            upper_bound = float(
                shapely_distance(
                    segment,
                    shapely_points(
                        mesh_xy[seed_position, 0], mesh_xy[seed_position, 1]
                    ),
                )
            )
            half_length = float(np.linalg.norm(end - start)) / 2.0
            candidates = np.asarray(
                tree.query_ball_point(midpoint, upper_bound + half_length + 1.0e-9),
                dtype=np.int64,
            )
            points = shapely_points(mesh_xy[candidates, 0], mesh_xy[candidates, 1])
            distances = np.asarray(shapely_distance(segment, points), dtype=float)
            local = int(np.argmin(distances))
            distance = float(distances[local])
            position = int(candidates[local])
            if (distance, position) < (best_distance, best_position):
                best_position, best_distance = position, distance
    if best_position < 0:
        raise ValueError("Hydrofabric COMID has no non-empty linear geometry")
    return best_position, best_distance


def _mapping_record(
    *,
    comid: object,
    granularity: MappingGranularity,
    part_idx: int | None,
    vertex_idx: int | None,
    mesh_position: int,
    hgrid: HgridNodes,
    distance_m: float,
    match_method: str,
    candidate_count: int,
    target_xy: tuple[float, float] | None,
) -> dict[str, object]:
    return {
        "COMID": comid,
        "granularity": granularity,
        "part_idx": part_idx,
        "vertex_idx": vertex_idx,
        "mesh_position": mesh_position,
        "mesh_node_id": int(hgrid.node_id[mesh_position]),
        "mesh_x": float(hgrid.x[mesh_position]),
        "mesh_y": float(hgrid.y[mesh_position]),
        "mesh_dp": float(hgrid.dp[mesh_position]),
        "distance_m": distance_m,
        "match_method": match_method,
        "is_fallback": match_method == "nearest_mesh_fallback",
        "candidate_count": candidate_count,
        "target_x_m": target_xy[0] if target_xy is not None else np.nan,
        "target_y_m": target_xy[1] if target_xy is not None else np.nan,
    }


def map_hydrofabric_to_mesh(
    hydrofabric: gpd.GeoDataFrame,
    hgrid: HgridNodes,
    in_river_mask: np.ndarray,
    search_radius_m: float,
    granularity: MappingGranularity = "comid",
    query_workers: int = -1,
    node_hydrofabric_matches: pd.DataFrame | None = None,
) -> pd.DataFrame:
    """Map each COMID, or each of its vertices, back to SCHISM mesh nodes.

    ``comid`` mode returns every in-river node within ``search_radius_m`` of
    the complete grouped COMID geometry. ``vertex`` mode returns the nearest
    in-river node within the radius for each coordinate. If no viable in-river
    candidate exists, both modes return the nearest unrestricted mesh node and
    mark it as ``nearest_mesh_fallback``. If ``node_hydrofabric_matches`` is
    supplied, it must contain ``mesh_node_id`` and ``COMID``; only in-river
    nodes already associated with the current COMID are viable. This optional
    filter makes the operation a strict inverse of an existing node-to-reach
    mapping and avoids borrowing a neighboring reach's nodes at confluences.
    """
    if hydrofabric.crs is None:
        raise ValueError("Hydrofabric has no CRS")
    if "COMID" not in hydrofabric.columns:
        raise ValueError("Hydrofabric is missing required column: COMID")
    if search_radius_m <= 0:
        raise ValueError("search_radius_m must be positive")
    if granularity not in {"comid", "vertex"}:
        raise ValueError("granularity must be 'comid' or 'vertex'")
    if len(hgrid.node_id) == 0:
        raise ValueError("No hgrid nodes are available")
    mask = np.asarray(in_river_mask, dtype=bool)
    if len(mask) != len(hgrid.node_id):
        raise ValueError("in_river_mask must match the hgrid node count")

    hydro_metric = hydrofabric.to_crs(MATCH_CRS)
    mesh_xy = project_lonlat(hgrid.x, hgrid.y)
    mesh_tree = cKDTree(mesh_xy)
    in_river_positions = np.flatnonzero(mask)
    in_river_xy = mesh_xy[in_river_positions]
    in_river_tree = cKDTree(in_river_xy) if len(in_river_xy) else None
    matched_positions_by_comid: dict[object, set[int]] | None = None
    if node_hydrofabric_matches is not None:
        required = {"mesh_node_id", "COMID"}
        if missing := sorted(required - set(node_hydrofabric_matches.columns)):
            raise ValueError(
                f"node_hydrofabric_matches is missing columns: {missing}"
            )
        position_by_node_id = {
            int(node_id): position
            for position, node_id in enumerate(hgrid.node_id)
            if mask[position]
        }
        matched_positions_by_comid = {}
        for row in node_hydrofabric_matches.itertuples(index=False):
            if pd.isna(row.COMID) or pd.isna(row.mesh_node_id):
                continue
            position = position_by_node_id.get(int(row.mesh_node_id))
            if position is not None:
                matched_positions_by_comid.setdefault(row.COMID, set()).add(position)

    def corresponding_positions(comid: object, positions: np.ndarray) -> np.ndarray:
        if matched_positions_by_comid is None:
            return positions
        allowed = matched_positions_by_comid.get(comid, set())
        return np.asarray(
            [position for position in positions if int(position) in allowed],
            dtype=np.int64,
        )

    records: list[dict[str, object]] = []

    for comid, geometry in _group_comid_geometries(hydro_metric):
        if granularity == "comid":
            local_positions = (
                _segment_radius_candidates(
                    geometry, in_river_xy, in_river_tree, search_radius_m
                )
                if in_river_tree is not None
                else np.empty(0, dtype=np.int64)
            )
            if len(local_positions):
                positions = in_river_positions[local_positions]
                positions = corresponding_positions(comid, positions)
            else:
                positions = np.empty(0, dtype=np.int64)
            if len(positions):
                points = shapely_points(mesh_xy[positions, 0], mesh_xy[positions, 1])
                distances = np.asarray(shapely_distance(geometry, points), dtype=float)
                order = np.lexsort((hgrid.node_id[positions], distances))
                for position, distance in zip(positions[order], distances[order]):
                    records.append(
                        _mapping_record(
                            comid=comid,
                            granularity=granularity,
                            part_idx=None,
                            vertex_idx=None,
                            mesh_position=int(position),
                            hgrid=hgrid,
                            distance_m=float(distance),
                            match_method="in_river_radius",
                            candidate_count=len(positions),
                            target_xy=None,
                        )
                    )
            else:
                position, distance = _nearest_mesh_position_to_line(
                    geometry, mesh_xy, mesh_tree
                )
                records.append(
                    _mapping_record(
                        comid=comid,
                        granularity=granularity,
                        part_idx=None,
                        vertex_idx=None,
                        mesh_position=position,
                        hgrid=hgrid,
                        distance_m=distance,
                        match_method="nearest_mesh_fallback",
                        candidate_count=0,
                        target_xy=None,
                    )
                )
            continue

        for part_idx, line in enumerate(_iter_lines(geometry)):
            coordinates = np.asarray(line.coords, dtype=float)[:, :2]
            for vertex_idx, coordinate in enumerate(coordinates):
                local_candidates = (
                    in_river_tree.query_ball_point(
                        coordinate, search_radius_m, workers=query_workers
                    )
                    if in_river_tree is not None
                    else []
                )
                if local_candidates:
                    local = np.asarray(local_candidates, dtype=np.int64)
                    positions = corresponding_positions(
                        comid, in_river_positions[local]
                    )
                    local = np.searchsorted(in_river_positions, positions)
                if len(local_candidates) and len(local):
                    distances = np.linalg.norm(
                        in_river_xy[local] - coordinate, axis=1
                    )
                    node_ids = hgrid.node_id[positions]
                    selected = int(np.lexsort((node_ids, distances))[0])
                    position = int(positions[selected])
                    distance = float(distances[selected])
                    method = "in_river_radius"
                    candidate_count = len(positions)
                else:
                    distance, position = mesh_tree.query(
                        coordinate, k=1, workers=query_workers
                    )
                    position, distance = int(position), float(distance)
                    method = "nearest_mesh_fallback"
                    candidate_count = 0
                records.append(
                    _mapping_record(
                        comid=comid,
                        granularity=granularity,
                        part_idx=part_idx,
                        vertex_idx=vertex_idx,
                        mesh_position=position,
                        hgrid=hgrid,
                        distance_m=distance,
                        match_method=method,
                        candidate_count=candidate_count,
                        target_xy=(float(coordinate[0]), float(coordinate[1])),
                    )
                )

    result = pd.DataFrame.from_records(records, columns=INVERSE_MAPPING_COLUMNS)
    if not result.empty:
        result["COMID"] = pd.array(result.COMID, dtype="Int64")
        result["part_idx"] = pd.array(result.part_idx, dtype="Int64")
        result["vertex_idx"] = pd.array(result.vertex_idx, dtype="Int64")
    return result
