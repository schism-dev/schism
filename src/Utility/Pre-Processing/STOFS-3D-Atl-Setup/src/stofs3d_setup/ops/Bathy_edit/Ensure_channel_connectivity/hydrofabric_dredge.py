#!/usr/bin/env python3
"""Apply hydrofabric bankfull depths through RiverMapper vertices and a KD-tree.

``selected_intervals`` is the authoritative COMID/depth mapping along each
RiverMapper centerline.  The mapping is materialized at RiverMapper stations,
propagated across the parallel arcs, and associated with nearest hgrid nodes.
Bank vertices provide the reference depth; only inner arc vertices request
dredging.  Multiple requests for one mesh node are reduced with a maximum so
the deepest request wins deterministically.  A reverse mesh-node search also
recovers multi-river junction nodes, screens them against local half-width and
signed channel-polygon/bank distance, and merges accepted targets with the
forward requests.

This module provides the serial production implementation used by the
``Ensure_channel_connectivity`` Bathy-edit task and a command-line driver for
standalone validation runs.
"""

from __future__ import annotations

import argparse
import copy
import json
import math
import os
import shutil
from dataclasses import asdict, dataclass, field
from pathlib import Path
from typing import Sequence

import geopandas as gpd
import numpy as np
import pandas as pd
from pyproj import CRS, Transformer
from scipy.spatial import cKDTree
from shapely import distance as shapely_distance
from shapely import intersects_xy, points as shapely_points, prepare
from shapely.geometry import LineString, Polygon, box
from shapely.validation import make_valid

try:
    from stofs3d_setup.ops.Bathy_edit.Ensure_channel_connectivity.hydrofabric_match import (
        ANALYSIS_GEOGRAPHIC_CRS,
        MATCH_CRS,
        normalized_bbox,
    )
except ModuleNotFoundError:  # Permit direct execution from this directory.
    from hydrofabric_match import (  # type: ignore[no-redef]
        ANALYSIS_GEOGRAPHIC_CRS,
        MATCH_CRS,
        normalized_bbox,
    )


DEFAULT_CONNECTIVITY_DIR = Path(
    "/sciclone/schism10/feiye/CIROH/Channel_Geometry/"
    "Connectivity_test/full_domain"
)
DEFAULT_HGRID = Path(
    "/sciclone/schism10/Hgrid_projects/STOFS3D-v7.4/v32e/"
    "Bathy_edit_channel_variations/Channel_variations_verified_pre_channel/"
    "pre_channel/hgrid_pre_channel.gr3"
)
DEFAULT_ARCS = DEFAULT_CONNECTIVITY_DIR / "input_cache/river_arcs.parquet"
DEFAULT_CENTERLINES = (
    DEFAULT_CONNECTIVITY_DIR / "input_cache/river_centerlines.parquet"
)
DEFAULT_MATCHES = DEFAULT_CONNECTIVITY_DIR / "hydrofabric_river_matches.gpkg"
DEFAULT_OUTPUT_DIR = Path("/tmp/hydrofabric_dredge_test")
DEFAULT_WATERSHED = Path(
    "/sciclone/schism10/Hgrid_projects/STOFS3D-v7.4/v32e/"
    "Clip/outputs/watershed.shp"
)
DEFAULT_EXCLUDE_REGIONS = (
    Path(
        "/sciclone/schism10/feiye/STOFS3D-v8/I15_v7/Bathy_edit/"
        "RiverArc_Dredge/watershed_ME.shp"
    ),
)


@dataclass(frozen=True)
class DredgeSettings:
    min_channel_depth_m: float = 1.0
    channel_depth_source: str = "hydrofabric"
    measured_from_high_bank: bool = True
    max_nearest_distance_m: float = 500.0
    intersection_search_radius_m: float = 200.0
    intersection_width_tolerance_m: float = 10.0
    intersection_bank_exclusion_fraction: float = 0.05
    # Rounded Test 1 P99 (6.156 m). Capping preserves connectivity across
    # upper-tail requests, including a small number of major-river nodes.
    max_dredging_delta_m: float = 6.0
    intersection_recovery: bool = True
    unmatched_policy: str = "baseline"
    query_workers: int = -1

    def __post_init__(self) -> None:
        if self.min_channel_depth_m < 0:
            raise ValueError("min_channel_depth_m must be non-negative")
        if self.channel_depth_source not in {"hydrofabric", "constant"}:
            raise ValueError(
                "channel_depth_source must be 'hydrofabric' or 'constant'"
            )
        if self.max_nearest_distance_m <= 0:
            raise ValueError("max_nearest_distance_m must be positive")
        if self.intersection_search_radius_m < 0:
            raise ValueError("intersection_search_radius_m must be non-negative")
        if self.intersection_width_tolerance_m < 0:
            raise ValueError("intersection_width_tolerance_m must be non-negative")
        if not 0 <= self.intersection_bank_exclusion_fraction < 0.5:
            raise ValueError(
                "intersection_bank_exclusion_fraction must be in [0, 0.5)"
            )
        if self.max_dredging_delta_m <= 0:
            raise ValueError("max_dredging_delta_m must be positive")
        if self.unmatched_policy not in {"baseline", "skip"}:
            raise ValueError("unmatched_policy must be 'baseline' or 'skip'")


@dataclass
class HgridNodes:
    node_id: np.ndarray
    x: np.ndarray
    y: np.ndarray
    dp: np.ndarray
    total_node_count: int


@dataclass
class DredgeResult:
    vertices: pd.DataFrame
    node_changes: pd.DataFrame
    summary: dict[str, object]
    requested_nodes: pd.DataFrame = field(default_factory=pd.DataFrame)
    intersection_candidates: pd.DataFrame = field(default_factory=pd.DataFrame)
    intersection_targets: pd.DataFrame = field(default_factory=pd.DataFrame)


def _to_analysis_geographic(frame: gpd.GeoDataFrame, context: str) -> gpd.GeoDataFrame:
    """Project a mask to the RiverMapper geographic CRS.

    The production Maine exclusion is WGS84 while RiverMapper coordinates are
    treated as NAD83.  Some installations lack the optional WGS84-to-NAD83
    datum grid and return non-finite coordinates.  In that specific
    geographic-to-geographic case, retain the finite longitude/latitude
    coordinates and apply the same NAD83 approximation used elsewhere in this
    driver.
    """
    if frame.crs is None:
        raise ValueError(f"{context} has no CRS")
    source_crs = CRS.from_user_input(frame.crs)
    target_crs = CRS.from_user_input(ANALYSIS_GEOGRAPHIC_CRS)
    if (
        source_crs.is_geographic
        and target_crs.is_geographic
        and source_crs != target_crs
    ):
        if not np.all(np.isfinite(frame.total_bounds)):
            raise ValueError(f"{context} contains non-finite coordinates")
        return frame.set_crs(ANALYSIS_GEOGRAPHIC_CRS, allow_override=True)
    try:
        projected = frame.to_crs(ANALYSIS_GEOGRAPHIC_CRS)
        if np.all(np.isfinite(projected.total_bounds)):
            return projected
    except Exception:
        projected = None
    if source_crs.is_geographic and target_crs.is_geographic:
        if not np.all(np.isfinite(frame.total_bounds)):
            raise ValueError(f"{context} contains non-finite coordinates")
        return frame.set_crs(ANALYSIS_GEOGRAPHIC_CRS, allow_override=True)
    raise ValueError(f"Could not project {context} to {ANALYSIS_GEOGRAPHIC_CRS}")


def load_effective_watershed(
    watershed_path: Path | Sequence[Path],
    exclude_region_paths: Sequence[Path] = (),
) -> gpd.GeoDataFrame:
    """Load one or more connectivity regions and subtract their exclusions."""
    if isinstance(watershed_path, (str, Path)):
        watershed_paths = (Path(watershed_path),)
    else:
        watershed_paths = tuple(Path(path) for path in watershed_path)
    if not watershed_paths:
        raise ValueError("At least one watershed path is required")

    effective_geometry = None
    for path in watershed_paths:
        watershed = _to_analysis_geographic(
            gpd.read_file(path), f"watershed {path}"
        )
        watershed_geometry = watershed.geometry.unary_union
        effective_geometry = (
            watershed_geometry
            if effective_geometry is None
            else effective_geometry.union(watershed_geometry)
        )
    for exclusion_path in exclude_region_paths:
        exclusion = _to_analysis_geographic(
            gpd.read_file(exclusion_path), f"watershed exclusion {exclusion_path}"
        )
        effective_geometry = effective_geometry.difference(
            exclusion.geometry.unary_union
        )
    if effective_geometry.is_empty:
        raise ValueError("The effective watershed is empty after exclusions")
    return gpd.GeoDataFrame(
        {"name": ["effective_watershed"]},
        geometry=[effective_geometry],
        crs=ANALYSIS_GEOGRAPHIC_CRS,
    )


def points_intersect_region(
    x: np.ndarray | pd.Series,
    y: np.ndarray | pd.Series,
    region: gpd.GeoDataFrame,
) -> np.ndarray:
    """Return a boundary-inclusive point mask for a polygonal region."""
    region_lonlat = _to_analysis_geographic(region, "point-selection region")
    geometry = region_lonlat.geometry.unary_union
    prepare(geometry)
    return np.asarray(
        intersects_xy(
            geometry,
            np.asarray(x, dtype=float),
            np.asarray(y, dtype=float),
        ),
        dtype=bool,
    )


def line_cumulative_measures(line) -> np.ndarray:
    coordinates = np.asarray(line.coords, dtype=float)[:, :2]
    if len(coordinates) == 0:
        return np.empty(0, dtype=float)
    return np.r_[0.0, np.cumsum(np.linalg.norm(np.diff(coordinates, axis=0), axis=1))]


def _validate_intervals(intervals: pd.DataFrame, tolerance_m: float = 1.0e-6) -> None:
    if intervals.empty:
        return
    ordered = intervals.sort_values(["start_m", "end_m"])
    starts = ordered.start_m.to_numpy(dtype=float)
    ends = ordered.end_m.to_numpy(dtype=float)
    if np.any(ends <= starts):
        raise ValueError("Selected intervals must have positive length")
    if len(ordered) > 1 and np.any(starts[1:] < ends[:-1] - tolerance_m):
        raise ValueError("Selected intervals overlap within one river")


def assign_intervals_to_stations(
    centerlines: gpd.GeoDataFrame,
    selected_intervals: gpd.GeoDataFrame,
) -> gpd.GeoDataFrame:
    """Materialize piecewise selected intervals at RiverMapper stations.

    Shared interval boundaries belong to the interval starting at that measure
    (half-open ``[start_m, end_m)`` behavior).  The final river endpoint is
    included in the final interval.
    """
    required_centerline = {"river_idx", "geometry"}
    required_interval = {
        "river_idx",
        "start_m",
        "end_m",
        "COMID",
        "bnk_depth",
        "geometry",
    }
    if missing := sorted(required_centerline - set(centerlines.columns)):
        raise ValueError(f"Centerlines are missing columns: {missing}")
    if missing := sorted(required_interval - set(selected_intervals.columns)):
        raise ValueError(f"Selected intervals are missing columns: {missing}")

    center_metric = centerlines.to_crs(MATCH_CRS)
    selected = selected_intervals.drop(columns="geometry").copy()
    interval_groups = {
        int(river_idx): group.sort_values(["start_m", "end_m"]).reset_index(drop=True)
        for river_idx, group in selected.groupby("river_idx", sort=False)
    }
    frames: list[gpd.GeoDataFrame] = []
    optional_fields = ["bnk_width", "quality_flag", "continuity_basis"]
    for river in center_metric.itertuples(index=False):
        river_idx = int(river.river_idx)
        coordinates = np.asarray(river.geometry.coords, dtype=float)[:, :2]
        measures = line_cumulative_measures(river.geometry)
        count = len(measures)
        data: dict[str, object] = {
            "river_idx": np.full(count, river_idx, dtype=np.int64),
            "station_idx": np.arange(count, dtype=np.int64),
            "centerline_m": measures,
            "COMID": pd.array([pd.NA] * count, dtype="Int64"),
            "bnk_depth": np.full(count, np.nan),
            "match_status": np.full(count, "unmatched", dtype=object),
            "match_method": np.full(count, None, dtype=object),
            "bnk_width": np.full(count, np.nan),
            "quality_flag": np.full(count, None, dtype=object),
            "continuity_basis": np.full(count, None, dtype=object),
        }
        intervals = interval_groups.get(river_idx)
        if intervals is not None and not intervals.empty:
            if "quality_flag" in intervals:
                continuity = intervals.loc[
                    intervals.quality_flag == "continuity_fill"
                ].copy()
                primary = intervals.loc[
                    intervals.quality_flag != "continuity_fill"
                ].copy()
            else:
                primary = intervals
                continuity = intervals.iloc[0:0]
            _validate_intervals(primary)

            # Continuity intervals can overlap the edge of a primary interval
            # because station-neighbor fills use station midpoints. Preserve
            # the original precedence: primary line matches first, then use
            # continuity only for stations that remain unmatched.
            for interval_pass, unmatched_only in (
                (primary, False),
                (continuity, True),
            ):
                if interval_pass.empty:
                    continue
                interval_pass = interval_pass.sort_values(
                    ["start_m", "end_m"]
                ).reset_index(drop=True)
                starts = interval_pass.start_m.to_numpy(dtype=float)
                ends = interval_pass.end_m.to_numpy(dtype=float)
                interval_positions = np.searchsorted(
                    starts, measures, side="right"
                ) - 1
                candidate = interval_positions >= 0
                clipped_positions = np.clip(
                    interval_positions, 0, len(interval_pass) - 1
                )
                candidate &= measures < ends[clipped_positions] + 1.0e-8
                final_endpoint = np.isclose(
                    measures, measures[-1], atol=1.0e-6
                )
                candidate |= final_endpoint & (
                    clipped_positions == len(interval_pass) - 1
                ) & (measures <= ends[clipped_positions] + 1.0e-6)
                if unmatched_only:
                    candidate &= data["match_status"] == "unmatched"
                station_positions = np.flatnonzero(candidate)
                if len(station_positions) == 0:
                    continue
                matched_intervals = interval_pass.iloc[
                    clipped_positions[candidate]
                ]
                data["COMID"][station_positions] = pd.array(
                    matched_intervals.COMID.astype(np.int64), dtype="Int64"
                )
                data["bnk_depth"][station_positions] = (
                    matched_intervals.bnk_depth.to_numpy(dtype=float)
                )
                data["match_status"][station_positions] = "matched"
                quality = (
                    matched_intervals.quality_flag.to_numpy(dtype=object)
                    if "quality_flag" in matched_intervals
                    else np.full(
                        len(matched_intervals), "nominal", dtype=object
                    )
                )
                data["match_method"][station_positions] = np.where(
                    quality == "continuity_fill",
                    "continuity_fill",
                    "line_overlap",
                )
                for field in optional_fields:
                    if field in matched_intervals:
                        data[field][station_positions] = matched_intervals[
                            field
                        ].to_numpy()
        frame = gpd.GeoDataFrame(
            data,
            geometry=gpd.points_from_xy(coordinates[:, 0], coordinates[:, 1]),
            crs=MATCH_CRS,
        )
        frames.append(frame)
    if not frames:
        return gpd.GeoDataFrame(
            columns=[
                "river_idx",
                "station_idx",
                "centerline_m",
                "COMID",
                "bnk_depth",
                "match_status",
                "match_method",
                "bnk_width",
                "quality_flag",
                "continuity_basis",
                "geometry",
            ],
            geometry="geometry",
            crs=MATCH_CRS,
        )
    return gpd.GeoDataFrame(
        pd.concat(frames, ignore_index=True), geometry="geometry", crs=MATCH_CRS
    )


def expand_stations_to_arc_vertices(
    arcs: gpd.GeoDataFrame, stations: gpd.GeoDataFrame
) -> pd.DataFrame:
    """Propagate centerline station attributes across every parallel arc."""
    required = {"river_idx", "local_arc_idx", "n_arcs", "geometry"}
    if missing := sorted(required - set(arcs.columns)):
        raise ValueError(f"River arcs are missing columns: {missing}")
    station_groups = {
        int(river_idx): group.sort_values("station_idx").reset_index(drop=True)
        for river_idx, group in stations.drop(columns="geometry").groupby(
            "river_idx", sort=False
        )
    }
    frames: list[pd.DataFrame] = []
    station_fields = [
        "centerline_m",
        "COMID",
        "bnk_depth",
        "bnk_width",
        "match_status",
        "match_method",
        "quality_flag",
        "continuity_basis",
    ]
    for region_field in (
        "station_inside_test_region",
        "station_inside_watershed",
        "station_inside_dredge_region",
    ):
        if region_field in stations:
            station_fields.append(region_field)
    for arc in arcs.itertuples(index=False):
        river_idx = int(arc.river_idx)
        station_group = station_groups.get(river_idx)
        if station_group is None:
            continue
        coordinates = np.asarray(arc.geometry.coords, dtype=float)[:, :2]
        if len(coordinates) != len(station_group):
            raise ValueError(
                f"River {river_idx}, arc {arc.local_arc_idx} has {len(coordinates)} "
                f"vertices but its centerline has {len(station_group)} stations"
            )
        count = len(coordinates)
        frame = pd.DataFrame(
            {
                "river_idx": np.full(count, river_idx, dtype=np.int64),
                "local_arc_idx": np.full(count, int(arc.local_arc_idx), dtype=np.int64),
                "station_idx": np.arange(count, dtype=np.int64),
                "n_arcs": np.full(count, int(arc.n_arcs), dtype=np.int64),
                "is_bank": np.full(
                    count,
                    int(arc.local_arc_idx) in {0, int(arc.n_arcs) - 1},
                    dtype=bool,
                ),
                "x": coordinates[:, 0],
                "y": coordinates[:, 1],
            }
        )
        for field in station_fields:
            frame[field] = station_group[field].to_numpy()
        frames.append(frame)
    if not frames:
        return pd.DataFrame()
    vertices = pd.concat(frames, ignore_index=True)
    vertices["COMID"] = vertices.COMID.astype("Int64")
    return vertices


def read_hgrid_nodes(
    path: Path,
    bbox_lonlat: tuple[float, float, float, float] | None = None,
    bbox_buffer_degrees: float = 0.05,
) -> HgridNodes:
    """Read hgrid node coordinates/depths, optionally limiting a test region."""
    if bbox_lonlat is not None:
        xmin, ymin, xmax, ymax = bbox_lonlat
        xmin -= bbox_buffer_degrees
        ymin -= bbox_buffer_degrees
        xmax += bbox_buffer_degrees
        ymax += bbox_buffer_degrees
    node_ids: list[int] = []
    x_values: list[float] = []
    y_values: list[float] = []
    depths: list[float] = []
    with path.open("r", encoding="utf-8") as stream:
        stream.readline()
        fields = stream.readline().split()
        if len(fields) < 2:
            raise ValueError(f"Invalid hgrid header: {path}")
        total_nodes = int(fields[1])
        for _ in range(total_nodes):
            row = stream.readline().split()
            node_id = int(row[0])
            x, y, dp = map(float, row[1:4])
            if bbox_lonlat is None or (xmin <= x <= xmax and ymin <= y <= ymax):
                node_ids.append(node_id)
                x_values.append(x)
                y_values.append(y)
                depths.append(dp)
    return HgridNodes(
        node_id=np.asarray(node_ids, dtype=np.int64),
        x=np.asarray(x_values),
        y=np.asarray(y_values),
        dp=np.asarray(depths),
        total_node_count=total_nodes,
    )


def project_lonlat(x: np.ndarray, y: np.ndarray) -> np.ndarray:
    """Project coordinates using the same NAD83 approximation as matching."""
    transformer = Transformer.from_crs(
        ANALYSIS_GEOGRAPHIC_CRS, MATCH_CRS, always_xy=True
    )
    projected_x, projected_y = transformer.transform(x, y)
    return np.column_stack((projected_x, projected_y))


def map_vertices_to_mesh(
    vertices: pd.DataFrame,
    hgrid: HgridNodes,
    query_workers: int = -1,
    mesh_candidate_mask: np.ndarray | None = None,
) -> pd.DataFrame:
    """Associate every RiverMapper vertex with its nearest eligible hgrid node."""
    if len(hgrid.node_id) == 0:
        raise ValueError("No hgrid nodes are available for the KD-tree")
    if mesh_candidate_mask is None:
        candidate_positions = np.arange(len(hgrid.node_id), dtype=np.int64)
    else:
        candidate_mask = np.asarray(mesh_candidate_mask, dtype=bool)
        if len(candidate_mask) != len(hgrid.node_id):
            raise ValueError("mesh_candidate_mask must match the hgrid node count")
        candidate_positions = np.flatnonzero(candidate_mask)
    if len(candidate_positions) == 0:
        raise ValueError("No eligible hgrid nodes are available for the KD-tree")
    mapped = vertices.copy()
    mesh_xy = project_lonlat(
        hgrid.x[candidate_positions], hgrid.y[candidate_positions]
    )
    vertex_xy = project_lonlat(
        mapped.x.to_numpy(dtype=float), mapped.y.to_numpy(dtype=float)
    )
    tree = cKDTree(mesh_xy)
    distances, tree_positions = tree.query(vertex_xy, k=1, workers=query_workers)
    positions = candidate_positions[tree_positions.astype(np.int64)]
    mapped["mesh_position"] = positions
    mapped["mesh_node_id"] = hgrid.node_id[positions]
    mapped["nearest_distance_m"] = distances
    mapped["mesh_dp"] = hgrid.dp[positions]
    return mapped


INTERSECTION_TARGET_COLUMNS = [
    "mesh_position",
    "mesh_node_id",
    "x",
    "y",
    "mesh_dp",
    "intersection_candidate_200m",
    "intersection_target",
    "forward_requested",
    "additional_target",
    "river_count",
    "strict_qualifying_river_count",
    "half_width_qualifying_river_count",
    "qualifying_river_count",
    "passes_strict_half_width",
    "passes_half_width_with_tolerance",
    "passes_bank_proximity_screen",
    "width_tolerance_m",
    "bank_exclusion_fraction",
    "nearby_vertex_count",
    "river_idx_list",
    "half_width_qualifying_river_idx_list",
    "qualifying_river_idx_list",
    "nearest_station_idx_list",
    "nearest_local_arc_idx_list",
    "river_distance_m_list",
    "centerline_distance_m_list",
    "local_rivermapper_width_m_list",
    "half_width_m_list",
    "normalized_centerline_distance_list",
    "passes_strict_half_width_list",
    "passes_tolerant_half_width_list",
    "nearest_bank_vertex_distance_m_list",
    "nearest_bank_line_distance_m_list",
    "signed_bank_distance_m_list",
    "inside_channel_polygon_list",
    "bank_vertex_proximity_ratio_list",
    "bank_proximity_ratio_list",
    "passes_bank_proximity_list",
    "passes_final_screen_list",
    "nearest_river_distance_m",
    "farthest_river_distance_m",
    "nearest_centerline_distance_m",
    "farthest_centerline_distance_m",
    "max_rivermapper_width_m",
    "max_hydrofabric_width_m",
    "proposed_target_dp",
    "target_dp_available",
    "would_deepen",
    "proposed_dredging_delta_m",
]


def _empty_intersection_targets() -> pd.DataFrame:
    return pd.DataFrame(columns=INTERSECTION_TARGET_COLUMNS)


def _station_rivermapper_widths(vertices: pd.DataFrame) -> pd.Series:
    """Return bank-to-bank width at each complete RiverMapper station."""
    banks = vertices.loc[vertices.is_bank, [
        "river_idx", "station_idx", "x", "y"
    ]].copy()
    if banks.empty:
        return pd.Series(dtype=float, name="rivermapper_width_m")
    bank_xy = project_lonlat(
        banks.x.to_numpy(dtype=float), banks.y.to_numpy(dtype=float)
    )
    banks["metric_x"] = bank_xy[:, 0]
    banks["metric_y"] = bank_xy[:, 1]

    def bank_span(group: pd.DataFrame) -> float:
        coordinates = group[["metric_x", "metric_y"]].to_numpy(dtype=float)
        if len(coordinates) < 2:
            return math.nan
        delta = coordinates[:, None, :] - coordinates[None, :, :]
        return float(np.sqrt(np.sum(delta * delta, axis=2)).max())

    widths = banks.groupby(
        ["river_idx", "station_idx"], sort=False
    )[["metric_x", "metric_y"]].apply(bank_span)
    widths.name = "rivermapper_width_m"
    return widths


def screen_intersection_mesh_nodes(
    vertices: pd.DataFrame,
    hgrid: HgridNodes,
    effective_watershed: gpd.GeoDataFrame,
    radius_m: float = 200.0,
    width_tolerance_m: float = 10.0,
    bank_exclusion_fraction: float = 0.05,
    bbox_lonlat: tuple[float, float, float, float] | None = None,
    query_workers: int = -1,
    mesh_candidate_mask: np.ndarray | None = None,
) -> pd.DataFrame:
    """Label mesh nodes lying within ``radius_m`` of multiple rivers.

    The spatial search retains one nearest RiverMapper vertex per
    ``(mesh node, river_idx)`` pair. A candidate becomes a target only when at
    least two rivers pass the local ``centerline distance <= half width +
    tolerance`` test. At least one of those rivers must lie farther than
    ``bank_exclusion_fraction * width`` from its continuous bank lines.
    Watershed-filtered vertices seed the broad search, while complete river
    geometry supplies centerlines, bank lines, and channel polygons so clipping
    cannot collapse a bank segment to its last in-watershed vertex. The
    returned targets are merged with forward vertex requests by the serial and
    MPI drivers; the deepest finite request wins. The source hgrid is changed
    only when the caller subsequently writes the resulting node changes.
    """
    if radius_m <= 0 or vertices.empty or len(hgrid.node_id) == 0:
        return _empty_intersection_targets()

    if mesh_candidate_mask is None:
        mesh_eligible = np.ones(len(hgrid.node_id), dtype=bool)
    else:
        mesh_eligible = np.asarray(mesh_candidate_mask, dtype=bool).copy()
        if len(mesh_eligible) != len(hgrid.node_id):
            raise ValueError("mesh_candidate_mask must match the hgrid node count")
    candidate_positions = np.flatnonzero(mesh_eligible)
    inside_region = np.zeros(len(hgrid.node_id), dtype=bool)
    inside_region[candidate_positions] = points_intersect_region(
        hgrid.x[candidate_positions],
        hgrid.y[candidate_positions],
        effective_watershed,
    )
    mesh_eligible &= inside_region
    if bbox_lonlat is not None:
        xmin, ymin, xmax, ymax = bbox_lonlat
        mesh_eligible &= (
            (hgrid.x >= xmin)
            & (hgrid.x <= xmax)
            & (hgrid.y >= ymin)
            & (hgrid.y <= ymax)
        )
    eligible_positions = np.flatnonzero(mesh_eligible)
    source = vertices.loc[vertices.inside_watershed].copy()
    if len(eligible_positions) == 0 or source.empty:
        return _empty_intersection_targets()

    mesh_xy = project_lonlat(
        hgrid.x[eligible_positions], hgrid.y[eligible_positions]
    )
    mesh_tree = cKDTree(mesh_xy)
    source_xy = project_lonlat(
        source.x.to_numpy(dtype=float), source.y.to_numpy(dtype=float)
    )
    source["metric_x"] = source_xy[:, 0]
    source["metric_y"] = source_xy[:, 1]
    widths = _station_rivermapper_widths(vertices)

    def adjacent_max(values: np.ndarray, positions: np.ndarray) -> np.ndarray:
        output = np.full(len(positions), np.nan)
        for output_idx, station_position in enumerate(positions):
            neighborhood = values[
                max(0, station_position - 1):min(len(values), station_position + 2)
            ]
            finite = neighborhood[np.isfinite(neighborhood)]
            if len(finite):
                output[output_idx] = float(np.max(finite))
        return output

    river_node_frames: list[pd.DataFrame] = []
    complete_rivers = vertices.groupby("river_idx", sort=False)
    for river_idx, search_group in source.groupby("river_idx", sort=False):
        group_xy = search_group[["metric_x", "metric_y"]].to_numpy(dtype=float)
        neighborhoods = mesh_tree.query_ball_point(
            group_xy, r=radius_m, workers=query_workers
        )
        counts = np.fromiter((len(items) for items in neighborhoods), dtype=np.int64)
        if counts.sum() == 0:
            continue
        local_mesh_position = np.concatenate(neighborhoods).astype(np.int64)
        group_position = np.repeat(
            np.arange(len(search_group), dtype=np.int64), counts
        )
        distance = np.linalg.norm(
            mesh_xy[local_mesh_position] - group_xy[group_position], axis=1
        )
        pairs = pd.DataFrame(
            {
                "local_mesh_position": local_mesh_position,
                "group_position": group_position,
                "distance_m": distance,
            }
        )
        vertex_counts = pairs.groupby("local_mesh_position", sort=False).size()
        nearest = pairs.loc[
            pairs.groupby("local_mesh_position", sort=False).distance_m.idxmin()
        ].copy()
        nearest["nearby_vertex_count"] = nearest.local_mesh_position.map(
            vertex_counts
        ).to_numpy(dtype=np.int64)
        nearest_vertices = search_group.iloc[
            nearest.group_position.to_numpy(dtype=np.int64)
        ]
        nearest["river_idx"] = int(river_idx)
        nearest["local_arc_idx"] = nearest_vertices.local_arc_idx.to_numpy(
            dtype=np.int64
        )

        geometry_group = complete_rivers.get_group(river_idx).copy()
        geometry_xy = project_lonlat(
            geometry_group.x.to_numpy(dtype=float),
            geometry_group.y.to_numpy(dtype=float),
        )
        geometry_group["metric_x"] = geometry_xy[:, 0]
        geometry_group["metric_y"] = geometry_xy[:, 1]

        station_data = geometry_group.groupby("station_idx", sort=True).agg(
            center_x=("metric_x", "mean"),
            center_y=("metric_y", "mean"),
            hydrofabric_width_m=("bnk_width", "first"),
            target_thalweg_dp=("target_thalweg_dp", "first"),
            station_bank_valid=("station_bank_valid", "first"),
        )
        station_ids = station_data.index.to_numpy(dtype=np.int64)
        station_xy = station_data[["center_x", "center_y"]].to_numpy(dtype=float)
        candidate_xy = mesh_xy[
            nearest.local_mesh_position.to_numpy(dtype=np.int64)
        ]
        nearest_station_positions = cKDTree(station_xy).query(
            candidate_xy, k=1, workers=query_workers
        )[1].astype(np.int64)
        nearest["station_idx"] = station_ids[nearest_station_positions]
        if len(station_xy) == 1:
            centerline_distance = np.linalg.norm(
                candidate_xy - station_xy[0], axis=1
            )
        else:
            centerline = LineString(station_xy)
            centerline_distance = np.asarray(
                shapely_distance(
                    centerline,
                    shapely_points(candidate_xy[:, 0], candidate_xy[:, 1]),
                ),
                dtype=float,
            )
        nearest["centerline_distance_m"] = centerline_distance

        station_keys = pd.MultiIndex.from_arrays(
            [np.full(len(station_ids), int(river_idx)), station_ids],
            names=["river_idx", "station_idx"],
        )
        station_widths = widths.reindex(station_keys).to_numpy(dtype=float)
        nearest["rivermapper_width_m"] = adjacent_max(
            station_widths, nearest_station_positions
        )
        hydrofabric_widths = station_data.hydrofabric_width_m.to_numpy(dtype=float)
        nearest["hydrofabric_width_m"] = adjacent_max(
            hydrofabric_widths, nearest_station_positions
        )
        proposed = station_data.target_thalweg_dp.to_numpy(dtype=float)[
            nearest_station_positions
        ]
        # Keep geometric targets for review but withhold their depth when a
        # very wide/boundary channel cannot map both banks within the 500 m
        # guard. Such channels are intentionally outside the small-watershed
        # scope of this workflow.
        valid_target = station_data.station_bank_valid.to_numpy(dtype=bool)[
            nearest_station_positions
        ]
        nearest["target_thalweg_dp"] = np.where(valid_target, proposed, np.nan)
        local_width = nearest.rivermapper_width_m.to_numpy(dtype=float)
        half_width = local_width / 2.0
        width_valid = np.isfinite(local_width) & (local_width > 0)
        nearest["half_width_m"] = half_width
        nearest["normalized_centerline_distance"] = np.where(
            width_valid, 2.0 * centerline_distance / local_width, np.nan
        )
        nearest["passes_strict_half_width"] = (
            width_valid & (centerline_distance <= half_width)
        )
        nearest["passes_tolerant_half_width"] = (
            width_valid
            & (centerline_distance <= half_width + width_tolerance_m)
        )
        nearest_bank_distance = np.full(len(nearest), np.nan)
        bank_vertices = geometry_group.loc[geometry_group.is_bank]
        for candidate_idx, station_position in enumerate(nearest_station_positions):
            local_station_ids = station_ids[
                max(0, station_position - 1):min(
                    len(station_ids), station_position + 2
                )
            ]
            local_banks = bank_vertices.loc[
                bank_vertices.station_idx.isin(local_station_ids),
                ["metric_x", "metric_y"],
            ].to_numpy(dtype=float)
            if len(local_banks):
                nearest_bank_distance[candidate_idx] = float(
                    np.linalg.norm(
                        local_banks - candidate_xy[candidate_idx], axis=1
                    ).min()
                )
        bank_ratio = np.where(
            width_valid, nearest_bank_distance / local_width, np.nan
        )
        nearest["nearest_bank_vertex_distance_m"] = nearest_bank_distance
        nearest["bank_vertex_proximity_ratio"] = bank_ratio

        candidate_points = shapely_points(candidate_xy[:, 0], candidate_xy[:, 1])
        bank_line_distance = np.full(len(nearest), np.inf)
        bank_coordinate_sets: list[np.ndarray] = []
        for _, bank_arc in bank_vertices.groupby("local_arc_idx", sort=False):
            bank_coordinates = (
                bank_arc.sort_values("station_idx")[["metric_x", "metric_y"]]
                .drop_duplicates()
                .to_numpy(dtype=float)
            )
            if len(bank_coordinates) == 0:
                continue
            bank_coordinate_sets.append(bank_coordinates)
            bank_geometry = (
                LineString(bank_coordinates)
                if len(bank_coordinates) > 1
                else shapely_points(bank_coordinates[0, 0], bank_coordinates[0, 1])
            )
            bank_line_distance = np.minimum(
                bank_line_distance,
                np.asarray(
                    shapely_distance(bank_geometry, candidate_points), dtype=float
                ),
            )
        bank_line_distance[~np.isfinite(bank_line_distance)] = np.nan
        inside_channel_polygon = centerline_distance <= half_width
        if (
            len(bank_coordinate_sets) >= 2
            and len(bank_coordinate_sets[0]) > 1
            and len(bank_coordinate_sets[-1]) > 1
        ):
            polygon_coordinates = np.vstack(
                [
                    bank_coordinate_sets[0],
                    bank_coordinate_sets[-1][::-1],
                    bank_coordinate_sets[0][:1],
                ]
            )
            channel_polygon = Polygon(polygon_coordinates)
            if not channel_polygon.is_valid:
                channel_polygon = make_valid(channel_polygon)
            inside_channel_polygon = np.asarray(
                intersects_xy(
                    channel_polygon, candidate_xy[:, 0], candidate_xy[:, 1]
                ),
                dtype=bool,
            )
        signed_bank_distance = np.where(
            inside_channel_polygon, bank_line_distance, -bank_line_distance
        )
        signed_bank_ratio = np.where(
            width_valid, signed_bank_distance / local_width, np.nan
        )
        nearest["nearest_bank_line_distance_m"] = bank_line_distance
        nearest["signed_bank_distance_m"] = signed_bank_distance
        nearest["inside_channel_polygon"] = inside_channel_polygon
        nearest["bank_proximity_ratio"] = signed_bank_ratio
        nearest["passes_bank_proximity"] = (
            width_valid
            & inside_channel_polygon
            & np.isfinite(signed_bank_distance)
            & (signed_bank_ratio > bank_exclusion_fraction)
        )
        nearest["passes_final_screen"] = (
            nearest.passes_tolerant_half_width
            & nearest.passes_bank_proximity
        )
        river_node_frames.append(nearest)

    if not river_node_frames:
        return _empty_intersection_targets()
    river_nodes = pd.concat(river_node_frames, ignore_index=True)
    grouped = river_nodes.groupby("local_mesh_position", sort=False)
    candidate_local_positions = grouped.river_idx.nunique()
    candidate_local_positions = candidate_local_positions.loc[
        candidate_local_positions >= 2
    ].index.to_numpy(dtype=np.int64)
    if len(candidate_local_positions) == 0:
        return _empty_intersection_targets()

    records: list[dict[str, object]] = []
    forward_requested_positions = set(
        vertices.loc[vertices.dredge_request, "mesh_position"].astype(int)
    )
    for local_position in candidate_local_positions:
        contributing = river_nodes.loc[
            river_nodes.local_mesh_position == local_position
        ].sort_values("river_idx")
        hgrid_position = int(eligible_positions[local_position])
        strict_qualifying = contributing.passes_strict_half_width.to_numpy(dtype=bool)
        half_width_qualifying = contributing.passes_tolerant_half_width.to_numpy(
            dtype=bool
        )
        strict_qualifying_count = int(strict_qualifying.sum())
        half_width_qualifying_count = int(half_width_qualifying.sum())
        passes_strict = strict_qualifying_count >= 2
        passes_half_width = half_width_qualifying_count >= 2
        final_qualifying = (
            half_width_qualifying
            & contributing.passes_bank_proximity.to_numpy(dtype=bool)
        )
        qualifying_count = int(final_qualifying.sum())
        any_qualifying_river_passes_bank_proximity = qualifying_count >= 1
        intersection_target = (
            passes_half_width and any_qualifying_river_passes_bank_proximity
        )
        half_width_contributors = contributing.loc[half_width_qualifying]
        target_contributors = contributing.loc[final_qualifying]
        proposed_values = target_contributors.target_thalweg_dp.to_numpy(dtype=float)
        finite_target = np.isfinite(proposed_values)
        proposed_target = (
            float(np.max(proposed_values[finite_target]))
            if finite_target.any()
            else math.nan
        )
        mesh_dp = float(hgrid.dp[hgrid_position])
        would_deepen = bool(
            intersection_target
            and np.isfinite(proposed_target)
            and proposed_target > mesh_dp
        )
        forward_requested = hgrid_position in forward_requested_positions
        records.append(
            {
                "mesh_position": hgrid_position,
                "mesh_node_id": int(hgrid.node_id[hgrid_position]),
                "x": float(hgrid.x[hgrid_position]),
                "y": float(hgrid.y[hgrid_position]),
                "mesh_dp": mesh_dp,
                "intersection_candidate_200m": True,
                "intersection_target": intersection_target,
                "forward_requested": forward_requested,
                "additional_target": intersection_target and not forward_requested,
                "river_count": int(contributing.river_idx.nunique()),
                "strict_qualifying_river_count": strict_qualifying_count,
                "half_width_qualifying_river_count": half_width_qualifying_count,
                "qualifying_river_count": qualifying_count,
                "passes_strict_half_width": passes_strict,
                "passes_half_width_with_tolerance": passes_half_width,
                "passes_bank_proximity_screen": (
                    any_qualifying_river_passes_bank_proximity
                ),
                "width_tolerance_m": float(width_tolerance_m),
                "bank_exclusion_fraction": float(bank_exclusion_fraction),
                "nearby_vertex_count": int(contributing.nearby_vertex_count.sum()),
                "river_idx_list": ",".join(
                    contributing.river_idx.astype(str).tolist()
                ),
                "half_width_qualifying_river_idx_list": ",".join(
                    half_width_contributors.river_idx.astype(str).tolist()
                ),
                "qualifying_river_idx_list": ",".join(
                    target_contributors.river_idx.astype(str).tolist()
                ),
                "nearest_station_idx_list": ",".join(
                    contributing.station_idx.astype(str).tolist()
                ),
                "nearest_local_arc_idx_list": ",".join(
                    contributing.local_arc_idx.astype(str).tolist()
                ),
                "river_distance_m_list": ",".join(
                    f"{value:.3f}" for value in contributing.distance_m
                ),
                "centerline_distance_m_list": ",".join(
                    f"{value:.3f}" for value in contributing.centerline_distance_m
                ),
                "local_rivermapper_width_m_list": ",".join(
                    f"{value:.3f}" for value in contributing.rivermapper_width_m
                ),
                "half_width_m_list": ",".join(
                    f"{value:.3f}" for value in contributing.half_width_m
                ),
                "normalized_centerline_distance_list": ",".join(
                    f"{value:.4f}"
                    for value in contributing.normalized_centerline_distance
                ),
                "passes_strict_half_width_list": ",".join(
                    "1" if value else "0" for value in strict_qualifying
                ),
                "passes_tolerant_half_width_list": ",".join(
                    "1" if value else "0" for value in half_width_qualifying
                ),
                "nearest_bank_vertex_distance_m_list": ",".join(
                    f"{value:.3f}"
                    for value in contributing.nearest_bank_vertex_distance_m
                ),
                "nearest_bank_line_distance_m_list": ",".join(
                    f"{value:.3f}"
                    for value in contributing.nearest_bank_line_distance_m
                ),
                "signed_bank_distance_m_list": ",".join(
                    f"{value:.3f}"
                    for value in contributing.signed_bank_distance_m
                ),
                "inside_channel_polygon_list": ",".join(
                    "1" if value else "0"
                    for value in contributing.inside_channel_polygon
                ),
                "bank_vertex_proximity_ratio_list": ",".join(
                    f"{value:.5f}"
                    for value in contributing.bank_vertex_proximity_ratio
                ),
                "bank_proximity_ratio_list": ",".join(
                    f"{value:.5f}" for value in contributing.bank_proximity_ratio
                ),
                "passes_bank_proximity_list": ",".join(
                    "1" if value else "0"
                    for value in contributing.passes_bank_proximity
                ),
                "passes_final_screen_list": ",".join(
                    "1" if value else "0" for value in final_qualifying
                ),
                "nearest_river_distance_m": float(contributing.distance_m.min()),
                "farthest_river_distance_m": float(contributing.distance_m.max()),
                "nearest_centerline_distance_m": float(
                    contributing.centerline_distance_m.min()
                ),
                "farthest_centerline_distance_m": float(
                    contributing.centerline_distance_m.max()
                ),
                "max_rivermapper_width_m": float(
                    contributing.rivermapper_width_m.max()
                ),
                "max_hydrofabric_width_m": float(
                    contributing.hydrofabric_width_m.max()
                ),
                "proposed_target_dp": proposed_target,
                "target_dp_available": bool(intersection_target and finite_target.any()),
                "would_deepen": would_deepen,
                "proposed_dredging_delta_m": (
                    max(0.0, proposed_target - mesh_dp)
                    if np.isfinite(proposed_target)
                    else math.nan
                ),
            }
        )
    return pd.DataFrame.from_records(records, columns=INTERSECTION_TARGET_COLUMNS)


REQUESTED_MESH_NODE_COLUMNS = [
    "mesh_position",
    "mesh_node_id",
    "x",
    "y",
    "original_dp",
    "requested_dp",
    "requested_dredging_delta_m",
    "final_dp",
    "dredging_delta_m",
    "would_deepen",
    "passes_max_dredging_delta",
    "dredging_delta_capped",
    "depth_screen_status",
    "forward_requested",
    "intersection_target",
    "request_source",
    "forward_vertex_requests",
    "intersection_qualifying_river_count",
    "request_count",
    "vertex_requests",
]


def finalize_mesh_requests(
    requests: pd.DataFrame,
    max_dredging_delta_m: float = 6.0,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Apply the shared serial/MPI depth cap and request classifications."""
    if max_dredging_delta_m <= 0:
        raise ValueError("max_dredging_delta_m must be positive")
    if requests.empty:
        empty = pd.DataFrame(columns=REQUESTED_MESH_NODE_COLUMNS)
        return empty, empty.copy()
    if requests.mesh_position.duplicated().any():
        raise ValueError("Mesh requests must be reduced to one row per mesh position")

    finalized = requests.copy()
    finalized["requested_dredging_delta_m"] = np.maximum(
        0.0, finalized.requested_dp - finalized.original_dp
    )
    requests_deepening = finalized.requested_dredging_delta_m > 1.0e-12
    finalized["passes_max_dredging_delta"] = (
        finalized.requested_dredging_delta_m
        <= max_dredging_delta_m + 1.0e-12
    )
    finalized["dredging_delta_capped"] = (
        requests_deepening & ~finalized.passes_max_dredging_delta
    )
    finalized["dredging_delta_m"] = np.minimum(
        finalized.requested_dredging_delta_m, max_dredging_delta_m
    )
    finalized["would_deepen"] = finalized.dredging_delta_m > 1.0e-12
    finalized["final_dp"] = finalized.original_dp + finalized.dredging_delta_m
    finalized["depth_screen_status"] = np.select(
        [~requests_deepening, finalized.dredging_delta_capped],
        ["no_deepening_needed", "capped_at_max_dredging_delta"],
        default="accepted",
    )
    finalized["request_source"] = np.select(
        [
            finalized.forward_requested & finalized.intersection_target,
            finalized.intersection_target,
        ],
        ["forward+intersection", "intersection"],
        default="forward",
    )
    finalized["request_count"] = (
        finalized.forward_vertex_requests
        + finalized.intersection_target.astype(np.int64)
    )
    # Preserve the existing changed-node field for downstream users.
    finalized["vertex_requests"] = finalized.forward_vertex_requests
    finalized = finalized[REQUESTED_MESH_NODE_COLUMNS]
    changed = finalized.loc[finalized.would_deepen].copy()
    return finalized, changed


def build_dredge_requested_mesh_nodes(
    vertices: pd.DataFrame,
    hgrid: HgridNodes,
    intersection_targets: pd.DataFrame | None = None,
    max_dredging_delta_m: float = 6.0,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Combine requests and cap implausibly large mesh-depth changes."""
    node_count = len(hgrid.node_id)
    requested_dp = np.full(node_count, -np.inf)
    forward_counts = np.zeros(node_count, dtype=np.int64)
    intersection_requested = np.zeros(node_count, dtype=bool)
    intersection_river_counts = np.zeros(node_count, dtype=np.int64)

    forward = vertices.loc[vertices.dredge_request]
    if not forward.empty:
        positions = forward.mesh_position.to_numpy(dtype=np.int64)
        np.maximum.at(
            requested_dp,
            positions,
            forward.target_thalweg_dp.to_numpy(dtype=float),
        )
        forward_counts = np.bincount(positions, minlength=node_count)

    if intersection_targets is not None and not intersection_targets.empty:
        valid = (
            intersection_targets.intersection_target.to_numpy(dtype=bool)
            & np.isfinite(
                intersection_targets.proposed_target_dp.to_numpy(dtype=float)
            )
        )
        intersections = intersection_targets.loc[valid]
        if not intersections.empty:
            positions = intersections.mesh_position.to_numpy(dtype=np.int64)
            np.maximum.at(
                requested_dp,
                positions,
                intersections.proposed_target_dp.to_numpy(dtype=float),
            )
            intersection_requested[positions] = True
            np.maximum.at(
                intersection_river_counts,
                positions,
                intersections.qualifying_river_count.to_numpy(dtype=np.int64),
            )

    positions = np.flatnonzero(np.isfinite(requested_dp))
    forward_requested = forward_counts[positions] > 0
    intersection_target = intersection_requested[positions]
    requests = pd.DataFrame(
        {
            "mesh_position": positions,
            "mesh_node_id": hgrid.node_id[positions],
            "x": hgrid.x[positions],
            "y": hgrid.y[positions],
            "original_dp": hgrid.dp[positions],
            "requested_dp": requested_dp[positions],
            "forward_requested": forward_requested,
            "intersection_target": intersection_target,
            "forward_vertex_requests": forward_counts[positions],
            "intersection_qualifying_river_count": (
                intersection_river_counts[positions]
            ),
        }
    )
    return finalize_mesh_requests(requests, max_dredging_delta_m)


def refresh_result_mesh_requests(
    result: DredgeResult,
    hgrid: HgridNodes,
    max_dredging_delta_m: float = 6.0,
) -> None:
    """Rebuild requested and changed nodes after intersection screening."""
    requested_nodes, node_changes = build_dredge_requested_mesh_nodes(
        result.vertices,
        hgrid,
        result.intersection_targets,
        max_dredging_delta_m=max_dredging_delta_m,
    )
    result.requested_nodes = requested_nodes
    result.node_changes = node_changes
    final_by_position = pd.Series(
        node_changes.final_dp.to_numpy(dtype=float),
        index=node_changes.mesh_position.to_numpy(dtype=np.int64),
    )
    result.vertices["final_node_dp"] = result.vertices.mesh_position.map(
        final_by_position
    )
    result.vertices["final_node_dp"] = result.vertices.final_node_dp.fillna(
        result.vertices.mesh_dp
    )
    result.vertices["node_dredging_delta_m"] = (
        result.vertices.final_node_dp - result.vertices.mesh_dp
    )
    changed_delta = node_changes.dredging_delta_m.to_numpy(dtype=float)
    result.summary.update(
        {
            "dredge_requested_mesh_nodes": len(requested_nodes),
            "changed_mesh_nodes": len(node_changes),
            "capped_dredging_mesh_nodes": int(
                requested_nodes.dredging_delta_capped.sum()
            ),
            "collision_mesh_nodes": int((node_changes.vertex_requests > 1).sum())
            if not node_changes.empty
            else 0,
            "dredging_delta_m": {
                "min": float(np.min(changed_delta)) if len(changed_delta) else None,
                "median": float(np.median(changed_delta))
                if len(changed_delta)
                else None,
                "p95": float(np.percentile(changed_delta, 95))
                if len(changed_delta)
                else None,
                "max": float(np.max(changed_delta)) if len(changed_delta) else None,
            },
        }
    )


def apply_station_depths(
    mapped_vertices: pd.DataFrame,
    hgrid: HgridNodes,
    settings: DredgeSettings,
) -> DredgeResult:
    """Compute bankfull targets and collision-safe mesh-node updates."""
    vertices = mapped_vertices.copy()
    keys = ["river_idx", "station_idx"]
    banks = vertices.loc[vertices.is_bank].copy()
    banks["bank_mapping_valid"] = (
        banks.nearest_distance_m <= settings.max_nearest_distance_m
    )
    bank_counts = banks.groupby(keys).bank_mapping_valid.sum().rename("valid_bank_count")
    if settings.measured_from_high_bank:
        bank_depth = banks.groupby(keys).mesh_dp.min().rename("bank_dp")
    else:
        bank_depth = banks.groupby(keys).mesh_dp.max().rename("bank_dp")
    station_bank = pd.concat([bank_depth, bank_counts], axis=1).reset_index()
    vertices = vertices.merge(station_bank, on=keys, how="left", validate="many_to_one")

    matched_depth = np.isfinite(vertices.bnk_depth.to_numpy(dtype=float)) & (
        vertices.match_status.to_numpy() == "matched"
    )
    if settings.channel_depth_source == "constant":
        effective_depth = np.full(len(vertices), settings.min_channel_depth_m)
    else:
        effective_depth = np.full(len(vertices), np.nan)
        effective_depth[matched_depth] = np.maximum(
            settings.min_channel_depth_m,
            vertices.loc[matched_depth, "bnk_depth"].to_numpy(dtype=float),
        )
        if settings.unmatched_policy == "baseline":
            effective_depth[~matched_depth] = settings.min_channel_depth_m
    vertices["effective_channel_depth_m"] = effective_depth
    vertices["target_thalweg_dp"] = vertices.bank_dp + effective_depth
    vertices["mapping_valid"] = (
        vertices.nearest_distance_m <= settings.max_nearest_distance_m
    )
    if "inside_test_region" not in vertices:
        vertices["inside_test_region"] = True
    if "inside_watershed" not in vertices:
        vertices["inside_watershed"] = True
    vertices["inside_test_region"] = vertices.inside_test_region.astype(bool)
    vertices["inside_watershed"] = vertices.inside_watershed.astype(bool)
    vertices["inside_dredge_region"] = (
        vertices.inside_test_region & vertices.inside_watershed
    )
    vertices["region_status"] = np.select(
        [
            ~vertices.inside_test_region,
            vertices.inside_test_region & ~vertices.inside_watershed,
        ],
        ["outside_test_region", "outside_watershed"],
        default="inside_dredge_region",
    )
    vertices["station_bank_valid"] = vertices.valid_bank_count >= 2
    vertices["dredge_request"] = (
        (~vertices.is_bank)
        & vertices.inside_dredge_region
        & vertices.mapping_valid
        & vertices.station_bank_valid
        & np.isfinite(vertices.target_thalweg_dp)
    )
    vertices["requested_dp"] = np.where(
        vertices.dredge_request,
        np.maximum(vertices.mesh_dp, vertices.target_thalweg_dp),
        vertices.mesh_dp,
    )

    requested_nodes, node_changes = build_dredge_requested_mesh_nodes(
        vertices,
        hgrid,
        max_dredging_delta_m=settings.max_dredging_delta_m,
    )
    final_by_position = pd.Series(
        node_changes.final_dp.to_numpy(dtype=float),
        index=node_changes.mesh_position.to_numpy(dtype=np.int64),
    )
    vertices["final_node_dp"] = vertices.mesh_position.map(final_by_position)
    vertices["final_node_dp"] = vertices.final_node_dp.fillna(vertices.mesh_dp)
    vertices["node_dredging_delta_m"] = vertices.final_node_dp - vertices.mesh_dp

    requested_count = int(vertices.dredge_request.sum())
    matched_request_count = int(
        (vertices.dredge_request & (vertices.match_status == "matched")).sum()
    )
    distance = vertices.loc[
        vertices.inside_dredge_region, "nearest_distance_m"
    ].to_numpy(dtype=float)
    changed_delta = node_changes.dredging_delta_m.to_numpy(dtype=float)
    station_values = vertices.sort_values(keys).drop_duplicates(keys)
    matched_station_depths = station_values.loc[
        station_values.match_status == "matched", "bnk_depth"
    ].to_numpy(dtype=float)
    requested_targets = vertices.loc[
        vertices.dredge_request, "target_thalweg_dp"
    ].to_numpy(dtype=float)
    inside_station_keys = vertices.loc[
        vertices.inside_dredge_region, keys
    ].drop_duplicates()
    matched_inside_station_keys = vertices.loc[
        vertices.inside_dredge_region & (vertices.match_status == "matched"),
        keys,
    ].drop_duplicates()
    summary: dict[str, object] = {
        "settings": asdict(settings),
        "river_count": int(vertices.river_idx.nunique()),
        "arc_vertex_count": len(vertices),
        "station_count": int(vertices[keys].drop_duplicates().shape[0]),
        "inside_dredge_region_station_count": len(inside_station_keys),
        "matched_station_count": int(
            vertices.loc[vertices.match_status == "matched", keys]
            .drop_duplicates()
            .shape[0]
        ),
        "matched_inside_dredge_region_station_count": len(
            matched_inside_station_keys
        ),
        "dredge_request_vertices": requested_count,
        "matched_dredge_request_vertices": matched_request_count,
        "fallback_dredge_request_vertices": requested_count - matched_request_count,
        "outside_dredge_region_vertices": int(
            (~vertices.inside_dredge_region).sum()
        ),
        "outside_test_region_vertices": int((~vertices.inside_test_region).sum()),
        "outside_watershed_vertices_inside_test_region": int(
            (vertices.inside_test_region & ~vertices.inside_watershed).sum()
        ),
        "invalid_distance_vertices_inside_dredge_region": int(
            (vertices.inside_dredge_region & ~vertices.mapping_valid).sum()
        ),
        "invalid_bank_station_count": int(
            vertices.loc[
                vertices.inside_dredge_region & ~vertices.station_bank_valid,
                keys,
            ].drop_duplicates().shape[0]
        ),
        "nearest_distance_m": {
            "min": float(np.min(distance)) if len(distance) else None,
            "median": float(np.median(distance)) if len(distance) else None,
            "p95": float(np.percentile(distance, 95)) if len(distance) else None,
            "p99": float(np.percentile(distance, 99)) if len(distance) else None,
            "max": float(np.max(distance)) if len(distance) else None,
        },
        "matched_bankfull_depth_m": {
            "min": float(np.min(matched_station_depths))
            if len(matched_station_depths)
            else None,
            "median": float(np.median(matched_station_depths))
            if len(matched_station_depths)
            else None,
            "p95": float(np.percentile(matched_station_depths, 95))
            if len(matched_station_depths)
            else None,
            "max": float(np.max(matched_station_depths))
            if len(matched_station_depths)
            else None,
        },
        "requested_target_thalweg_dp": {
            "min": float(np.min(requested_targets)) if len(requested_targets) else None,
            "median": float(np.median(requested_targets))
            if len(requested_targets)
            else None,
            "p95": float(np.percentile(requested_targets, 95))
            if len(requested_targets)
            else None,
            "max": float(np.max(requested_targets)) if len(requested_targets) else None,
        },
        "unique_requested_mesh_nodes": int(
            vertices.loc[vertices.dredge_request, "mesh_position"].nunique()
        ),
        "changed_mesh_nodes": len(node_changes),
        "capped_dredging_mesh_nodes": int(
            requested_nodes.dredging_delta_capped.sum()
        ),
        "collision_mesh_nodes": int((node_changes.vertex_requests > 1).sum())
        if not node_changes.empty
        else 0,
        "dredging_delta_m": {
            "min": float(np.min(changed_delta)) if len(changed_delta) else None,
            "median": float(np.median(changed_delta)) if len(changed_delta) else None,
            "p95": float(np.percentile(changed_delta, 95)) if len(changed_delta) else None,
            "max": float(np.max(changed_delta)) if len(changed_delta) else None,
        },
    }
    summary["dredge_requested_mesh_nodes"] = len(requested_nodes)
    return DredgeResult(
        vertices=vertices,
        node_changes=node_changes,
        summary=summary,
        requested_nodes=requested_nodes,
    )


def write_hgrid_with_updates(
    source: Path, output: Path, node_changes: pd.DataFrame
) -> None:
    """Stream-copy an hgrid while replacing only explicitly changed depths."""
    updates = dict(
        zip(
            node_changes.mesh_node_id.astype(int),
            node_changes.final_dp.astype(float),
        )
    )
    output.parent.mkdir(parents=True, exist_ok=True)
    temporary = output.with_suffix(output.suffix + ".tmp")
    with source.open("r", encoding="utf-8") as src, temporary.open(
        "w", encoding="utf-8"
    ) as dst:
        dst.write(src.readline())
        counts_line = src.readline()
        dst.write(counts_line)
        counts = counts_line.split()
        if len(counts) < 2:
            raise ValueError(f"Invalid hgrid header: {source}")
        total_nodes = int(counts[1])
        for _ in range(total_nodes):
            line = src.readline()
            fields = line.split()
            node_id = int(fields[0])
            if node_id in updates:
                dst.write(
                    f"{fields[0]} {fields[1]} {fields[2]} {updates[node_id]:.10f}\n"
                )
            else:
                dst.write(line)
        shutil.copyfileobj(src, dst, length=16 * 1024 * 1024)
    os.replace(temporary, output)


def write_2dm_with_updates(
    source: Path, output: Path, node_changes: pd.DataFrame
) -> None:
    """Stream-copy an SMS 2DM mesh while replacing selected ND depths."""
    updates = dict(
        zip(
            node_changes.mesh_node_id.astype(int),
            node_changes.final_dp.astype(float),
        )
    )
    output.parent.mkdir(parents=True, exist_ok=True)
    temporary = output.with_suffix(output.suffix + ".tmp")
    with source.open("r", encoding="utf-8") as src, temporary.open(
        "w", encoding="utf-8"
    ) as dst:
        for line in src:
            if not line.startswith("ND "):
                dst.write(line)
                continue
            fields = line.split()
            node_id = int(fields[1])
            if node_id in updates:
                fields[4] = f"{updates[node_id]:.10f}"
                dst.write(" ".join(fields) + "\n")
            else:
                dst.write(line)
    os.replace(temporary, output)


def _atomic_parquet(frame: pd.DataFrame, path: Path) -> None:
    temporary = path.with_suffix(".tmp.parquet")
    frame.to_parquet(temporary, index=False, compression="zstd")
    os.replace(temporary, path)


def write_diagnostics(
    output_dir: Path,
    result: DredgeResult,
    bbox_lonlat: tuple[float, float, float, float] | None,
    write_gpkg: bool,
    effective_watershed: gpd.GeoDataFrame | None = None,
) -> dict[str, str]:
    output_dir.mkdir(parents=True, exist_ok=True)
    vertex_path = output_dir / "river_vertex_mesh_mapping.parquet"
    requested_node_path = output_dir / "dredge_requested_mesh_nodes.parquet"
    node_path = output_dir / "dredged_mesh_nodes.parquet"
    intersection_candidate_path = output_dir / "intersection_candidates_200m.parquet"
    intersection_path = output_dir / "intersection_target_mesh_nodes.parquet"
    _atomic_parquet(result.vertices, vertex_path)
    _atomic_parquet(result.requested_nodes, requested_node_path)
    _atomic_parquet(result.node_changes, node_path)
    _atomic_parquet(result.intersection_candidates, intersection_candidate_path)
    _atomic_parquet(result.intersection_targets, intersection_path)
    products = {
        "vertex_mapping": str(vertex_path),
        "requested_nodes": str(requested_node_path),
        "node_changes": str(node_path),
        "intersection_candidates": str(intersection_candidate_path),
        "intersection_targets": str(intersection_path),
    }
    if write_gpkg:
        gpkg_path = output_dir / "hydrofabric_dredge_diagnostics.gpkg"
        if gpkg_path.exists():
            gpkg_path.unlink()
        vertices = result.vertices
        changed_vertices = vertices.loc[vertices.node_dredging_delta_m > 0].copy()
        dredge_requests = vertices.loc[vertices.dredge_request].copy()
        flagged_distance = vertices.loc[
            vertices.inside_dredge_region & ~vertices.mapping_valid
        ].copy()
        outside_watershed = vertices.loc[
            vertices.inside_test_region & ~vertices.inside_watershed
        ].copy()
        layers: dict[str, gpd.GeoDataFrame] = {}
        if not dredge_requests.empty:
            layers["dredge_request_vertices"] = gpd.GeoDataFrame(
                dredge_requests,
                geometry=gpd.points_from_xy(dredge_requests.x, dredge_requests.y),
                crs=ANALYSIS_GEOGRAPHIC_CRS,
            )
        if not changed_vertices.empty:
            layers["changed_vertices"] = gpd.GeoDataFrame(
                changed_vertices,
                geometry=gpd.points_from_xy(changed_vertices.x, changed_vertices.y),
                crs=ANALYSIS_GEOGRAPHIC_CRS,
            )
        if not result.requested_nodes.empty:
            requested_nodes = result.requested_nodes
            layers["dredge_requested_mesh_nodes"] = gpd.GeoDataFrame(
                requested_nodes,
                geometry=gpd.points_from_xy(requested_nodes.x, requested_nodes.y),
                crs=ANALYSIS_GEOGRAPHIC_CRS,
            )
        if not result.node_changes.empty:
            changes = result.node_changes
            layers["changed_mesh_nodes"] = gpd.GeoDataFrame(
                changes,
                geometry=gpd.points_from_xy(changes.x, changes.y),
                crs=ANALYSIS_GEOGRAPHIC_CRS,
            )
        if not result.requested_nodes.empty:
            capped = result.requested_nodes.loc[
                result.requested_nodes.dredging_delta_capped
            ].copy()
            if not capped.empty:
                layers["capped_dredging_delta"] = gpd.GeoDataFrame(
                    capped,
                    geometry=gpd.points_from_xy(capped.x, capped.y),
                    crs=ANALYSIS_GEOGRAPHIC_CRS,
                )
                raw_gt10 = capped.loc[
                    capped.requested_dredging_delta_m > 10.0
                ].copy()
                if not raw_gt10.empty:
                    layers["capped_raw_delta_gt10m"] = gpd.GeoDataFrame(
                        raw_gt10,
                        geometry=gpd.points_from_xy(raw_gt10.x, raw_gt10.y),
                        crs=ANALYSIS_GEOGRAPHIC_CRS,
                    )
        if not result.intersection_candidates.empty:
            candidates = result.intersection_candidates
            layers["intersection_candidates_200m"] = gpd.GeoDataFrame(
                candidates,
                geometry=gpd.points_from_xy(candidates.x, candidates.y),
                crs=ANALYSIS_GEOGRAPHIC_CRS,
            )
        if not result.intersection_targets.empty:
            targets = result.intersection_targets
            layers["intersection_target_mesh_nodes"] = gpd.GeoDataFrame(
                targets,
                geometry=gpd.points_from_xy(targets.x, targets.y),
                crs=ANALYSIS_GEOGRAPHIC_CRS,
            )
        if not flagged_distance.empty:
            layers["flagged_nearest_distance"] = gpd.GeoDataFrame(
                flagged_distance,
                geometry=gpd.points_from_xy(flagged_distance.x, flagged_distance.y),
                crs=ANALYSIS_GEOGRAPHIC_CRS,
            )
        if not outside_watershed.empty:
            layers["outside_watershed"] = gpd.GeoDataFrame(
                outside_watershed,
                geometry=gpd.points_from_xy(outside_watershed.x, outside_watershed.y),
                crs=ANALYSIS_GEOGRAPHIC_CRS,
            )
        if bbox_lonlat is not None:
            layers["test_region"] = gpd.GeoDataFrame(
                {"name": ["test_region"]},
                geometry=[box(*bbox_lonlat)],
                crs=ANALYSIS_GEOGRAPHIC_CRS,
            )
        if effective_watershed is not None:
            watershed_layer = effective_watershed.copy()
            if bbox_lonlat is not None:
                watershed_layer.geometry = watershed_layer.geometry.intersection(
                    box(*bbox_lonlat)
                )
                watershed_layer = watershed_layer.loc[
                    ~watershed_layer.geometry.is_empty
                ].copy()
            if not watershed_layer.empty:
                layers["effective_watershed"] = watershed_layer
        for layer_name, frame in layers.items():
            frame.to_file(gpkg_path, layer=layer_name, driver="GPKG", mode="w")
        products["diagnostics_gpkg"] = str(gpkg_path)
    return products


def run_hydrofabric_dredge(
    hgrid: HgridNodes,
    river_arcs_file: str | Path,
    river_centerlines_file: str | Path,
    matches_gpkg_file: str | Path,
    effective_watershed: gpd.GeoDataFrame,
    settings: DredgeSettings,
    bbox_lonlat: tuple[float, float, float, float] | None = None,
) -> DredgeResult:
    """Calculate hydrofabric-informed connectivity changes for one hgrid."""
    centerlines = gpd.read_parquet(river_centerlines_file)
    if bbox_lonlat is not None:
        selection = gpd.GeoDataFrame(
            geometry=[box(*bbox_lonlat)], crs=ANALYSIS_GEOGRAPHIC_CRS
        ).to_crs(centerlines.crs).geometry.iloc[0]
        centerlines = centerlines.loc[
            centerlines.geometry.intersects(selection)
        ].copy()
    if centerlines.empty:
        raise ValueError("No RiverMapper centerlines intersect the dredging region")

    river_ids = centerlines.river_idx.astype(np.int64).to_numpy()
    arcs = gpd.read_parquet(river_arcs_file)
    arcs = arcs.loc[arcs.river_idx.isin(river_ids)].copy()
    selected = gpd.read_file(matches_gpkg_file, layer="selected_intervals")
    selected = selected.loc[selected.river_idx.isin(river_ids)].copy()
    print(
        f"Loaded {len(centerlines)} rivers, {len(arcs)} arcs, and "
        f"{len(selected)} selected intervals",
        flush=True,
    )

    stations = assign_intervals_to_stations(centerlines, selected)
    if bbox_lonlat is not None:
        test_region = box(*bbox_lonlat)
        stations["station_inside_test_region"] = (
            stations.to_crs(ANALYSIS_GEOGRAPHIC_CRS)
            .geometry.intersects(test_region)
            .to_numpy()
        )
    else:
        stations["station_inside_test_region"] = True
    station_lonlat = stations.to_crs(ANALYSIS_GEOGRAPHIC_CRS)
    stations["station_inside_watershed"] = points_intersect_region(
        station_lonlat.geometry.x,
        station_lonlat.geometry.y,
        effective_watershed,
    )
    stations["station_inside_dredge_region"] = (
        stations.station_inside_test_region & stations.station_inside_watershed
    )
    vertices = expand_stations_to_arc_vertices(arcs, stations)
    if vertices.empty:
        raise ValueError("No RiverMapper arc vertices are available for dredging")
    if bbox_lonlat is not None:
        xmin, ymin, xmax, ymax = bbox_lonlat
        vertices["inside_test_region"] = (
            (vertices.x >= xmin)
            & (vertices.x <= xmax)
            & (vertices.y >= ymin)
            & (vertices.y <= ymax)
        )
    else:
        vertices["inside_test_region"] = True
    vertices["inside_watershed"] = points_intersect_region(
        vertices.x, vertices.y, effective_watershed
    )
    mesh_inside_watershed = points_intersect_region(
        hgrid.x, hgrid.y, effective_watershed
    )
    print(
        f"Built test set with {len(vertices)} arc vertices and "
        f"{len(hgrid.node_id)} hgrid nodes",
        flush=True,
    )
    mapped = map_vertices_to_mesh(
        vertices,
        hgrid,
        settings.query_workers,
        mesh_candidate_mask=mesh_inside_watershed,
    )
    result = apply_station_depths(mapped, hgrid, settings)
    if settings.intersection_recovery:
        result.intersection_candidates = screen_intersection_mesh_nodes(
            result.vertices,
            hgrid,
            effective_watershed,
            radius_m=settings.intersection_search_radius_m,
            width_tolerance_m=settings.intersection_width_tolerance_m,
            bank_exclusion_fraction=settings.intersection_bank_exclusion_fraction,
            bbox_lonlat=bbox_lonlat,
            query_workers=settings.query_workers,
        )
    else:
        result.intersection_candidates = _empty_intersection_targets()
    result.intersection_targets = result.intersection_candidates.loc[
        result.intersection_candidates.intersection_target
    ].copy()
    refresh_result_mesh_requests(
        result, hgrid, max_dredging_delta_m=settings.max_dredging_delta_m
    )

    candidates = result.intersection_candidates
    targets = result.intersection_targets
    result.summary.update(
        {
            "intersection_candidate_mesh_nodes_200m": len(candidates),
            "strict_half_width_intersection_target_mesh_nodes": int(
                candidates.passes_strict_half_width.sum()
            ) if not candidates.empty else 0,
            "tolerant_half_width_intersection_target_mesh_nodes": int(
                candidates.passes_half_width_with_tolerance.sum()
            ) if not candidates.empty else 0,
            "bank_proximity_rejected_intersection_mesh_nodes": int(
                (
                    candidates.passes_half_width_with_tolerance
                    & ~candidates.intersection_target
                ).sum()
            ) if not candidates.empty else 0,
            "intersection_target_mesh_nodes": len(targets),
            "additional_intersection_target_mesh_nodes": int(
                targets.additional_target.sum()
            ) if not targets.empty else 0,
            "intersection_targets_that_would_deepen": int(
                targets.would_deepen.sum()
            ) if not targets.empty else 0,
            "additional_intersection_targets_that_would_deepen": int(
                (targets.additional_target & targets.would_deepen).sum()
            ) if not targets.empty else 0,
        }
    )
    return result


def ensure_channel_connectivity(
    hgrid_obj,
    *,
    river_arcs_file: str | Path = DEFAULT_ARCS,
    river_centerlines_file: str | Path = DEFAULT_CENTERLINES,
    matches_gpkg_file: str | Path = DEFAULT_MATCHES,
    region_gdf_file_list: Sequence[str | Path] | None = None,
    exclude_region_gdf_file_list: Sequence[str | Path] | None = None,
    output_dir: str | Path,
    min_channel_depth: float = 1.0,
    channel_depth_source: str = "hydrofabric",
    measured_from_high_bank: bool = True,
    max_nearest_distance_m: float = 500.0,
    max_dredging_delta_m: float = 6.0,
    intersection_search_radius_m: float = 200.0,
    intersection_width_tolerance_m: float = 10.0,
    intersection_bank_exclusion_fraction: float = 0.05,
    intersection_recovery: bool = True,
    unmatched_policy: str = "baseline",
    query_workers: int = -1,
    write_gpkg: bool = False,
):
    """Apply the hydrofabric connectivity workflow to an in-memory SCHISM grid."""
    if output_dir is None:
        raise ValueError("output_dir is required")
    required_files = {
        "river_arcs_file": Path(river_arcs_file),
        "river_centerlines_file": Path(river_centerlines_file),
        "matches_gpkg_file": Path(matches_gpkg_file),
    }
    for label, path in required_files.items():
        if not path.is_file():
            raise FileNotFoundError(f"{label} does not exist: {path}")

    if region_gdf_file_list is None:
        region_paths = (DEFAULT_WATERSHED,)
    elif isinstance(region_gdf_file_list, (str, Path)):
        raise TypeError("region_gdf_file_list must be a sequence of paths")
    else:
        region_paths = tuple(Path(path) for path in region_gdf_file_list)
    if not region_paths:
        raise ValueError("region_gdf_file_list cannot be empty")
    for path in region_paths:
        if not path.is_file():
            raise FileNotFoundError(f"Region file does not exist: {path}")

    if exclude_region_gdf_file_list is None:
        exclude_paths = DEFAULT_EXCLUDE_REGIONS
    elif isinstance(exclude_region_gdf_file_list, (str, Path)):
        raise TypeError("exclude_region_gdf_file_list must be a sequence of paths")
    else:
        exclude_paths = tuple(Path(path) for path in exclude_region_gdf_file_list)
    for path in exclude_paths:
        if not path.is_file():
            raise FileNotFoundError(f"Exclusion file does not exist: {path}")

    x = np.asarray(hgrid_obj.x, dtype=float)
    y = np.asarray(hgrid_obj.y, dtype=float)
    dp = np.asarray(hgrid_obj.dp, dtype=float)
    if x.ndim != 1 or y.ndim != 1 or dp.ndim != 1:
        raise ValueError("hgrid x, y, and dp must be one-dimensional arrays")
    if not (len(x) == len(y) == len(dp)):
        raise ValueError("hgrid x, y, and dp must have equal lengths")
    if len(x) == 0:
        raise ValueError("hgrid contains no nodes")
    if not np.isfinite(x).all() or not np.isfinite(y).all():
        raise ValueError("hgrid coordinates must be finite")
    if not np.isfinite(dp).all():
        raise ValueError("hgrid depths must be finite")

    settings = DredgeSettings(
        min_channel_depth_m=min_channel_depth,
        channel_depth_source=channel_depth_source,
        measured_from_high_bank=measured_from_high_bank,
        max_nearest_distance_m=max_nearest_distance_m,
        max_dredging_delta_m=max_dredging_delta_m,
        intersection_search_radius_m=intersection_search_radius_m,
        intersection_width_tolerance_m=intersection_width_tolerance_m,
        intersection_bank_exclusion_fraction=(
            intersection_bank_exclusion_fraction
        ),
        intersection_recovery=intersection_recovery,
        unmatched_policy=unmatched_policy,
        query_workers=query_workers,
    )
    effective_watershed = load_effective_watershed(region_paths, exclude_paths)
    hgrid = HgridNodes(
        node_id=np.arange(1, len(x) + 1, dtype=np.int64),
        x=x.copy(),
        y=y.copy(),
        dp=dp.copy(),
        total_node_count=len(x),
    )
    result = run_hydrofabric_dredge(
        hgrid,
        required_files["river_arcs_file"],
        required_files["river_centerlines_file"],
        required_files["matches_gpkg_file"],
        effective_watershed,
        settings,
    )

    output_path = Path(output_dir)
    output_path.mkdir(parents=True, exist_ok=True)
    products = write_diagnostics(
        output_path,
        result,
        bbox_lonlat=None,
        write_gpkg=write_gpkg,
        effective_watershed=effective_watershed,
    )
    result.summary.update(
        {
            "bbox_lonlat": None,
            "inputs": {
                **{label: str(path) for label, path in required_files.items()},
                "regions": [str(path) for path in region_paths],
                "exclude_regions": [str(path) for path in exclude_paths],
            },
            "products": products,
        }
    )
    summary_path = output_path / "summary.json"
    summary_path.write_text(json.dumps(result.summary, indent=2) + "\n")

    dredged = copy.deepcopy(hgrid_obj)
    if not result.node_changes.empty:
        positions = result.node_changes.mesh_position.to_numpy(dtype=np.int64)
        dredged.dp[positions] = result.node_changes.final_dp.to_numpy(dtype=float)
    print(json.dumps(result.summary, indent=2), flush=True)
    return dredged


def get_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--hgrid", type=Path, default=DEFAULT_HGRID)
    parser.add_argument("--river-arcs", type=Path, default=DEFAULT_ARCS)
    parser.add_argument("--river-centerlines", type=Path, default=DEFAULT_CENTERLINES)
    parser.add_argument("--matches-gpkg", type=Path, default=DEFAULT_MATCHES)
    parser.add_argument(
        "--watershed",
        action="append",
        type=Path,
        default=None,
        help=(
            "Region included in dredging; repeat for multiple files. The "
            "default watershed is used when this option is omitted."
        ),
    )
    parser.add_argument(
        "--exclude-region",
        action="append",
        type=Path,
        default=None,
        help=(
            "Polygon subtracted from the watershed; repeat for multiple files. "
            "The legacy Maine exclusion is used when this option is omitted."
        ),
    )
    parser.add_argument("--output-dir", type=Path, default=DEFAULT_OUTPUT_DIR)
    parser.add_argument(
        "--bbox",
        nargs=4,
        type=float,
        default=None,
        metavar=("MIN_LON", "MIN_LAT", "MAX_LON", "MAX_LAT"),
        help="Optional local-test bounds. Omit only for a guarded full-domain run.",
    )
    parser.add_argument("--min-channel-depth-m", type=float, default=1.0)
    parser.add_argument(
        "--channel-depth-source",
        choices=("hydrofabric", "constant"),
        default="hydrofabric",
        help=(
            "Use matched hydrofabric bankfull depths with the configured "
            "minimum/fallback, or use the minimum depth as a constant for "
            "every station."
        ),
    )
    parser.add_argument("--max-nearest-distance-m", type=float, default=500.0)
    parser.add_argument(
        "--max-dredging-delta-m",
        type=float,
        default=6.0,
        help=(
            "Cap rather than reject a mesh-node request when it would deepen "
            "the existing bathymetry by more than this amount (default: 6 m)."
        ),
    )
    parser.add_argument("--intersection-search-radius-m", type=float, default=200.0)
    parser.add_argument("--intersection-width-tolerance-m", type=float, default=10.0)
    parser.add_argument(
        "--intersection-bank-exclusion-fraction", type=float, default=0.05
    )
    parser.add_argument(
        "--unmatched-policy", choices=("baseline", "skip"), default="baseline"
    )
    parser.add_argument("--measured-from-lower-bank", action="store_true")
    parser.add_argument(
        "--disable-intersection-recovery",
        action="store_true",
        help="Use forward RiverMapper-vertex requests only.",
    )
    parser.add_argument("--query-workers", type=int, default=-1)
    parser.add_argument("--write-hgrid", action="store_true")
    parser.add_argument("--write-gpkg", action="store_true")
    return parser


def dredge_settings_from_namespace(
    args: argparse.Namespace,
    *,
    query_workers: int | None = None,
) -> DredgeSettings:
    """Build settings identically for the serial and MPI command-line paths."""
    workers = args.query_workers if query_workers is None else query_workers
    return DredgeSettings(
        min_channel_depth_m=args.min_channel_depth_m,
        channel_depth_source=args.channel_depth_source,
        measured_from_high_bank=not args.measured_from_lower_bank,
        max_nearest_distance_m=args.max_nearest_distance_m,
        max_dredging_delta_m=args.max_dredging_delta_m,
        intersection_search_radius_m=args.intersection_search_radius_m,
        intersection_width_tolerance_m=args.intersection_width_tolerance_m,
        intersection_bank_exclusion_fraction=(
            args.intersection_bank_exclusion_fraction
        ),
        intersection_recovery=not args.disable_intersection_recovery,
        unmatched_policy=args.unmatched_policy,
        query_workers=workers,
    )


def main(argv: Sequence[str] | None = None) -> int:
    args = get_parser().parse_args(argv)
    bbox_lonlat = normalized_bbox(args.bbox) if args.bbox is not None else None
    settings = dredge_settings_from_namespace(args)
    watershed_regions = (
        tuple(args.watershed)
        if args.watershed is not None
        else (DEFAULT_WATERSHED,)
    )
    exclude_regions = (
        tuple(args.exclude_region)
        if args.exclude_region is not None
        else DEFAULT_EXCLUDE_REGIONS
    )
    effective_watershed = load_effective_watershed(
        watershed_regions, exclude_regions
    )
    hgrid = read_hgrid_nodes(args.hgrid, bbox_lonlat=bbox_lonlat)
    result = run_hydrofabric_dredge(
        hgrid,
        args.river_arcs,
        args.river_centerlines,
        args.matches_gpkg,
        effective_watershed,
        settings,
        bbox_lonlat=bbox_lonlat,
    )
    products = write_diagnostics(
        args.output_dir,
        result,
        bbox_lonlat,
        write_gpkg=args.write_gpkg,
        effective_watershed=effective_watershed,
    )
    if args.write_hgrid:
        output_hgrid = args.output_dir / "hgrid_hydrofabric_dredged.gr3"
        write_hgrid_with_updates(args.hgrid, output_hgrid, result.node_changes)
        products["hgrid"] = str(output_hgrid)
    result.summary.update(
        {
            "bbox_lonlat": list(bbox_lonlat) if bbox_lonlat is not None else None,
            "inputs": {
                "hgrid": str(args.hgrid),
                "river_arcs": str(args.river_arcs),
                "river_centerlines": str(args.river_centerlines),
                "matches_gpkg": str(args.matches_gpkg),
                "watershed": (
                    str(watershed_regions[0])
                    if len(watershed_regions) == 1
                    else None
                ),
                "regions": [str(path) for path in watershed_regions],
                "exclude_regions": [str(path) for path in exclude_regions],
            },
            "products": products,
        }
    )
    summary_path = args.output_dir / "summary.json"
    summary_path.write_text(json.dumps(result.summary, indent=2) + "\n")
    print(json.dumps(result.summary, indent=2), flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
