#!/usr/bin/env python3
"""Prototype RiverMapper-to-hydrofabric matching on a bounded region.

The matching is deliberately line based.  A hydrofabric reach first has to
overlap a contiguous portion of a RiverMapper centerline within the distance
and direction tolerances.  RiverMapper stations inherit COMIDs only after the
accepted line intervals have been resolved.

The input files used by this workflow are large.  Regional derivatives are
cached with source fingerprints and the requested bounds so repeated tuning
does not reread the national shapefile, full RiverMapper map, or full hgrid.
"""

from __future__ import annotations

import argparse
import json
import math
import os
import tempfile
from dataclasses import dataclass
from pathlib import Path
from typing import Iterator, Sequence

import geopandas as gpd
import numpy as np
import pandas as pd
from pyproj import CRS
from shapely.geometry import LineString, Point, box
from shapely.ops import substring


DEFAULT_HYDROFABRIC = Path(
    "/sciclone/schism10/feiye/CIROH/Channel_Geometry/Dataset/"
    "Bankfull_Meanflow_CONUS_Stream_Reaches/"
    "Bankfull_Meanflow_CONUS_Stream_Reaches.shp"
)
DEFAULT_RIVER_MAP = Path(
    "/sciclone/schism10/Hgrid_projects/STOFS3D-v7/v19_RiverMapper/Outputs/"
    "bora_v19.1.v19_ie_v18_3_nwm_clipped_in_cudem_missing_tiles_20-core/"
    "total_river_arcs_extra.map"
)
DEFAULT_HGRID = Path("/sciclone/schism10/feiye/STOFS3D-v7.4/RUN100j/hgrid.gr3")
DEFAULT_BBOX = (-81.37073, 31.77, -81.0104, 32.10649286)
DEFAULT_CACHE_DIR = Path(
    "/sciclone/schism10/feiye/CIROH/Channel_Geometry/Connectivity_test/savannah"
)
# The RiverMapper/hgrid coordinates are WGS 84 and the hydrofabric is NAD83.
# This environment cannot construct its WGS84<->NAD83 datum transformation
# (PROJ returns infinite coordinates).  For this regional matching operation,
# treat longitude/latitude input coordinates as NAD83 before projecting to
# CONUS Albers.  The resulting sub-metre datum approximation is negligible
# relative to the 500 m candidate radius and avoids silently invalid geometry.
INPUT_COORDINATE_CRS = "EPSG:4326"
ANALYSIS_GEOGRAPHIC_CRS = "EPSG:4269"
MATCH_CRS = "EPSG:5070"
HYDRO_COLUMNS = [
    "COMID",
    "bnk_depth",
    "bnk_width",
    "StreamOrde",
    "TotDASqKM",
    "geometry",
]

CANDIDATE_COLUMNS = [
    "river_idx",
    "hydro_feature_idx",
    "COMID",
    "bnk_depth",
    "bnk_width",
    "start_m",
    "end_m",
    "overlap_m",
    "river_coverage",
    "short_line_coverage",
    "mean_distance_m",
    "max_distance_m",
    "mean_angle_deg",
    "geometry",
]
SELECTED_COLUMNS = [
    "river_idx",
    "COMID",
    "bnk_depth",
    "bnk_width",
    "start_m",
    "end_m",
    "mean_distance_m",
    "mean_angle_deg",
    "length_m",
    "geometry",
]
STATION_COLUMNS = [
    "river_idx",
    "station_idx",
    "centerline_m",
    "COMID",
    "bnk_depth",
    "match_status",
    "match_method",
    "continuity_basis",
    "geometry",
]


@dataclass(frozen=True)
class MatchSettings:
    search_radius_m: float = 500.0
    sample_spacing_m: float = 25.0
    max_angle_deg: float = 75.0
    min_overlap_m: float = 25.0
    min_short_reach_fraction: float = 0.10
    review_max_angle_deg: float = 30.0
    review_min_overlap_m: float = 100.0
    review_max_distance_m: float = 250.0
    continuity_max_interval_gap_m: float = 25.0


def normalized_bbox(values: Sequence[float]) -> tuple[float, float, float, float]:
    """Return bounds ordered as min longitude, min latitude, max longitude, max latitude."""
    if len(values) != 4:
        raise ValueError("bbox must contain four values: min_lon min_lat max_lon max_lat")
    x1, y1, x2, y2 = (float(value) for value in values)
    return min(x1, x2), min(y1, y2), max(x1, x2), max(y1, y2)


def buffered_region(
    bbox_lonlat: tuple[float, float, float, float], buffer_m: float
) -> gpd.GeoDataFrame:
    """Create a metric buffer around a longitude/latitude bounding box."""
    region = gpd.GeoDataFrame(
        {"name": ["region"]},
        geometry=[box(*bbox_lonlat)],
        crs=ANALYSIS_GEOGRAPHIC_CRS,
    ).to_crs(MATCH_CRS)
    region.geometry = region.geometry.buffer(buffer_m)
    return region


def source_fingerprint(path: Path) -> dict[str, object]:
    resolved = path.resolve()
    stat = resolved.stat()
    return {
        "path": str(resolved),
        "size": stat.st_size,
        "mtime_ns": stat.st_mtime_ns,
    }


def shapefile_fingerprint(path: Path) -> dict[str, object]:
    """Fingerprint geometry, attributes, index, and projection sidecars."""
    return {
        suffix: source_fingerprint(path.with_suffix(suffix))
        for suffix in (".shp", ".shx", ".dbf", ".prj")
    }


def cache_is_current(
    metadata_path: Path,
    expected: dict[str, object],
    required_files: Sequence[Path],
) -> bool:
    if not metadata_path.is_file() or not all(path.is_file() for path in required_files):
        return False
    try:
        actual = json.loads(metadata_path.read_text())
    except (OSError, json.JSONDecodeError):
        return False
    return actual == expected


def write_metadata(path: Path, metadata: dict[str, object]) -> None:
    temporary = path.with_suffix(path.suffix + ".tmp")
    temporary.write_text(json.dumps(metadata, indent=2, sort_keys=True) + "\n")
    os.replace(temporary, path)


def _atomic_geoparquet(gdf: gpd.GeoDataFrame, path: Path) -> None:
    temporary = path.with_suffix(".tmp.parquet")
    gdf.to_parquet(temporary, index=False, compression="zstd")
    os.replace(temporary, path)


def records_geodataframe(
    records: list[dict[str, object]], columns: Sequence[str], crs: str
) -> gpd.GeoDataFrame:
    """Construct a GeoDataFrame with a stable schema even when it is empty."""
    if records:
        return gpd.GeoDataFrame(records, geometry="geometry", crs=crs)
    return gpd.GeoDataFrame(
        {column: pd.Series(dtype="object") for column in columns},
        geometry="geometry",
        crs=crs,
    )


def build_hydrofabric_cache(
    source: Path,
    cache_dir: Path,
    bbox_lonlat: tuple[float, float, float, float],
    buffer_m: float,
    force: bool = False,
) -> Path:
    """Cache only hydrofabric reaches near the requested region in a metric CRS."""
    cache_path = cache_dir / "hydrofabric.parquet"
    metadata_path = cache_dir / "hydrofabric.metadata.json"
    metadata = {
        "kind": "hydrofabric",
        "source": shapefile_fingerprint(source),
        "bbox_lonlat": list(bbox_lonlat),
        "buffer_m": buffer_m,
        "cache_crs": MATCH_CRS,
        "columns": HYDRO_COLUMNS,
    }
    if not force and cache_is_current(metadata_path, metadata, [cache_path]):
        print(f"Reusing hydrofabric cache: {cache_path}", flush=True)
        return cache_path

    source_crs = CRS.from_user_input(gpd.read_file(source, rows=1).crs)
    region_metric = buffered_region(bbox_lonlat, buffer_m)
    source_bounds = tuple(region_metric.to_crs(source_crs).total_bounds)
    print(f"Reading hydrofabric features in source-CRS bounds {source_bounds}", flush=True)
    hydro = gpd.read_file(source, bbox=source_bounds)
    missing = sorted(set(HYDRO_COLUMNS) - set(hydro.columns))
    if missing:
        raise ValueError(f"Hydrofabric is missing required columns: {missing}")
    hydro = hydro[HYDRO_COLUMNS].copy().to_crs(MATCH_CRS)
    region_shape = region_metric.geometry.iloc[0]
    hydro = hydro.loc[hydro.geometry.intersects(region_shape)].reset_index(drop=True)
    hydro.insert(0, "hydro_feature_idx", np.arange(len(hydro), dtype=np.int64))
    _atomic_geoparquet(hydro, cache_path)
    write_metadata(metadata_path, metadata)
    print(f"Cached {len(hydro)} hydrofabric reaches: {cache_path}", flush=True)
    return cache_path


def iter_sms_arcs_fast(path: Path) -> Iterator[np.ndarray]:
    """Read SMS arcs without SMS_MAP's repeated NumPy-array appends.

    The standard SMS_MAP reader repeatedly appends global nodes to an ndarray,
    which is prohibitively slow for the 150 MB diagnostic map used here.
    """
    nodes: list[tuple[float, float, float]] = []
    with path.open("r", encoding="utf-8") as stream:
        while line := stream.readline():
            keyword = line.strip().split(" ", 1)[0]
            if keyword == "NODE":
                xyz = stream.readline().split()
                if not xyz or xyz[0] != "XY":
                    raise ValueError("Expected XY immediately after NODE in SMS map")
                nodes.append((float(xyz[1]), float(xyz[2]), float(xyz[3])))
                continue
            if keyword != "ARC":
                continue

            endpoint_indices: tuple[int, int] | None = None
            while arc_line := stream.readline():
                fields = arc_line.split()
                if not fields:
                    continue
                if fields[0] == "NODES":
                    endpoint_indices = (int(fields[1]) - 1, int(fields[2]) - 1)
                elif fields[0] == "ARCVERTICES":
                    if endpoint_indices is None:
                        raise ValueError("ARCVERTICES encountered before NODES")
                    count = int(fields[1])
                    vertices = np.empty((count, 3), dtype=float)
                    for index in range(count):
                        vertices[index] = np.asarray(stream.readline().split()[:3], dtype=float)
                    start, end = endpoint_indices
                    yield np.vstack((nodes[start], vertices, nodes[end]))
                    break


def decode_nrow(encoded_z: float) -> int:
    """Decode RiverMapper's first two-digit field (the number of parallel arcs)."""
    field_count = int(encoded_z)
    if field_count <= 0 or field_count > 6:
        raise ValueError(f"Invalid RiverMapper encoded z value: {encoded_z}")
    scaled = int(round(encoded_z * 100**field_count))
    digits = str(scaled)[1:]
    if len(digits) < 2:
        raise ValueError(f"Cannot decode the arc count from z={encoded_z}")
    return int(digits[:2])


def build_river_cache(
    source: Path,
    cache_dir: Path,
    bbox_lonlat: tuple[float, float, float, float] | None,
    buffer_m: float,
    force: bool = False,
) -> tuple[Path, Path]:
    """Cache RiverMapper arcs and their assembled centerlines.

    If ``bbox_lonlat`` is ``None``, retain the complete RiverMapper map. This
    mode is used by the MPI driver so the expensive SMS parsing is performed
    once and reused by subsequent full-domain runs.
    """
    arcs_path = cache_dir / "river_arcs.parquet"
    centerlines_path = cache_dir / "river_centerlines.parquet"
    metadata_path = cache_dir / "river_map.metadata.json"
    metadata = {
        "kind": "river_map",
        "source": source_fingerprint(source),
        "bbox_lonlat": list(bbox_lonlat) if bbox_lonlat is not None else None,
        "buffer_m": buffer_m if bbox_lonlat is not None else None,
        "cache_crs": ANALYSIS_GEOGRAPHIC_CRS,
    }
    if not force and cache_is_current(
        metadata_path, metadata, [arcs_path, centerlines_path]
    ):
        print(f"Reusing RiverMapper cache: {centerlines_path}", flush=True)
        return arcs_path, centerlines_path

    selection_shape = None
    if bbox_lonlat is not None:
        selection_shape = (
            buffered_region(bbox_lonlat, buffer_m)
            .to_crs(ANALYSIS_GEOGRAPHIC_CRS)
            .geometry.iloc[0]
        )
    arc_records: list[dict[str, object]] = []
    centerline_records: list[dict[str, object]] = []
    arc_iterator = iter(iter_sms_arcs_fast(source))
    river_idx = 0
    while True:
        try:
            first_arc = next(arc_iterator)
        except StopIteration:
            break
        n_arcs = decode_nrow(float(first_arc[0, 2]))
        group = [first_arc]
        try:
            group.extend(next(arc_iterator) for _ in range(n_arcs - 1))
        except StopIteration as exc:
            raise ValueError("Incomplete RiverMapper arc group at end of SMS map") from exc
        station_counts = {len(arc) for arc in group}
        if len(station_counts) != 1:
            raise ValueError(
                f"River {river_idx} has unaligned arc sizes: {sorted(station_counts)}"
            )
        center_xy = np.mean(np.stack([arc[:, :2] for arc in group]), axis=0)
        centerline = LineString(center_xy)
        if centerline.is_empty or centerline.length == 0 or (
            selection_shape is not None and not centerline.intersects(selection_shape)
        ):
            river_idx += 1
            continue
        for local_arc_idx, arc in enumerate(group):
            arc_records.append(
                {
                    "river_idx": river_idx,
                    "local_arc_idx": local_arc_idx,
                    "n_arcs": n_arcs,
                    "n_stations": len(arc),
                    "geometry": LineString(arc),
                }
            )
        centerline_records.append(
            {
                "river_idx": river_idx,
                "n_arcs": n_arcs,
                "n_stations": len(center_xy),
                "geometry": centerline,
            }
        )
        river_idx += 1

    arcs = gpd.GeoDataFrame(
        arc_records, geometry="geometry", crs=ANALYSIS_GEOGRAPHIC_CRS
    )
    centerlines = gpd.GeoDataFrame(
        centerline_records, geometry="geometry", crs=ANALYSIS_GEOGRAPHIC_CRS
    )
    _atomic_geoparquet(arcs, arcs_path)
    _atomic_geoparquet(centerlines, centerlines_path)
    write_metadata(metadata_path, metadata)
    print(
        f"Cached {len(centerlines)} river centerlines and {len(arcs)} arcs",
        flush=True,
    )
    return arcs_path, centerlines_path


def build_hgrid_cache(
    source: Path,
    cache_dir: Path,
    bbox_lonlat: tuple[float, float, float, float],
    buffer_m: float,
    force: bool = False,
) -> Path:
    """Cache hgrid nodes near the region without loading full element connectivity."""
    cache_path = cache_dir / "hgrid_nodes.npz"
    metadata_path = cache_dir / "hgrid.metadata.json"
    metadata = {
        "kind": "hgrid_nodes",
        "source": source_fingerprint(source),
        "bbox_lonlat": list(bbox_lonlat),
        "buffer_m": buffer_m,
        "crs": INPUT_COORDINATE_CRS,
    }
    if not force and cache_is_current(metadata_path, metadata, [cache_path]):
        print(f"Reusing hgrid-node cache: {cache_path}", flush=True)
        return cache_path

    selection_bounds = (
        buffered_region(bbox_lonlat, buffer_m)
        .to_crs(ANALYSIS_GEOGRAPHIC_CRS)
        .total_bounds
    )
    xmin, ymin, xmax, ymax = selection_bounds
    node_ids: list[int] = []
    x_values: list[float] = []
    y_values: list[float] = []
    depths: list[float] = []
    with source.open("r", encoding="utf-8") as stream:
        stream.readline()
        counts = stream.readline().split()
        if len(counts) < 2:
            raise ValueError(f"Invalid hgrid header in {source}")
        n_nodes = int(counts[1])
        for _ in range(n_nodes):
            fields = stream.readline().split()
            node_id, x, y, depth = int(fields[0]), *map(float, fields[1:4])
            if xmin <= x <= xmax and ymin <= y <= ymax:
                node_ids.append(node_id)
                x_values.append(x)
                y_values.append(y)
                depths.append(depth)
    temporary = cache_path.with_suffix(".tmp.npz")
    np.savez_compressed(
        temporary,
        node_id=np.asarray(node_ids, dtype=np.int64),
        x=np.asarray(x_values),
        y=np.asarray(y_values),
        dp=np.asarray(depths),
    )
    os.replace(temporary, cache_path)
    write_metadata(metadata_path, metadata)
    print(f"Cached {len(node_ids)} hgrid nodes: {cache_path}", flush=True)
    return cache_path


def sample_measures(length: float, spacing: float) -> np.ndarray:
    if length <= 0:
        return np.array([0.0])
    count = max(1, int(math.ceil(length / spacing)))
    return np.linspace(0.0, length, count + 1)


def line_tangent(line: LineString, measure: float, window: float) -> np.ndarray:
    start = max(0.0, measure - window)
    end = min(line.length, measure + window)
    if end <= start:
        return np.array([0.0, 0.0])
    a = line.interpolate(start)
    b = line.interpolate(end)
    vector = np.array([b.x - a.x, b.y - a.y])
    norm = np.linalg.norm(vector)
    return vector / norm if norm > 0 else vector


def local_angle_degrees(
    river: LineString,
    river_measure: float,
    hydro: LineString,
    point: Point,
    window: float,
) -> float:
    river_vector = line_tangent(river, river_measure, window)
    hydro_vector = line_tangent(hydro, hydro.project(point), window)
    if not river_vector.any() or not hydro_vector.any():
        return 90.0
    cosine = float(np.clip(abs(np.dot(river_vector, hydro_vector)), 0.0, 1.0))
    return math.degrees(math.acos(cosine))


def contiguous_true_runs(mask: np.ndarray) -> list[tuple[int, int]]:
    padded = np.r_[False, mask, False].astype(np.int8)
    transitions = np.diff(padded)
    starts = np.flatnonzero(transitions == 1)
    ends = np.flatnonzero(transitions == -1) - 1
    return list(zip(starts.tolist(), ends.tolist()))


def candidate_overlap_intervals(
    river_idx: int,
    river: LineString,
    hydro_row: pd.Series,
    settings: MatchSettings,
) -> list[dict[str, object]]:
    """Find contiguous, direction-compatible partial overlaps for one line pair."""
    hydro = hydro_row.geometry
    measures = sample_measures(river.length, settings.sample_spacing_m)
    points = [river.interpolate(measure) for measure in measures]
    distances = np.asarray([point.distance(hydro) for point in points])
    tangent_window = max(settings.sample_spacing_m, 10.0)
    angles = np.asarray(
        [
            local_angle_degrees(river, measure, hydro, point, tangent_window)
            for measure, point in zip(measures, points)
        ]
    )
    accepted = (distances <= settings.search_radius_m) & (
        angles <= settings.max_angle_deg
    )
    records: list[dict[str, object]] = []
    for start_idx, end_idx in contiguous_true_runs(accepted):
        start_measure = max(0.0, measures[start_idx] - settings.sample_spacing_m / 2)
        end_measure = min(
            river.length, measures[end_idx] + settings.sample_spacing_m / 2
        )
        overlap_length = end_measure - start_measure
        short_reference = min(river.length, hydro.length)
        short_fraction = overlap_length / short_reference if short_reference > 0 else 0.0
        if (
            overlap_length < settings.min_overlap_m
            and short_fraction < settings.min_short_reach_fraction
        ):
            continue
        interval_slice = slice(start_idx, end_idx + 1)
        records.append(
            {
                "river_idx": river_idx,
                "hydro_feature_idx": int(hydro_row.hydro_feature_idx),
                "COMID": int(hydro_row.COMID),
                "bnk_depth": float(hydro_row.bnk_depth),
                "bnk_width": float(hydro_row.bnk_width),
                "start_m": start_measure,
                "end_m": end_measure,
                "overlap_m": overlap_length,
                "river_coverage": overlap_length / river.length,
                "short_line_coverage": short_fraction,
                "mean_distance_m": float(np.mean(distances[interval_slice])),
                "max_distance_m": float(np.max(distances[interval_slice])),
                "mean_angle_deg": float(np.mean(angles[interval_slice])),
                "geometry": substring(river, start_measure, end_measure),
            }
        )
    return records


def resolve_candidate_intervals(
    river_idx: int,
    river: LineString,
    candidates: pd.DataFrame,
    hydro_by_feature: pd.DataFrame,
    settings: MatchSettings,
) -> list[dict[str, object]]:
    """Resolve overlapping candidate intervals using local distance and direction."""
    if candidates.empty:
        return []
    boundaries = sorted(
        set([0.0, river.length, *candidates.start_m.tolist(), *candidates.end_m.tolist()])
    )
    cells: list[dict[str, object]] = []
    tangent_window = max(settings.sample_spacing_m, 10.0)
    for start_measure, end_measure in zip(boundaries[:-1], boundaries[1:]):
        if end_measure - start_measure <= 1.0e-6:
            continue
        midpoint = 0.5 * (start_measure + end_measure)
        eligible = candidates.loc[
            (candidates.start_m <= midpoint) & (candidates.end_m >= midpoint)
        ]
        if eligible.empty:
            continue
        point = river.interpolate(midpoint)
        best: tuple[float, pd.Series, float, float] | None = None
        for _, candidate in eligible.iterrows():
            hydro = hydro_by_feature.loc[int(candidate.hydro_feature_idx)].geometry
            distance = point.distance(hydro)
            angle = local_angle_degrees(
                river, midpoint, hydro, point, tangent_window
            )
            score = distance + settings.search_radius_m * angle / settings.max_angle_deg
            if best is None or score < best[0]:
                best = (score, candidate, distance, angle)
        assert best is not None
        _, candidate, distance, angle = best
        cells.append(
            {
                "river_idx": river_idx,
                "COMID": int(candidate.COMID),
                "bnk_depth": float(candidate.bnk_depth),
                "bnk_width": float(candidate.bnk_width),
                "start_m": start_measure,
                "end_m": end_measure,
                "mean_distance_m": distance,
                "mean_angle_deg": angle,
            }
        )

    merged: list[dict[str, object]] = []
    for cell in cells:
        if (
            merged
            and merged[-1]["COMID"] == cell["COMID"]
            and math.isclose(float(merged[-1]["end_m"]), float(cell["start_m"]))
        ):
            old_length = float(merged[-1]["end_m"]) - float(merged[-1]["start_m"])
            new_length = float(cell["end_m"]) - float(cell["start_m"])
            total_length = old_length + new_length
            for field in ("mean_distance_m", "mean_angle_deg"):
                merged[-1][field] = (
                    float(merged[-1][field]) * old_length
                    + float(cell[field]) * new_length
                ) / total_length
            merged[-1]["end_m"] = cell["end_m"]
        else:
            merged.append(cell.copy())
    for record in merged:
        record["length_m"] = float(record["end_m"]) - float(record["start_m"])
        record["geometry"] = substring(
            river, float(record["start_m"]), float(record["end_m"])
        )
    return merged


def match_centerlines(
    rivers: gpd.GeoDataFrame,
    hydro: gpd.GeoDataFrame,
    settings: MatchSettings,
) -> tuple[gpd.GeoDataFrame, gpd.GeoDataFrame, gpd.GeoDataFrame]:
    """Match regional RiverMapper centerlines to partial hydrofabric intervals."""
    rivers_metric = rivers.to_crs(MATCH_CRS).reset_index(drop=True)
    hydro_metric = hydro.to_crs(MATCH_CRS).reset_index(drop=True)
    hydro_by_feature = hydro_metric.set_index("hydro_feature_idx", drop=False)
    candidate_records: list[dict[str, object]] = []
    selected_records: list[dict[str, object]] = []
    station_records: list[dict[str, object]] = []

    for _, river_row in rivers_metric.iterrows():
        river_idx = int(river_row.river_idx)
        river = river_row.geometry
        possible_positions = hydro_metric.sindex.query(
            river.buffer(settings.search_radius_m), predicate="intersects"
        )
        local_candidates: list[dict[str, object]] = []
        for position in possible_positions:
            records = candidate_overlap_intervals(
                river_idx, river, hydro_metric.iloc[int(position)], settings
            )
            local_candidates.extend(records)
            candidate_records.extend(records)
        candidate_frame = pd.DataFrame(local_candidates)
        selected = resolve_candidate_intervals(
            river_idx, river, candidate_frame, hydro_by_feature, settings
        )
        selected_records.extend(selected)

        coordinates = list(river.coords)
        cumulative = np.r_[
            0.0,
            np.cumsum(
                [
                    math.hypot(b[0] - a[0], b[1] - a[1])
                    for a, b in zip(coordinates[:-1], coordinates[1:])
                ]
            ),
        ]
        for station_idx, (coordinate, measure) in enumerate(zip(coordinates, cumulative)):
            containing = [
                record
                for record in selected
                if float(record["start_m"]) - 1.0e-6
                <= measure
                <= float(record["end_m"]) + 1.0e-6
            ]
            match = containing[0] if containing else None
            station_records.append(
                {
                    "river_idx": river_idx,
                    "station_idx": station_idx,
                    "centerline_m": measure,
                    "COMID": int(match["COMID"]) if match else None,
                    "bnk_depth": float(match["bnk_depth"]) if match else np.nan,
                    "match_status": "matched" if match else "unmatched",
                    "match_method": "line_overlap" if match else None,
                    "continuity_basis": None,
                    "geometry": Point(coordinate),
                }
            )

    candidates = records_geodataframe(
        candidate_records, CANDIDATE_COLUMNS, MATCH_CRS
    )
    selected = records_geodataframe(selected_records, SELECTED_COLUMNS, MATCH_CRS)
    stations = records_geodataframe(station_records, STATION_COLUMNS, MATCH_CRS)
    return candidates, selected, stations


def write_diagnostics(
    output_path: Path,
    hydro: gpd.GeoDataFrame,
    rivers: gpd.GeoDataFrame,
    candidates: gpd.GeoDataFrame,
    selected: gpd.GeoDataFrame,
    stations: gpd.GeoDataFrame,
) -> None:
    if output_path.exists():
        output_path.unlink()
    flagged = selected.loc[selected.quality_flag != "nominal"].copy()
    flagged["diagnostic_group"] = flagged.quality_flag
    unmatched = stations.loc[
        (stations.match_status == "unmatched") & stations.inside_test_region
    ].copy()
    continuity_intervals = selected.loc[
        selected.quality_flag == "continuity_fill"
    ].copy()
    continuity_stations = stations.loc[
        stations.match_method == "continuity_fill"
    ].copy()
    layers = {
        "hydrofabric": hydro,
        "river_centerlines": rivers,
        "candidate_overlaps": candidates,
        "selected_intervals": selected,
        "station_matches": stations,
        "flagged_for_review": flagged,
        "unmatched_stations": unmatched,
        "unmatched_no_valid_reach": unmatched.loc[
            unmatched.diagnostic_group == "no_valid_reach_within_radius"
        ],
        "unmatched_direction_mismatch": unmatched.loc[
            unmatched.diagnostic_group == "direction_mismatch"
        ],
        "unmatched_insufficient_overlap": unmatched.loc[
            unmatched.diagnostic_group == "insufficient_contiguous_overlap"
        ],
        "continuity_fill_intervals": continuity_intervals,
        "continuity_filled_stations": continuity_stations,
        "continuity_by_station_neighbors": continuity_stations.loc[
            continuity_stations.continuity_basis == "station_neighbors"
        ],
        "continuity_by_interval_neighbors": continuity_stations.loc[
            continuity_stations.continuity_basis == "interval_neighbors"
        ],
    }
    for layer, frame in layers.items():
        if frame.empty:
            continue
        frame.to_crs(ANALYSIS_GEOGRAPHIC_CRS).to_file(
            output_path, layer=layer, driver="GPKG", mode="w"
        )


def classify_unmatched_stations(
    stations: gpd.GeoDataFrame,
    rivers: gpd.GeoDataFrame,
    hydro: gpd.GeoDataFrame,
    settings: MatchSettings,
) -> None:
    """Classify unmatched stations by coverage, direction, or overlap failure."""
    reasons = np.full(len(stations), "matched", dtype=object)
    unmatched_positions = np.flatnonzero(
        stations.match_status.to_numpy() == "unmatched"
    )
    rivers_metric = rivers.to_crs(MATCH_CRS).set_index("river_idx")
    hydro_metric = hydro.to_crs(MATCH_CRS).reset_index(drop=True)
    hydro_index = hydro_metric.sindex
    tangent_window = max(settings.sample_spacing_m, 10.0)
    for position in unmatched_positions:
        point = stations.geometry.iloc[position]
        nearby = hydro_index.query(
            point.buffer(settings.search_radius_m), predicate="intersects"
        )
        if not len(nearby):
            reasons[position] = "no_valid_reach_within_radius"
            continue

        station = stations.iloc[position]
        river = rivers_metric.loc[int(station.river_idx)].geometry
        direction_compatible = False
        for hydro_position in nearby:
            angle = local_angle_degrees(
                river,
                float(station.centerline_m),
                hydro_metric.geometry.iloc[int(hydro_position)],
                point,
                tangent_window,
            )
            if angle <= settings.max_angle_deg:
                direction_compatible = True
                break
        reasons[position] = (
            "insufficient_contiguous_overlap"
            if direction_compatible
            else "direction_mismatch"
        )
    stations["diagnostic_group"] = reasons


def fill_same_comid_gaps(
    stations: gpd.GeoDataFrame,
    rivers: gpd.GeoDataFrame,
    hydro: gpd.GeoDataFrame,
    selected: gpd.GeoDataFrame,
    settings: MatchSettings,
) -> gpd.GeoDataFrame:
    """Fill same-COMID gaps using station and selected-interval continuity."""
    rivers_metric = rivers.to_crs(MATCH_CRS).set_index("river_idx")
    hydro_by_comid = hydro.drop_duplicates("COMID").set_index("COMID")
    interval_records: list[dict[str, object]] = []

    for river_idx, group in stations.groupby("river_idx", sort=False):
        group = group.sort_values("station_idx")
        row_indices = group.index.to_numpy()
        fillable = (
            (group.match_status.to_numpy() == "unmatched")
            & group.inside_test_region.to_numpy(dtype=bool)
        )
        for start_position, end_position in contiguous_true_runs(fillable):
            if start_position == 0 or end_position == len(group) - 1:
                continue
            left = group.iloc[start_position - 1]
            right = group.iloc[end_position + 1]
            if (
                left.match_status != "matched"
                or right.match_status != "matched"
                or not bool(left.inside_test_region)
                or not bool(right.inside_test_region)
                or pd.isna(left.COMID)
                or pd.isna(right.COMID)
                or int(left.COMID) != int(right.COMID)
            ):
                continue

            comid = int(left.COMID)
            run_indices = row_indices[start_position : end_position + 1]
            stations.loc[run_indices, "COMID"] = comid
            stations.loc[run_indices, "bnk_depth"] = float(left.bnk_depth)
            stations.loc[run_indices, "match_status"] = "matched"
            stations.loc[run_indices, "match_method"] = "continuity_fill"
            stations.loc[run_indices, "continuity_basis"] = "station_neighbors"
            stations.loc[run_indices, "diagnostic_group"] = "continuity_fill"

            first = group.iloc[start_position]
            last = group.iloc[end_position]
            start_measure = 0.5 * (
                float(left.centerline_m) + float(first.centerline_m)
            )
            end_measure = 0.5 * (
                float(last.centerline_m) + float(right.centerline_m)
            )
            river = rivers_metric.loc[int(river_idx)].geometry
            hydro_row = hydro_by_comid.loc[comid]
            interval_records.append(
                {
                    "river_idx": int(river_idx),
                    "COMID": comid,
                    "bnk_depth": float(left.bnk_depth),
                    "bnk_width": float(hydro_row.bnk_width),
                    "start_m": start_measure,
                    "end_m": end_measure,
                    "mean_distance_m": np.nan,
                    "mean_angle_deg": np.nan,
                    "length_m": end_measure - start_measure,
                    "quality_flag": "continuity_fill",
                    "continuity_basis": "station_neighbors",
                    "geometry": substring(river, start_measure, end_measure),
                }
            )

    # A short failed sample can split one COMID into two selected intervals even
    # when no RiverMapper station falls in the left interval.  Bridge such a gap
    # when the immediately adjacent intervals agree on COMID.
    for river_idx, river_intervals in selected.groupby("river_idx", sort=False):
        river_intervals = river_intervals.sort_values(["start_m", "end_m"])
        records = list(river_intervals.itertuples(index=False))
        if len(records) < 2:
            continue
        river = rivers_metric.loc[int(river_idx)].geometry
        river_station_mask = stations.river_idx == river_idx
        for left, right in zip(records[:-1], records[1:]):
            gap_start = float(left.end_m)
            gap_end = float(right.start_m)
            gap_length = gap_end - gap_start
            if (
                gap_length <= 1.0e-6
                or gap_length > settings.continuity_max_interval_gap_m
                or int(left.COMID) != int(right.COMID)
            ):
                continue
            gap_station_mask = (
                river_station_mask
                & (stations.match_status == "unmatched")
                & stations.inside_test_region
                & (stations.centerline_m >= gap_start - 1.0e-6)
                & (stations.centerline_m <= gap_end + 1.0e-6)
            )
            gap_station_indices = stations.index[gap_station_mask]
            if len(gap_station_indices) == 0:
                continue

            comid = int(left.COMID)
            hydro_row = hydro_by_comid.loc[comid]
            stations.loc[gap_station_indices, "COMID"] = comid
            stations.loc[gap_station_indices, "bnk_depth"] = float(left.bnk_depth)
            stations.loc[gap_station_indices, "match_status"] = "matched"
            stations.loc[gap_station_indices, "match_method"] = "continuity_fill"
            stations.loc[
                gap_station_indices, "continuity_basis"
            ] = "interval_neighbors"
            stations.loc[gap_station_indices, "diagnostic_group"] = "continuity_fill"
            interval_records.append(
                {
                    "river_idx": int(river_idx),
                    "COMID": comid,
                    "bnk_depth": float(left.bnk_depth),
                    "bnk_width": float(hydro_row.bnk_width),
                    "start_m": gap_start,
                    "end_m": gap_end,
                    "mean_distance_m": np.nan,
                    "mean_angle_deg": np.nan,
                    "length_m": gap_length,
                    "quality_flag": "continuity_fill",
                    "continuity_basis": "interval_neighbors",
                    "geometry": substring(river, gap_start, gap_end),
                }
            )

    if not interval_records:
        return gpd.GeoDataFrame(
            columns=[
                "river_idx",
                "COMID",
                "bnk_depth",
                "bnk_width",
                "start_m",
                "end_m",
                "mean_distance_m",
                "mean_angle_deg",
                "length_m",
                "quality_flag",
                "continuity_basis",
                "geometry",
            ],
            geometry="geometry",
            crs=MATCH_CRS,
        )
    return gpd.GeoDataFrame(interval_records, geometry="geometry", crs=MATCH_CRS)


def run_matching_pipeline(
    rivers: gpd.GeoDataFrame,
    hydro_for_matching: gpd.GeoDataFrame,
    settings: MatchSettings,
    region_geometry_metric: object | None = None,
) -> tuple[
    gpd.GeoDataFrame,
    gpd.GeoDataFrame,
    gpd.GeoDataFrame,
    gpd.GeoDataFrame,
]:
    """Run matching, diagnostics, and continuity recovery for complete rivers.

    ``region_geometry_metric`` limits which stations are reported as being in
    the test region. Passing ``None`` marks all supplied stations as in-domain,
    which is appropriate for MPI work units containing complete river groups.
    """
    candidates, selected, stations = match_centerlines(
        rivers, hydro_for_matching, settings
    )
    classify_unmatched_stations(stations, rivers, hydro_for_matching, settings)
    if region_geometry_metric is None:
        stations["inside_test_region"] = True
    else:
        station_metric = stations.to_crs(MATCH_CRS)
        stations["inside_test_region"] = station_metric.geometry.intersects(
            region_geometry_metric
        ).to_numpy()

    if selected.empty:
        selected["quality_flag"] = pd.Series(dtype="object")
    else:
        selected["quality_flag"] = np.select(
            [
                selected.length_m < settings.review_min_overlap_m,
                selected.mean_distance_m > settings.review_max_distance_m,
                selected.mean_angle_deg > settings.review_max_angle_deg,
            ],
            ["short_interval", "large_distance", "large_angle"],
            default="nominal",
        )
    selected["continuity_basis"] = None
    continuity_intervals = fill_same_comid_gaps(
        stations, rivers, hydro_for_matching, selected, settings
    )
    if not continuity_intervals.empty:
        selected = gpd.GeoDataFrame(
            pd.concat([selected, continuity_intervals], ignore_index=True),
            geometry="geometry",
            crs=MATCH_CRS,
        )
    selected["diagnostic_group"] = selected.quality_flag
    return candidates, selected, stations, continuity_intervals


def write_overview_png(
    output_path: Path,
    bbox_lonlat: tuple[float, float, float, float],
    hydro: gpd.GeoDataFrame,
    rivers: gpd.GeoDataFrame,
    selected: gpd.GeoDataFrame,
    stations: gpd.GeoDataFrame,
    search_radius_m: float,
) -> None:
    """Plot nominal matches, individual review flags, and unmatched groups."""
    matplotlib_cache = Path(tempfile.gettempdir()) / "stofs3d-matplotlib-cache"
    matplotlib_cache.mkdir(parents=True, exist_ok=True)
    os.environ.setdefault("MPLCONFIGDIR", str(matplotlib_cache))
    import matplotlib.pyplot as plt
    from matplotlib.lines import Line2D

    hydro_plot = hydro.to_crs(ANALYSIS_GEOGRAPHIC_CRS)
    rivers_plot = rivers.to_crs(ANALYSIS_GEOGRAPHIC_CRS)
    selected_plot = selected.to_crs(ANALYSIS_GEOGRAPHIC_CRS)
    stations_plot = stations.to_crs(ANALYSIS_GEOGRAPHIC_CRS)

    figure, axis = plt.subplots(figsize=(12, 11), dpi=160)
    hydro_plot.plot(ax=axis, color="#b7b7b7", linewidth=0.45, zorder=1)
    rivers_plot.plot(ax=axis, color="black", linewidth=0.65, zorder=2)

    line_groups = {
        "nominal": ("#15803d", "Nominal selected interval"),
        "short_interval": ("#f28e2b", "Review: short interval"),
        "large_distance": ("#9467bd", "Review: large distance"),
        "large_angle": ("#17becf", "Review: large angle"),
        "continuity_fill": ("#8c564b", "Retained: same-COMID continuity fill"),
    }
    legend_handles = [
        Line2D([0], [0], color="#b7b7b7", label="Hydrofabric"),
        Line2D([0], [0], color="black", label="RiverMapper centerline"),
    ]
    for group, (color, label) in line_groups.items():
        subset = selected_plot.loc[selected_plot.quality_flag == group]
        if subset.empty:
            continue
        subset.plot(ax=axis, color=color, linewidth=1.7, zorder=3)
        legend_handles.append(
            Line2D([0], [0], color=color, label=f"{label} ({len(subset)})")
        )

    point_groups = {
        "no_valid_reach_within_radius": (
            "#d62728",
            "o",
            f"Unmatched: no valid-depth reach within {search_radius_m:g} m",
        ),
        "direction_mismatch": (
            "#e377c2",
            "x",
            "Unmatched: direction mismatch",
        ),
        "insufficient_contiguous_overlap": (
            "#ff7f0e",
            "+",
            "Unmatched: insufficient contiguous overlap",
        ),
    }
    regional_stations = stations_plot.loc[stations_plot.inside_test_region]
    for group, (color, marker, label) in point_groups.items():
        subset = regional_stations.loc[
            regional_stations.diagnostic_group == group
        ]
        if subset.empty:
            continue
        subset.plot(
            ax=axis,
            color=color,
            marker=marker,
            markersize=6 if marker == "o" else 12,
            linewidth=0.7,
            zorder=5,
        )
        legend_handles.append(
            Line2D(
                [0],
                [0],
                color=color,
                marker=marker,
                linestyle="",
                label=f"{label} ({len(subset)})",
            )
        )

    xmin, ymin, xmax, ymax = bbox_lonlat
    axis.set_xlim(xmin, xmax)
    axis.set_ylim(ymin, ymax)
    axis.set_aspect("equal")
    axis.set_title("RiverMapper–hydrofabric partial-overlap diagnostics")
    axis.legend(handles=legend_handles, loc="upper left", fontsize=8)
    figure.tight_layout()
    figure.savefig(output_path)
    plt.close(figure)


def get_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--hydrofabric", type=Path, default=DEFAULT_HYDROFABRIC)
    parser.add_argument("--river-map", type=Path, default=DEFAULT_RIVER_MAP)
    parser.add_argument("--hgrid", type=Path, default=DEFAULT_HGRID)
    parser.add_argument(
        "--bbox",
        nargs=4,
        type=float,
        default=DEFAULT_BBOX,
        metavar=("MIN_LON", "MIN_LAT", "MAX_LON", "MAX_LAT"),
    )
    parser.add_argument("--cache-dir", type=Path, default=DEFAULT_CACHE_DIR)
    parser.add_argument("--search-radius-m", type=float, default=500.0)
    parser.add_argument("--sample-spacing-m", type=float, default=25.0)
    parser.add_argument("--max-angle-deg", type=float, default=75.0)
    parser.add_argument("--min-overlap-m", type=float, default=25.0)
    parser.add_argument("--min-short-reach-fraction", type=float, default=0.10)
    parser.add_argument("--review-max-angle-deg", type=float, default=30.0)
    parser.add_argument("--review-min-overlap-m", type=float, default=100.0)
    parser.add_argument("--review-max-distance-m", type=float, default=250.0)
    parser.add_argument("--continuity-max-interval-gap-m", type=float, default=25.0)
    parser.add_argument("--force-cache", action="store_true")
    return parser


def main(argv: Sequence[str] | None = None) -> int:
    args = get_parser().parse_args(argv)
    bbox_lonlat = normalized_bbox(args.bbox)
    settings = MatchSettings(
        search_radius_m=args.search_radius_m,
        sample_spacing_m=args.sample_spacing_m,
        max_angle_deg=args.max_angle_deg,
        min_overlap_m=args.min_overlap_m,
        min_short_reach_fraction=args.min_short_reach_fraction,
        review_max_angle_deg=args.review_max_angle_deg,
        review_min_overlap_m=args.review_min_overlap_m,
        review_max_distance_m=args.review_max_distance_m,
        continuity_max_interval_gap_m=args.continuity_max_interval_gap_m,
    )
    args.cache_dir.mkdir(parents=True, exist_ok=True)
    hydro_path = build_hydrofabric_cache(
        args.hydrofabric,
        args.cache_dir,
        bbox_lonlat,
        settings.search_radius_m,
        force=args.force_cache,
    )
    _, centerlines_path = build_river_cache(
        args.river_map,
        args.cache_dir,
        bbox_lonlat,
        settings.search_radius_m,
        force=args.force_cache,
    )
    hgrid_cache = build_hgrid_cache(
        args.hgrid,
        args.cache_dir,
        bbox_lonlat,
        settings.search_radius_m,
        force=args.force_cache,
    )

    hydro = gpd.read_parquet(hydro_path)
    hydro["depth_valid"] = np.isfinite(hydro.bnk_depth) & (hydro.bnk_depth > 0)
    hydro_for_matching = hydro.loc[hydro.depth_valid].copy()
    rivers = gpd.read_parquet(centerlines_path)
    test_region_metric = buffered_region(bbox_lonlat, 0.0).geometry.iloc[0]
    candidates, selected, stations, continuity_intervals = run_matching_pipeline(
        rivers,
        hydro_for_matching,
        settings,
        region_geometry_metric=test_region_metric,
    )
    output_path = args.cache_dir / "hydrofabric_river_matches.gpkg"
    write_diagnostics(output_path, hydro, rivers, candidates, selected, stations)
    overview_path = args.cache_dir / "match_overview.png"
    write_overview_png(
        overview_path,
        bbox_lonlat,
        hydro,
        rivers,
        selected,
        stations,
        settings.search_radius_m,
    )
    matched_stations = int((stations.match_status == "matched").sum())
    regional_stations = stations.loc[stations.inside_test_region]
    regional_matched = int((regional_stations.match_status == "matched").sum())
    summary = {
        "bbox_lonlat": list(bbox_lonlat),
        "settings": settings.__dict__,
        "hydrofabric_reaches": len(hydro),
        "valid_depth_hydrofabric_reaches": len(hydro_for_matching),
        "river_centerlines": len(rivers),
        "candidate_intervals": len(candidates),
        "selected_intervals": len(selected),
        "continuity_fill_intervals": len(continuity_intervals),
        "continuity_filled_stations": int(
            (stations.match_method == "continuity_fill").sum()
        ),
        "continuity_fill_basis": {
            str(group): int(count)
            for group, count in stations.loc[
                stations.match_method == "continuity_fill", "continuity_basis"
            ].value_counts().items()
        },
        "stations": len(stations),
        "matched_stations": matched_stations,
        "matched_station_fraction": matched_stations / len(stations) if len(stations) else 0.0,
        "stations_inside_test_region": len(regional_stations),
        "matched_stations_inside_test_region": regional_matched,
        "matched_station_fraction_inside_test_region": (
            regional_matched / len(regional_stations) if len(regional_stations) else 0.0
        ),
        "selected_interval_diagnostic_groups": {
            str(group): int(count)
            for group, count in selected.quality_flag.value_counts().items()
        },
        "regional_unmatched_diagnostic_groups": {
            str(group): int(count)
            for group, count in regional_stations.loc[
                regional_stations.match_status == "unmatched",
                "diagnostic_group",
            ].value_counts().items()
        },
        "hgrid_cache": str(hgrid_cache),
        "diagnostics": str(output_path),
        "overview_png": str(overview_path),
    }
    summary_path = args.cache_dir / "summary.json"
    summary_path.write_text(json.dumps(summary, indent=2) + "\n")
    print(json.dumps(summary, indent=2), flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
