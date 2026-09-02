#!/usr/bin/env python3
"""MPI full-domain hydrofabric dredging with spatial mesh-tile ownership.

Each mesh tile has one owner rank.  Ranks load complete RiverMapper rivers
whose bounds touch their tiles, so reverse intersection screening sees every
nearby river without truncating bank polygons.  Forward and reverse requests
are written independently by rank and reduced deterministically by rank 0.
"""

from __future__ import annotations

import argparse
import json
import os
import sys
import traceback
from pathlib import Path
from typing import Sequence

os.environ.setdefault("OMP_NUM_THREADS", "1")
os.environ.setdefault("OPENBLAS_NUM_THREADS", "1")
os.environ.setdefault("MKL_NUM_THREADS", "1")

import geopandas as gpd
import numpy as np
import pandas as pd
from mpi4py import MPI
from scipy.spatial import cKDTree
from shapely.geometry import box

try:
    from stofs3d_setup.ops.Bathy_edit.Ensure_channel_connectivity.hydrofabric_dredge import (
        ANALYSIS_GEOGRAPHIC_CRS,
        DEFAULT_ARCS,
        DEFAULT_CENTERLINES,
        DEFAULT_EXCLUDE_REGIONS,
        DEFAULT_MATCHES,
        DEFAULT_WATERSHED,
        HgridNodes,
        MATCH_CRS,
        _atomic_parquet,
        _empty_intersection_targets,
        apply_station_depths,
        assign_intervals_to_stations,
        dredge_settings_from_namespace,
        expand_stations_to_arc_vertices,
        finalize_mesh_requests,
        load_effective_watershed,
        points_intersect_region,
        project_lonlat,
        read_hgrid_nodes,
        refresh_result_mesh_requests,
        screen_intersection_mesh_nodes,
        write_2dm_with_updates,
        write_hgrid_with_updates,
    )
    from stofs3d_setup.ops.Bathy_edit.Ensure_channel_connectivity.hydrofabric_match import (
        normalized_bbox,
    )
except ModuleNotFoundError:  # Permit direct execution from this directory.
    from hydrofabric_dredge import (  # type: ignore[no-redef]
        ANALYSIS_GEOGRAPHIC_CRS,
        DEFAULT_ARCS,
        DEFAULT_CENTERLINES,
        DEFAULT_EXCLUDE_REGIONS,
        DEFAULT_MATCHES,
        DEFAULT_WATERSHED,
        HgridNodes,
        MATCH_CRS,
        _atomic_parquet,
        _empty_intersection_targets,
        apply_station_depths,
        assign_intervals_to_stations,
        dredge_settings_from_namespace,
        expand_stations_to_arc_vertices,
        finalize_mesh_requests,
        load_effective_watershed,
        points_intersect_region,
        project_lonlat,
        read_hgrid_nodes,
        refresh_result_mesh_requests,
        screen_intersection_mesh_nodes,
        write_2dm_with_updates,
        write_hgrid_with_updates,
    )
    from hydrofabric_match import normalized_bbox  # type: ignore[no-redef]


DEFAULT_HGRID = Path(
    "/sciclone/schism10/Hgrid_projects/STOFS3D-v7.4/v32e/"
    "Bathy_edit_channel_variations/Channel_variations_verified_pre_channel/"
    "pre_channel/hgrid_pre_channel.gr3"
)
DEFAULT_OUTPUT_DIR = Path(
    "/sciclone/schism10/feiye/CIROH/Channel_Geometry/"
    "Connectivity_test/dredge_full_v32e"
)
HGRID_ARRAY_NAMES = (
    "node_id",
    "x",
    "y",
    "dp",
    "metric_xy",
    "inside_test_region",
    "inside_watershed",
)


def log(rank: int, message: str) -> None:
    print(f"[rank {rank}] {message}", flush=True)


def file_fingerprint(path: Path) -> dict[str, object]:
    stat = path.stat()
    return {
        "path": str(path.resolve()),
        "size": stat.st_size,
        "mtime_ns": stat.st_mtime_ns,
    }


def vector_fingerprint(path: Path) -> list[dict[str, object]]:
    files = sorted(path.parent.glob(f"{path.stem}.*")) if path.suffix != ".gpkg" else [path]
    return [file_fingerprint(item) for item in files]


def atomic_save_array(path: Path, values: np.ndarray) -> None:
    temporary = path.with_suffix(path.suffix + ".tmp")
    with temporary.open("wb") as stream:
        np.save(stream, values, allow_pickle=False)
    os.replace(temporary, path)


def hgrid_total_node_count(path: Path) -> int:
    with path.open("r", encoding="utf-8") as stream:
        stream.readline()
        fields = stream.readline().split()
    if len(fields) < 2:
        raise ValueError(f"Invalid hgrid header: {path}")
    return int(fields[1])


def prepare_hgrid_input(source: Path, cache_dir: Path, force: bool = False) -> Path:
    """Return a GR3 input, converting a canonical SMS 2DM checkpoint once.

    Bathy_edit's cmvd stage normally writes a 2DM file, while the dredging
    reader and streaming GR3 writer require SCHISM ordering.  The conversion
    is cached by the source file fingerprint so an unchanged checkpoint is
    never converted twice.
    """
    suffix = source.suffix.lower()
    if suffix in {".gr3", ".ll"}:
        return source
    if suffix != ".2dm":
        raise ValueError(
            f"Unsupported hgrid format {source.suffix!r}: expected .2dm, .gr3, or .ll"
        )

    cache_dir.mkdir(parents=True, exist_ok=True)
    output = cache_dir / "source_hgrid_from_2dm.gr3"
    metadata_path = cache_dir / "source_hgrid_from_2dm.metadata.json"
    metadata = {
        "kind": "hydrofabric_dredge_2dm_to_gr3_v1",
        "source_2dm": file_fingerprint(source),
    }
    if not force and cache_is_current(metadata_path, metadata, [output]):
        return output

    # Import lazily: only workflows whose canonical checkpoint is 2DM need
    # pylib's SMS reader and the associated full-mesh memory allocation.
    from pylib import sms2grd

    grid = sms2grd(str(source))
    temporary = output.with_suffix(".tmp.gr3")
    grid.write_hgrid(str(temporary), fmt=1)
    if hgrid_total_node_count(temporary) != int(grid.np):
        temporary.unlink(missing_ok=True)
        raise RuntimeError(f"Converted GR3 node count is invalid: {temporary}")
    os.replace(temporary, output)
    metadata_path.write_text(json.dumps(metadata, indent=2) + "\n")
    return output


def cache_is_current(
    metadata_path: Path,
    expected: dict[str, object],
    products: Sequence[Path],
) -> bool:
    if not metadata_path.exists() or not all(path.exists() for path in products):
        return False
    try:
        existing = json.loads(metadata_path.read_text())
    except (OSError, json.JSONDecodeError):
        return False
    return existing == expected


def build_input_cache(
    cache_dir: Path,
    hgrid_path: Path,
    matches_path: Path,
    effective_watershed: gpd.GeoDataFrame,
    watershed_paths: Sequence[Path],
    exclude_regions: Sequence[Path],
    bbox_lonlat: tuple[float, float, float, float] | None,
    force: bool,
) -> tuple[Path, Path]:
    """Build reusable hgrid arrays and a selected-interval Parquet cache."""
    cache_dir.mkdir(parents=True, exist_ok=True)
    hgrid_cache = cache_dir / "hgrid_nodes"
    hgrid_cache.mkdir(parents=True, exist_ok=True)
    array_paths = {name: hgrid_cache / f"{name}.npy" for name in HGRID_ARRAY_NAMES}
    hgrid_metadata_path = hgrid_cache / "metadata.json"
    hgrid_metadata = {
        "kind": "hydrofabric_dredge_hgrid_nodes_v1",
        "hgrid": file_fingerprint(hgrid_path),
        "watersheds": [vector_fingerprint(path) for path in watershed_paths],
        "exclude_regions": [vector_fingerprint(path) for path in exclude_regions],
        "bbox_lonlat": list(bbox_lonlat) if bbox_lonlat is not None else None,
        "metric_crs": MATCH_CRS,
        "source_total_node_count": hgrid_total_node_count(hgrid_path),
    }
    if force or not cache_is_current(
        hgrid_metadata_path, hgrid_metadata, list(array_paths.values())
    ):
        nodes = read_hgrid_nodes(hgrid_path, bbox_lonlat=bbox_lonlat)
        metric_xy = project_lonlat(nodes.x, nodes.y)
        inside_test = np.ones(len(nodes.node_id), dtype=bool)
        if bbox_lonlat is not None:
            xmin, ymin, xmax, ymax = bbox_lonlat
            inside_test = (
                (nodes.x >= xmin)
                & (nodes.x <= xmax)
                & (nodes.y >= ymin)
                & (nodes.y <= ymax)
            )
        inside_watershed = points_intersect_region(
            nodes.x, nodes.y, effective_watershed
        )
        values = {
            "node_id": nodes.node_id,
            "x": nodes.x,
            "y": nodes.y,
            "dp": nodes.dp,
            "metric_xy": metric_xy,
            "inside_test_region": inside_test,
            "inside_watershed": inside_watershed,
        }
        for name, array in values.items():
            atomic_save_array(array_paths[name], np.asarray(array))
        hgrid_metadata_path.write_text(json.dumps(hgrid_metadata, indent=2) + "\n")

    interval_path = cache_dir / "selected_intervals.parquet"
    interval_metadata_path = cache_dir / "selected_intervals.metadata.json"
    interval_metadata = {
        "kind": "hydrofabric_selected_intervals_v1",
        "matches": file_fingerprint(matches_path),
        "layer": "selected_intervals",
    }
    if force or not cache_is_current(
        interval_metadata_path, interval_metadata, [interval_path]
    ):
        selected = gpd.read_file(matches_path, layer="selected_intervals")
        temporary = interval_path.with_suffix(".tmp.parquet")
        selected.to_parquet(temporary, index=False, compression="zstd")
        os.replace(temporary, interval_path)
        interval_metadata_path.write_text(
            json.dumps(interval_metadata, indent=2) + "\n"
        )
    return hgrid_cache, interval_path


def load_hgrid_cache(cache_dir: Path) -> tuple[HgridNodes, np.ndarray, np.ndarray, np.ndarray]:
    metadata = json.loads((cache_dir / "metadata.json").read_text())
    arrays = {
        name: np.load(cache_dir / f"{name}.npy", mmap_mode="r")
        for name in HGRID_ARRAY_NAMES
    }
    hgrid = HgridNodes(
        node_id=arrays["node_id"],
        x=arrays["x"],
        y=arrays["y"],
        dp=arrays["dp"],
        total_node_count=int(metadata["source_total_node_count"]),
    )
    return (
        hgrid,
        arrays["metric_xy"],
        arrays["inside_test_region"],
        arrays["inside_watershed"],
    )


def tile_codes(metric_xy: np.ndarray, tile_size_m: float) -> np.ndarray:
    indices = np.floor(metric_xy / tile_size_m).astype(np.int64)
    return (indices[:, 0] << np.int64(32)) | (
        indices[:, 1] & np.int64(0xFFFFFFFF)
    )


def decode_tile(code: int) -> tuple[int, int]:
    value = int(code)
    tile_x = value >> 32
    tile_y = value & 0xFFFFFFFF
    if tile_y >= 0x80000000:
        tile_y -= 0x100000000
    return tile_x, tile_y


def tile_polygon(code: int, tile_size_m: float, halo_m: float = 0.0):
    tile_x, tile_y = decode_tile(code)
    xmin = tile_x * tile_size_m - halo_m
    ymin = tile_y * tile_size_m - halo_m
    return box(
        xmin,
        ymin,
        xmin + tile_size_m + 2 * halo_m,
        ymin + tile_size_m + 2 * halo_m,
    )


def _geometry_vertex_count(geometry) -> int:
    """Count vertices in line or multipart line geometry."""
    if geometry is None or geometry.is_empty:
        return 0
    if hasattr(geometry, "geoms"):
        return sum(_geometry_vertex_count(part) for part in geometry.geoms)
    return len(geometry.coords)


def river_geometry_inventory(
    centerlines: gpd.GeoDataFrame,
    arcs: gpd.GeoDataFrame,
) -> pd.DataFrame:
    """Record geometry counts used to prove that rank inputs are complete."""
    centerline_rows = pd.DataFrame(
        {
            "river_idx": centerlines.river_idx.to_numpy(dtype=np.int64),
            "centerline_vertex_count": centerlines.geometry.map(
                _geometry_vertex_count
            ).to_numpy(dtype=np.int64),
        }
    ).groupby("river_idx", as_index=False).agg(
        centerline_count=("river_idx", "size"),
        centerline_vertex_count=("centerline_vertex_count", "sum"),
    )
    arc_rows = pd.DataFrame(
        {
            "river_idx": arcs.river_idx.to_numpy(dtype=np.int64),
            "arc_vertex_count": arcs.geometry.map(
                _geometry_vertex_count
            ).to_numpy(dtype=np.int64),
        }
    ).groupby("river_idx", as_index=False).agg(
        arc_count=("river_idx", "size"),
        arc_vertex_count=("arc_vertex_count", "sum"),
    )
    inventory = centerline_rows.merge(
        arc_rows, on="river_idx", how="outer", validate="one_to_one"
    ).fillna(0)
    count_columns = [
        "river_idx",
        "centerline_count",
        "centerline_vertex_count",
        "arc_count",
        "arc_vertex_count",
    ]
    inventory[count_columns] = inventory[count_columns].astype(np.int64)
    return inventory[count_columns].sort_values("river_idx").reset_index(drop=True)


def validate_complete_river_geometry(
    centerlines: gpd.GeoDataFrame,
    arcs: gpd.GeoDataFrame,
    expected: pd.DataFrame,
) -> None:
    """Raise when a rank received clipped or incomplete river geometry."""
    actual = river_geometry_inventory(centerlines, arcs)
    columns = [
        "centerline_count",
        "centerline_vertex_count",
        "arc_count",
        "arc_vertex_count",
    ]
    expected_indexed = expected.set_index("river_idx")[columns].sort_index()
    actual_indexed = actual.set_index("river_idx")[columns].sort_index()
    if expected_indexed.equals(actual_indexed):
        return
    comparison = expected_indexed.join(
        actual_indexed,
        how="outer",
        lsuffix="_expected",
        rsuffix="_actual",
    )
    differing = comparison.loc[
        comparison.isna().any(axis=1)
        | np.any(
            comparison.filter(like="_expected").to_numpy()
            != comparison.filter(like="_actual").to_numpy(),
            axis=1,
        )
    ]
    raise RuntimeError(
        "MPI rank received incomplete river geometry for river_idx values "
        f"{differing.index.astype(int).tolist()}: {differing.to_dict('index')}"
    )


def plan_tiles(
    metric_xy: np.ndarray,
    inside_test_region: np.ndarray,
    arcs: gpd.GeoDataFrame,
    tile_size_m: float,
    halo_m: float,
    mpi_size: int,
) -> tuple[
    list[list[int]],
    list[list[int]],
    list[int],
    pd.DataFrame,
    list[list[int]],
    list[int],
]:
    """Balance core mesh tiles and return complete rivers needed by each rank."""
    core_xy = np.asarray(metric_xy[np.asarray(inside_test_region, dtype=bool)])
    core_codes = tile_codes(core_xy, tile_size_m)
    unique_tiles, node_counts = np.unique(core_codes, return_counts=True)
    if mpi_size > len(unique_tiles):
        raise ValueError(
            f"Requested {mpi_size} ranks for only {len(unique_tiles)} mesh tiles; "
            "use fewer ranks or a smaller --tile-size-m"
        )

    arcs_metric = arcs.to_crs(MATCH_CRS)
    bounds = arcs_metric.geometry.bounds
    river_meta = pd.DataFrame(
        {
            "river_idx": arcs_metric.river_idx.to_numpy(dtype=np.int64),
            "minx": bounds.minx.to_numpy(),
            "miny": bounds.miny.to_numpy(),
            "maxx": bounds.maxx.to_numpy(),
            "maxy": bounds.maxy.to_numpy(),
            "vertex_count": arcs_metric.geometry.map(
                lambda geometry: len(geometry.coords)
            ).to_numpy(dtype=np.int64),
        }
    ).groupby("river_idx", as_index=False).agg(
        minx=("minx", "min"),
        miny=("miny", "min"),
        maxx=("maxx", "max"),
        maxy=("maxy", "max"),
        vertex_count=("vertex_count", "sum"),
    )
    river_boxes = gpd.GeoDataFrame(
        river_meta[["river_idx", "vertex_count"]],
        geometry=[
            box(row.minx, row.miny, row.maxx, row.maxy)
            for row in river_meta.itertuples(index=False)
        ],
        crs=MATCH_CRS,
    )
    tile_frame = gpd.GeoDataFrame(
        {"tile": unique_tiles.astype(np.int64)},
        geometry=[tile_polygon(int(code), tile_size_m, halo_m) for code in unique_tiles],
        crs=MATCH_CRS,
    )
    joined = gpd.sjoin(tile_frame, river_boxes, predicate="intersects", how="left")
    tile_rivers = {
        int(tile): sorted(group.river_idx.dropna().astype(np.int64).unique().tolist())
        for tile, group in joined.groupby("tile", sort=False)
        if group.river_idx.notna().any()
    }
    vertex_weights = (
        joined.dropna(subset=["river_idx"])
        .drop_duplicates(["tile", "river_idx"])
        .groupby("tile").vertex_count.sum()
        .to_dict()
    )
    weights = {
        int(tile): int(nodes) + int(vertex_weights.get(int(tile), 0))
        for tile, nodes in zip(unique_tiles, node_counts)
        if int(tile) in tile_rivers
    }
    if mpi_size > len(weights):
        raise ValueError(
            f"Requested {mpi_size} ranks for only {len(weights)} river-active "
            "mesh tiles; use fewer ranks or a smaller --tile-size-m"
        )
    assignments: list[list[int]] = [[] for _ in range(mpi_size)]
    loads = [0] * mpi_size
    for tile, weight in sorted(weights.items(), key=lambda item: item[1], reverse=True):
        owner = int(np.argmin(loads))
        assignments[owner].append(tile)
        loads[owner] += weight
    spatial_rank_rivers = [
        sorted({river for tile in owned for river in tile_rivers.get(int(tile), [])})
        for owned in assignments
    ]
    forward_rank_rivers: list[list[int]] = [[] for _ in range(mpi_size)]
    forward_loads = [0] * mpi_size
    for row in river_meta.sort_values("vertex_count", ascending=False).itertuples(
        index=False
    ):
        owner = int(np.argmin(forward_loads))
        forward_rank_rivers[owner].append(int(row.river_idx))
        forward_loads[owner] += int(row.vertex_count)
    rank_rivers = [
        sorted(set(spatial) | set(forward))
        for spatial, forward in zip(spatial_rank_rivers, forward_rank_rivers)
    ]
    return (
        assignments,
        rank_rivers,
        loads,
        river_meta,
        forward_rank_rivers,
        forward_loads,
    )


def map_vertices_with_tree(
    vertices: pd.DataFrame,
    hgrid: HgridNodes,
    mesh_tree: cKDTree,
    tree_mesh_positions: np.ndarray | None = None,
) -> pd.DataFrame:
    if tree_mesh_positions is None:
        tree_mesh_positions = np.arange(len(hgrid.node_id), dtype=np.int64)
    else:
        tree_mesh_positions = np.asarray(tree_mesh_positions, dtype=np.int64)
    if len(tree_mesh_positions) != mesh_tree.n:
        raise ValueError("tree_mesh_positions must match the KD-tree point count")
    mapped = vertices.copy()
    vertex_xy = project_lonlat(
        mapped.x.to_numpy(dtype=float), mapped.y.to_numpy(dtype=float)
    )
    distances, tree_positions = mesh_tree.query(vertex_xy, k=1, workers=1)
    positions = tree_mesh_positions[tree_positions.astype(np.int64)]
    mapped["mesh_position"] = positions
    mapped["mesh_node_id"] = hgrid.node_id[positions]
    mapped["nearest_distance_m"] = distances
    mapped["mesh_dp"] = hgrid.dp[positions]
    return mapped


def rank_part_paths(parts_dir: Path, rank: int) -> dict[str, Path]:
    prefix = parts_dir / f"rank_{rank:04d}"
    return {
        "vertices": prefix.with_name(prefix.name + "_vertices.parquet"),
        "requests": prefix.with_name(prefix.name + "_requests.parquet"),
        "candidates": prefix.with_name(prefix.name + "_candidates.parquet"),
    }


def reduce_requested_nodes(
    requests: pd.DataFrame,
    max_dredging_delta_m: float = 6.0,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    if requests.empty:
        return finalize_mesh_requests(requests, max_dredging_delta_m)
    grouped = requests.groupby("mesh_position", sort=True)
    identity_columns = ["mesh_node_id", "x", "y", "original_dp"]
    inconsistent = (
        grouped[identity_columns].nunique(dropna=False).gt(1).any(axis=1)
    )
    if inconsistent.any():
        positions = inconsistent.index[inconsistent].astype(int).tolist()
        raise RuntimeError(
            "MPI ranks disagree on hgrid identity/depth for mesh positions "
            f"{positions[:20]}"
        )
    reduced = grouped.agg(
        mesh_node_id=("mesh_node_id", "first"),
        x=("x", "first"),
        y=("y", "first"),
        original_dp=("original_dp", "first"),
        requested_dp=("requested_dp", "max"),
        forward_requested=("forward_requested", "max"),
        intersection_target=("intersection_target", "max"),
        forward_vertex_requests=("forward_vertex_requests", "sum"),
        intersection_qualifying_river_count=(
            "intersection_qualifying_river_count",
            "max",
        ),
    ).reset_index()
    return finalize_mesh_requests(reduced, max_dredging_delta_m)


def write_gpkg(
    path: Path,
    requested: pd.DataFrame,
    changed: pd.DataFrame,
    candidates: pd.DataFrame,
    targets: pd.DataFrame,
    effective_watershed: gpd.GeoDataFrame,
) -> None:
    if path.exists():
        path.unlink()
    layers = {
        "dredge_requested_mesh_nodes": requested,
        "changed_mesh_nodes": changed,
        "intersection_candidates_200m": candidates,
        "intersection_target_mesh_nodes": targets,
    }
    if not requested.empty:
        capped = requested.loc[requested.dredging_delta_capped].copy()
        layers["capped_dredging_delta"] = capped
        layers["capped_raw_delta_gt10m"] = capped.loc[
            capped.requested_dredging_delta_m > 10.0
        ].copy()
    for layer, frame in layers.items():
        if frame.empty:
            continue
        gpd.GeoDataFrame(
            frame,
            geometry=gpd.points_from_xy(frame.x, frame.y),
            crs=ANALYSIS_GEOGRAPHIC_CRS,
        ).to_file(path, layer=layer, driver="GPKG", mode="w")
    effective_watershed.to_file(
        path, layer="effective_watershed", driver="GPKG", mode="w"
    )


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
    parser.add_argument("--exclude-region", action="append", type=Path, default=None)
    parser.add_argument("--output-dir", type=Path, default=DEFAULT_OUTPUT_DIR)
    parser.add_argument(
        "--input-cache-dir",
        type=Path,
        default=None,
        help="Reusable cache directory; defaults to OUTPUT_DIR/input_cache.",
    )
    parser.add_argument("--bbox", nargs=4, type=float, default=None)
    parser.add_argument("--tile-size-m", type=float, default=50_000.0)
    parser.add_argument("--min-channel-depth-m", type=float, default=1.0)
    parser.add_argument(
        "--channel-depth-source",
        choices=("hydrofabric", "constant"),
        default="hydrofabric",
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
    parser.add_argument("--unmatched-policy", choices=("baseline", "skip"), default="baseline")
    parser.add_argument("--measured-from-lower-bank", action="store_true")
    parser.add_argument("--disable-intersection-recovery", action="store_true")
    parser.add_argument("--force-input-cache", action="store_true")
    parser.add_argument("--overwrite", action="store_true")
    parser.add_argument("--keep-parts", action="store_true")
    parser.add_argument("--write-gpkg", action="store_true")
    parser.add_argument("--write-hgrid", action="store_true")
    return parser


def main(argv: Sequence[str] | None = None) -> int:
    comm = MPI.COMM_WORLD
    rank = comm.Get_rank()
    mpi_size = comm.Get_size()
    args = get_parser().parse_args(argv)
    if args.tile_size_m <= 0:
        raise ValueError("--tile-size-m must be positive")
    bbox_lonlat = normalized_bbox(args.bbox) if args.bbox is not None else None
    settings = dredge_settings_from_namespace(args, query_workers=1)
    watershed_regions = (
        tuple(args.watershed)
        if args.watershed is not None
        else (DEFAULT_WATERSHED,)
    )
    exclude_regions = tuple(args.exclude_region) if args.exclude_region else DEFAULT_EXCLUDE_REGIONS
    effective_watershed = load_effective_watershed(
        watershed_regions, exclude_regions
    )
    cache_dir = (
        args.input_cache_dir
        if args.input_cache_dir is not None
        else args.output_dir / "input_cache"
    )
    parts_dir = args.output_dir / "mpi_parts"

    if rank == 0:
        if args.output_dir.exists() and (args.output_dir / "summary.json").exists() and not args.overwrite:
            raise FileExistsError(
                f"Output exists: {args.output_dir}; pass --overwrite to replace products"
            )
        args.output_dir.mkdir(parents=True, exist_ok=True)
        parts_dir.mkdir(parents=True, exist_ok=True)
        prepared_hgrid = prepare_hgrid_input(
            args.hgrid, cache_dir, force=args.force_input_cache
        )
        hgrid_cache, interval_cache = build_input_cache(
            cache_dir,
            prepared_hgrid,
            args.matches_gpkg,
            effective_watershed,
            watershed_regions,
            exclude_regions,
            bbox_lonlat,
            args.force_input_cache,
        )
        hgrid_root, metric_xy_root, inside_test_root, _ = load_hgrid_cache(hgrid_cache)
        centerlines_root = gpd.read_parquet(args.river_centerlines)
        if bbox_lonlat is not None:
            selection = gpd.GeoDataFrame(
                geometry=[box(*bbox_lonlat)], crs=ANALYSIS_GEOGRAPHIC_CRS
            ).to_crs(centerlines_root.crs).geometry.iloc[0]
            centerlines_root = centerlines_root.loc[
                centerlines_root.geometry.intersects(selection)
            ].copy()
        river_ids = centerlines_root.river_idx.to_numpy(dtype=np.int64)
        arcs_root = gpd.read_parquet(args.river_arcs)
        arcs_root = arcs_root.loc[arcs_root.river_idx.isin(river_ids)].copy()
        (
            assignments,
            rank_rivers,
            planned_loads,
            river_meta,
            forward_rank_rivers,
            forward_loads,
        ) = plan_tiles(
            metric_xy_root,
            inside_test_root,
            arcs_root,
            args.tile_size_m,
            max(settings.max_nearest_distance_m, settings.intersection_search_radius_m),
            mpi_size,
        )
        geometry_inventory = river_geometry_inventory(centerlines_root, arcs_root)
        payload = {
            "hgrid_cache": str(hgrid_cache),
            "prepared_hgrid": str(prepared_hgrid),
            "interval_cache": str(interval_cache),
            "assignments": assignments,
            "rank_rivers": rank_rivers,
            "forward_rank_rivers": forward_rank_rivers,
            "planned_loads": planned_loads,
            "planned_forward_vertex_loads": forward_loads,
            "river_geometry_inventory": geometry_inventory.to_dict("records"),
            "rivers": len(centerlines_root),
            "source_total_node_count": hgrid_root.total_node_count,
        }
        log(rank, f"planned {sum(map(len, assignments))} tiles; loads={planned_loads}")
    else:
        payload = None
    payload = comm.bcast(payload, root=0)

    hgrid, metric_xy, inside_test, inside_watershed = load_hgrid_cache(
        Path(payload["hgrid_cache"])
    )
    forward_mesh_positions = np.flatnonzero(
        np.asarray(inside_watershed, dtype=bool)
    )
    if len(forward_mesh_positions) == 0:
        raise ValueError("No in-watershed hgrid nodes are available for the KD-tree")
    mesh_tree = cKDTree(
        np.asarray(metric_xy)[forward_mesh_positions],
        compact_nodes=False,
        balanced_tree=False,
    )
    all_tile_codes = tile_codes(np.asarray(metric_xy), args.tile_size_m)
    own_tiles = set(payload["assignments"][rank])
    own_node_mask = np.isin(all_tile_codes, np.asarray(sorted(own_tiles), dtype=np.int64))
    core_node_mask = own_node_mask & np.asarray(inside_test, dtype=bool)
    own_river_ids = np.asarray(payload["rank_rivers"][rank], dtype=np.int64)

    centerlines = gpd.read_parquet(args.river_centerlines)
    centerlines = centerlines.loc[centerlines.river_idx.isin(own_river_ids)].copy()
    arcs = gpd.read_parquet(args.river_arcs)
    arcs = arcs.loc[arcs.river_idx.isin(own_river_ids)].copy()
    expected_geometry = pd.DataFrame.from_records(
        payload["river_geometry_inventory"]
    )
    expected_geometry = expected_geometry.loc[
        expected_geometry.river_idx.isin(own_river_ids)
    ].copy()
    validate_complete_river_geometry(centerlines, arcs, expected_geometry)
    selected = gpd.read_parquet(payload["interval_cache"])
    selected = selected.loc[selected.river_idx.isin(own_river_ids)].copy()
    stations = assign_intervals_to_stations(centerlines, selected)
    station_lonlat = stations.to_crs(ANALYSIS_GEOGRAPHIC_CRS)
    # This field describes the requested geographic test region, not MPI
    # ownership. Vertex ownership is assigned later from the mapped mesh tile.
    stations["station_inside_test_region"] = True
    if bbox_lonlat is not None:
        xmin, ymin, xmax, ymax = bbox_lonlat
        stations["station_inside_test_region"] &= (
            (station_lonlat.geometry.x >= xmin)
            & (station_lonlat.geometry.x <= xmax)
            & (station_lonlat.geometry.y >= ymin)
            & (station_lonlat.geometry.y <= ymax)
        ).to_numpy()
    stations["station_inside_watershed"] = points_intersect_region(
        station_lonlat.geometry.x,
        station_lonlat.geometry.y,
        effective_watershed,
    )
    stations["station_inside_dredge_region"] = (
        stations.station_inside_test_region & stations.station_inside_watershed
    )
    vertices = expand_stations_to_arc_vertices(arcs, stations)
    vertex_in_test_region = np.ones(len(vertices), dtype=bool)
    if bbox_lonlat is not None:
        xmin, ymin, xmax, ymax = bbox_lonlat
        vertex_in_test_region &= (
            (vertices.x >= xmin) & (vertices.x <= xmax)
            & (vertices.y >= ymin) & (vertices.y <= ymax)
        )
    vertices["inside_watershed"] = points_intersect_region(
        vertices.x, vertices.y, effective_watershed
    )
    mapped = map_vertices_with_tree(
        vertices, hgrid, mesh_tree, forward_mesh_positions
    )
    # Complete rivers have one balanced forward/diagnostic owner. Intersection
    # candidates remain spatially owned by mesh tile.
    mapped["inside_test_region"] = vertex_in_test_region & mapped.river_idx.isin(
        payload["forward_rank_rivers"][rank]
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
            query_workers=1,
            mesh_candidate_mask=core_node_mask,
        )
    else:
        result.intersection_candidates = _empty_intersection_targets()
    result.intersection_targets = result.intersection_candidates.loc[
        result.intersection_candidates.intersection_target
    ].copy()
    refresh_result_mesh_requests(
        result, hgrid, max_dredging_delta_m=settings.max_dredging_delta_m
    )

    part_paths = rank_part_paths(parts_dir, rank)
    owned_vertices = result.vertices.loc[result.vertices.inside_test_region].copy()
    _atomic_parquet(owned_vertices, part_paths["vertices"])
    _atomic_parquet(result.requested_nodes, part_paths["requests"])
    _atomic_parquet(result.intersection_candidates, part_paths["candidates"])
    local_stats = {
        "rank": rank,
        "tiles": len(own_tiles),
        "rivers_loaded": len(own_river_ids),
        "owned_vertices": len(owned_vertices),
        "requests": len(result.requested_nodes),
        "intersection_candidates": len(result.intersection_candidates),
        "intersection_targets": len(result.intersection_targets),
    }
    log(rank, json.dumps(local_stats, sort_keys=True))
    all_stats = comm.gather(local_stats, root=0)
    comm.barrier()

    if rank == 0:
        assert all_stats is not None
        paths = [rank_part_paths(parts_dir, index) for index in range(mpi_size)]
        vertices_all = pd.concat(
            [pd.read_parquet(item["vertices"]) for item in paths], ignore_index=True
        )
        requests_all = pd.concat(
            [pd.read_parquet(item["requests"]) for item in paths], ignore_index=True
        )
        candidate_frames = [pd.read_parquet(item["candidates"]) for item in paths]
        candidates = pd.concat(candidate_frames, ignore_index=True)
        if not candidates.empty:
            if candidates.mesh_position.duplicated().any():
                raise RuntimeError("Intersection candidate mesh ownership is not unique")
            candidates = candidates.sort_values("mesh_position").reset_index(drop=True)
        requested, changed = reduce_requested_nodes(
            requests_all,
            max_dredging_delta_m=settings.max_dredging_delta_m,
        )
        if not candidates.empty:
            forward_positions = requested.loc[
                requested.forward_requested, "mesh_position"
            ].to_numpy(dtype=np.int64)
            candidates["forward_requested"] = candidates.mesh_position.isin(
                forward_positions
            )
            candidates["additional_target"] = (
                candidates.intersection_target & ~candidates.forward_requested
            )
        targets = candidates.loc[candidates.intersection_target].copy()
        final_by_node = pd.Series(changed.final_dp.to_numpy(), index=changed.mesh_node_id)
        vertices_all["final_node_dp"] = vertices_all.mesh_node_id.map(final_by_node)
        vertices_all["final_node_dp"] = vertices_all.final_node_dp.fillna(vertices_all.mesh_dp)
        vertices_all["node_dredging_delta_m"] = vertices_all.final_node_dp - vertices_all.mesh_dp

        products = {
            "vertex_mapping": args.output_dir / "river_vertex_mesh_mapping.parquet",
            "requested_nodes": args.output_dir / "dredge_requested_mesh_nodes.parquet",
            "node_changes": args.output_dir / "dredged_mesh_nodes.parquet",
            "intersection_candidates": args.output_dir / "intersection_candidates_200m.parquet",
            "intersection_targets": args.output_dir / "intersection_target_mesh_nodes.parquet",
        }
        _atomic_parquet(vertices_all, products["vertex_mapping"])
        _atomic_parquet(requested, products["requested_nodes"])
        _atomic_parquet(changed, products["node_changes"])
        _atomic_parquet(candidates, products["intersection_candidates"])
        _atomic_parquet(targets, products["intersection_targets"])
        if args.write_gpkg:
            gpkg_path = args.output_dir / "hydrofabric_dredge_diagnostics.gpkg"
            write_gpkg(gpkg_path, requested, changed, candidates, targets, effective_watershed)
            products["diagnostics_gpkg"] = gpkg_path
        if args.write_hgrid:
            output_hgrid = args.output_dir / "hgrid_hydrofabric_dredged.gr3"
            write_hgrid_with_updates(
                Path(payload["prepared_hgrid"]), output_hgrid, changed
            )
            products["hgrid"] = output_hgrid
            products["hgrid_gr3"] = output_hgrid
            if args.hgrid.suffix.lower() == ".2dm":
                output_2dm = args.output_dir / "hgrid_hydrofabric_dredged.2dm"
                write_2dm_with_updates(args.hgrid, output_2dm, changed)
                products["hgrid_2dm"] = output_2dm

        changed_delta = changed.dredging_delta_m.to_numpy(dtype=float)
        summary = {
            "bbox_lonlat": list(bbox_lonlat) if bbox_lonlat is not None else None,
            "settings": settings.__dict__,
            "mpi_ranks": mpi_size,
            "tile_size_m": args.tile_size_m,
            "work_tiles": sum(item["tiles"] for item in all_stats),
            "planned_loads": payload["planned_loads"],
            "planned_forward_vertex_loads": payload[
                "planned_forward_vertex_loads"
            ],
            "rank_statistics": all_stats,
            "source_total_hgrid_nodes": payload["source_total_node_count"],
            "cached_hgrid_nodes": len(hgrid.node_id),
            "river_count": payload["rivers"],
            "arc_vertex_count": len(vertices_all),
            "dredge_request_vertices": int(vertices_all.dredge_request.sum()),
            "unique_forward_requested_mesh_nodes": int(requested.forward_requested.sum()),
            "intersection_candidate_mesh_nodes_200m": len(candidates),
            "intersection_target_mesh_nodes": len(targets),
            "dredge_requested_mesh_nodes": len(requested),
            "changed_mesh_nodes": len(changed),
            "capped_dredging_mesh_nodes": int(
                requested.dredging_delta_capped.sum()
            ),
            "dredging_delta_m": {
                "min": float(np.min(changed_delta)) if len(changed_delta) else None,
                "median": float(np.median(changed_delta)) if len(changed_delta) else None,
                "p95": float(np.percentile(changed_delta, 95)) if len(changed_delta) else None,
                "max": float(np.max(changed_delta)) if len(changed_delta) else None,
            },
            "inputs": {
                "hgrid": str(args.hgrid),
                "prepared_hgrid": str(payload["prepared_hgrid"]),
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
            "products": {key: str(value) for key, value in products.items()},
        }
        (args.output_dir / "summary.json").write_text(json.dumps(summary, indent=2) + "\n")
        print(json.dumps(summary, indent=2), flush=True)
        if not args.keep_parts:
            for item in paths:
                for path in item.values():
                    path.unlink(missing_ok=True)
            try:
                parts_dir.rmdir()
            except OSError:
                pass
    comm.barrier()
    return 0


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except SystemExit:
        raise
    except Exception:
        traceback.print_exc()
        try:
            MPI.COMM_WORLD.Abort(1)
        finally:
            sys.exit(1)
