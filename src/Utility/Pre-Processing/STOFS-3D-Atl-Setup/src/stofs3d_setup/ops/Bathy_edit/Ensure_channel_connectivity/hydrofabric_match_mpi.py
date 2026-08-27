#!/usr/bin/env python3
"""MPI full-domain RiverMapper-to-hydrofabric matching.

Complete RiverMapper rivers are grouped into spatial work tiles and assigned
to MPI ranks with a greedy station-count balance. Each rank reads only the
indexed hydrofabric features near its assigned rivers. Rank-local GeoPackages
are merged by rank 0, so concurrent processes never write the same file.
"""

from __future__ import annotations

import argparse
import json
import os
import subprocess
import sys
import traceback
import warnings
from collections import Counter
from pathlib import Path
from typing import Sequence

import fiona
import geopandas as gpd
import numpy as np
import pandas as pd
from mpi4py import MPI
from shapely.geometry import box

try:
    from stofs3d_setup.ops.Bathy_edit.Ensure_channel_connectivity.hydrofabric_match import (
        ANALYSIS_GEOGRAPHIC_CRS,
        DEFAULT_HYDROFABRIC,
        DEFAULT_RIVER_MAP,
        MATCH_CRS,
        MatchSettings,
        build_river_cache,
        cache_is_current,
        normalized_bbox,
        run_matching_pipeline,
        shapefile_fingerprint,
        write_metadata,
    )
except ModuleNotFoundError:  # Permit direct execution from this source directory.
    from hydrofabric_match import (  # type: ignore[no-redef]
        ANALYSIS_GEOGRAPHIC_CRS,
        DEFAULT_HYDROFABRIC,
        DEFAULT_RIVER_MAP,
        MATCH_CRS,
        MatchSettings,
        build_river_cache,
        cache_is_current,
        normalized_bbox,
        run_matching_pipeline,
        shapefile_fingerprint,
        write_metadata,
    )


DEFAULT_OUTPUT_DIR = Path(
    "/sciclone/schism10/feiye/CIROH/Channel_Geometry/"
    "Connectivity_test/full_domain"
)
HYDRO_LAYER = "hydrofabric"
HYDRO_CACHE_COLUMNS = [
    "OBJECTID",
    "COMID",
    "bnk_depth",
    "bnk_width",
    "StreamOrde",
    "TotDASqKM",
]


def log(rank: int, message: str) -> None:
    print(f"[rank {rank}] {message}", flush=True)


def run_command(command: list[str]) -> None:
    subprocess.run(command, check=True)


def build_indexed_hydrofabric_cache(
    source: Path,
    cache_path: Path,
    source_bounds: tuple[float, float, float, float],
    force: bool = False,
) -> Path:
    """Create a spatially indexed, valid-depth, metric hydrofabric cache."""
    metadata_path = cache_path.with_suffix(".metadata.json")
    metadata = {
        "kind": "indexed_hydrofabric",
        "source": shapefile_fingerprint(source),
        "source_bounds": list(source_bounds),
        "cache_crs": MATCH_CRS,
        "layer": HYDRO_LAYER,
        "columns": HYDRO_CACHE_COLUMNS,
        "depth_filter": "bnk_depth IS NOT NULL AND bnk_depth > 0",
    }
    if not force and cache_is_current(metadata_path, metadata, [cache_path]):
        print(f"Reusing indexed hydrofabric cache: {cache_path}", flush=True)
        return cache_path

    cache_path.parent.mkdir(parents=True, exist_ok=True)
    temporary = cache_path.with_name(f"{cache_path.stem}.tmp.gpkg")
    if temporary.exists():
        temporary.unlink()
    xmin, ymin, xmax, ymax = source_bounds
    command = [
        "ogr2ogr",
        "-f",
        "GPKG",
        str(temporary),
        str(source),
        "-nln",
        HYDRO_LAYER,
        "-dim",
        "XY",
        "-t_srs",
        MATCH_CRS,
        "-spat",
        str(xmin),
        str(ymin),
        str(xmax),
        str(ymax),
        "-where",
        "bnk_depth IS NOT NULL AND bnk_depth > 0",
        "-select",
        ",".join(HYDRO_CACHE_COLUMNS),
        "-lco",
        "SPATIAL_INDEX=YES",
    ]
    print("Building indexed hydrofabric cache (one streamed source scan)", flush=True)
    run_command(command)
    os.replace(temporary, cache_path)
    write_metadata(metadata_path, metadata)
    print(f"Created indexed hydrofabric cache: {cache_path}", flush=True)
    return cache_path


def buffered_source_bounds(
    rivers_metric: gpd.GeoDataFrame, search_radius_m: float
) -> tuple[float, float, float, float]:
    footprint = gpd.GeoDataFrame(
        geometry=[box(*rivers_metric.total_bounds).buffer(search_radius_m)],
        crs=MATCH_CRS,
    ).to_crs(ANALYSIS_GEOGRAPHIC_CRS)
    return tuple(float(value) for value in footprint.total_bounds)


def assign_tiles(
    rivers_metric: gpd.GeoDataFrame, tile_size_m: float, mpi_size: int
) -> tuple[gpd.GeoDataFrame, list[list[str]], list[int]]:
    """Assign complete rivers to tiles, then balance tiles across ranks."""
    centroids = rivers_metric.geometry.centroid
    tile_x = np.floor(centroids.x.to_numpy() / tile_size_m).astype(np.int64)
    tile_y = np.floor(centroids.y.to_numpy() / tile_size_m).astype(np.int64)
    tiled = rivers_metric.copy()
    tiled["work_tile"] = [f"{x}_{y}" for x, y in zip(tile_x, tile_y)]
    weights = (
        tiled.groupby("work_tile", sort=False).n_stations.sum().astype(int).to_dict()
    )
    assignments: list[list[str]] = [[] for _ in range(mpi_size)]
    loads = [0] * mpi_size
    for tile, weight in sorted(weights.items(), key=lambda item: item[1], reverse=True):
        rank = int(np.argmin(loads))
        assignments[rank].append(str(tile))
        loads[rank] += int(weight)
    return tiled, assignments, loads


def read_hydrofabric_slice(
    cache_path: Path, rivers_metric: gpd.GeoDataFrame, search_radius_m: float
) -> gpd.GeoDataFrame:
    xmin, ymin, xmax, ymax = rivers_metric.total_bounds
    bounds = (
        xmin - search_radius_m,
        ymin - search_radius_m,
        xmax + search_radius_m,
        ymax + search_radius_m,
    )
    hydro = gpd.read_file(cache_path, layer=HYDRO_LAYER, bbox=bounds)
    missing = sorted(set(HYDRO_CACHE_COLUMNS) - set(hydro.columns))
    if missing:
        raise ValueError(f"Indexed hydrofabric cache is missing columns: {missing}")
    hydro["hydro_feature_idx"] = hydro.OBJECTID.astype(np.int64)
    return hydro


def concat_geodataframes(frames: list[gpd.GeoDataFrame]) -> gpd.GeoDataFrame:
    if not frames:
        raise ValueError("At least one GeoDataFrame is required")
    nonempty = [frame for frame in frames if not frame.empty]
    if not nonempty:
        return frames[0].copy()
    with warnings.catch_warnings():
        warnings.filterwarnings(
            "ignore",
            message="The behavior of DataFrame concatenation with empty or all-NA",
            category=FutureWarning,
        )
        combined = pd.concat(nonempty, ignore_index=True)
    return gpd.GeoDataFrame(
        combined, geometry="geometry", crs=nonempty[0].crs
    )


def diagnostic_layers(
    candidates: gpd.GeoDataFrame,
    selected: gpd.GeoDataFrame,
    stations: gpd.GeoDataFrame,
) -> dict[str, gpd.GeoDataFrame]:
    flagged = selected.loc[selected.quality_flag != "nominal"].copy()
    unmatched = stations.loc[
        (stations.match_status == "unmatched") & stations.inside_test_region
    ].copy()
    continuity_intervals = selected.loc[
        selected.quality_flag == "continuity_fill"
    ].copy()
    continuity_stations = stations.loc[
        stations.match_method == "continuity_fill"
    ].copy()
    return {
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


def write_rank_part(
    path: Path, layers: dict[str, gpd.GeoDataFrame]
) -> dict[str, int]:
    if path.exists():
        path.unlink()
    counts: dict[str, int] = {}
    for layer_name, frame in layers.items():
        counts[layer_name] = len(frame)
        if frame.empty:
            continue
        frame.to_crs(ANALYSIS_GEOGRAPHIC_CRS).to_file(
            path, layer=layer_name, driver="GPKG", mode="w"
        )
    return counts


def merge_rank_parts(
    output_path: Path,
    rivers: gpd.GeoDataFrame,
    part_paths: list[Path],
    rank_counts: list[dict[str, int]],
    hydro_cache: Path,
    include_hydrofabric: bool,
) -> None:
    if output_path.exists():
        output_path.unlink()
    rivers.to_crs(ANALYSIS_GEOGRAPHIC_CRS).to_file(
        output_path, layer="river_centerlines", driver="GPKG"
    )
    if include_hydrofabric:
        run_command(
            [
                "ogr2ogr",
                "-update",
                "-overwrite",
                str(output_path),
                str(hydro_cache),
                HYDRO_LAYER,
                "-nln",
                HYDRO_LAYER,
                "-t_srs",
                ANALYSIS_GEOGRAPHIC_CRS,
            ]
        )

    layer_names = sorted({name for counts in rank_counts for name in counts})
    for layer_name in layer_names:
        sources = [
            path
            for path, counts in zip(part_paths, rank_counts)
            if counts.get(layer_name, 0) > 0
        ]
        for source_index, source in enumerate(sources):
            action = "-overwrite" if source_index == 0 else "-append"
            run_command(
                [
                    "ogr2ogr",
                    "-update",
                    action,
                    str(output_path),
                    str(source),
                    layer_name,
                    "-nln",
                    layer_name,
                ]
            )


def layer_count(path: Path, layer: str) -> int:
    with fiona.open(path, layer=layer) as collection:
        return len(collection)


def get_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--hydrofabric", type=Path, default=DEFAULT_HYDROFABRIC)
    parser.add_argument("--river-map", type=Path, default=DEFAULT_RIVER_MAP)
    parser.add_argument("--output-dir", type=Path, default=DEFAULT_OUTPUT_DIR)
    parser.add_argument(
        "--bbox",
        nargs=4,
        type=float,
        default=None,
        metavar=("MIN_LON", "MIN_LAT", "MAX_LON", "MAX_LAT"),
        help="Optional validation subset; omit for the full RiverMapper domain.",
    )
    parser.add_argument("--tile-size-m", type=float, default=100_000.0)
    parser.add_argument("--search-radius-m", type=float, default=500.0)
    parser.add_argument("--sample-spacing-m", type=float, default=25.0)
    parser.add_argument("--max-angle-deg", type=float, default=75.0)
    parser.add_argument("--min-overlap-m", type=float, default=25.0)
    parser.add_argument("--min-short-reach-fraction", type=float, default=0.10)
    parser.add_argument("--review-max-angle-deg", type=float, default=30.0)
    parser.add_argument("--review-min-overlap-m", type=float, default=100.0)
    parser.add_argument("--review-max-distance-m", type=float, default=250.0)
    parser.add_argument("--continuity-max-interval-gap-m", type=float, default=25.0)
    parser.add_argument("--force-input-cache", action="store_true")
    parser.add_argument("--overwrite", action="store_true")
    parser.add_argument("--keep-parts", action="store_true")
    parser.add_argument(
        "--include-hydrofabric-layer",
        action="store_true",
        help="Duplicate the indexed hydrofabric into the final diagnostic GPKG.",
    )
    return parser


def main(argv: Sequence[str] | None = None) -> int:
    comm = MPI.COMM_WORLD
    rank = comm.Get_rank()
    mpi_size = comm.Get_size()
    args = get_parser().parse_args(argv)
    if args.tile_size_m <= 0:
        raise ValueError("--tile-size-m must be positive")
    bbox_lonlat = normalized_bbox(args.bbox) if args.bbox is not None else None
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
    output_path = args.output_dir / "hydrofabric_river_matches.gpkg"
    cache_dir = args.output_dir / "input_cache"
    hydro_cache = cache_dir / "hydrofabric_indexed.gpkg"
    parts_dir = args.output_dir / "mpi_parts"

    if rank == 0:
        if output_path.exists() and not args.overwrite:
            raise FileExistsError(
                f"Output exists: {output_path}; pass --overwrite to replace it"
            )
        args.output_dir.mkdir(parents=True, exist_ok=True)
        cache_dir.mkdir(parents=True, exist_ok=True)
        parts_dir.mkdir(parents=True, exist_ok=True)
        _, centerlines_path = build_river_cache(
            args.river_map,
            cache_dir,
            bbox_lonlat,
            settings.search_radius_m,
            force=args.force_input_cache,
        )
        rivers = gpd.read_parquet(centerlines_path)
        rivers_metric = rivers.to_crs(MATCH_CRS)
        if bbox_lonlat is None:
            with fiona.open(args.hydrofabric) as source_collection:
                source_bounds = tuple(
                    float(value) for value in source_collection.bounds
                )
        else:
            source_bounds = buffered_source_bounds(
                rivers_metric, settings.search_radius_m
            )
        build_indexed_hydrofabric_cache(
            args.hydrofabric,
            hydro_cache,
            source_bounds,
            force=args.force_input_cache,
        )
        tiled_rivers, assignments, planned_loads = assign_tiles(
            rivers_metric, args.tile_size_m, mpi_size
        )
        tile_count = int(tiled_rivers.work_tile.nunique())
        if mpi_size > tile_count:
            raise ValueError(
                f"Requested {mpi_size} MPI ranks for only {tile_count} work tiles; "
                "use fewer ranks or a smaller --tile-size-m"
            )
        log(
            rank,
            f"planned {tile_count} tiles and "
            f"{len(tiled_rivers)} complete rivers; station loads={planned_loads}",
        )
        payload = (tiled_rivers, assignments, planned_loads)
    else:
        rivers = None
        payload = None

    tiled_rivers, assignments, planned_loads = comm.bcast(payload, root=0)
    own_tiles = assignments[rank]
    if bbox_lonlat is None:
        region_geometry_metric = None
    else:
        region_geometry_metric = gpd.GeoDataFrame(
            geometry=[box(*bbox_lonlat)], crs=ANALYSIS_GEOGRAPHIC_CRS
        ).to_crs(MATCH_CRS).geometry.iloc[0]

    candidates_parts: list[gpd.GeoDataFrame] = []
    selected_parts: list[gpd.GeoDataFrame] = []
    station_parts: list[gpd.GeoDataFrame] = []
    for tile_number, tile in enumerate(own_tiles, start=1):
        tile_rivers = tiled_rivers.loc[tiled_rivers.work_tile == tile].drop(
            columns="work_tile"
        )
        hydro = read_hydrofabric_slice(
            hydro_cache, tile_rivers, settings.search_radius_m
        )
        candidates, selected, stations, _ = run_matching_pipeline(
            tile_rivers,
            hydro,
            settings,
            region_geometry_metric=region_geometry_metric,
        )
        candidates_parts.append(candidates)
        selected_parts.append(selected)
        station_parts.append(stations)
        log(
            rank,
            f"tile {tile_number}/{len(own_tiles)} {tile}: "
            f"{len(tile_rivers)} rivers, {len(hydro)} reaches, "
            f"{len(stations)} stations",
        )

    if candidates_parts:
        candidates = concat_geodataframes(candidates_parts)
        selected = concat_geodataframes(selected_parts)
        stations = concat_geodataframes(station_parts)
    else:
        raise RuntimeError(f"Rank {rank} received no work tiles; use fewer MPI ranks")
    layers = diagnostic_layers(candidates, selected, stations)
    part_path = parts_dir / f"rank_{rank:04d}.gpkg"
    layer_counts = write_rank_part(part_path, layers)
    local_stats = {
        "rank": rank,
        "tiles": len(own_tiles),
        "rivers": len(tiled_rivers.loc[tiled_rivers.work_tile.isin(own_tiles)]),
        "stations": len(stations),
        "matched_stations": int((stations.match_status == "matched").sum()),
        "inside_stations": int(stations.inside_test_region.sum()),
        "inside_matched_stations": int(
            ((stations.match_status == "matched") & stations.inside_test_region).sum()
        ),
        "continuity_filled_stations": int(
            (stations.match_method == "continuity_fill").sum()
        ),
        "continuity_fill_basis": stations.loc[
            stations.match_method == "continuity_fill", "continuity_basis"
        ].value_counts().astype(int).to_dict(),
        "unmatched_diagnostic_groups": stations.loc[
            (stations.match_status == "unmatched") & stations.inside_test_region,
            "diagnostic_group",
        ].value_counts().astype(int).to_dict(),
        "selected_diagnostic_groups": selected.quality_flag.value_counts()
        .astype(int)
        .to_dict(),
        "layer_counts": layer_counts,
    }
    all_stats = comm.gather(local_stats, root=0)
    comm.barrier()

    if rank == 0:
        assert rivers is not None and all_stats is not None
        part_paths = [parts_dir / f"rank_{index:04d}.gpkg" for index in range(mpi_size)]
        merge_rank_parts(
            output_path,
            rivers,
            part_paths,
            [stats["layer_counts"] for stats in all_stats],
            hydro_cache,
            args.include_hydrofabric_layer,
        )
        total_stations = sum(stats["stations"] for stats in all_stats)
        matched_stations = sum(stats["matched_stations"] for stats in all_stats)
        inside_stations = sum(stats["inside_stations"] for stats in all_stats)
        inside_matched = sum(stats["inside_matched_stations"] for stats in all_stats)
        continuity_basis: Counter[str] = Counter()
        unmatched_groups: Counter[str] = Counter()
        selected_groups: Counter[str] = Counter()
        for stats in all_stats:
            continuity_basis.update(stats["continuity_fill_basis"])
            unmatched_groups.update(stats["unmatched_diagnostic_groups"])
            selected_groups.update(stats["selected_diagnostic_groups"])
        total_layer_counts = {
            layer: sum(stats["layer_counts"].get(layer, 0) for stats in all_stats)
            for layer in all_stats[0]["layer_counts"]
        }
        summary = {
            "bbox_lonlat": list(bbox_lonlat) if bbox_lonlat is not None else None,
            "settings": settings.__dict__,
            "mpi_ranks": mpi_size,
            "tile_size_m": args.tile_size_m,
            "work_tiles": sum(stats["tiles"] for stats in all_stats),
            "planned_station_loads": planned_loads,
            "rank_statistics": all_stats,
            "hydrofabric_cache": str(hydro_cache),
            "valid_depth_hydrofabric_reaches": layer_count(
                hydro_cache, HYDRO_LAYER
            ),
            "river_centerlines": len(rivers),
            "candidate_intervals": total_layer_counts["candidate_overlaps"],
            "selected_intervals": total_layer_counts["selected_intervals"],
            "continuity_fill_intervals": total_layer_counts[
                "continuity_fill_intervals"
            ],
            "stations": total_stations,
            "matched_stations": matched_stations,
            "matched_station_fraction": (
                matched_stations / total_stations if total_stations else 0.0
            ),
            "stations_inside_domain": inside_stations,
            "matched_stations_inside_domain": inside_matched,
            "matched_station_fraction_inside_domain": (
                inside_matched / inside_stations if inside_stations else 0.0
            ),
            "continuity_filled_stations": sum(
                stats["continuity_filled_stations"] for stats in all_stats
            ),
            "continuity_fill_basis": dict(continuity_basis),
            "selected_interval_diagnostic_groups": dict(selected_groups),
            "unmatched_diagnostic_groups": dict(unmatched_groups),
            "diagnostics": str(output_path),
            "hydrofabric_embedded_in_diagnostics": args.include_hydrofabric_layer,
        }
        summary_path = args.output_dir / "summary.json"
        summary_path.write_text(json.dumps(summary, indent=2) + "\n")
        print(json.dumps(summary, indent=2), flush=True)
        if not args.keep_parts:
            for part_path in part_paths:
                part_path.unlink(missing_ok=True)
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
