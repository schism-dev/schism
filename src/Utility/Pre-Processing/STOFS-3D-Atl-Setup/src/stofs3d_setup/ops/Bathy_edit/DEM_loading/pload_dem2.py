#!/usr/bin/env python3
'''
Load bathymetry for SCHISM hgrid

Recommended usage:
Run this script in the SCHISM source code folder where it originally resides.
Select a suitable YAML recipe with --version or provide a custom recipe with
--config. Version-specific source lists, working directories, and merge rules
are kept in configs/*.yaml; Python implements the allowed operations.

Needs to prepare the following files/folder:
   - hgrid.ll: hgrid in lon/lat

   - DEM tiles: pre-processed from *.tif and *.asc to python's *.npz format for faster reading.

   - DEM*.json: specifying DEM tiles used, and arranged from lower priority to higher priority.
         See samples in this folder.

   - region (optional): A polygon shapefile specifiying regions where special treatment is applied.
        Sample *.shp is provided in this folder.


Example:
   mpirun -np 48 python pload_dem2.py --version stofs3d_v7p4
'''


import argparse
import errno
import json
import os
import re
import shutil
import socket
import sys
import time
from copy import deepcopy
from glob import glob
from pathlib import Path
from typing import Any

try:
    import yaml
except ImportError:  # pragma: no cover - depends on runtime environment
    yaml = None

import numpy as np
from mpi4py import MPI
import geopandas as gpd
from pylib import load_bathymetry, zdata, convert_dem_format

import rasterio
from rasterio.warp import calculate_default_transform, reproject, Resampling

# Attempt to import experimental grid reader for speed-up
if 'gulf' in socket.gethostname():
    from pylib_experimental.schism_file import cread_schism_hgrid as schism_read
    print('Using c++ function to accelerate hgrid reading')
else:
    from pylib import schism_grid as schism_read
    print('Using python function to read hgrid')

script_dir = os.path.dirname(os.path.realpath(__file__))
SCRIPT_DIR = Path(__file__).resolve().parent
CONFIG_DIR = SCRIPT_DIR / "configs"
DEFAULT_VERSION = "stofs3d_v7p4"


# ---------------------- MPI Utilities ----------------------
def initialize_mpi():
    """
    Initialize MPI communicator and return size and rank.
    """
    comm = MPI.COMM_WORLD
    return comm, comm.Get_size(), comm.Get_rank()


# ---------------------- File Handling ----------------------
def check_dem_files(dem_file_headers, dem_dir):
    """
    Check for the existence of DEM files in the specified directory.
    Raises FileNotFoundError if any files are missing.
    """
    for header in dem_file_headers:
        fnames = [i for i in os.listdir(dem_dir) if i.startswith(header)]
        if not fnames:
            raise FileNotFoundError(errno.ENOENT, os.strerror(errno.ENOENT), f'{header}*.npz')


def sort_dem_files(dem_file_headers, dem_dir):
    """
    Sort DEM files by ID numbers.
    """
    fnames0 = np.array([i for i in os.listdir(dem_dir) if i.endswith('.npz')])
    fnames_sort = []

    for header in dem_file_headers:
        fnames_sub = np.array([fname for fname in fnames0 if fname.startswith(header)])
        if len(fnames_sub) == 1:
            fnames_sort.extend(fnames_sub)
            continue

        # Sort by ID number, the last number before the suffix ".*"
        fid = np.array([int(re.findall(r'(\d+)', i)[-1]) for i in fnames_sub])
        sind = np.argsort(fid)
        fnames_sort.extend(fnames_sub[sind])

    return np.array(fnames_sort), fnames0


# ---------------------- Bathymetry Loading ----------------------
def initialize_grid_data():
    """
    Initialize the zdata object for saving grid data.
    """
    grid_data = zdata()
    grid_data.dp = dict()
    grid_data.sind = dict()
    return grid_data


def distribute_files(fnames_sort, nproc, myrank):
    """
    Distribute DEM files among MPI processes.
    """
    fnames, inum = [], []
    for m, fname in enumerate(fnames_sort):
        if m % nproc == myrank:
            fnames.append(fname)
            inum.append(m)
    return np.array(fnames), inum


def load_bathymetry_on_core(gd, dem_dir, fnames, reverse_sign, myrank):
    """
    Load bathymetry on each core and return a zdata object with results.
    """
    grid_data = initialize_grid_data()
    for fname in fnames:
        bname = fname.split('.')[0]
        dpi, sindi = load_bathymetry(gd.x, gd.y, os.path.join(dem_dir, fname), fmt=1)
        if reverse_sign:
            dpi = -dpi
        grid_data.dp[bname] = dpi
        grid_data.sind[bname] = sindi
        print(f'Finished reading {fname}, myrank={myrank}')
        sys.stdout.flush()
    return grid_data


def gather_and_combine_results(comm, myrank, grid_data):
    """
    Gather and combine bathymetry results from all cores.
    """
    gathered_data = comm.gather(grid_data, root=0)

    if myrank == 0:
        combined_data = initialize_grid_data()
        for this_data in gathered_data:
            combined_data.dp.update(this_data.dp)
            combined_data.sind.update(this_data.sind)
        return combined_data
    else:
        return None


def set_final_bathymetry(grid_data, gd, fnames_sort, fnames0):
    """
    Set the final bathymetry in the grid object and return diagnostic information.
    """
    did = np.ones(gd.np, dtype=int)
    dname = []

    for i, fname in enumerate(fnames_sort):
        bname = fname.split('.')[0]
        sind = grid_data.sind[bname]
        dp = grid_data.dp[bname]
        gd.dp[sind] = dp
        did[sind] = i + 1
        dname.append([k for k in fnames0 if k.startswith(fname)][0])

    return did, dname


def prep_dir(wdir, dem_json_list, region_shpfile_list):
    """
    Prepare directory for output files.
    """
    os.system(f'cp {__file__} {wdir}')
    for dem_json_file in dem_json_list:
        os.system(f'cp {script_dir}/{dem_json_file} {wdir}')
    for region_shpfile in region_shpfile_list:
        os.system(f'cp {script_dir}/{Path(region_shpfile).stem}.* {wdir}')


def write_diagnostics(grdout, did, dname):
    """
    Write diagnostic files with DEM IDs and names.
    """
    with open(f'{grdout}_dem_id', 'w', encoding='utf-8') as fid:
        fid.writelines(f'{i}\n' for i in did)

    with open(f'{grdout}_dem_name', 'w', encoding='utf-8') as fid:
        fid.writelines(f'{i + 1}: {name}\n' for i, name in enumerate(dname))


# ---------------------- core routine ----------------------
def pload_dem(grd, grdout, dem_json, dem_dir, reverse_sign=True):
    """
    Core routine to load bathymetry using mpi.
    """
    comm, nproc, myrank = initialize_mpi()

    # read inputs
    if myrank == 0:
        gd = schism_read(grd)
        gd.dp[:] = -9999  # Reset depth to -9999
        dem_info_dict = json.load(open(dem_json, encoding='utf-8'))
        check_dem_files(dem_info_dict.keys(), dem_dir)
        start_time = time.time()
    else:
        gd, dem_info_dict = None, None

    # Broadcast grid and DEM info to all cores
    gd = comm.bcast(gd, 0)
    dem_info_dict = comm.bcast(dem_info_dict, 0)

    # Sort DEM files
    fnames_sort, fnames0 = sort_dem_files(dem_info_dict.keys(), dem_dir)

    # Distribute DEM files among cores
    fnames, _ = distribute_files(fnames_sort, nproc, myrank)

    # Load bathymetry on each core
    grid_data = load_bathymetry_on_core(gd, dem_dir, fnames, reverse_sign, myrank)

    # Gather and combine results from all cores
    comm.Barrier()
    combined_data = gather_and_combine_results(comm, myrank, grid_data)  # returns None for non-root processes

    # Save combined results and write output if root process
    if myrank == 0:
        did, dname = set_final_bathymetry(combined_data, gd, fnames_sort, fnames0)
        if grdout is not None:
            gd.write_hgrid(grdout)
            write_diagnostics(grdout, did, dname)

        print(f'Total time used: {time.time() - start_time:.2f} s')
        return gd
    else:
        return None


def max_dp_in_region(grid_list: list, region_file: str, primary_grid_idx: int = 0, reverse_sign: bool = False):
    """
    In a specified region, copy the depth of the first grid in a list
    and change it to the maximum depth of all grids.
    If region_file is None, the operation applies to the entire grid.

    if reverse_sign is True, the minimum depth is taken instead of the maximum.
    """

    dp = deepcopy(grid_list[primary_grid_idx].dp)
    if region_file is not None:
        region_gdf = gpd.read_file(region_file)

        points_gdf = gpd.GeoDataFrame(
            geometry=gpd.points_from_xy(x=grid_list[0].x, y=grid_list[0].y), crs='epsg:4326')

        points_in_region = gpd.sjoin(
            points_gdf, region_gdf.to_crs('epsg:4326'), how='inner', predicate='within'
        ).index.values
        print(f'{len(points_in_region)} nodes in region.')

        if reverse_sign:
            dp_processed = np.min([gd.dp for gd in grid_list], axis=0)
        else:
            dp_processed = np.max([gd.dp for gd in grid_list], axis=0)
        dp[points_in_region] = dp_processed[points_in_region]
    else:
        dp = np.max([gd.dp for gd in grid_list], axis=0)

    return dp


def direct_replace_dp_in_region(base_grid, replacement_grid, region_file: str):
    """
    Directly replace depth values (dp) in a specified region.

    Within the region defined by region_file, the depth values from
    replacement_grid are used to overwrite those in base_grid.

    However, invalid values in replacement_grid (e.g., -9999) are ignored,
    and the original base_grid values are preserved at those locations.
    """

    dp = deepcopy(base_grid.dp)

    # Read region shapefile
    region_gdf = gpd.read_file(region_file)

    # Create GeoDataFrame of grid points
    points_gdf = gpd.GeoDataFrame(
        geometry=gpd.points_from_xy(x=base_grid.x, y=base_grid.y),
        crs='epsg:4326'
    )

    # Find indices of points inside the region
    points_in_region = gpd.sjoin(
        points_gdf,
        region_gdf.to_crs('epsg:4326'),
        how='inner',
        predicate='within'
    ).index.values

    print(f'{len(points_in_region)} nodes found inside the replacement region.')

    # Define valid replacement values (exclude invalid sentinel values)
    valid_mask = (replacement_grid.dp != -9999) & np.isfinite(replacement_grid.dp)

    # Only replace where both conditions are satisfied:
    # (1) inside the region
    # (2) replacement value is valid
    replace_idx = points_in_region[valid_mask[points_in_region]]

    print(f'{len(replace_idx)} nodes actually replaced (valid values only).')

    # Perform replacement
    dp[replace_idx] = replacement_grid.dp[replace_idx]

    return dp

# ---------------------- Config-driven workflow ----------------------
def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description='Load SCHISM hgrid bathymetry from a YAML/JSON DEM recipe.'
    )
    parser.add_argument(
        '--config',
        type=Path,
        help='YAML or JSON workflow config. Overrides --version.',
    )
    parser.add_argument(
        '--version',
        default=DEFAULT_VERSION,
        help=f'Named YAML recipe in {CONFIG_DIR}. Default: {DEFAULT_VERSION}',
    )
    parser.add_argument(
        '--write-template',
        type=Path,
        help='Write the selected named recipe to a YAML/JSON template and exit.',
    )
    return parser.parse_args()


def named_config_file(version: str) -> Path:
    config_file = CONFIG_DIR / f'{version}.yaml'
    if not config_file.exists():
        versions = ', '.join(sorted(path.stem for path in CONFIG_DIR.glob('*.yaml')))
        raise FileNotFoundError(f"Unknown DEM recipe '{version}'. Available recipes: {versions}")
    return config_file


def load_config(config_file: Path | None, version: str = DEFAULT_VERSION) -> dict[str, Any]:
    if config_file is None:
        config_file = named_config_file(version)

    suffix = config_file.suffix.lower()
    with config_file.open(encoding='utf-8') as stream:
        if suffix == '.json':
            config = json.load(stream)
        elif suffix in {'.yaml', '.yml'}:
            if yaml is None:
                raise ImportError('PyYAML is required to read YAML configs. Use JSON or install yaml.')
            config = yaml.safe_load(stream)
        else:
            raise ValueError(f'Unsupported config format: {config_file}')

    config['_config_file'] = str(config_file)
    return config


def write_template(output_file: Path, config: dict[str, Any]) -> None:
    output_file.parent.mkdir(parents=True, exist_ok=True)
    config = {key: value for key, value in config.items() if not key.startswith('_')}
    if yaml is not None:
        text = yaml.safe_dump(config, sort_keys=False)
    else:
        text = json.dumps(config, indent=2)
    output_file.write_text(text, encoding='utf-8')
    print(f'Wrote template config: {output_file}')


def as_path(path_value: str | Path, base_dir: Path) -> Path:
    path = Path(path_value)
    if path.is_absolute():
        return path
    return base_dir / path


def copy_if_different(source_file: Path, target_file: Path) -> None:
    source_file = source_file.resolve()
    target_file = target_file.resolve()
    if source_file == target_file:
        return
    shutil.copy2(source_file, target_file)


def copy_provenance_files(config: dict[str, Any], workdir: Path) -> None:
    copy_if_different(Path(__file__), workdir / Path(__file__).name)
    if config.get('_config_file') is not None:
        config_file = Path(config['_config_file'])
        copy_if_different(config_file, workdir / config_file.name)

    for source in config['sources']:
        source_file = as_path(source['dem_json'], SCRIPT_DIR)
        copy_if_different(source_file, workdir / source_file.name)

    for rule in config.get('merge_rules', []):
        region = rule.get('region')
        if region is None:
            continue
        region_stem = as_path(region, SCRIPT_DIR).with_suffix('')
        for region_file in SCRIPT_DIR.glob(f'{region_stem.name}.*'):
            copy_if_different(region_file, workdir / region_file.name)


def source_grid(rule_source: str, loaded_grids: dict[str, Any], current_grid: Any) -> Any:
    if rule_source == 'current':
        return current_grid
    return loaded_grids[rule_source]


def run_config(config: dict[str, Any]) -> None:
    comm, _, myrank = initialize_mpi()

    workdir = Path(config['workdir'])
    input_grid = as_path(config.get('input_grid', 'hgrid.ll'), workdir)
    output_grid = as_path(config['output_grid'], workdir)
    dem_dir = Path(config['dem_dir'])
    reverse_sign = bool(config.get('reverse_sign', True))

    if myrank == 0:
        workdir.mkdir(parents=True, exist_ok=True)
        if config.get('copy_provenance', True):
            copy_provenance_files(config, workdir)
        print(f"Running DEM recipe: {config.get('name', 'unnamed')}")
        print(f'workdir: {workdir}')
        print(f'DEM dir: {dem_dir}')
        print(f"sources: {[source['name'] for source in config['sources']]}")

    comm.Barrier()

    loaded_grids = {}
    for source in config['sources']:
        source_name = source['name']
        dem_json = as_path(source['dem_json'], workdir)
        loaded_grids[source_name] = pload_dem(
            grd=str(input_grid),
            grdout=None,
            dem_json=str(dem_json),
            dem_dir=str(dem_dir),
            reverse_sign=reverse_sign,
        )
        comm.Barrier()
        if myrank == 0:
            print(f'---------Loaded grid from {source_name}: {dem_json.name}----------\n')

    if myrank == 0:
        first_source = config['sources'][0]['name']
        current_grid = deepcopy(loaded_grids[first_source])

        for rule in config.get('merge_rules', []):
            method = rule['method']
            region = rule.get('region')
            region_file = None if region is None else str(as_path(region, workdir))
            print(f"Applying merge rule: {rule.get('name', method)} ({method})")

            if method == 'max_depth':
                grids = [source_grid(name, loaded_grids, current_grid) for name in rule['sources']]
                current_grid.dp = max_dp_in_region(grids, region_file=region_file)
            elif method == 'min_depth':
                grids = [source_grid(name, loaded_grids, current_grid) for name in rule['sources']]
                current_grid.dp = max_dp_in_region(grids, region_file=region_file, reverse_sign=True)
            elif method == 'direct_replace':
                replacement = source_grid(rule['replacement'], loaded_grids, current_grid)
                current_grid.dp = direct_replace_dp_in_region(
                    current_grid,
                    replacement,
                    region_file=region_file,
                )
            else:
                raise ValueError(f'Unsupported merge method: {method}')

        current_grid.save(str(output_grid))
        print(f'Wrote DEM-loaded grid: {output_grid}')

    comm.Barrier()


def convert_tif_crs(tif_path, output_path="converted_to_lonlat.tif"):
    """Prints the original CRS of a GeoTIFF file and converts to WGS84 (lon-lat) if needed."""
    with rasterio.open(tif_path) as src:
        original_crs = src.crs
        print("Original CRS:", original_crs)

        # Check if already in lon-lat (WGS84)
        if original_crs == "EPSG:4326":
            print("The file is already in lon-lat (WGS84). No conversion needed.")
            return

        # Convert to WGS84 (EPSG:4326)
        transform, width, height = calculate_default_transform(
            original_crs, "EPSG:4326", src.width, src.height, *src.bounds
        )

        # Create a new dataset with the transformed data
        kwargs = src.meta.copy()
        kwargs.update({
            "crs": "EPSG:4326",
            "transform": transform,
            "width": width,
            "height": height
        })

        with rasterio.open(output_path, "w", **kwargs) as dst:
            for i in range(1, src.count + 1):  # Loop through bands
                reproject(
                    source=rasterio.band(src, i),
                    destination=rasterio.band(dst, i),
                    src_transform=src.transform,
                    src_crs=original_crs,
                    dst_transform=transform,
                    dst_crs="EPSG:4326",
                    resampling=Resampling.nearest
                )

        print(f"Converted file saved as: {output_path}")


def sample_convert_dem():
    """
    Sample converting *.tif to *.npz,
    may need:
        pip install tifffile
        pip install imagecodecs
    """
    # input_tifs = glob('/sciclone/schism10/Hgrid_projects/DEMs/USGS_1M_Stationwide/*.tif')
    # tif_lonlat_files = []
    # for tif_file in input_tifs:
    #     tif_outname = tif_file.replace('.tif', '_lonlat.tif')
    #     convert_tif_crs(tif_file, tif_outname)
    #     tif_lonlat_files.append(tif_outname)

    tif_lonlat_files = glob('/sciclone/schism10/Hgrid_projects/DEMs/USGS_1M_Stationwide/*_lonlat.tif')
    for tif_lonlat_file in tif_lonlat_files:
        convert_dem_format(tif_lonlat_file, sname=tif_lonlat_file.replace('.tif', '.npz'))


def main() -> None:
    args = parse_args()
    config = load_config(args.config, version=args.version)
    if args.write_template is not None:
        write_template(args.write_template, config)
    else:
        run_config(config)


if __name__ == '__main__':
    main()
