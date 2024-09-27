#!/usr/bin/env python3
'''
Load bathymetry for SCHISM hgrid

Needs to prepare the following files/folder:
   - hgrid.ll: hgrid in lon/lat

   - DEM tiles: pre-processed from *.tif and *.asc to python's *.npz format for faster reading.

   - DEM*.json: specifying DEM tiles used, and arranged from lower priority to higher priority.
         See samples in this folder.

   - region (optional): A polygon shapefile specifiying regions where special treatment is applied.
        See sample procedure for STOFS-3D v7 in "def stofs_v7" in "pload_dem2.py".
        Sample *.shp is provided in this folder.


 See sample procedures in "def stofs3d_v6" and "def stofs3d_v7".
'''


import os
import sys
import re
import time
import json
import errno
from copy import deepcopy

import numpy as np
from mpi4py import MPI
import geopandas as gpd
from pylib import load_bathymetry, zdata


# Attempt to import experimental grid reader for speed-up
try:
    from pylib_experimental.schism_file import cread_schism_hgrid as read_hgrid
except ImportError:
    from pylib import schism_grid as read_hgrid


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
        gd = read_hgrid(grd)
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


def max_dp_in_region(grid_list: list, region_file: str, primary_grid_idx: int = 0):
    """
    In a specified region, copy the depth of the first grid in a list
    and change it to the maximum depth of all grids.
    If region_file is None, the operation applies to the entire grid.
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

        dp_max = np.max([gd.dp for gd in grid_list], axis=0)
        dp[points_in_region] = dp_max[points_in_region]
    else:
        dp = np.max([gd.dp for gd in grid_list], axis=0)

    return dp


# ---------------------- Sample Usage ----------------------
def sample_max_dp_usage():
    '''Sample usage of the max_dp_in_region function.'''

    wdir = '/sciclone/schism10/feiye/STOFS3D-v7/Inputs/v7_test/Bathy_edit/DEM_loading/'
    gd1 = read_hgrid(f'{wdir}/Original/hgrid.ll.dem_loaded.mpi.gr3')
    gd2 = read_hgrid(f'{wdir}/BlueTopo/hgrid.ll.dem_loaded.mpi.gr3')

    dp = max_dp_in_region([gd1, gd2], region_file=f'{wdir}/v18_s2_v1_polys_dissolved.shp')

    gd1.save(f'{wdir}/hgrid_max_dp.gr3', value=dp)


def stofs3d_v6():
    """
    Sample usage of the pload_dem function, corresponding to STOFS-3D v6
    """
    wdir = '/sciclone/schism10/feiye/STOFS3D-v7/Inputs/v7_test/Bathy_edit/DEM_loading/'
    pload_dem(
        grd=f'{wdir}/hgrid.ll',
        grdout=f'{wdir}/hgrid.ll.dem_loaded.mpi.gr3',
        dem_json=f'{wdir}/DEM_info_original.json',
        dem_dir='/sciclone/schism10/Hgrid_projects/DEMs/npz2/',
        reverse_sign=True
    )


def stofs3d_v7_hercules():
    """
    Similar to stofs3d_v7 below, but with minor tweaks 
    to resolve some unknown issues on Hercules.
    """
    # ----------- inputs -------------------
    wdir = '/work/noaa/nosofs/feiye/STOFS-3D-Atl-Example-Setup/DEM_loading_example/'
    dem_dir = '/work2/noaa/nos-surge/feiye/npz2/'
    dem_json_list = [
        f'{wdir}/DEM_info_original.json',
        f'{wdir}/DEM_info_with_bluetopo.json',
    ]
    dem_region_shpfile = f'{wdir}/v18_s2_v1_polys_dissolved.shp'
    output_fname = f'{wdir}/hgrid.ll.dem_loaded.mpi.gr3'
    # ---------------------------------------

    comm, _, myrank = initialize_mpi()

    loaded_grid_fnames = []
    for dem_json in dem_json_list:
        grdout = f'{dem_json}.gr3'
        loaded_grid_fnames.append(grdout)
        # Load grids in parallel
        pload_dem(grd=f'{wdir}/hgrid.ll', grdout=grdout, dem_json=dem_json,
                  dem_dir=dem_dir, reverse_sign=True)

        comm.Barrier()
        if myrank == 0:
            print(f'---------Loaded grid from {dem_json}----------\n')

    # On root process, take the maximum depth from loaded_grids
    if myrank == 0:
        print(f'Taking the maximum depth from loaded grids: {loaded_grid_fnames}')
        loaded_grids = [read_hgrid(f) for f in loaded_grid_fnames]
        dp = max_dp_in_region(loaded_grids, region_file=dem_region_shpfile)
        loaded_grids[0].save(output_fname, value=dp)

    comm.Barrier()


def stofs3d_v7():
    """
    Load bathymetry for STOFS3D-v7.

    Two DEM sets are used to load bathymetry from different sources.
    The first one corresponds to v6, the second one includes BlueTopo.
    This leads to two bathymetry-loaded grids.

    A region (Louisiana) is defined where the v6 grid seems too shallow.
    The v7 grid depth takes the larger value from the two grids inside region,
    i.e., where BlueTopo leads to deeper channels than v6.
    And it takes the v6 value outside the region.
    """
    # ----------- inputs -------------------
    # wdir = '/work/noaa/nosofs/feiye/STOFS-3D-Atl-Example-Setup/DEM_loading_example/'
    # dem_dir = '/work2/noaa/nos-surge/feiye/npz2/'
    wdir = '/sciclone/schism10/feiye/STOFS3D-v7/Inputs/v7_test/Bathy_edit/DEM_loading/'
    dem_dir = '/sciclone/schism10/Hgrid_projects/DEMs/npz2/'
    dem_json_list = [
        f'{wdir}/DEM_info_original.json',
        f'{wdir}/DEM_info_with_bluetopo.json',
    ]
    dem_region_shpfile = f'{wdir}/v18_s2_v1_polys_dissolved.shp'
    output_fname = f'{wdir}/hgrid.ll.dem_loaded.mpi.gr3'
    # ---------------------------------------

    comm, _, myrank = initialize_mpi()

    loaded_grids = []
    for dem_json in dem_json_list:
        # Load grids in parallel
        loaded_grids.append(
            pload_dem(grd=f'{wdir}/hgrid.ll', grdout=None, dem_json=dem_json,
                      dem_dir=dem_dir, reverse_sign=True)  # returns None for non-root
        )  # On non-root processes, loaded_grids only contains None
        comm.Barrier()  # wait for all cores to finish populating loaded_grids
        if myrank == 0:
            print(f'---------Loaded grid from {dem_json}----------\n')

    # On root process, take the maximum depth from loaded_grids
    if myrank == 0:
        print(f'Taking the maximum depth from loaded grids: {loaded_grids}')
        dp = max_dp_in_region(loaded_grids, region_file=dem_region_shpfile)
        loaded_grids[0].save(output_fname, value=dp)

    comm.Barrier()


if __name__ == '__main__':
    stofs3d_v7_hercules()
