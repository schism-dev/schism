#!/usr/bin/env python3
'''
Load bathymetry for SCHISM hgrid

Recommended usage:
Run this script in the SCHISM source code folder where it originally resides.
Select a suitable pre-configuration (e.g., stofs3d_v8()) or start a new one,
set parameters at the beginning of the function, then call it in main()

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
import socket
import re
import time
import json
import errno
from pathlib import Path
from copy import deepcopy
from glob import glob

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


# ---------------------- Sample Usage ----------------------
def sample_max_dp_usage():
    '''Sample usage of the max_dp_in_region function.'''

    wdir = '/sciclone/schism10/feiye/STOFS3D-v7/Inputs/v7_test/Bathy_edit/DEM_loading/'
    gd1 = schism_read(f'{wdir}/Original/hgrid.ll.dem_loaded.mpi.gr3')
    gd2 = schism_read(f'{wdir}/BlueTopo/hgrid.ll.dem_loaded.mpi.gr3')

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
        loaded_grids = [schism_read(f) for f in loaded_grid_fnames]
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
    wdir = '/sciclone/schism10/feiye/STOFS3D-v8/I13b/Bathy_edit/DEM_loading/'
    dem_dir = '/sciclone/schism10/Hgrid_projects/DEMs/npz2/'
    dem_json_list = [
        'DEM_info_original.json',
        'DEM_info_with_bluetopo.json',
    ]
    dem_region_shpfile = 'v18_s2_v1_polys_dissolved.shp'
    output_fname = f'{wdir}/hgrid.ll.dem_loaded.mpi.gr3'
    # ---------------------------------------

    comm, _, myrank = initialize_mpi()
    if myrank == 0:  # copy this script to the working directory to keep a record
        os.system(f'cp {__file__} {wdir}')
        for dem_json_file in dem_json_list:
            os.system(f'cp {dem_json_file} {wdir}')
        os.system(f'cp {dem_region_shpfile} {wdir}')

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


def stofs3d_v8_LA():
    """
    Load bathymetry for STOFS3D-v8 Louisiana.

    Two DEM sets are used to load bathymetry from different sources: v2a and v3a.
        DEM_info.json v2a: v1 (DEM_info_original.json, used for STOFS-3D-Atl v6) +
                           Bluetopo, same as DEM_info_with_BlueTopo.json

        DEM_info.json v3a: v1 + NGOM_CoNED_2022 (mostly LA, but including Mobile Bay, AL)¬

        Take the maximum depth from the two grids

        Note: since the first grid is the primary grid by default,
        Bluetopo prevails in the entire domain
    """
    # dem_dir = '/work2/noaa/nos-surge/feiye/npz2/'
    wdir = '/sciclone/schism10/feiye/STOFS3D-v8/I20/Bathy_edit/DEM_loading/'
    dem_dir = '/sciclone/schism10/Hgrid_projects/DEMs/npz2/'
    dem_json_list = ['DEM_info_v2a.json', 'DEM_info_v3a.json']
    dem_region_shpfile = 'v18_s2_v1_polys_dissolved.shp'
    output_fname = f'{wdir}/hgrid.ll.dem_loaded.mpi.gr3'
    # ---------------------------------------

    comm, _, myrank = initialize_mpi()
    if myrank == 0:  # copy this script to the working directory to keep a record
        os.system(f'cp {__file__} {wdir}')
        for dem_json_file in dem_json_list:
            os.system(f'cp {dem_json_file} {wdir}')
        os.system(f'cp {dem_region_shpfile} {wdir}')

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


def stofs3d_v8():
    """
    Load bathymetry for STOFS3D-v8.

    Needs hgrid.ll in the working directory.

    Three DEM sets are used to load bathymetry from different sources.
    The first one corresponds to v6, the second one includes BlueTopo.
    The third one includes NGOM_CoNED_2022.
    This leads to three bathymetry-loaded grids.

    A region (Louisiana) is defined where the v6 grid seems too shallow.
    The v8 grid depth takes the larger value from the three grids inside region,
    i.e., where BlueTopo and NGOM leads to deeper channels than v6.
    And it takes the v6 value outside the region.
    """
    # ----------- inputs -------------------
    wdir = '/sciclone/schism10/feiye/STOFS3D-v8/I15a_v7/Bathy_edit/DEM_loading/'
    dem_dir = '/sciclone/schism10/Hgrid_projects/DEMs/npz2/'
    dem_json_list = [
        'DEM_info_original_patched.json',
        'DEM_info_with_bluetopo.json',
        'DEM_info_v3a.json',  # v3a has CoNED 2022 NGOM for LA
        'DEM_info_Statewide.json',  # , v4 added USGS 1M Statewide for CT and RI
    ]
    max_dem_region_shpfile = 'bluetopo_regions4.shp'
    min_dem_region_shpfile = 'breakwaters_poly.shp'
    output_fname = f'{wdir}/hgrid.ll.dem_loaded.mpi.gr3'
    # ---------------------------------------

    comm, _, myrank = initialize_mpi()
    if myrank == 0:  # copy this script to the working directory to keep a record
        prep_dir(wdir, dem_json_list, [max_dem_region_shpfile, min_dem_region_shpfile])
        print(f'preparing files in {wdir}, DEM list: {dem_json_list}')

    comm.Barrier()

    loaded_grids = []
    for dem_json in dem_json_list:
        # Load grids in parallel
        loaded_grids.append(
            pload_dem(grd=f'{wdir}/hgrid.ll', grdout=None, dem_json=f'{wdir}/{dem_json}',
                      dem_dir=dem_dir, reverse_sign=True)  # returns None for non-root
        )  # On non-root processes, loaded_grids only contains None
        comm.Barrier()  # wait for all cores to finish populating loaded_grids
        if myrank == 0:
            print(f'---------Loaded grid from {dem_json}----------\n')

    # On root process, take the maximum depth from loaded_grids
    if myrank == 0:
        print(f'processing loaded grids: {loaded_grids}')
        hgrid_final = deepcopy(loaded_grids[0])
        hgrid_final.dp = max_dp_in_region(
            loaded_grids[:-1], region_file=f'{wdir}/{max_dem_region_shpfile}')
        hgrid_final.dp = max_dp_in_region(
            [hgrid_final, loaded_grids[-1]],
            region_file=f'{wdir}/{min_dem_region_shpfile}',
            reverse_sign=True)

        hgrid_final.save(output_fname)

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


if __name__ == '__main__':
    # sample_convert_dem()
    stofs3d_v8()
