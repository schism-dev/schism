"""
This script provides a driver for grouping thalwegs based on their parent DEM tiles,
then allocate groups to mpi cores,
and finally calls the function "make_river_map" to sequentially process each group on each core

Usage: 
Import and call the function "river_map_mpi_driver",
see sample_parallel.py in the installation directory.
"""


import os
import time
from mpi4py import MPI
from glob import glob
import numpy as np
import geopandas as gpd
from shapely.ops import polygonize
from RiverMapper.river_map_tif_preproc import find_thalweg_tile, Tif2XYZ
from RiverMapper.make_river_map import make_river_map, clean_intersections, geos2SmsArcList
from RiverMapper.SMS import merge_maps, SMS_MAP
from RiverMapper.util import silentremove


def my_mpi_idx(N, size, rank):
    my_idx = np.zeros((N, ), dtype=bool)
    n_per_rank, _ = divmod(N, size)
    n_per_rank = n_per_rank + 1
    my_idx[rank*n_per_rank:min((rank+1)*n_per_rank, N)] = True
    return my_idx

def river_map_mpi_driver(
    dems_json_file = './dems.json',  # files for all DEM tiles
    thalweg_shp_fname='',
    output_dir = './',
    thalweg_buffer = 1000,
    i_DEM_cache = True,
    comm = MPI.COMM_WORLD
):
    '''
    Driver for the parallel execution of make_river_map.py

    Thalwegs are grouped based on the DEM tiles associated with each thalweg.
    For each thalweg, its associated DEM tiles are those needed for determining
    the elevations on all thalweg points, as well as
    the elevations within a buffer zone of the thalweg (within which left and right banks will be sought)

    One core can be responsible for one or more thalweg groups,
    which are fed to make_river_map.py one at a time

    Summary of the input parameters:
    thalweg_buffer: unit in meters. This is the search range on either side of the thalweg.
                    Because banks will be searched within this range,
                    its value is needed now to associate DEM tiles with each thalweg
    i_DEM_cache : Whether or not to read DEM info from cache.
                  Reading from original *.tif files can be slow, so the default option is True
    '''

    # deprecated (fast enough without caching)
    i_thalweg_cache = False  # Whether or not to read thalweg info from cache.
                             # The cache file saves coordinates, index, curvature, and direction at all thalweg points

    rank = comm.Get_rank()
    size = comm.Get_size()

    if rank == 0: print('\n---------------------------------grouping thalwegs---------------------------------\n')
    comm.Barrier()

    if rank == 0:
        print(f'A total of {size} core(s) used.')
        silentremove(output_dir)
        os.makedirs(output_dir, exist_ok=True)

        # group thalwegs (with the option of cache)
        thalwegs2tile_groups, tile_groups_files, tile_groups2thalwegs = find_thalweg_tile(
            dems_json_file=dems_json_file,
            thalweg_shp_fname=thalweg_shp_fname,
            thalweg_buffer = thalweg_buffer,
            iNoPrint=bool(rank), # only rank 0 prints to screen
            i_thalweg_cache=i_thalweg_cache
        )
    else:
        thalwegs2tile_groups, tile_groups_files, tile_groups2thalwegs = None, None, None

    thalwegs2tile_groups = comm.bcast(thalwegs2tile_groups, root=0)
    tile_groups_files = comm.bcast(tile_groups_files, root=0)
    tile_groups2thalwegs = comm.bcast(tile_groups2thalwegs, root=0)

    if rank == 0:
        print(f'Thalwegs are divided into {len(tile_groups2thalwegs)} groups.')
        for i, tile_group2thalwegs in enumerate(tile_groups2thalwegs):
            print(f'[ Group {i+1} ]-----------------------------------------------------------------------\n' + \
                  f'Group {i+1} includes the following thalwegs (idx starts from 0): {tile_group2thalwegs}\n' + \
                  f'Group {i+1} needs the following DEMs: {tile_groups_files[i]}\n')

    comm.barrier()
    if rank == 0: print('\n---------------------------------caching DEM tiles---------------------------------\n')
    comm.barrier()

    if i_DEM_cache:
        unique_tile_files = []
        for group in tile_groups_files:
            for file in group:
                if (file not in unique_tile_files) and (file is not None):
                    unique_tile_files.append(file)
        unique_tile_files = np.array(unique_tile_files)

        for tif_fname in unique_tile_files[my_mpi_idx(len(unique_tile_files), size, rank)]:
            _, is_new_cache = Tif2XYZ(tif_fname=tif_fname)
            if is_new_cache:
                print(f'[Rank: {rank} cached DEM {tif_fname}.')
            else:
                print(f'[Rank: {rank} read DEM from cache {tif_fname}.')

    comm.Barrier()
    if rank == 0: print('\n---------------------------------assign groups to each core---------------------------------\n')
    comm.Barrier()

    my_idx = my_mpi_idx(N=len(tile_groups_files), size=size, rank=rank)
    my_tile_groups = tile_groups_files[my_idx]
    my_tile_groups_thalwegs = tile_groups2thalwegs[my_idx]
    print(f'Rank {rank} handles Group {np.squeeze(np.argwhere(my_idx))}\n')

    comm.Barrier()
    if rank == 0: print('\n---------------------------------beginning map generation---------------------------------\n')
    comm.Barrier()
    time_all_groups_start = time.time()

    # each core handles its assigned groups sequentially
    for i, (my_tile_group, my_tile_group_thalwegs) in enumerate(zip(my_tile_groups, my_tile_groups_thalwegs)):
        time_this_group_start = time.time()
        make_river_map(
            tif_fnames = my_tile_group,
            thalweg_shp_fname = thalweg_shp_fname,
            selected_thalweg = my_tile_group_thalwegs,
            output_dir = output_dir,
            output_prefix = f'{rank}_{i}_',
            mpi_print_prefix = f'[Rank {rank}, Group {i+1} of {len(my_tile_groups)}] ',
        )

        print(f'Rank {rank}: Group {i+1} run time: {time.time()-time_this_group_start} seconds.')

    print(f'Rank {rank}: total run time: {time.time()-time_all_groups_start} seconds.')

    comm.Barrier()

    # merge outputs on each core (for testing purposes)
    merge_maps(f'{output_dir}/{rank}_*_total_arcs.map', merged_fname=f'{output_dir}/Rank{rank}_total_arcs.map')
    comm.Barrier()

    # merge outputs from all ranks
    if rank == 0:
        print(f'\n------------------ merging outputs from all cores --------------\n')
        # sms maps
        merge_maps(f'{output_dir}/Rank*_total_arcs.map', merged_fname=f'{output_dir}/total_arcs.map')
        merge_maps(f'{output_dir}/*intersection_joints*.map', merged_fname=f'{output_dir}/total_intersection_joints.map')
        merge_maps(f'{output_dir}/*river_arcs.map', merged_fname=f'{output_dir}/total_river_arcs.map')
        merge_maps(f'{output_dir}/*centerlines.map', merged_fname=f'{output_dir}/total_centerlines.map')
        merge_maps(f'{output_dir}/*bank_final*.map', merged_fname=f'{output_dir}/total_banks_final.map')
        # shapefiles
        gpd.pd.concat([gpd.read_file(x).to_crs('epsg:4326') for x in glob(f'{output_dir}/*river_outline*.shp')]).to_file(f'{output_dir}/total_river_outline.shp')
        gpd.pd.concat([gpd.read_file(x).to_crs('epsg:4326') for x in glob(f'{output_dir}/*bomb*.shp')]).to_file(f'{output_dir}/total_bomb_polygons.shp')
        
        # 
        print(f'\n--------------- final clean-ups on intersections near inter-subdomain interfaces ----\n')
        total_arcs_cleaned = clean_intersections(
            lines=SMS_MAP(filename=f'{output_dir}/total_arcs.map').to_GeoDataFrame(),
            target_polygons=gpd.read_file(f'{output_dir}/total_bomb_polygons.shp'),
            snap_points=SMS_MAP(filename=f'{output_dir}/total_intersection_joints.map').detached_nodes
        )
        SMS_MAP(arcs=geos2SmsArcList(total_arcs_cleaned)).writer(filename=f'{output_dir}/total_arcs.map')

        total_arcs_cleaned_polys = [poly for poly in polygonize(gpd.GeoSeries(total_arcs_cleaned))]
        gpd.GeoDataFrame(
            index=range(len(total_arcs_cleaned_polys)), crs='epsg:4326', geometry=total_arcs_cleaned_polys
        ).to_file(filename=f'{output_dir}/total_river_arc_polygons.shp', driver="ESRI Shapefile")

        print(f'>>>>>>>> Total run time: {time.time()-time_all_groups_start} seconds >>>>>>>>')
