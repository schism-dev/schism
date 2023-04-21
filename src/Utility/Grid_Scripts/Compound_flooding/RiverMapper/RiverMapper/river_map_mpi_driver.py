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
import pickle
import geopandas as gpd
from shapely.ops import polygonize
from RiverMapper.river_map_tif_preproc import find_thalweg_tile, Tif2XYZ
from RiverMapper.make_river_map import make_river_map, clean_intersections, geos2SmsArcList, Geoms_XY, clean_arcs
from RiverMapper.SMS import merge_maps, SMS_MAP
from RiverMapper.util import silentremove


def my_mpi_idx(N, size, rank):
    '''
    Distribute N tasks to {size} ranks.
    The return value is a bool vector of the shape (N, ),
    with True indices indicating tasks for the current rank.
    '''
    i_my_groups = np.zeros((N, ), dtype=bool)
    groups = np.array_split(range(N), size)  # n_per_rank, _ = divmod(N, size)
    my_group_ids = groups[rank]
    i_my_groups[my_group_ids] = True
    return my_group_ids, i_my_groups

def merge_outputs(output_dir):
    print(f'\n------------------ merging outputs from all cores --------------\n')
    time_merge_start = time.time()

    # sms maps
    total_arcs_map = merge_maps(f'{output_dir}/*_total_arcs.map', merged_fname=f'{output_dir}/total_arcs.map')

    total_intersection_joints = merge_maps(f'{output_dir}/*intersection_joints*.map', merged_fname=f'{output_dir}/total_intersection_joints.map').detached_nodes

    total_river_map = merge_maps(f'{output_dir}/*river_arcs.map', merged_fname=f'{output_dir}/total_river_arcs.map')
    if total_river_map is not None:
        total_river_arcs = total_river_map.arcs
    else:
        total_river_arcs = None

    total_centerlines = merge_maps(f'{output_dir}/*centerlines.map', merged_fname=f'{output_dir}/total_centerlines.map')
    merge_maps(f'{output_dir}/*bank_final*.map', merged_fname=f'{output_dir}/total_banks_final.map')

    # # shapefiles
    river_outline_files = glob(f'{output_dir}/*_river_outline.shp')
    if len(river_outline_files) > 0:
        gpd.pd.concat([gpd.read_file(x).to_crs('epsg:4326') for x in river_outline_files]).to_file(f'{output_dir}/total_river_outline.shp')

    bomb_polygon_files = glob(f'{output_dir}/*_bomb_polygons.shp')
    if len(bomb_polygon_files) > 0:
        gpd.pd.concat([gpd.read_file(x).to_crs('epsg:4326') for x in bomb_polygon_files]).to_file(f'{output_dir}/total_bomb_polygons.shp')

    print(f'Merging outputs took: {time.time()-time_merge_start} seconds.')
    return [total_arcs_map, total_intersection_joints, total_river_arcs, total_centerlines]


def final_clean_up(output_dir, total_arcs_map, snap_points, i_blast_intersection=False, total_river_arcs=None):
    pass

def river_map_mpi_driver(
    dems_json_file = './dems.json',  # files for all DEM tiles
    thalweg_shp_fname='',
    output_dir = './',
    thalweg_buffer = 1000,
    i_DEM_cache = True,
    i_blast_intersection = False,
    comm = MPI.COMM_WORLD
):
    '''
    Driver for the parallel execution of make_river_map.py

    Thalwegs are grouped based on the DEM tiles associated with each thalweg.
    For each thalweg, its associated DEM tiles are those needed for determining
    the elevations on all thalweg points, as well as
    the elevations within a buffer zone along the thalweg
    (within which the positions of left and right banks will be searched)

    One core can be responsible for one or more thalweg groups,
    which are fed to make_river_map.py one at a time

    Summary of the input parameters:
    thalweg_buffer: in meters. This is the search range on either side of the thalweg.
                    Because banks will be searched within this range,
                    its value is needed now to identify parent DEM tiles of each thalweg
    i_DEM_cache : Whether or not to read DEM info from cache.
                  Reading from original *.tif files can be a little slower than reading from the *.pkl cache,
                  so the default option is True
    '''

    # deprecated (fast enough without caching)
    i_thalweg_cache = False  # Whether or not to read thalweg info from cache.
                             # The cache file saves coordinates, index, curvature, and direction at all thalweg points
    # i_grouping_cache: Whether or not to read grouping info from cache,
    #                   which is useful when the same DEMs and thalweg_shp_fname are used.
    #                   A cache file named "dems_json_file + thalweg_shp_fname_grouping.cache" will be saved regardless of the option value.
    #                   This is usually fast even without reading cache.
    i_grouping_cache = True; iValidateCache = False
    cache_folder = './Cache/'

    rank = comm.Get_rank()
    size = comm.Get_size()

    if rank == 0:
        time_start = time_grouping_start = time.time()
        print('\n---------------------------------grouping thalwegs---------------------------------\n')

    comm.Barrier()

    thalwegs2tile_groups, tile_groups_files, tile_groups2thalwegs = None, None, None
    if rank == 0:
        print(f'A total of {size} core(s) used.')
        silentremove(output_dir)
        os.makedirs(output_dir, exist_ok=True)

        if i_grouping_cache:
            os.makedirs(cache_folder, exist_ok=True)
            cache_name = cache_folder + \
                os.path.basename(dems_json_file) + '_' + \
                os.path.basename(thalweg_shp_fname) + '_grouping.cache'
            try:
                with open(cache_name, 'rb') as file:
                    print(f'Reading grouping info from cache ...')
                    thalwegs2tile_groups, tile_groups_files, tile_groups2thalwegs = pickle.load(file)
            except FileNotFoundError:
                print(f"Grouping cache does not exist at {cache_folder}. Cache will be generated after grouping.")

        if thalwegs2tile_groups is None:
            thalwegs2tile_groups, tile_groups_files, tile_groups2thalwegs = find_thalweg_tile(
                dems_json_file=dems_json_file,
                thalweg_shp_fname=thalweg_shp_fname,
                thalweg_buffer = thalweg_buffer,
                iNoPrint=bool(rank), # only rank 0 prints to screen
                i_thalweg_cache=i_thalweg_cache
            )
            if i_grouping_cache:
                with open(cache_name, 'wb') as file:
                    pickle.dump([thalwegs2tile_groups, tile_groups_files, tile_groups2thalwegs], file)

    thalwegs2tile_groups = comm.bcast(thalwegs2tile_groups, root=0)
    tile_groups_files = comm.bcast(tile_groups_files, root=0)
    tile_groups2thalwegs = comm.bcast(tile_groups2thalwegs, root=0)

    if rank == 0:
        print(f'Thalwegs are divided into {len(tile_groups2thalwegs)} groups.')
        for i, tile_group2thalwegs in enumerate(tile_groups2thalwegs):
            print(f'[ Group {i+1} ]-----------------------------------------------------------------------\n' + \
                  f'Group {i+1} includes the following thalwegs (idx starts from 0): {tile_group2thalwegs}\n' + \
                  f'Group {i+1} needs the following DEMs: {tile_groups_files[i]}\n')
        print(f'Grouping took: {time.time()-time_grouping_start} seconds')

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

        if iValidateCache:
            for tif_fname in unique_tile_files[my_mpi_idx(len(unique_tile_files), size, rank)]:
                _, is_new_cache = Tif2XYZ(tif_fname=tif_fname)
                if is_new_cache:
                    print(f'[Rank: {rank} cached DEM {tif_fname}')
                else:
                    print(f'[Rank: {rank} validated existing cache for {tif_fname}')

    comm.Barrier()
    if rank == 0: print('\n---------------------------------assign groups to each core---------------------------------\n')
    comm.Barrier()

    my_group_ids, i_my_groups = my_mpi_idx(N=len(tile_groups_files), size=size, rank=rank)
    my_tile_groups = tile_groups_files[i_my_groups]
    my_tile_groups_thalwegs = tile_groups2thalwegs[i_my_groups]
    print(f'Rank {rank} handles Group {np.squeeze(np.argwhere(i_my_groups))}\n')

    comm.Barrier()
    if rank == 0: print('\n---------------------------------beginning map generation---------------------------------\n')
    comm.Barrier()
    time_all_groups_start = time.time()

    for i, (my_group_id, my_tile_group, my_tile_group_thalwegs) in enumerate(zip(my_group_ids, my_tile_groups, my_tile_groups_thalwegs)):
        time_this_group_start = time.time()
        print(f'Rank {rank}: Group {i+1} (global: {my_group_id}) started ...')
        if True:  # my_group_id in range(0, 200):
            make_river_map(
                tif_fnames = my_tile_group,
                thalweg_shp_fname = thalweg_shp_fname,
                selected_thalweg = my_tile_group_thalwegs,
                output_dir = output_dir,
                output_prefix = f'Group_{my_group_id}_{rank}_{i}_',
                mpi_print_prefix = f'[Rank {rank}, Group {i+1} of {len(my_tile_groups)}, global: {my_group_id}] ',
                i_blast_intersection=i_blast_intersection
            )
        else:
            pass  # print(f'Rank {rank}: Group {my_group_id} failed')

        print(f'Rank {rank}: Group {i+1} (global: {my_group_id}) run time: {time.time()-time_this_group_start} seconds.')

    print(f'Rank {rank}: total run time: {time.time()-time_all_groups_start} seconds.')

    comm.Barrier()

    # finalize
    if rank == 0:
        # merge outputs from all ranks
        total_arcs_map, total_intersection_joints, total_river_arcs, total_centerlines = merge_outputs(output_dir)

        print(f'\n--------------- final clean-ups --------------------------------------------------------\n')
        time_final_cleanup_start = time.time()

        total_arcs_cleaned = [arc for arc in total_arcs_map.to_GeoDataFrame().geometry.unary_union.geoms]
        if not i_blast_intersection:
            bomb_polygons = gpd.read_file(f'{output_dir}/total_bomb_polygons.shp')
            total_arcs_cleaned = clean_intersections(arcs=total_arcs_cleaned, target_polygons=bomb_polygons, snap_points=total_intersection_joints)
        total_arcs_cleaned = clean_arcs(total_arcs_cleaned)
        SMS_MAP(arcs=geos2SmsArcList(total_arcs_cleaned)).writer(filename=f'{output_dir}/total_arcs.map')

        total_arcs_cleaned_polys = [poly for poly in polygonize(gpd.GeoSeries(total_arcs_cleaned))]
        gpd.GeoDataFrame(
            index=range(len(total_arcs_cleaned_polys)), crs='epsg:4326', geometry=total_arcs_cleaned_polys
        ).to_file(filename=f'{output_dir}/total_arcs.shp', driver="ESRI Shapefile")

        # river_arcs_cleaned = clean_river_arcs(total_river_arcs, total_arcs_cleaned)
        # total_river_outline_polys = generate_river_outline_polys(river_arcs_cleaned)
        # if len(total_river_outline_polys) > 0:
        #     gpd.GeoDataFrame(
        #         index=range(len(total_river_outline_polys)), crs='epsg:4326', geometry=total_river_outline_polys
        #     ).to_file(filename=f'{output_dir}/total_river_outline.shp', driver="ESRI Shapefile")
        # else:
        #     print(f'Warning: total_river_outline_polys empty')

        print(f'Final clean-ups took: {time.time()-time_final_cleanup_start} seconds.')
    
        # delete per-core outputs
        silentremove(glob(f'{output_dir}/Group*'))
        print(f'>>>>>>>> Total run time: {time.time()-time_start} seconds >>>>>>>>')
