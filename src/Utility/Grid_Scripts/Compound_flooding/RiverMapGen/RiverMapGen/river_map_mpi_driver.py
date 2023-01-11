import os
import time
from mpi4py import MPI
import glob
import numpy as np
import pickle
from RiverMapGen.river_map_tif_preproc import find_thalweg_tile
from RiverMapGen.make_river_map import make_river_map
from RiverMapGen.SMS import merge_maps


comm = MPI.COMM_WORLD
rank = comm.Get_rank()
size = comm.Get_size()

def my_mpi_idx(N, size, rank):
    my_idx = np.zeros((N, ), dtype=bool)
    n_per_rank, _ = divmod(N, size)
    n_per_rank = n_per_rank + 1
    # return slice())
    my_idx[rank*n_per_rank:min((rank+1)*n_per_rank, N)] = True
    return my_idx

def river_map_mpi_driver(
    dems_json_file = './dems.json',  # files for all DEM tiles
    thalweg_shp_fname='',
    output_dir = './',
    thalweg_buffer = 1000,
    cache_folder = './Cache/',
    i_DEM_cache = True,
    i_thalweg_cache = False,
    i_grouping_cache = False
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
    cache_folder: Used to store temporary variables, so that the next run with similar inputs can be executed faster
    i_DEM_cache : Whether or not to read DEM info from cache.
                  The cache file saves box info of all DEM tiles.
                  Reading from original tiles can be slow, so the default option is True
    i_thalweg_cache: Whether or not to read thalweg info from cache.
                     The cache file saves coordinates, index, curvature, and direction at all thalweg points
                     This is usually fast even without reading cache.
    i_grouping_cache: Whether or not to read grouping info from cache,
                      which is useful when the same DEMs and thalweg_shp_fname are used.
                      A cache file named "dems_json_file + thalweg_shp_fname_grouping.cache" will be saved regardless of the option value.
                      This is usually fast even without reading cache.
    '''
    if rank == 0: print(f'A total of {size} core(s) used.')
    comm.barrier()

    if rank == 0:
        if os.path.exists(output_dir):
            os.system(f"rm -r {output_dir}")
        else:
            os.makedirs(output_dir, exist_ok=True)
    comm.barrier()

    # group thalwegs (with the option of cache)
    cache_name = cache_folder + \
        os.path.basename(dems_json_file) + '_' + \
        os.path.basename(thalweg_shp_fname) + '_grouping.cache'
    # Remove the existing cache if i_grouping_cache is False,
    # this assumes the cache file needs to be updated
    if (not i_grouping_cache) and os.path.exists(cache_name): os.remove(cache_name)
    # core 1 calculates the grouping, then saves a cache if the cache file does not exist
    if rank == 0:
        if not os.path.exists(cache_name):
            thalwegs2tile_groups, tile_groups_files, tile_groups2thalwegs = \
                find_thalweg_tile(
                    dems_json_file=dems_json_file,
                    thalweg_shp_fname=thalweg_shp_fname,
                    thalweg_buffer = 1000,
                    iNoPrint=bool(rank), # only rank 0 prints to screen
                    cache_folder=cache_folder,
                    i_DEM_cache=i_DEM_cache, i_thalweg_cache=i_thalweg_cache
                )
            with open(cache_name, 'wb') as file:
                pickle.dump([thalwegs2tile_groups, tile_groups_files, tile_groups2thalwegs], file)
    comm.barrier()
    # all cores read from cache
    with open(cache_name, 'rb') as file:
        print(f'Reading grouping info from cache ...')
        thalwegs2tile_groups, tile_groups_files, tile_groups2thalwegs = pickle.load(file)
    if rank == 0:
        print(f'Thalwegs are divided into {len(tile_groups2thalwegs)} groups.')
        for i, tile_group2thalwegs in enumerate(tile_groups2thalwegs):
            print(f'[ Group {i+1} ]-----------------------------------------------------------------------\n' + \
                  f'Group {i+1} includes the following thalwegs (idx starts from 0): {tile_group2thalwegs}\n' + \
                  f'Group {i+1} needs the following DEMs: {tile_groups_files[i]}\n')
    comm.barrier()

    # assign groups to each core
    my_idx = my_mpi_idx(N=len(tile_groups_files), size=size, rank=rank)
    # my_idx.fill(False); my_idx[[50]] = True
    my_tile_groups = tile_groups_files[my_idx]
    my_tile_groups_thalwegs = tile_groups2thalwegs[my_idx]
    print(f'Rank {rank} handles Group {np.squeeze(np.argwhere(my_idx))}\n')
    comm.Barrier()

    time_all_groups_start = time.time()

    # each core handles its assigned groups sequentially
    for i, (my_tile_group, my_tile_group_thalwegs) in enumerate(zip(my_tile_groups, my_tile_groups_thalwegs)):
        time_this_group_start = time.time()
        make_river_map(
            tif_fnames = my_tile_group,
            thalweg_shp_fname = thalweg_shp_fname,
            thalweg_smooth_shp_fname = None,  # '/GA_riverstreams_cleaned_corrected_utm17N.shp'
            selected_thalweg = my_tile_group_thalwegs,
            cache_folder=cache_folder,
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
        merge_maps(f'{output_dir}/Rank*_total_arcs.map', merged_fname=f'{output_dir}/total_arcs.map')
        merge_maps(f'{output_dir}/*intersection_joints*.map', merged_fname=f'{output_dir}/total_intersection_joints.map')
        merge_maps(f'{output_dir}/*inner_arcs.map', merged_fname=f'{output_dir}/total_inner_arcs.map')
        merge_maps(f'{output_dir}/*centerlines.map', merged_fname=f'{output_dir}/total_centerlines.map')
        merge_maps(f'{output_dir}/*bank_final*.map', merged_fname=f'{output_dir}/total_banks_final.map')

        xyz_files = glob.glob(f'{output_dir}/*.xyz')
        os.system(f'cat {" ".join(xyz_files)} > {output_dir}/intersection_res.xyz')

        print(f'Total run time: {time.time()-time_all_groups_start} seconds.')

if __name__ == "__main__":

    # ------------------------- sample input ---------------------------
    dems_json_file = './Inputs/DEMs/dems.json'  # specifying files for all DEM tiles
    thalweg_shp_fname='./Inputs/Shapefiles/LA_local.shp'
    output_dir = './Outputs/' +  f'{os.path.basename(thalweg_shp_fname).split(".")[0]}_{size}cores/'
    # ------------------------- end input section ---------------------------

    river_map_mpi_driver(
        dems_json_file=dems_json_file,
        thalweg_shp_fname=thalweg_shp_fname,
        output_dir=output_dir,
    )
