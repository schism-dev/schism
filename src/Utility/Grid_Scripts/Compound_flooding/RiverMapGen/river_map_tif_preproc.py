# %%
from difflib import ndiff
from functools import cache
import json
from statistics import median
from osgeo import gdal
from glob import glob
import copy
from schism_py_pre_post.Rivermap.make_river_map import Tif2XYZ, get_all_points_from_shp
import numpy as np
import os
import matplotlib.pyplot as plt
import shapefile
from schism_py_pre_post.Grid.SMS import lonlat2cpp, cpp2lonlat
from pylib import schism_grid
import pickle
import math

def parse_dem_tiles(dem_code, dem_tile_digits):
    if dem_code == 0:
        return [-1]  # no DEM found

    dem_tile_ids = []
    n_tiles = int(math.log10(dem_code)/dem_tile_digits) + 1
    if n_tiles > 4:
        raise Exception("Some thalweg points belong to more than 4 tiles from one DEM source, clean up the DEM tiles first.")
    for digit in reversed(range(n_tiles)):
        x, dem_code = divmod(dem_code, 10**(digit*dem_tile_digits))
        dem_tile_ids.append(int(x-1))
    return dem_tile_ids

def get_tif_boxes(tif_files:list):
    tif_box = []
    for i, tif_file in enumerate(tif_files):
        print(f'reading tifs: {i+1} of {len(tif_files)}, {tif_file}')
        tif = Tif2XYZ(tif_file)
        tif_box.append([min(tif.x), min(tif.y), max(tif.x), max(tif.y)])
    return tif_box

def reproject_tifs(tif_files:list, srcSRS='EPSG:4326', dstSRS='EPSG:26917', outdir='./'):
    '''
    The function failed on CRM tiles, for which use the cmd line tool gdalwarp, e.g.:
    gdalwarp -s_srs EPSG:4326 -t_srs EPSG:26917 -of GTiff ../Lonlat/crm_vol5.tif crm_vol5.26917.tif
    '''
    for i, tif_file in enumerate(tif_files):
        print(f'reprojecting tifs: {i+1} of {len(tif_files)}, {tif_file}')
        epsg = dstSRS.split(':')[1]
        tif_outfile = outdir + os.path.basename(tif_file).split('.')[0] + '.' + epsg + '.tif'
        if not os.path.exists(tif_outfile):
            g = gdal.Warp(tif_outfile, tif_file, srcSRS=srcSRS, dstSRS=dstSRS)
            g = None


def pts_in_box(pts, box):
    in_box = (pts[:, 0] >  box[0]) * (pts[:, 0] <= box[2]) * \
             (pts[:, 1] >  box[1]) * (pts[:, 1] <= box[3])
    return in_box

def find_parent_box(pts, boxes, i_overlap=False):
    ndigits = int(math.log10(len(boxes))) + 1  # number of digits needed for representing tile id, e.g., CuDEM (819 tiles) needs 3 digits
    parent = np.zeros((len(pts), 1), dtype='int')
    digits = np.zeros((len(pts), 1), dtype='int')
    for j, box in enumerate(boxes):
        in_box = pts_in_box(pts[:,:2], box)
        parent[in_box] += ((j+1) * 10.0 ** (digits[in_box])).astype(int)  # save multiple tiles in an integer, e.g.,
                                                                          # 100101 of CuDEM (819 tiles) means tile 100 and tile 101;
                                                                          # 12 of CRM (6 tiles) means tile 1 and tile 2;
        digits[in_box] += ndigits
    # plt.hist(parent, bins=len(np.unique(parent)))
    # np.savetxt('thalweg_parent.xyz', np.c_[pts[:,:2], parent])
    return parent

def tile2dem_file(dem_dict, dem_order, tile_code):
    DEM_id, tile_id = int(tile_code.real), int(tile_code.imag)
    if tile_id != -1:
        return dem_dict[dem_order[DEM_id]]['file_list'][tile_id]
    else:
        return None

def find_thalweg_tile(
    dems_json_file='dems.json',
    thalweg_shp_fname='/sciclone/schism10/feiye/STOFS3D-v5/Inputs/v14/GA_riverstreams_cleaned_utm17N.shp',
    thalweg_buffer=1000,
    cache_folder=None,
    iNoPrint=True, i_DEM_cache=True, i_thalweg_cache=True
):
    '''
    Assign thalwegs to DEM tiles
    '''
    # read DEMs
    with open(dems_json_file) as d:
        dem_dict = json.load(d)

    # get the box of each tile of each DEM
    dem_order = []
    for k, v in dem_dict.items():
        if not iNoPrint: print(f"reading dem: {dem_dict[k]['name']}")
        dem_order.append(k)
        if cache_folder is None:
            cache_folder = os.path.dirname(os.path.abspath(dem_dict[k]['glob_pattern']))  # same as *.shp's folder
        cache_name =  cache_folder + \
            '/' + dem_dict[k]['name'] + '.cache'
        # Remove the existing cache if i_DEM_cache is False,
        # this assumes the cache file needs to be updated
        if (not i_DEM_cache) and os.path.exists(cache_name): os.remove(cache_name)

        if i_DEM_cache and os.path.exists(cache_name):
            if not iNoPrint: print(f"cache read for dem: {dem_dict[k]['name']}")
            with open(cache_name, 'rb') as file:
                tmp_dict = pickle.load(file)
                dem_dict[k]['file_list'] = tmp_dict['file_list']
                dem_dict[k]['boxes'] = tmp_dict['boxes']
        else:
            dem_dict[k]['file_list'] = glob(dem_dict[k]['glob_pattern'])
            dem_dict[k]['boxes'] = get_tif_boxes(dem_dict[k]['file_list'])
            tmp_dict = {
                'file_list': dem_dict[k]['file_list'],
                'boxes': dem_dict[k]['boxes']
            }
            with open(cache_name, 'wb') as file:
                pickle.dump(tmp_dict, file)

    # read thalwegs
    print(f'Reading thalwegs from {thalweg_shp_fname} ...')
    thalweg_shp_fname = thalweg_shp_fname
    xyz, l2g, curv, perp = get_all_points_from_shp(thalweg_shp_fname, iCache=i_thalweg_cache, cache_folder=cache_folder)

    # find DEM tiles for all thalwegs' points
    print(f'finding DEM tiles for each thalweg ...')
    x_cpp, y_cpp = lonlat2cpp(xyz[:, 0], xyz[:, 1])
    xt_right = x_cpp + thalweg_buffer * np.cos(perp)
    yt_right = y_cpp + thalweg_buffer * np.sin(perp)
    xt_left = x_cpp + thalweg_buffer * np.cos(perp + np.pi)
    yt_left = y_cpp + thalweg_buffer * np.sin(perp + np.pi)
    # find thalweg itself and two search boundaries (one on each side)
    thalwegs2dems = [find_parent_box(xyz[:,:2], dem_dict[k]['boxes']) for k in dem_dict.keys()]
    thalwegs_right2dems = [find_parent_box(np.array(cpp2lonlat(xt_right, yt_right)).T, dem_dict[k]['boxes']) for k in dem_dict.keys()]
    thalwegs_left2dems = [find_parent_box(np.array(cpp2lonlat(xt_left, yt_left)).T, dem_dict[k]['boxes']) for k in dem_dict.keys()]

    # how many digits in tile numbers
    dems_tile_digits = [int(math.log10(len(dem_dict[k]['boxes'])))+1 for k in dem_dict.keys()]

    # use cpp projection hereafter, because meter unit is easier for parameterization
    xyz[:, 0], xyz[:, 1] = x_cpp, y_cpp

    thalwegs = []
    thalwegs_parents = []
    for i, idx in enumerate(l2g):  # enumerate thalwegs
        line = xyz[idx,:]  # one segment of a thalweg
        thalwegs.append(line)

        thalweg_parents = []  # one thalweg can have parent tiles from all DEM sources
        for i_dem, [thalwegs2dem, thalwegs_left2dem, thalwegs_right2dem] in enumerate(zip(thalwegs2dems, thalwegs_right2dems, thalwegs_left2dems)):
            # find all DEM tiles that a thalweg (including its left and right search boundaries) touches
            thalweg2dem = np.unique(np.r_[thalwegs2dem[idx], thalwegs_left2dem[idx], thalwegs_right2dem[idx]]).tolist()
            for dem_code in thalweg2dem:
                thalweg_parents += [complex(i_dem, x) for x in parse_dem_tiles(dem_code, dems_tile_digits[i_dem])]
                pass
            # thalweg_parents += [complex(i_dem, x) for x in thalweg2dem]  # real part is DEM id; complex part is tile id
        thalwegs_parents.append(thalweg_parents)

    # Group thalwegs: thalwegs from the same group have the same parent tiles
    print(f'grouping thalwegs ...')
    groups = []
    group_id = 0
    thalweg2group = -np.ones((len(thalwegs)), dtype=int)
    for i, thalweg_parents in enumerate(thalwegs_parents):
        if thalweg_parents not in groups:
            groups.append(thalweg_parents)
            thalweg2group[i] = group_id
            group_id += 1
        else:
            for j, x in enumerate(groups):
                if x == thalweg_parents:
                    thalweg2group[i] = j

    # reduce groups: merge smaller groups into larger groups
    groups = np.array(groups, dtype=object)
    ngroup = len(groups)
    grp2large_grp = ngroup * np.ones((ngroup+1,), dtype=int)  # add a dummy mapping at the end
    for i1, group1 in enumerate(groups):
        for i2, group2 in enumerate(groups):
            if len(group1) < len(group2) and all(elem in group2 for elem in group1):
                # print(f'{group1} is contained in {group2}')
                grp2large_grp[i1] = i2
                break

    # But some large groups are still contained in larger groups,
    # get to the bottom of the family tree (e.g., parent's parent's parent ...)
    parents = np.squeeze(grp2large_grp[grp2large_grp])  # parent's parent
    idx = parents != len(groups)  # where parent's parent exists
    while any(idx):
        grp2large_grp[idx] = parents[idx]  # reset parent to parent's parent
        parents = parents[parents]  # advance family tree
        idx = parents != len(groups)  # get the idx where parent's parent still exists

    idx = grp2large_grp==len(groups)  # where parent's parent is no-existent
    grp2large_grp[idx] = np.arange(len(groups)+1)[idx]  # parent group is self
    grp2large_grp = grp2large_grp[:-1]  # remove the dummy group at the end

    large_groups = groups[np.unique(grp2large_grp)]
    if not iNoPrint: print(f'number of groups after reduction: {len(large_groups)}')
    group_lens = [len(x) for x in large_groups]
    if not iNoPrint: print(f'group lengths: min {min(group_lens)}; max {max(group_lens)}; mean {np.mean(group_lens)}')

    thalweg2group = grp2large_grp[thalweg2group]
    map_grp = dict(zip(np.unique(grp2large_grp), np.arange(len(np.unique(grp2large_grp)))))
    thalweg2large_group = np.array([map_grp[x] for x in thalweg2group])

    large_group2thalwegs = [[] for _ in range(len(large_groups))]
    for i, x in enumerate(thalweg2large_group):
        large_group2thalwegs[x].append(i)

    large_groups_files = copy.deepcopy(large_groups)
    for i, group in enumerate(large_groups):
        for j, tile_code in enumerate(group):
            large_groups_files[i][j] = tile2dem_file(dem_dict=dem_dict, dem_order=dem_order, tile_code=tile_code)

    # histogram
    # plt.hist(thalweg2large_group, bins=len(np.unique(thalweg2large_group)))
    # plt.show()

    return thalweg2large_group, large_groups_files, np.array(large_group2thalwegs, dtype=object)

if __name__ == "__main__":
    # find_thalweg_tile()
    # %%
    # Reproject
    # tif_files = glob(f'/sciclone/schism10/feiye/STOFS3D-v5/Inputs/v14/GA_parallel/CRM/Lonlat/*.tif')
    # reproject_tifs(tif_files, 'EPSG:26917', outdir='/sciclone/schism10/feiye/STOFS3D-v5/Inputs/v14/GA_parallel/CRM/UTM17/')

    # Merge small coned tiles into larger ones (similar to CuDEM's tile size)
    # cudem_files = glob(f'/sciclone/schism10/feiye/STOFS3D-v5/Inputs/v14/GA_parallel/CuDEM/*.tif')
    # cudem_boxes = get_tif_boxes(cudem_files)

    # coned_files = glob(f'/sciclone/schism10/feiye/STOFS3D-v5/Inputs/v14/GA_parallel/CoNED/Original/*.tif')
    # coned_boxes = get_tif_boxes(coned_files)

    # coned_centers = np.c_[
    #     (np.array(coned_boxes)[:, 0]+np.array(coned_boxes)[:, 2])/2,
    #     (np.array(coned_boxes)[:, 1]+np.array(coned_boxes)[:, 3])/2,
    # ]

    # for i, cudem_box in enumerate(cudem_boxes):
    #     in_box = (coned_centers[:, 0] >  cudem_box[0]) * \
    #              (coned_centers[:, 0] <= cudem_box[2]) * \
    #              (coned_centers[:, 1] >  cudem_box[1]) * \
    #              (coned_centers[:, 1] <= cudem_box[3])
    #     in_box_files = (np.array(coned_files)[in_box]).tolist()
    #     g = gdal.Warp(f"/sciclone/schism10/feiye/STOFS3D-v5/Inputs/v14/GA_parallel/CoNED/Combined/GA_CoNED_merged_{i}.tif",
    #                   in_box_files, format="GTiff", options=["COMPRESS=LZW", "TILED=YES"])
    #     g = None # Close file and flush to disk

    pass

