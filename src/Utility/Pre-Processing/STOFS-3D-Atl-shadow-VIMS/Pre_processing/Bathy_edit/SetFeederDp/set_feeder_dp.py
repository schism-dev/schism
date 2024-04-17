from email.mime import base
from schism_py_pre_post.Grid.SourceSinkIn import source_sink, SourceSinkIn
from schism_py_pre_post.Grid.SMS import lonlat2cpp
from schism_py_pre_post.Timeseries.TimeHistory import TimeHistory
from schism_py_pre_post.Grid.Hgrid_extended import read_schism_hgrid_cached
import os
import numpy as np
from scipy import spatial
import pickle
import copy
from pylib_essentials.schism_file import sms2grd


def nearest_neighbour(points_a, points_b):
    tree = spatial.cKDTree(points_b)
    return np.array(tree.query(points_a)[1]).reshape(-1,), np.array(tree.query(points_a)[0]).reshape(-1,)

def dist(points_group_A, points_group_B):
    points_A = np.squeeze(points_group_A.view(np.complex128))
    points_B = np.squeeze(points_group_B.view(np.complex128))
    return np.absolute(points_A-points_B)

def set_feeder_dp(feeder_info_dir='', hgrid_obj=None, hgrid_obj_no_feeder=None):
    # Read feeder channel info
    with open(f'{feeder_info_dir}/feeder.pkl', 'rb') as file:
        [feeder_l2g, feeder_points, feeder_heads, feeder_bases] = pickle.load(file)

    # pair feeder channel points to grid points
    f2g, _ = nearest_neighbour(feeder_points, np.c_[hgrid_obj.x, hgrid_obj.y])
    # pair feeder channel bases to grid points
    f2g_base, _ = nearest_neighbour(feeder_bases[:, :2], np.c_[hgrid_obj.x, hgrid_obj.y])

    # find outside grid points
    feeder_in_grid = hgrid_obj_no_feeder.inside_grid(feeder_points).astype(bool)

    for i, id in enumerate(feeder_l2g):
        gd_points_in_feeder = f2g[id]
        gd_points_in_feeder_in_grid = feeder_in_grid[id]
        gd_points_in_external_feeder = gd_points_in_feeder[~gd_points_in_feeder_in_grid]

        # Option 1, use min depth (highest z)
        # hgrid_obj.dp[gd_points_in_external_feeder] = np.min(hgrid_obj.dp[gd_points_in_feeder])

        # Option 2, use the depth of feeder bases (near the interface between land boundary and feeder channel)
        base_point_in_grid = f2g_base[i]
        hgrid_obj.dp[gd_points_in_external_feeder] = hgrid_obj.dp[base_point_in_grid]

    return hgrid_obj

if __name__ == "__main__":
    # sample usage
    wdir ='./'

    # the current grid
    gd0 = read_schism_hgrid_cached(f'{wdir}/hgrid.before_feeder_dp.ll', overwrite_cache=False)
    # gd.save(f'{wdir}/hgrid.copy.ll', fmt=1)

    # a grid without feeder is needed to identify which feeder points are outside and should be deepened
    # Only the boundary matters, the interior of the grid doesn't matter,
    # so if you don't have a grid without feeders, you can just generate a simplified grid with the lbnd_ocean map
    gd_no_feeder = sms2grd('/sciclone/schism10/Hgrid_projects/STOFS3D-v7/v20.0/no_feeder.2dm')
    gd_no_feeder.proj(prj0='esri:102008', prj1='epsg:4326')

    gd = set_feeder_dp(
        feeder_info_dir='/sciclone/schism10/Hgrid_projects/STOFS3D-v7/v20.0/Feeder/',
        hgrid_obj=gd0, hgrid_obj_no_feeder=gd_no_feeder
    )
    dp_diff = gd.dp - gd0.dp
    print(f'min deepening: {min(dp_diff)}; max deepening: {max(dp_diff)} \n')
    gd.save(f'{wdir}/hgrid.ll', fmt=0)
