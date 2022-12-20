#!/usr/bin/env python3

import shapefile
import matplotlib.pyplot as plt
from pylib import schism_grid, inside_polygon
import numpy as np
from time import time
import os
# import inpoly_operations
from schism_py_pre_post.Grid.SourceSinkIn import source_sink
from schism_py_pre_post.Grid.Hgrid_ported import read_schism_hgrid_cached
from schism_py_pre_post.Timeseries.TimeHistory import TimeHistory
from pathlib import Path
import STOFS3D_scripts


def find_ele_in_shpfile(shapefile_name, hgrid_name: str = ""):
    '''
    Find element index within polygons of a list of shapefiles

    hgrid_name:
        input hgrid to be modified
    '''
    # hgrid
    t = time()
    gd = read_schism_hgrid_cached(hgrid_name, overwrite_cache=False)
    print(f'Reading hgrid took {time()-t} seconds')

    gd.compute_ctr()

    # set default value outside polygons
    t = time()
    # shapefile # records = sf.records() # sf.shapeType # len(sf) # s = sf.shape(idx)
    sf = shapefile.Reader(shapefile_name)
    shapes = sf.shapes()

    ind_list = []
    for shp in shapes:
        poly_xy = np.array(shp.points).T
        ind = inside_polygon(np.c_[gd.xctr, gd.yctr], poly_xy[0], poly_xy[1])  # 1: true; 0: false
        ind = ind.astype('bool')
        ind_list.append(ind)

    print(f'Processing {shapefile_name} took {time()-t} seconds')

    return [ind_list, gd]


def set_constant_sink(wdir='./', shapefile_name='levee_4_pump_polys.shp'):
    # copy datafiles
    mydir = os.path.dirname(STOFS3D_scripts.__file__)
    shp_basename = Path(shapefile_name).stem
    print(f'copying shapefiles from {mydir}/Datafiles/Levee_shapefiles/\n')
    os.system(f"cp {mydir}/Datafiles/Levee_shapefiles/{shp_basename}.* {wdir}/")

    # set pump capacities
    [ele_idx_list, gd] = find_ele_in_shpfile(
        # order matters, latter prevails
        shapefile_name=f'{wdir}/{shapefile_name}',
        hgrid_name=f"{wdir}/hgrid.utm.gr3",
    )
    gd.compute_area()

    pump_sinks = np.zeros((gd.ne, ), dtype=float)
    pump_sink_eles = np.zeros((gd.ne, ), dtype=bool)
    for k, ele_idx in enumerate(ele_idx_list):
        # uniform 0.5 inch/h = 0.00000352777 m/s
        # uni_sinks = np.maximum(-0.00000352777 * gd.area[ele_idx], sinks)
        pump_sinks[ele_idx] = -0.00000352777 * gd.area[ele_idx]
        pump_sink_eles[ele_idx] = True
    
    # set background sink on land
    land = gd.dpe < -2
    background_sink = np.zeros((gd.ne, ), dtype=float)
    background_sink[land] = -0.5/100/3600 * gd.area[land]  # 0.5 cm/h
    const_sinks = np.minimum(background_sink, pump_sinks)

    total_sink_eles = (land + pump_sink_eles).astype(bool)
    # better be a list, for instantiating added_ss
    total_sink_ele_ids = (np.argwhere(total_sink_eles).reshape(-1, ) + 1).tolist()

    # write sink files only
    vsink = TimeHistory(
        data_array=np.r_[np.c_[np.array([0.0]), const_sinks[total_sink_eles].reshape(1, -1)],
                         np.c_[np.array([1e10]), const_sinks[total_sink_eles].reshape(1, -1)]]
    )
    vsink.writer(f'{wdir}/vsink.th')
    np.savetxt(
        f'{wdir}/sink.in',
        np.array([[len(total_sink_ele_ids)] + total_sink_ele_ids]), fmt='%i'
    )

    # assemble source_sink instance
    added_ss = source_sink(
        source_dir=None,
        source_eles=[total_sink_ele_ids[0]],  # dummy source point for constructing added_ss
        sink_eles=total_sink_ele_ids
    )
    for i, col in enumerate(added_ss.vsource.df.columns[1:]):
        added_ss.vsource.df[col] = 0
    const_sinks_actual = const_sinks[total_sink_eles]
    for i, col in enumerate(added_ss.vsink.df.columns[1:]):
        added_ss.vsink.df[col] = const_sinks_actual[i]

    return added_ss

if __name__ == "__main__":
    background_ss = set_constant_sink(wdir='/sciclone/schism10/feiye/STOFS3D-v6/Inputs/I23q/Const_sinks/')
    original_ss = source_sink(source_dir='/sciclone/schism10/feiye/STOFS3D-v6/Inputs/I23q/NWM/')
    total_ss = original_ss + background_ss
    total_ss.writer('/sciclone/schism10/feiye/STOFS3D-v6/Inputs/I23q/Total_ss/')
    pass
    total_ss.writer('/sciclone/schism10/feiye/STOFS3D-v6/Inputs/I23q/Total_ss/')
    

