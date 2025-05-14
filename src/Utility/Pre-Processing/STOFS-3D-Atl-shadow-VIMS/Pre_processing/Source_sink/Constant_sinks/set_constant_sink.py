#!/usr/bin/env python3
import shapefile
import numpy as np
from time import time
from pathlib import Path

from pylib_essentials.utility_functions import inside_polygon
from pylib_essentials.schism_file import schism_grid
from pylib_essentials.schism_file import source_sink, TimeHistory


def find_ele_in_shpfile(shapefile_name, hgrid:schism_grid):
    '''
    Find element index within polygons of a list of shapefiles
    '''
    hgrid.compute_ctr()

    # set default value outside polygons
    t = time()
    sf = shapefile.Reader(shapefile_name)
    shapes = sf.shapes()

    ind_list = []
    for shp in shapes:
        poly_xy = np.array(shp.points).T
        ind = inside_polygon(np.c_[hgrid.xctr, hgrid.yctr], poly_xy[0], poly_xy[1])  # 1: true; 0: false
        ind = ind.astype('bool')
        ind_list.append(ind)

    print(f'Processing {shapefile_name} took {time()-t} seconds')

    return ind_list


def set_constant_sink(wdir='./', shapefile_name='levee_4_pump_polys.shp', hgrid_utm:schism_grid=None):
    # copy datafiles
    gd = hgrid_utm

    # set pump capacities
    ele_idx_list = find_ele_in_shpfile(
        # order matters, latter prevails
        shapefile_name=f'{wdir}/{shapefile_name}',
        hgrid=gd,
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
    # index starts from 1
    total_sink_ele_ids = (np.argwhere(total_sink_eles).reshape(-1, ) + 1).tolist()
    # index starts from 0
    total_sink_ele_idx = np.argwhere(total_sink_eles).reshape(-1, )

    # build vsink.th
    vsink = TimeHistory(
        data_array=np.r_[np.c_[np.array([0.0]), const_sinks[total_sink_eles].reshape(1, -1)],
                         np.c_[np.array([100*365*86400]), const_sinks[total_sink_eles].reshape(1, -1)]],
        columns=total_sink_ele_ids,
    )

    # save sink ele ids and sink values for operational use
    np.savetxt(f'{wdir}/pump_sinks.txt', np.c_[total_sink_ele_ids, const_sinks[total_sink_eles]], fmt='%i %10.5f')
    # save as xyz format
    gd.compute_ctr()
    my_ele_xyz = np.c_[gd.xctr[total_sink_ele_idx], gd.yctr[total_sink_ele_idx], const_sinks[total_sink_eles]]
    np.savetxt(f'{wdir}/sinks.xyz', my_ele_xyz)

    # build a source_sink object
    const_source_sink = source_sink(vsource=None, vsink=vsink, msource=None)
    # write source sink files, actually only vsink.th and source_sink.in
    const_source_sink.writer(output_dir=wdir)

    return const_source_sink

if __name__ == "__main__":
    my_ss = source_sink.from_files('/sciclone/schism10/feiye/STOFS3D-v8/I202503/Source_sink/')
    eles = np.array(my_ss.source_eles.astype(str))
    # read a json file
    import json
    json_file = Path('/sciclone/schism10/feiye/STOFS3D-v8/I202503/Source_sink/sources.json')
    with open(json_file, 'r') as f:
        ele2fid = json.load(f)
    eles1 = np.array(list(ele2fid.keys()))

    np.all(np.equal(eles, eles1))

    print('Done')


