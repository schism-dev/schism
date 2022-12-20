#!/usr/bin/env python3

import shapefile
import matplotlib.pyplot as plt
from pylib import schism_grid, inside_polygon
import numpy as np
from time import time
import os
# import inpoly_operations
from schism_py_pre_post.Grid.SourceSinkIn import source_sink
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
    gr3_pkl = f'{os.path.splitext(hgrid_name)[0]}.pkl'
    if os.path.exists(gr3_pkl):
        gd = schism_grid(gr3_pkl)
    else:
        gd = schism_grid(hgrid_name)
        gd.save(gr3_pkl)
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


def set_pump(wdir='./', original_sourc_sink_dir='./original_ss/', shapefile_name='levee_4_pump_polys.shp'):
    # copy datafiles
    mydir = os.path.dirname(STOFS3D_scripts.__file__)
    shp_basename = Path(shapefile_name).stem
    print(f'copying shapefiles from {mydir}/Datafiles/Levee_shapefiles/\n')
    os.system(f"cp {mydir}/Datafiles/Levee_shapefiles/{shp_basename}.* {wdir}/")

    # read pump capacities
    sf = shapefile.Reader(shapefile_name)
    pump_capacities = np.array([x[-1] for x in sf.records()]) * 0.028316847

    [ele_idx_list, gd] = find_ele_in_shpfile(
        # order matters, latter prevails
        shapefile_name=shapefile_name,
        hgrid_name=f"{wdir}/hgrid.utm.gr3",
    )
    gd.compute_area()

    sink_eles_pump = np.array([])
    sink_pump = np.array([])
    total_ss = source_sink(source_dir=original_sourc_sink_dir)
    for k, [ele_idx, pump_capacity] in enumerate(zip(ele_idx_list, pump_capacities)):
        area_weight = gd.area[ele_idx]/sum(gd.area[ele_idx])
        sinks = - pump_capacity * area_weight
        added_sink_eles = np.where(ele_idx)[0] + 1  # ele_idx starts from 0

        # uniform 0.5 inch/h = 0.00000352777 m/s
        # uni_sinks = np.maximum(-0.00000352777 * gd.area[ele_idx], sinks)
        uni_sinks = -0.00000352777 * gd.area[ele_idx]
        plt.plot(sinks, label='sinks')
        plt.plot(uni_sinks, label='uni')
        plt.legend()
        # plt.show()
        sinks = uni_sinks

        added_ss = source_sink(
            source_dir=None,
            source_eles=[added_sink_eles[0]],  # dummy source point for constructing added_ss
            sink_eles=added_sink_eles.tolist()
        )
        sink_eles_pump = np.r_[sink_eles_pump, added_sink_eles]
        sink_pump = np.r_[sink_pump, sinks]

        for i, col in enumerate(added_ss.vsource.df.columns[1:]):
            added_ss.vsource.df[col] = 0
        for i, col in enumerate(added_ss.vsink.df.columns[1:]):
            added_ss.vsink.df[col] = sinks[i]
        # added_ss.writer(dirname=f"{wdir}/Test{k}/")

        total_ss = total_ss + added_ss
        if np.isnan(total_ss.vsink.data).any():
            raise Exception('nan found in sink')
    total_ss.writer(dirname=wdir)
    np.savetxt(f'{wdir}/pump_sinks.txt', np.c_[sink_eles_pump, sink_pump], fmt='%i %10.5f')
    pass

if __name__ == "__main__":
    set_pump()

#!/usr/bin/env python3

import shapefile
import matplotlib.pyplot as plt
from pylib import schism_grid, inside_polygon
import numpy as np
from time import time
import os
# import inpoly_operations
from schism_py_pre_post.Grid.SourceSinkIn import source_sink
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
    gr3_pkl = f'{os.path.splitext(hgrid_name)[0]}.pkl'
    if os.path.exists(gr3_pkl):
        gd = schism_grid(gr3_pkl)
    else:
        gd = schism_grid(hgrid_name)
        gd.save(gr3_pkl)
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


def set_pump(wdir='./', original_sourc_sink_dir='./original_ss/', shapefile_name='levee_4_pump_polys.shp'):
    # copy datafiles
    mydir = os.path.dirname(STOFS3D_scripts.__file__)
    shp_basename = Path(shapefile_name).stem
    print(f'copying shapefiles from {mydir}/Datafiles/Levee_shapefiles/\n')
    os.system(f"cp {mydir}/Datafiles/Levee_shapefiles/{shp_basename}.* {wdir}/")

    # read pump capacities
    sf = shapefile.Reader(shapefile_name)
    pump_capacities = np.array([x[-1] for x in sf.records()]) * 0.028316847

    [ele_idx_list, gd] = find_ele_in_shpfile(
        # order matters, latter prevails
        shapefile_name=shapefile_name,
        hgrid_name=f"{wdir}/hgrid.utm.gr3",
    )
    gd.compute_area()

    sink_eles_pump = np.array([])
    sink_pump = np.array([])
    total_ss = source_sink(source_dir=original_sourc_sink_dir)
    for k, [ele_idx, pump_capacity] in enumerate(zip(ele_idx_list, pump_capacities)):
        area_weight = gd.area[ele_idx]/sum(gd.area[ele_idx])
        sinks = - pump_capacity * area_weight
        added_sink_eles = np.where(ele_idx)[0] + 1  # ele_idx starts from 0

        # uniform 0.5 inch/h = 0.00000352777 m/s
        # uni_sinks = np.maximum(-0.00000352777 * gd.area[ele_idx], sinks)
        uni_sinks = -0.00000352777 * gd.area[ele_idx]
        plt.plot(sinks, label='sinks')
        plt.plot(uni_sinks, label='uni')
        plt.legend()
        # plt.show()
        sinks = uni_sinks

        added_ss = source_sink(
            source_dir=None,
            source_eles=[added_sink_eles[0]],  # dummy source point for constructing added_ss
            sink_eles=added_sink_eles.tolist()
        )
        sink_eles_pump = np.r_[sink_eles_pump, added_sink_eles]
        sink_pump = np.r_[sink_pump, sinks]

        for i, col in enumerate(added_ss.vsource.df.columns[1:]):
            added_ss.vsource.df[col] = 0
        for i, col in enumerate(added_ss.vsink.df.columns[1:]):
            added_ss.vsink.df[col] = sinks[i]
        # added_ss.writer(dirname=f"{wdir}/Test{k}/")

        total_ss = total_ss + added_ss
        if np.isnan(total_ss.vsink.data).any():
            raise Exception('nan found in sink')
    total_ss.writer(dirname=wdir)
    np.savetxt(f'{wdir}/pump_sinks.txt', np.c_[sink_eles_pump, sink_pump], fmt='%i %10.5f')
    pass

if __name__ == "__main__":
    set_pump()
