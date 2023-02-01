# %%
from pylib import schism_grid, grd2sms
import os
from schism_py_pre_post.Grid.Hgrid_extended import read_schism_hgrid_cached
from schism_py_pre_post.Grid.SMS import lonlat2cpp
import pathlib
import copy
import pickle


def set_ocean_bnd(hgrid_name='', gd:schism_grid=None, south_east=[-60.04, 8.56], north_east=[-60.04, 45.82]):
    if gd is None:
        gd = read_schism_hgrid_cached(hgrid_name, overwrite_cache=True)
    gd.compute_bnd(bxy=[south_east[0], north_east[0], south_east[1], north_east[1]])

    if hgrid_name != '':
        gd.save(hgrid_name, fmt=1)

    return gd

def gen_hgrid_formats(hgrid_name='', gd:schism_grid=None, write_bnd=True):
    if gd is None:
        gd = read_schism_hgrid_cached(hgrid_name, overwrite_cache=True)
    else:
        hgrid_name = gd.source_file

    gd.lon, gd.lat = gd.x, gd.y

    dirname = os.path.dirname(hgrid_name)
    file_basename = os.path.basename(hgrid_name)
    file_extension = pathlib.Path(hgrid_name).suffix

    print('outputing hgrid.cpp')
    gd_cpp = copy.deepcopy(gd)
    gd_cpp.x, gd_cpp.y = lonlat2cpp(lon=gd.x, lat=gd.y, lon0=-77.07, lat0=24.0)
    gd_cpp.save(f'{dirname}/hgrid.cpp.gr3', fmt=int(write_bnd))
    os.system(f'mv {dirname}/hgrid.cpp.gr3 {dirname}/hgrid.cpp')
    gd.cpp_x, gd.cpp_y = gd_cpp.x, gd_cpp.y

    print('outputing hgrid in UTM')
    gd.proj(prj0='epsg:4326', prj1='epsg:26918')
    gd.write_hgrid(f'{dirname}/hgrid.utm.26918.gr3', fmt=int(write_bnd))
    gd.utm_x, gd.utm_y = gd.x, gd.y
    print('outputing *.2dm')
    grd2sms(gd, f'{dirname}/hgrid.utm.2dm')

    print('outputing hgrid.102008.gr3')
    gd.x, gd.y = gd.lon, gd.lat
    gd.proj(prj0='epsg:4326', prj1='esri:102008')
    gd.write_hgrid(f'{dirname}/hgrid.102008.gr3', fmt=int(write_bnd))
    gd.x_102008, gd.y_102008 = gd.x, gd.y

    print('saving *.pkl, which has x, y of all needed projections')
    with open(f'{dirname}/hgrids.pkl', 'wb') as file:
        pickle.dump(gd, file)

    # with open(f'{dirname}/hgrids.pkl', 'rb') as file:
    #    gd_test = pickle.load(file)
    print('finish generating hgrids in different projections/formats')
    pass

if __name__ == "__main__":
    # Sample usage
    wdir = '/sciclone/schism10/feiye/STOFS3D-v5/Inputs/v14/Parallel/SMS_proj/v14.42_post_proc2/'
    # Add boundary and generate different hgrid formats
    gd = set_ocean_bnd(f'{wdir}/hgrid.ll')
    gen_hgrid_formats(gd=gd)  # put bnd in the wdir before this step
    pass