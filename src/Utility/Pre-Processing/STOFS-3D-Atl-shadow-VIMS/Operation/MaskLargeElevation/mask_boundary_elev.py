import xarray
import numpy as np
from pylib import schism_grid, grd2sms
from schism_py_pre_post.Geometry.inpoly import find_node_in_shpfiles
import os
from netCDF4 import Dataset

def make_mask(hg: schism_grid, high_ground_thres=-19.0, additional_mask_file=None):
    
    hg.compute_bnd()
    hg.compute_node_ball()
    # hg.save('/sciclone/schism10/feiye/STOFS3D-v4/Inputs/v4_20220715_update/hgrid.pkl')

    # find nodes near lbnd
    lbnd_nd_idx = np.hstack(hg.ilbn)
    ilbnd_nd = np.zeros((hg.np), dtype=bool)
    ilbnd_nd[lbnd_nd_idx] = True

    lbnd_side_idx = np.sum(ilbnd_nd[hg.isidenode], axis=1).astype(bool)

    ilbnd_ele = np.unique(np.reshape(hg.isdel[lbnd_side_idx], (-1, )))
    ilbnd_ele = ilbnd_ele[ilbnd_ele>=0]

    lbnd_nd_2_tier_idx = np.unique(np.reshape(hg.elnode[ilbnd_ele], (-1, )))
    lbnd_nd_2_tier_idx = lbnd_nd_2_tier_idx[lbnd_nd_2_tier_idx>=0].astype(int)
    ilbnd_nd_2_tier = np.zeros((hg.np), dtype=bool)
    ilbnd_nd_2_tier[lbnd_nd_2_tier_idx] = True

    # find high ground
    i_high_ground = hg.dp < high_ground_thres

    # get additional mask points
    if additional_mask_file is None:
        i_additional = np.zeros((hg.np), dtype=bool)
    else:
        i_additional = find_node_in_shpfiles(shapefile_names=[additional_mask_file], gd=hg)

    imask = i_high_ground + ilbnd_nd_2_tier + i_additional

    return imask

def overwrite_mask_nc(imask, outputfile):
    ds = xarray.open_dataset(outputfile)
    ds['idmask'][:] = imask.astype(int)
    ds.to_netcdf('tmp.nc')
    os.system(f'mv tmp.nc {outputfile}')
    ds.close()

def write_mask_nc(imask, outputfile):
    np = len(imask)
    with Dataset(outputfile, "w", format="NETCDF4") as fout:
        #dimensions
        fout.createDimension('nSCHISM_hgrid_node', np)

        #variables
        fout.createVariable('idmask', 'i', ('nSCHISM_hgrid_node',))
        fout['idmask'].long_name="Node mask"
        fout['idmask'][:]=imask.astype(int)

def diagnose_mask(imask: np.ndarray):
    # diagnostic outputs
    maxelev = schism_grid('/sciclone/schism10/feiye/STOFS3D-v4/Shared_with_NOAA/Mask_lbnd/maxelev.gr3')
    maxelev.dp[imask] = np.nan
    grd2sms(grd=maxelev, sms='/sciclone/schism10/feiye/STOFS3D-v4/Shared_with_NOAA/Mask_lbnd/maxelev.mask.2dm')

    hg = schism_grid('/sciclone/schism10/feiye/STOFS3D-v4/Inputs/v4_20220715_update/hgrid.ll')
    ds = xarray.open_dataset('/sciclone/schism10/feiye/STOFS3D-v4/Shared_with_NOAA/Mask_lbnd/out2d_1.mask.nc')
    hg.dp[:] = ds['elevation'][-1, :]
    ds.close()
    grd2sms(grd=hg, sms='/sciclone/schism10/feiye/STOFS3D-v4/Shared_with_NOAA/Mask_lbnd/out2d_1.mask.2dm')


if __name__ == '__main__':
    hg = schism_grid('/sciclone/schism10/feiye/STOFS3D-v4/Inputs/v4_20220715_update/hgrid.pkl')

    imask = make_mask(hg, high_ground_thres=-19.0, additional_mask_file='./Shapefiles/additional_masks.shp')
    write_mask_nc(imask, outputfile='/sciclone/schism10/feiye/STOFS3D-v4/Shared_with_NOAA/Mask_lbnd/lbnd_node_mask1.nc')
    