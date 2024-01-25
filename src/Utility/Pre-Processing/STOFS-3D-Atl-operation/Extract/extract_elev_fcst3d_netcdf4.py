'''
Author: L. Cui (lcui@vims.edu)
Usage: python extract_slab_fcst_netcdf4.py YYYYMMDD

July 22, 2021:
    -Changed netcdf format from NETCDF4 to NETCDF3_CLASSIC to
     use MFDataset
'''

import sys
from datetime import datetime
from time import time 
import argparse

import numpy as np
from netCDF4 import Dataset

if __name__ == '__main__':

    t0=time()

    argparser = argparse.ArgumentParser()
#   argparser.add_argument('date', type=datetime.fromisoformat, help='input file date')
    argparser.add_argument('date',  help='input file date')
    args=argparser.parse_args()

    #date=datetime(2021, 6, 7)
    date=args.date 
#    fpath = "/sciclone/schism10/hyu05/NOAA_NWM/oper_3D/fcst"
    fpath = "./"
#   ds=Dataset(f"{fpath}/{date.strftime('%Y%m%d')}/schout_{date.strftime('%Y%m%d')}.nc")
    ds=Dataset(f"{fpath}/outputs/out2d_{date}.nc")
    #units=ds['time'].units
    #base_date=ds['time'].base_date
    #print(f"{fpath}/{date}/schout_{date}.nc")
    #ds=Dataset(f"{fpath}/{date}/schout_{date}.nc")

    #get coordinates
    x=ds['SCHISM_hgrid_node_x'][:]
    y=ds['SCHISM_hgrid_node_y'][:]
    depth=ds['depth'][:]
    elements=ds['SCHISM_hgrid_face_nodes'][:,:]

    NP=len(depth)
  
    #get elements and split quads into tris and make index start from 0
    elements=ds['SCHISM_hgrid_face_nodes'][:,:]
    tris = []
    for ele in elements:
        ele=np.ma.masked_values(ele, -1)
        ele=ele[~ele.mask]
        #breakpoint()

        if len(ele) == 3:
            tris.append([ele[0]-1, ele[1]-1, ele[2]-1])
        elif len(ele) == 4:
            tris.append([ele[0]-1, ele[1]-1, ele[3]-1])
            tris.append([ele[1]-1, ele[2]-1, ele[3]-1])
    NE=len(tris)
    print(f'NE is {NE}')
    NV=3

    #get wetdry nodes
    #wd_nodes=ds['wetdry_node'][:,:]
    wd_nodes=ds['dryFlagNode'][:,:]
    #elev=ds['elev'][:,:]
    elev=ds['elevation'][:,:]

    #get times
    times=ds['time'][:]
    print(times)
    ntimes=len(times)

    #elev[np.where(elev>10000)]=-99999

    #outdir= '/sciclone/pscr/lcui01/ICOGS3D_dev'
    outdir = './extract'
    #with Dataset(f"{outdir}/schout_elev_fcst3D_{date}.nc", "w", format="NETCDF3_CLASSIC") as fout:
    with Dataset(f"{outdir}/schout_elev_fcst3D_{date}.nc", "w", format="NETCDF4") as fout:
        #dimensions
        fout.createDimension('nSCHISM_hgrid_node', NP)
        fout.createDimension('nSCHISM_hgrid_face', NE)
        fout.createDimension('nMaxSCHISM_hgrid_face_nodes', NV)
        fout.createDimension('time', None)

        #variables
        fout.createVariable('time', 'f', ('time',))
        fout['time'].long_name="Time"
        #fout['time'].units = units #f'seconds since {date.year}-{date.month}-{date.day} 00:00:00 +0000'
        #fout['time'].base_date= base_date #(date.year, date.month, date.day, 0)
        fout['time'].standard_name="time"
        fout['time'][:] = times

        fout.createVariable('SCHISM_hgrid_face_nodes', 'i', ('nSCHISM_hgrid_face', 'nMaxSCHISM_hgrid_face_nodes',))
        fout['SCHISM_hgrid_face_nodes'].long_name="Horizontal Element Table"
        fout['SCHISM_hgrid_face_nodes'].cf_role="face_node_connectivity"
        fout['SCHISM_hgrid_face_nodes'].start_index=1
        fout['SCHISM_hgrid_face_nodes'][:]=np.array(tris)

        fout.createVariable('SCHISM_hgrid_node_x', 'f', ('nSCHISM_hgrid_node',))
        fout['SCHISM_hgrid_node_x'].long_name="node x-coordinate"
        fout['SCHISM_hgrid_node_x'].standard_name="longitude"
        fout['SCHISM_hgrid_node_x'].units="degrees_east"
        fout['SCHISM_hgrid_node_x'].mesh="SCHISM_hgrid"
        fout['SCHISM_hgrid_node_x'][:]=x

        fout.createVariable('SCHISM_hgrid_node_y', 'f', ('nSCHISM_hgrid_node',))
        fout['SCHISM_hgrid_node_y'].long_name="node y-coordinate"
        fout['SCHISM_hgrid_node_y'].standard_name="latitude"
        fout['SCHISM_hgrid_node_y'].units="degrees_north"
        fout['SCHISM_hgrid_node_y'].mesh="SCHISM_hgrid"
        fout['SCHISM_hgrid_node_y'][:]=y

        #fout.createVariable('depth', 'f', ('nSCHISM_hgrid_node',), fill_value=-99999.)
        fout.createVariable('depth', 'f', ('nSCHISM_hgrid_node',))
        fout['depth'].long_name="bathymetry"
        fout['depth'].units="m"
        fout['depth'].mesh="SCHISM_hgrid"
        fout['depth'][:]=depth

        #fout.createVariable('elev', 'f', ('time', 'nSCHISM_hgrid_node',), fill_value=-99999.)
        fout.createVariable('elev', 'f', ('time', 'nSCHISM_hgrid_node',))
        fout['elev'].long_name="water elevation"
        fout['elev'].units="m"
        fout['elev'].mesh="SCHISM_hgrid"
        fout['elev'][:]=elev

        fout.title = 'SCHISM Model output'
        fout.source = 'SCHISM model output version v10'
        fout.references = 'http://ccrm.vims.edu/schismweb/'

    print(f'It took {time()-t0} to interpolate')
