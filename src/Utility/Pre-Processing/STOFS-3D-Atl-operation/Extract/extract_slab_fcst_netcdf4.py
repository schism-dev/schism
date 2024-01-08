import sys
from datetime import datetime
from time import time 
import argparse

import numpy as np
import numpy.ma as ma
from netCDF4 import Dataset

from generate_adcirc import split_quads  # modified by FY

def get_zcor_interp_coefficient(zcor, zinter, kbp):
    '''
    inputs:
        -zcor: zcor[np,nvrt] for each time record.
        -zinter: zinter[np] depth where currents will be interpolated at
        -kbp
    outputs:
        -k1[np]: integer, k-level at each node
        -coeff[np]: interpolation coefficient
    '''
    #surface 
    idxs=zinter>=zcor[:,-1]
    k1[idxs]=nvrt-2
    coeff[idxs]=1.0

    #bottom
    idxs=zinter<zcor[:,0]
    k1[idxs]=kbp[idxs]
    coeff[idxs]=0.0

    for k in np.arange(nvrt-1):
        idxs=(zinter>=zcor[:,k])*(zinter<zcor[:,k+1])
        k1[idxs]=k
        coeff[idxs]=(zinter[idxs]-zcor[idxs,k])/(zcor[idxs,k+1]-zcor[idxs,k])

    if sum(np.isnan(np.r_[k1,coeff])) != 0:
        sys.exit('Check vertical interpolation')
    return np.array(k1).astype('int'), np.array(coeff)


if __name__ == '__main__':

    '''
    Usage: python extract_slab_fcst_netcdf4.py stack
    Input: 
        1. Assign work directory path (for example: fpath='.')

        2. All netcdf files are under {fpath}/outputs/:
           {fpath}/outputs/out2d_*.nc
           {fpath}/outputs/horizontalVelX_*.nc
           {fpath}/outputs/horizontalVelY_*.nc
           {fpath}/outputs/salinity_*.nc
           {fpath}/outputs/temperature_*.nc

        3. output directory:
           outdir = '{fpath}/extract'

        4. vgrid:
           {fpath}/vgrid.in


    Output:
        schout_UV4.5m_{stack}.nc
    '''

    t0=time()

    argparser = argparse.ArgumentParser()
    argparser.add_argument('stack', help='input stack id')
    args = argparser.parse_args()
    sid = args.stack

    #1. work directory
    fpath = "."

    #2. netcdf files
    ds_2d = Dataset(f"{fpath}/outputs/out2d_{sid}.nc")
    ds_u = Dataset(f"{fpath}/outputs/horizontalVelX_{sid}.nc")
    ds_v = Dataset(f"{fpath}/outputs/horizontalVelY_{sid}.nc")
    ds_s = Dataset(f"{fpath}/outputs/salinity_{sid}.nc")
    ds_t = Dataset(f"{fpath}/outputs/temperature_{sid}.nc")
    #units=ds['time'].units
    #base_date=ds['time'].base_date

    #3. output directory
    outdir= './extract'

    #4. get kbp and sigma from vgrid.in
    fid=open(f'{fpath}/vgrid.in','r')
    lines=fid.readlines()
    fid.close()
    ivcor=int(lines[0].strip().split()[0])
    nvrt=int(lines[1].strip().split()[0])
    lines=lines[2:]
    sline=np.array(lines[0].split()).astype('float')
    if sline.min()<0:
        kbp=np.array([int(i.split()[1])-1 for i in lines])
        NP=len(kbp)
        print(f'NP is {NP}')
        sigma=-np.ones([NP,nvrt])
        for i, line in enumerate(lines):
            sigma[i,kbp[i]:]=np.array(line.strip().split()[2:]).astype('float')
    else:
        sline=sline.astype('int')
        kbp=sline-1
        NP=len(sline)
        sigma=np.array([line.split()[1:] for line in lines[1:]]).T.astype('float')
        fpm=sigma<-1
        sigma[fpm]=-1

    #print(np.unique(kbp))

    #get coordinates
    x=ds_2d['SCHISM_hgrid_node_x'][:]
    y=ds_2d['SCHISM_hgrid_node_y'][:]
    depth=ds_2d['depth'][:]

    #get wetdry nodes
    #wd_nodes=ds['wetdry_node'][:,:]
    elev2d=ds_2d['elevation'][:,:]
    #maxelevation
    maxelev=np.max(elev2d,axis=0)
    #get mask
    idry=np.where(maxelev+depth <= 1e-6)
    elev2d[:,idry]=-99999

    #get elements and split quads into tris
    elements=ds_2d['SCHISM_hgrid_face_nodes'][:,:]
    tris = split_quads(elements=elements)  # modified by FY
    NE=len(tris)
    NV=3
    print(f'NE is {NE}')

    #get times
    times = ds_2d['time'][:]
    #print(times)
    ntimes = len(times)
    #ntimes = 2

    #variables for surface salt/temp/u/v
    temp_sur = np.full((ntimes, NP), np.nan)
    salt_sur = np.full((ntimes, NP), np.nan)

    uvel_sur = np.full((ntimes, NP), np.nan)
    vvel_sur = np.full((ntimes, NP), np.nan)

    #variables for u/v at 4.5m below the surface
    uvel_inter = np.full((ntimes, NP), np.nan)
    vvel_inter = np.full((ntimes, NP), np.nan)

    #variables for bottom salt/temp/u/v
    temp_bot = np.full((ntimes, NP), np.nan)
    salt_bot = np.full((ntimes, NP), np.nan)
    uvel_bot = np.full((ntimes, NP), np.nan)
    vvel_bot = np.full((ntimes, NP), np.nan)

    for it in np.arange(ntimes):
        print(it)
        elev=ds_2d['elevation'][it,:]
        
        #surface
        temp_sur[it,:]=ds_t['temperature'][it,:,-1]
        salt_sur[it,:]=ds_s['salinity'][it,:,-1]
        uvel_sur[it,:]=np.squeeze(ds_u['horizontalVelX'][it,:,-1])
        vvel_sur[it,:]=np.squeeze(ds_v['horizontalVelY'][it,:,-1])

        #all levels
        salt_tmp = np.squeeze(ds_s['salinity'][it,:,:])
        temp_tmp = np.squeeze(ds_t['temperature'][it,:,:])

        uvel=np.squeeze(ds_u['horizontalVelX'][it,:,:])
        vvel=np.squeeze(ds_v['horizontalVelY'][it,:,:])

        #compute z#cor
        zcor=(depth[:,None]+elev[:,None])*sigma
    
        level=[-4.5]

        k1=np.full((NP), np.nan)
        coeff=np.full((NP), np.nan)
        zinter=np.ones(NP)*level+elev

        k1, coeff = get_zcor_interp_coefficient(zcor, zinter, kbp)

        #bottom salt/temp/u/v
        temp_bot[it, :]=temp_tmp[np.arange(NP), kbp]
        salt_bot[it, :]=salt_tmp[np.arange(NP), kbp]

        uvel_bot[it, :]=uvel[np.arange(NP), kbp-1]
        vvel_bot[it, :]=vvel[np.arange(NP), kbp-1]

        #tmp=np.array(salt[np.arange(NP),k1]*(1-coeff)+salt[np.arange(NP),k1+1]*coeff)
        #salt_inter[it, :]=np.squeeze(tmp)

        #interpolate at level
        tmp=np.array(uvel[np.arange(NP),k1]*(1-coeff)+uvel[np.arange(NP),k1+1]*coeff)
        uvel_inter[it, :]=np.squeeze(tmp)

        tmp=np.array(vvel[np.arange(NP),k1]*(1-coeff)+vvel[np.arange(NP),k1+1]*coeff)
        vvel_inter[it, :]=np.squeeze(tmp)
        #print(f'It took {time()-t0} to interpolate')

    #Mask dry nodes
    temp_sur[:,idry]=-99999
    salt_sur[:,idry]=-99999
    uvel_sur[:,idry]=-99999
    vvel_sur[:,idry]=-99999
    temp_bot[:,idry]=-99999
    salt_bot[:,idry]=-99999
    uvel_bot[:,idry]=-99999
    vvel_bot[:,idry]=-99999

    #u/v at 4.5m
    uvel_inter[:,idry]=-99999
    vvel_inter[:,idry]=-99999

    #change fill_values
    elev2d[np.where(elev2d>10000)]=-99999
    temp_sur[np.where(temp_sur>10000)]=-99999
    salt_sur[np.where(salt_sur>10000)]=-99999
    uvel_sur[np.where(uvel_sur>10000)]=-99999
    vvel_sur[np.where(vvel_sur>10000)]=-99999

    temp_bot[np.where(temp_bot>10000)]=-99999
    salt_bot[np.where(salt_bot>10000)]=-99999
    uvel_bot[np.where(uvel_bot>10000)]=-99999
    vvel_bot[np.where(vvel_bot>10000)]=-99999

    uvel_inter[np.where(uvel_inter>10000)]=-99999
    vvel_inter[np.where(vvel_inter>10000)]=-99999

#   outdir= '/home1/06923/hyu05/work/oper_3D/run/extract'
    with Dataset(f"{outdir}/schout_UV4.5m_{sid}.nc", "w", format="NETCDF4") as fout:
        #dimensions
        fout.createDimension('time', None)
        fout.createDimension('nSCHISM_hgrid_node', NP)
        fout.createDimension('nSCHISM_hgrid_face', NE)
        fout.createDimension('nMaxSCHISM_hgrid_face_nodes', NV)

        #variables
        fout.createVariable('time', 'f', ('time',))
        fout['time'].long_name="Time"
        #fout['time'].units = units #f'seconds since {date.year}-{date.month}-{date.day} 00:00:00 UTC'
        #fout['time'].base_date=base_date #(date.year, date.month, date.day, 0)
        fout['time'].standard_name="time"
        fout['time'][:] = times

        fout.createVariable('SCHISM_hgrid_node_x', 'f8', ('nSCHISM_hgrid_node',))
        fout['SCHISM_hgrid_node_x'].long_name="node x-coordinate"
        fout['SCHISM_hgrid_node_x'].standard_name="longitude"
        fout['SCHISM_hgrid_node_x'].units="degrees_east"
        fout['SCHISM_hgrid_node_x'].mesh="SCHISM_hgrid"
        fout['SCHISM_hgrid_node_x'][:]=x

        fout.createVariable('SCHISM_hgrid_node_y', 'f8', ('nSCHISM_hgrid_node',))
        fout['SCHISM_hgrid_node_y'].long_name="node y-coordinate"
        fout['SCHISM_hgrid_node_y'].standard_name="latitude"
        fout['SCHISM_hgrid_node_y'].units="degrees_north"
        fout['SCHISM_hgrid_node_y'].mesh="SCHISM_hgrid"
        fout['SCHISM_hgrid_node_y'][:]=y

        fout.createVariable('SCHISM_hgrid_face_nodes', 'i', ('nSCHISM_hgrid_face','nMaxSCHISM_hgrid_face_nodes',))
        fout['SCHISM_hgrid_face_nodes'].long_name="element"
        fout['SCHISM_hgrid_face_nodes'].standard_name="face_node_connectivity"
        fout['SCHISM_hgrid_face_nodes'].start_index=1
        fout['SCHISM_hgrid_face_nodes'].units="nondimensional"
        fout['SCHISM_hgrid_face_nodes'][:]=np.array(tris)

        fout.createVariable('depth', 'f', ('nSCHISM_hgrid_node',))
        fout['depth'].long_name="bathymetry"
        fout['depth'].units="m"
        fout['depth'].mesh="SCHISM_hgrid"
        fout['depth'][:]=depth

        fout.createVariable('elev', 'f8', ('time', 'nSCHISM_hgrid_node',), fill_value=-99999)
        fout['elev'].long_name="water elevation"
        fout['elev'].units="m"
        fout['elev'].mesh="SCHISM_hgrid"
        #fout['elev'].missing_value=np.nan
        fout['elev'][:,:]=elev2d

        fout.createVariable('temp_surface','f8', ('time', 'nSCHISM_hgrid_node',),fill_value=-99999)
        fout['temp_surface'].long_name="sea surface temperature"
        fout['temp_surface'].units="deg C"
        #fout['temp'].missing_value=np.nan
        fout['temp_surface'][:,:]=temp_sur

        fout.createVariable('temp_bottom','f8', ('time', 'nSCHISM_hgrid_node',),fill_value=-99999)
        fout['temp_bottom'].long_name="Bottom temperature"
        fout['temp_bottom'].units="deg C"
        #fout['temp'].missing_value=np.nan
        fout['temp_bottom'][:,:]=temp_bot

        fout.createVariable('salt_surface','f8', ('time', 'nSCHISM_hgrid_node',), fill_value=-99999)
        fout['salt_surface'].long_name="sea surface salinity"
        fout['salt_surface'].units="psu"
        #fout['salt'].missing_value=np.nan
        fout['salt_surface'][:,:]=salt_sur

        fout.createVariable('salt_bottom','f8', ('time', 'nSCHISM_hgrid_node',), fill_value=-99999)
        fout['salt_bottom'].long_name="Bottom salinity"
        fout['salt_bottom'].units="psu"
        #fout['salt'].missing_value=np.nan
        fout['salt_bottom'][:,:]=salt_bot

        fout.createVariable('uvel_surface','f8', ('time', 'nSCHISM_hgrid_node',), fill_value=-99999)
        fout['uvel_surface'].long_name="U-component at the surface"
        fout['uvel_surface'].units="m/s"
        #fout['uvel'].missing_value=np.nan
        fout['uvel_surface'][:,:]=uvel_sur

        fout.createVariable('vvel_surface','f8', ('time', 'nSCHISM_hgrid_node',), fill_value=-99999)
        fout['vvel_surface'].long_name="V-component at the surface"
        fout['vvel_surface'].units="m/s"
        #fout['vvel'].missing_value=np.nan
        fout['vvel_surface'][:,:]=vvel_sur

        fout.createVariable('uvel_bottom','f8', ('time', 'nSCHISM_hgrid_node',), fill_value=-99999)
        fout['uvel_bottom'].long_name="U-component at the bottom"
        fout['uvel_bottom'].units="m/s"
        #fout['uvel'].missing_value=np.nan
        fout['uvel_bottom'][:,:]=uvel_bot

        fout.createVariable('vvel_bottom','f8', ('time', 'nSCHISM_hgrid_node',), fill_value=-99999)
        fout['vvel_bottom'].long_name="V-component at the bottom"
        fout['vvel_bottom'].units="m/s"
        #fout['vvel'].missing_value=np.nan
        fout['vvel_bottom'][:,:]=vvel_bot

        fout.createVariable('uvel4.5','f8', ('time', 'nSCHISM_hgrid_node',), fill_value=-99999)
        fout['uvel4.5'].long_name="U-component at 4.5m below free surface"
        fout['uvel4.5'].units="m/s"
        #fout['uvel'].missing_value=np.nan
        fout['uvel4.5'][:,:]=uvel_inter

        fout.createVariable('vvel4.5','f8', ('time', 'nSCHISM_hgrid_node',), fill_value=-99999)
        fout['vvel4.5'].long_name="V-component at 4.5m below free surface"
        fout['vvel4.5'].units="m/s"
        #fout['vvel'].missing_value=np.nan
        fout['vvel4.5'][:,:]=vvel_inter
        
        fout.title = 'SCHISM Model output'
        fout.source = 'SCHISM model output version v10'
        fout.references = 'http://ccrm.vims.edu/schismweb/'

    print(f'It took {time()-t0} to interpolate')
