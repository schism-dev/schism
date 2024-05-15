import os
import sys
import gc
from datetime import datetime, timedelta
from time import time
from dateutil import parser
import argparse

import numpy as np
from netCDF4 import Dataset

from generate_adcirc import split_quads  # modified by FY

import psutil
process = psutil.Process(os.getpid())

def read_vgrid(vgrid_file):
    time_start = time()
    fid=open(vgrid_file,'r')
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
    print(f"read vgrid took {time()-time_start}")

    return nvrt, sigma

from dataclasses import dataclass, field
from typing import Any
@dataclass
class VerticalInterpInfo:
    it: int
    zinter: float
    level: Any = field(default_factory=np.array)
    weight: Any = field(default_factory=np.array)

def vertical_interp(var, zcor, zinter, bottom_index, vertical_interp_info=None):
    '''
    inputs:
        -var: var[np,nvrt] for each time record.
        -zcor: zcor[np,nvrt] for each time record.
        -zinter[np, ]: zinter depth where currents will be interpolated at
        -bottom_index[np, ]: bottom level, 0-based index
    '''

    if vertical_interp_info is None:
        print('calculating vertical interpolation weights')
        t1 = time()

        NP, nvrt = zcor.shape
        row_indices = np.arange(NP)
        target_level = np.zeros(NP, dtype=int)  # zinter should be between zcor[:, target_level] and zcor[:, target_level+1]
                                                # , i.e., target_level is the lower level
        target_level_weight = np.zeros(NP, dtype=float)  # linear interpolation, i.e., lower_level_weight + upper_level_weight = 1.0

        # get all points where zinter is above the surface
        idx = zinter >= zcor[:,-1]
        target_level[idx] = nvrt-2  # i.e., surface level (nvrt-1) is the upper level
        target_level_weight[idx] = 0.0  # weight of the lower level is 0.0, i.e., all from the upper level

        # get all points where zinter is below the bottom
        idx = zinter < zcor[row_indices,bottom_index]
        target_level[idx] = bottom_index[idx]
        target_level_weight[idx] = 1.0  # weight of the lower level is 1.0, i.e., all from the lower level

        #intermediate
        for k in np.arange(nvrt-1):
            idx=(zinter>=zcor[:,k])*(zinter<zcor[:,k+1])
            target_level[idx] = k
            target_level_weight[idx] = (zcor[idx,k+1] - zinter[idx]) /(zcor[idx,k+1] - zcor[idx,k])

        print(f"calculating vertical interpolation weights took {time()-t1}")

    else:
        target_level = vertical_interp_info.level
        target_level_weight = vertical_interp_info.weight

    if any(np.isnan(target_level_weight)) or any(np.isnan(target_level)):
        raise ValueError('NaN values in target_level or target_level_weight')

    row_indices = np.arange(len(target_level))
    var_interp = var[row_indices, target_level] * target_level_weight + var[row_indices, target_level + 1] * (1.0 - target_level_weight)

    return [var_interp, target_level, target_level_weight]

def mask_var(data, idry, threshold=10000, mask_value=-99999):
    '''
    Mask dry nodes and large values.
    data is a 2D array of shape (ntimes, NP).
    idry is a boolean array of shape (ntimes, NP) that indicates dry nodes.
    '''
    data[idry] = mask_value  # Mask dry nodes
    data[data > threshold] = mask_value  # mask large values, to be consistent with older version of the script

    return data

# a dictionary of the data to be extracted
native_var_dict = {
    'temp_surface': {
        'file_prefix': 'temperature',
        'var_name': 'temperature',
        'long name': 'sea surface temperature',
        'units': 'deg C',
        'level': 'surface'
    },
    'temp_bottom': {
        'file_prefix': 'temperature',
        'var_name': 'temperature',
        'long name': 'Bottom temperature',
        'units': 'deg C',
        'level': 'bottom'
    },
    'salt_surface': {
        'file_prefix': 'salinity',
        'var_name': 'salinity',
        'long name': 'sea surface salinity',
        'units': 'psu',
        'level': 'surface'
    },
    'salt_bottom': {
        'file_prefix': 'salinity',
        'var_name': 'salinity',
        'long name': 'Bottom salinity',
        'units': 'psu',
        'level': 'bottom'
    },
    'uvel_surface': {
        'file_prefix': 'horizontalVelX',
        'var_name': 'horizontalVelX',
        'long name': 'U-component at the surface',
        'units': 'm/s',
        'level': 'surface'
    },
    'uvel_bottom': {
        'file_prefix': 'horizontalVelX',
        'var_name': 'horizontalVelX',
        'long name': 'U-component at the bottom',
        'units': 'm/s',
        'level': 'near bottom'
    },
    'vvel_surface': {
        'file_prefix': 'horizontalVelY',
        'var_name': 'horizontalVelY',
        'long name': 'V-component at the surface',
        'units': 'm/s',
        'level': 'surface'
    },
    'vvel_bottom': {
        'file_prefix': 'horizontalVelY',
        'var_name': 'horizontalVelY',
        'long name': 'V-component at the bottom',
        'units': 'm/s',
        'level': 'near bottom'
    }
}

interpolated_var_dict = {
    'uvel4.5': {
        'file_prefix': 'horizontalVelX',
        'var_name': 'horizontalVelX',
        'long name': 'U-component at 4.5m below free surface',
        'units': 'm/s',
        'level': -4.5  # 4.5 m below the surface
    },
    'vvel4.5': {
        'file_prefix': 'horizontalVelY',
        'var_name': 'horizontalVelY',
        'long name': 'V-component at 4.5m below free surface',
        'units': 'm/s',
        'level': -4.5  # 4.5 m below the surface
    }
}

if __name__ == '__main__':

    '''
    Usage:
        python extract_slab_fcst_netcdf4.py --stack N
    , where N is the stack id

    The outputs of this script are saved under './extract/', or
    you can set it with an optional argument:
        python extract_slab_fcst_netcdf4.py --stack N  --output_dir some_dir/

    Run this script under a schism run directory that has an outputs/ folder.
    All SCHISM netcdf files should be under {current_dir}/outputs/, including:
     {fpath}/outputs/out2d_*.nc
     {fpath}/outputs/horizontalVelX_*.nc
     {fpath}/outputs/horizontalVelY_*.nc
     {fpath}/outputs/salinity_*.nc
     {fpath}/outputs/temperature_*.nc

    Output:
        {output_dir}/schout_UV4.5m_{stack}.nc
    '''

    t0=time()

    # --------------------- parse input arguments ---------------------
    # parse input arguments
    argparser = argparse.ArgumentParser()
    argparser.add_argument('--stack', required=True, help='input stack id')
    argparser.add_argument('--output_dir', default='./extract/', help='A SCHISM run folder that has outputs/')
    argparser.add_argument('--mem_save_mode', default=False, help='Setting the memory saving mode: true for less memory consumption but slightly slower')
    args = argparser.parse_args()

    sid = args.stack
    outdir= args.output_dir
    mem_save_mode = args.mem_save_mode
    os.makedirs(outdir, exist_ok=True)

    # --------------------- specify schism outputs ---------------------
    fpath = './'  # should be a SCHISM run directory that contains outputs/

    # open 3D schism outputs
    schism_output = {
        'temperature': Dataset(f"{fpath}/outputs/temperature_{sid}.nc")['temperature'],
        'salinity': Dataset(f"{fpath}/outputs/salinity_{sid}.nc")['salinity'],
        'horizontalVelX': Dataset(f"{fpath}/outputs/horizontalVelX_{sid}.nc")['horizontalVelX'],
        'horizontalVelY': Dataset(f"{fpath}/outputs/horizontalVelY_{sid}.nc")['horizontalVelY'],
        'zCoordinates': Dataset(f'{fpath}/outputs/zCoordinates_{sid}.nc')['zCoordinates']
    }
    print(f"Memory usage before reading files: {process.memory_info().rss / 1024 ** 2:.2f} MB")

    # --------------------- basic info, should be same for all input files ---------------------
    # process time information from out2d_*.nc; the time info is the same for all files
    ds = Dataset(f"{fpath}/outputs/out2d_{sid}.nc")
    base_date_str = ds['time'].base_date.split()
    base_datetime = datetime(int(base_date_str[0]), int(base_date_str[1]), int(base_date_str[2]), 0, 0, 0)
    base_date_str = base_datetime.strftime('%Y-%m-%d %H:%M:%S UTC')

    time_units_str = ds['time'].units.split("since")[1]
    time_units_datetime = parser.parse(time_units_str)

    # check time zone is UTC
    if time_units_datetime.tzinfo is None or time_units_datetime.tzinfo.utcoffset(time_units_datetime) != timedelta(0):
        raise ValueError("Time zone is not UTC")

    time_units_str = f"seconds since {time_units_datetime.strftime('%Y-%m-%d %H:%M:%S UTC')}"

    #get coordinates
    depth = np.array(ds['depth'])
    elev2d = np.array(ds['elevation'])

    # nvrt, sigma = read_vgrid(f'{fpath}/vgrid.in')

    x=ds['SCHISM_hgrid_node_x'][:]
    y=ds['SCHISM_hgrid_node_y'][:]
    NP = len(x)
    bottom_index_node = np.array(ds['bottom_index_node']) - 1  # minus 1 to convert to 0-based index

    # set dry mask
    idry = elev2d + depth.reshape(1, -1) <= 1e-6  # inundation <= 1e-6 m is considered dry
    elev2d[idry]=-99999  # mask dry nodes with -99999

    #get elements and split quads into tris
    elements=ds['SCHISM_hgrid_face_nodes'][:,:]
    tris = split_quads(elements=elements)  # modified by FY
    NE=len(tris)
    max_ele_nodes = 3  # max number of nodes per element
    print(f'number of elements increased from {len(elements)} to {NE} after splitting quads to triangles')

    # extract time
    times = ds['time'][:]
    ntimes = len(times)
    time_indices, node_indices = np.ogrid[0:ntimes, 0:NP]
    # Done processing basic information
    ds.close()

    print(f"Memory usage after extracting basic info: {process.memory_info().rss / 1024 ** 2:.2f} MB")

    # write the output file
    with Dataset(f"{outdir}/schout_UV4.5m_{sid}.nc", "w", format="NETCDF4") as fout:
        # basic information
        #dimensions
        fout.createDimension('time', None)
        fout.createDimension('nSCHISM_hgrid_node', NP)
        fout.createDimension('nSCHISM_hgrid_face', NE)
        fout.createDimension('nMaxSCHISM_hgrid_face_nodes', max_ele_nodes)

        #variables
        fout.createVariable('time', 'f', ('time',))
        fout['time'].long_name="Time"
        fout['time'].units = time_units_str #f'seconds since {date.year}-{date.month}-{date.day} 00:00:00 UTC'
        fout['time'].base_date=base_date_str #(date.year, date.month, date.day, 0)
        fout['time'].standard_name="time"
        fout['time'][:] = times
        del times; gc.collect()

        fout.createVariable('SCHISM_hgrid_node_x', 'f8', ('nSCHISM_hgrid_node',))
        fout['SCHISM_hgrid_node_x'].long_name="node x-coordinate"
        fout['SCHISM_hgrid_node_x'].standard_name="longitude"
        fout['SCHISM_hgrid_node_x'].units="degrees_east"
        fout['SCHISM_hgrid_node_x'].mesh="SCHISM_hgrid"
        fout['SCHISM_hgrid_node_x'][:]=x
        del x; gc.collect()

        fout.createVariable('SCHISM_hgrid_node_y', 'f8', ('nSCHISM_hgrid_node',))
        fout['SCHISM_hgrid_node_y'].long_name="node y-coordinate"
        fout['SCHISM_hgrid_node_y'].standard_name="latitude"
        fout['SCHISM_hgrid_node_y'].units="degrees_north"
        fout['SCHISM_hgrid_node_y'].mesh="SCHISM_hgrid"
        fout['SCHISM_hgrid_node_y'][:]=y
        del y; gc.collect()

        fout.createVariable('SCHISM_hgrid_face_nodes', 'i', ('nSCHISM_hgrid_face','nMaxSCHISM_hgrid_face_nodes',))
        fout['SCHISM_hgrid_face_nodes'].long_name="element"
        fout['SCHISM_hgrid_face_nodes'].standard_name="face_node_connectivity"
        fout['SCHISM_hgrid_face_nodes'].start_index=1
        fout['SCHISM_hgrid_face_nodes'].units="nondimensional"
        fout['SCHISM_hgrid_face_nodes'][:]=np.array(tris)
        del tris; gc.collect()

        fout.createVariable('depth', 'f', ('nSCHISM_hgrid_node',))
        fout['depth'].long_name="bathymetry"
        fout['depth'].units="m"
        fout['depth'].mesh="SCHISM_hgrid"
        fout['depth'][:]=depth
        # del depth; gc.collect()

        fout.createVariable('elev', 'f8', ('time', 'nSCHISM_hgrid_node',), fill_value=-99999)
        fout['elev'].long_name="water elevation"
        fout['elev'].units="m"
        fout['elev'].mesh="SCHISM_hgrid"
        #fout['elev'].missing_value=np.nan
        fout['elev'][:,:]=elev2d
        # elev2d is still needed

        # create 2D variables, both native and interpolated
        for var_name, var_info in {**native_var_dict, **interpolated_var_dict}.items():
            fout.createVariable(var_name, 'f8', ('time', 'nSCHISM_hgrid_node',), fill_value=-99999)
            fout[var_name].long_name=var_info['long name']
            fout[var_name].units=var_info['units']

        if not mem_save_mode:
            # native variables that don't need to be interpolated are read in at all time steps
            last_var = ''
            for var_name, var_info in native_var_dict.items():
                t1 = time()
                print(f"\nProcessing {var_name}")

                # some 3D data can be reused, only read when needed
                if last_var != var_info['var_name']:
                    data_3d = np.array(schism_output[var_info['var_name']])
                    last_var = var_info['var_name']

                if var_info['level'] == 'surface':
                    data = data_3d[:, :, -1]
                elif var_info['level'] == 'bottom':
                    data = data_3d[:, np.arange(NP), bottom_index_node]
                elif var_info['level'] == 'near bottom':
                    data = data_3d[:, np.arange(NP), bottom_index_node+1]

                data = mask_var(data, idry)
                print(f"Memory usage after extracting {var_name}: {process.memory_info().rss / 1024 ** 2:.2f} MB")
                fout[var_name][:] = data
                print(f"Processing {var_name} took {time()-t1} seconds\n")

        # do interpolation time-step by time-step to save memory
        if mem_save_mode:
            var_dict = {**native_var_dict, **interpolated_var_dict}
        else:
            var_dict = interpolated_var_dict

        for it in range(ntimes):
            print(f"\nProcessing time {it+1}/{ntimes}")
            t_it = time()
            vertical_interp_info = VerticalInterpInfo(-1, 0.0, None, None)  # force reset for a new time step
            for var_name, var_info in var_dict.items():
                print(f"Processing {var_name}")

                if var_info['level'] == 'surface':
                    data = schism_output[var_info['var_name']][it, :, -1]
                elif var_info['level'] == 'bottom':
                    data_3d = np.array(schism_output[var_info['var_name']][it, :, :])
                    data = data_3d[np.arange(NP), bottom_index_node]
                elif var_info['level'] == 'near bottom':
                    data_3d = np.array(schism_output[var_info['var_name']][it, :, :])
                    data = data_3d[np.arange(NP), bottom_index_node + 1]
                else:  # specified zinter, to be interpolated
                    if type(var_info['level']) not in [int, float]:
                        raise ValueError('Invalid level value')
                    zinter = var_info['level'] + elev2d[it, :]  # zinter is relative to the free surface
                    data_3d = np.array(schism_output[var_info['var_name']][it, :, :])
                    # inun_depth = elev2d[it, :] + depth.reshape(1,-1)
                    # zcor = sigma[np.newaxis,:] * inun_depth[:, np.newaxis] + elev2d[it, :, np.newaxis]
                    zcor = np.array(schism_output['zCoordinates'][it, :])

                    # the vertical interpolation is expensive, so we need to check if we can reuse the weights
                    # if the time index or the target z has changed, we need to recalculate the weights
                    if vertical_interp_info.it != it or vertical_interp_info.zinter != var_info['level']:
                        data, level, weight = vertical_interp(data_3d[:, :], zcor[:, :], zinter[:], bottom_index_node[:])
                        vertical_interp_info = VerticalInterpInfo(it, float(var_info['level']), level, weight)
                    else:
                        data, _, _ = vertical_interp(data_3d[:, :], zcor[:, :], zinter[:], bottom_index_node[:], vertical_interp_info)

                # mask dry nodes and large values
                data = mask_var(data, idry[it, :])
                fout[var_name][it, :] = data

                print(f"Memory usage after extracting {var_name}: {process.memory_info().rss / 1024 ** 2:.2f} MB")
            print(f"Processing time step {it+1} took {time()-t_it} seconds\n")

        fout.title = 'SCHISM Model output'
        fout.source = 'SCHISM model output version v10'
        fout.references = 'http://ccrm.vims.edu/schismweb/'

    print(f'Extraction took {time()-t0} seconds')
