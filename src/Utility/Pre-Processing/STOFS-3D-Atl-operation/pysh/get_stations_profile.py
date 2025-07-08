#!/usr/bin/env python3
from time import time
import os
import sys
from datetime import datetime
import argparse

import numpy as np
from netCDF4 import Dataset, stringtochar

from pylib import loadz, read_schism_hgrid, read_schism_vgrid, save_schism_grid, proj, \
                 zdata, read_schism_bpfile

def read_station_file(station_file_name):
    with open(station_file_name) as f:
        f.readline()
        f.readline()
        station_name=[]
        lon=[]
        lat=[]
        for line in f.read().splitlines():
            if '!' in line:
                station_name.append(line.split('!')[-1])
                lon.append(line.split()[1])
                lat.append(line.split()[2])
    return np.array(lon).astype('float'), np.array(lat).astype('float'), np.array(station_name)

if __name__ == '__main__':
    '''
    Example usage: python get_stations_profile.py --date 2000-01-01-12 --stack_start 1 --stack_end 2 --output_dir ./
    Run this script in the schism run directory, whichi should include:
      1. hgrid.gr3
      2. vgrid.in
      3. station.in
      4. outputs/
    '''
    argparser = argparse.ArgumentParser()
    argparser.add_argument('--date', type=datetime.fromisoformat, help='input file date')
    argparser.add_argument('--stack_start', type=int, required=True)
    argparser.add_argument('--stack_end', type=int, required=True)
    argparser.add_argument('--output_dir', type=str, required=True)
    args = argparser.parse_args()

    # -------------------------- input paramters  ----------------------------------
    # Input (1), command line input: date (e.g., 2000-01-01-12, yyyy-mm-dd-HH)
    startdate=args.date

    # Input (2), stack numbers
    stack_start=args.stack_start
    stack_end=args.stack_end

    # Input (3), the ouput dir of the extracted profile
    output_dir = args.output_dir.strip()
    # ----------------------- end input parameters------------------------------

    print('List of your inputs:\n')
    print(f'startdate: {startdate}\n')
    print(f'stack_start: {stack_start}\n')
    print(f'stack_end: {stack_end}\n')
    print(f'output dir: {output_dir}\n')

    if not os.path.exists(output_dir):
        os.mkdir(output_dir)

    #hgrid.gr3, vgrid.in, outputs/)
    #hgrid = './hgrid.gr3'
    #vgrid = './vgrid.in'
    schism_outputs_dir = './outputs'

    #station files
    stationID = './station.in'

    #Read grid
    grid = './grid.npz'
    if os.path.exists(grid):
        hgrid = loadz(grid).hgrid
        vgrid = loadz(grid).vgrid
    else:
        save_schism_grid()
        hgrid = loadz(grid).hgrid
        vgrid = loadz(grid).vgrid


    #generate stacks list
    stacks = [i for i in np.arange(stack_start, stack_end+1)]
    print(f'Stacks are {stacks}!')

    #Initialize variables
    times = []
    zeta = []
    uwind = []
    vwind = []
    salt = []
    temp = []
    uvel = []
    vvel = []
    zcor = []
    svars2d = {'elevation': zeta, 'windSpeedX': uwind, 'windSpeedY': vwind}
    svars3d = {'salinity': salt, 'temperature': temp, 'horizontalVelX': uvel, 'horizontalVelY': vvel, 'zCoordinates': zcor}
    #out2d: SCHISM_hgrid_node_x, SCHISM_hgrid_node_y, depth, SCHISM_hgrid_face_node

    #get station id
    # lon, lat, station_name = read_station_file(stationID)
    # nstation = len(station_name)

    #Read station bp file
    bp = read_schism_bpfile(stationID)
    nstation = len(bp.station)
    lon = bp.x
    lat = bp.y
    station_name = bp.station

    #Compute area coordinate for stations
    bp.ie, bp.ip, bp.acor = hgrid.compute_acor(np.c_[lon, lat], fmt=1, out=0)
    #check pts inside grid
    pts_outside_grid = np.nonzero(bp.ie == -1)[0]
    if len(pts_outside_grid) != 0:
        for pt in pts_outside_grid:
            print(f'Warning: Station {station_name[pt]} ({bp.x[pt]}, {bp.y[pt]}) is outside model domain, ')
            print('the value of the nearest node will be used.\n')

    bp.depth = hgrid.dp[bp.ip]
    bp.depth0 = (bp.depth*bp.acor).sum(axis=1)

    #get sigma
    if vgrid.ivcor == 1:
        bp.sigma = vgrid.sigma[bp.ip]
        bp.kbp = vgrid.kbp[bp.ip]
        vgrid.sigma = None

    for i, istack in enumerate(stacks):
        #get elevation
        ds2d = Dataset(f'{schism_outputs_dir}/out2d_{istack}.nc')

        #get times
        times2 = ds2d['time'][:]
        ntimes = len(times2)
        [times.append(it) for it in times2]

        eis = []
        k1s = []
        k2s = []
        rats = []
        for it in np.arange(ntimes):
            print(f'Processing time {it} from stack {istack} for 2D variables!')
            for var, value in svars2d.items():
                trii = ds2d.variables[var][it][bp.ip]
                tri = (trii*bp.acor).sum(axis=1)
                value.append(tri)

        ds2d.close()

        #svars3d = {'salinity': salt, 'temperature': temp, 'horizontalVelX': uvel, 'horizontalVelY': vvel, 'zCoordinates': zcor}
        for var, value in svars3d.items():
            ds3d = Dataset(f'{schism_outputs_dir}/{var}_{istack}.nc')
            ndims = ds3d.variables[var].ndim
            dimname = ds3d.variables[var].dimensions

            for it in np.arange(ntimes):

                print(f'Processing time {it} from stack {istack} for 3D variables {var}!')

                if ('nSCHISM_hgrid_node' in dimname):
                    trii = ds3d[var][it][bp.ip]
                elif ('nSCHISM_hgrid_face' in dimname):
                    trii = ds3d[var][it][bp.ie]
                else:
                    sys.exit(f'Unkown variable format: {var}!')

                #horizontal interpolation
                if ('nSCHISM_hgrid_node' in dimname):
                    if ndims == 2:
                        tri = (trii*bp.acor).sum(axis=1)
                    if ndims == 3:
                        tri = (trii*bp.acor[..., None]).sum(axis=1)
                    if ndims == 4:
                        tri = (trii*bp.acor[..., None, None]).sum(axis=1)
                        rat = rat[:, None]
                else:
                    tri = trii

                value.append(tri)
            ds3d.close()


    namelen = 100
    with Dataset(f'{output_dir}/stofs_stations_profile_{stack_start}_{stack_end}.nc', "w", format="NETCDF3_CLASSIC") as fout:
        fout.createDimension('station', nstation)
        fout.createDimension('namelen', namelen)
        fout.createDimension('siglay', vgrid.nvrt)
        fout.createDimension('time', None)

        #variables
        fout.createVariable('time', 'f4', ('time',))
        fout['time'].long_name="Time"
        fout['time'].units = f'seconds since {startdate.year}-{startdate.month:02d}-{startdate.day:02d} {startdate.hour:02d}:00:00 UTC'
        fout['time'].base_date=f'{startdate.year}-{startdate.month:02d}-{startdate.day:02d} {startdate.hour:02d}:00:00 UTC'
        fout['time'].standard_name="time"
        fout['time'][:] = times

        fout.createVariable('station_name', 'c', ('station','namelen',))
        fout['station_name'].long_name="station name"
        names=np.empty((nstation,), 'S'+repr(namelen))
        for i in np.arange(nstation):
            str_in=station_name[i]
            #strlen=len(str_in)
            #str_out=list(str_in)
            #tmp="".join(str_out[j] for j in range(strlen))
            #names[i]=str(station_id[i])+" "+tmp
            names[i]=str_in[:-1]
        namesc=stringtochar(names)
        fout['station_name'][:]=namesc

        fout.createVariable('lon', 'f4', ('station',))
        fout['lon'].long_name="longitude"
        fout['lon'].standard_name="longitude"
        fout['lon'].units="degrees_east"
        fout['lon'].positive="east"
        fout['lon'][:]=lon

        fout.createVariable('lat', 'f4', ('station',))
        fout['lat'].long_name="latitude"
        fout['lat'].standard_name="latitude"
        fout['lat'].units="degrees_north"
        fout['lat'].positive="north"
        fout['lat'][:]=lat

        fout.createVariable('depth', 'f4', ('station',))
        fout['depth'].long_name = "Bathymetry"
        fout['depth'].standard_name = "depth"
        fout['depth'].units = "meters"
        fout['depth'][:] = bp.depth0

        fout.createVariable('zeta', 'f4', ('time', 'station',), fill_value=-99999.)
        fout['zeta'].long_name="water surface elevation above navd88"
        fout['zeta'].standard_name="sea_surface_height_above_navd88"
        fout['zeta'].units="m"
        fout['zeta'][:,:]=np.array(zeta)

        fout.createVariable('zCoordinates', 'f4', ('time', 'station', 'siglay',), fill_value=-99999.)
        fout['zCoordinates'].long_name = "vertical coordinate, positive upward"
        fout['zCoordinates'].standard_name = "vertical coordinate"
        fout['zCoordinates'].units = "m"
        fout['zCoordinates'][:,:,:]=np.array(zcor)

        fout.createVariable('salinity', 'f4', ('time', 'station', 'siglay',), fill_value=-99999.)
        fout['salinity'].long_name = "salinity"
        fout['salinity'].standard_name = "sea_water_salinity"
        fout['salinity'].units = "psu"
        fout['salinity'][:,:,:]=np.array(salt)

        fout.createVariable('temperature', 'f4', ('time', 'station', 'siglay',), fill_value=-99999.)
        fout['temperature'].long_name = "temperature"
        fout['temperature'].standard_name = "sea_water_temperature"
        fout['temperature'].units = "degree_C"
        fout['temperature'][:,:,:]=np.array(temp)

        fout.createVariable('u', 'f4', ('time', 'station', 'siglay',), fill_value=-99999.)
        fout['u'].long_name = "Eastward Water Velocity"
        fout['u'].standard_name = "eastward_sea_water_velocity"
        fout['u'].units = "meters s-1"
        fout['u'][:,:,:]=np.array(uvel)

        fout.createVariable('v', 'f4', ('time', 'station', 'siglay',), fill_value=-99999.)
        fout['v'].long_name = "Northward Water Velocity"
        fout['v'].standard_name = "northward_sea_water_velocity"
        fout['v'].units = "meters s-1"
        fout['v'][:,:,:]=np.array(vvel)

        fout.createVariable('uwind_speed', 'f4', ('time', 'station',), fill_value=-99999.)
        fout['uwind_speed'].long_name = "Eastward Wind Velocity"
        fout['uwind_speed'].standard_name = "eastward_wind"
        fout['uwind_speed'].units = "meters s-1"
        fout['uwind_speed'][:,:]=np.array(uwind)

        fout.createVariable('vwind_speed', 'f4', ('time', 'station',), fill_value=-99999.)
        fout['vwind_speed'].long_name = "Northward Wind Velocity"
        fout['vwind_speed'].standard_name = "northward_wind"
        fout['vwind_speed'].units = "meters s-1"
        fout['vwind_speed'][:,:]=np.array(vwind)

        fout.title = 'SCHISM Model output'
        fout.references = 'http://ccrm.vims.edu/schismweb/'

