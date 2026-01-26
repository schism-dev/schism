import os
import sys
from datetime import datetime,timedelta
import logging
import pathlib
import tempfile
from typing import Union
from time import time

import numpy as np
import scipy as sp
from numba import jit, prange
import netCDF4 as nc
from netCDF4 import Dataset
from matplotlib.transforms import Bbox
import xarray as xr

logger = logging.getLogger(__name__)

def convert_longitude(ds):
    lon_name = 'lon'
    ds['_lon_adjusted'] = xr.where(
        ds[lon_name] > 180,
        ds[lon_name] - 360,
        ds[lon_name])
    ds = (
        ds.swap_dims({lon_name: '_lon_adjusted'})
        .sel(**{'_lon_adjusted': sorted(ds._lon_adjusted)})
        .drop(lon_name)
    )
    ds = ds.rename({'_lon_adjusted': lon_name})
    return ds

def get_database(date, Bbox=None):
    if date >= datetime(2018, 12, 4):
        database = f'GLBy0.08/expt_93.0'
    elif date >= datetime(2018, 1, 1) and date < datetime(2018, 12, 4):
        database = f'GLBv0.08/expt_93.0'
    elif date >= datetime(2017, 10, 1) and date < datetime(2018, 1, 1):
        database = f'GLBv0.08/expt_92.9'
    elif date >= datetime(2017, 6, 1) and date < datetime(2017, 10, 1):
        database = f'GLBv0.08/expt_57.7'
    elif date >= datetime(2017, 2, 1) and date < datetime(2017, 6, 1):
        database = f'GLBv0.08/expt_92.8'
    elif date >= datetime(2016, 5, 1) and date < datetime(2017, 2, 1):
        database = f'GLBv0.08/expt_57.2'
    elif date >= datetime(2016, 1, 1) and date < datetime(2016, 5, 1):
        database = f'GLBv0.08/expt_56.3'
    elif date >= datetime(1994, 1, 1) and date < datetime(2016, 1, 1):
        database = f'GLBv0.08/expt_53.X/data/{date.year}'
    else:
        logger.info('No data for {date}')
    return database

def get_idxs(date, database, bbox):

    if date.strftime("%Y-%m-%d") >= datetime.now().strftime("%Y-%m-%d"):
        date2 = datetime.now() - timedelta(days=1)
        baseurl = f'https://tds.hycom.org/thredds/dodsC/{database}/FMRC/runs/GLBy0.08_930_FMRC_RUN_{date2.strftime("%Y-%m-%dT12:00:00Z")}?depth[0:1:-1],lat[0:1:-1],lon[0:1:-1],time[0:1:-1]'
    else:
        baseurl=f'https://tds.hycom.org/thredds/dodsC/{database}?lat[0:1:-1],lon[0:1:-1],time[0:1:-1],depth[0:1:-1]'

    ds=Dataset(baseurl)
    time1=ds['time']
    times=nc.num2date(time1,units=time1.units,only_use_cftime_datetimes=False)
    
    lon=ds['lon'][:]
    lat=ds['lat'][:]
    dep=ds['depth'][:]
    lat_idxs=np.where((lat>=bbox.ymin-2.0)&(lat<=bbox.ymax+2.0))[0]
    lon_idxs=np.where((lon>=bbox.xmin-2.0) & (lon<=bbox.xmax+2.0))[0]
    lon=lon[lon_idxs]
    lat=lat[lat_idxs]
    #logger.info(lon_idxs)
    #logger.info(lat_idxs)
    lon_idx1=lon_idxs[0].item()
    lon_idx2=lon_idxs[-1].item()
    #logger.info(f'lon_idx1 is {lon_idx1}, lon_idx2 is {lon_idx2}')
    lat_idx1=lat_idxs[0].item()
    lat_idx2=lat_idxs[-1].item()
    #logger.info(f'lat_idx1 is {lat_idx1}, lat_idx2 is {lat_idx2}')
    
    for ilon in np.arange(len(lon)):
        if lon[ilon] > 180:
            lon[ilon] = lon[ilon]-360.
    #lonc=(np.max(lon)+np.min(lon))/2.0
    #logger.info(f'lonc is {lonc}')
    #latc=(np.max(lat)+np.min(lat))/2.0
    #logger.info(f'latc is {latc}')
    x2, y2=transform_ll_to_cpp(lon, lat)

    idxs=np.where( date == times)[0]
    #check if time_idx is empty
    if len(idxs) == 0:
        #If there is missing data, use the data from the next days, the maximum searching days is 3. Otherwise, stop.
        for i in np.arange(0,3):
            date_before=(date + timedelta(days=int(i)+1)) #.astype(datetime)
            logger.info(f'Try replacing the missing data from {date_before}')
            idxs=np.where(date_before == times)[0]
            if len(idxs) == 0:
                continue
            else:
                break
    if len(idxs) ==0:
        logger.info(f'No date for date {date}')
        sys.exit()
    time_idx=idxs.item()  

    ds.close()

    return time_idx, lon_idx1, lon_idx2, lat_idx1, lat_idx2, x2, y2

def transform_ll_to_cpp(lon, lat, lonc=-77.07, latc=24.0):
    #lonc=(np.max(lon)+np.min(lon))/2.0
    #logger.info(f'lonc is {lonc}')
    #latc=(np.max(lat)+np.min(lat))/2.0
    #logger.info(f'latc is {latc}')
    longitude=lon/180*np.pi
    latitude=lat/180*np.pi
    radius=6378206.4
    loncc=lonc/180*np.pi
    latcc=latc/180*np.pi
    lon_new=[radius*(longitude[i]-loncc)*np.cos(latcc) for i in np.arange(len(longitude))]
    lat_new=[radius*latitude[i] for i in np.arange(len(latitude))]

    return np.array(lon_new), np.array(lat_new)

def interp_to_points_3d(dep, y2, x2, bxyz, val):
    idxs = np.where(abs(val) > 10000)
    val[idxs] = float('nan')

    val_fd = sp.interpolate.RegularGridInterpolator((dep,y2,x2),np.squeeze(val),'linear', bounds_error=False, fill_value = float('nan'))
    val_int = val_fd(bxyz)
    idxs = np.isnan(val_int)
    if np.sum(idxs) != 0:
        val_int[idxs] = sp.interpolate.griddata(bxyz[~idxs,:], val_int[~idxs], bxyz[idxs,:],'nearest')
    idxs = np.isnan(val_int)
    if np.sum(idxs) != 0:
        logger.info(f'There is still missing value for {val}')
        sys.exit()
    return val_int

def interp_to_points_2d(y2, x2, bxy, val):
    idxs = np.where(abs(val) > 10000)
    val[idxs] = float('nan')

    val_fd = sp.interpolate.RegularGridInterpolator((y2,x2),np.squeeze(val),'linear', bounds_error=False, fill_value = float('nan'))
    val_int = val_fd(bxy)
    idxs = np.isnan(val_int)
    if np.sum(idxs) != 0:
        val_int[idxs] = sp.interpolate.griddata(bxy[~idxs,:], val_int[~idxs], bxy[idxs,:],'nearest')
    idxs = np.isnan(val_int)
    if np.sum(idxs) != 0:
        logger.info(f'There is still missing value for {val}')
        sys.exit()
    return val_int

class OpenBoundaryInventory:

    def __init__(self, hgrid):
        self.hgrid = hgrid

    def fetch_data(self, outdir: Union[str, os.PathLike], start_date, rnday, ocean_bnd_ids = [0], elev2D=True): 
        outdir = pathlib.Path(outdir)

        self.start_date = start_date
        self.rnday=rnday
        self.timevector=np.arange(
            self.start_date,
            self.start_date + timedelta(days=self.rnday+1),
            timedelta(days=1)).astype(datetime)

        #Get open boundary 
        gdf=self.hgrid.boundaries.open.copy()
        opbd=[]
        for ibnd in ocean_bnd_ids:
            opbd.extend(list(gdf.iloc[ibnd].indexes))
        blon = self.hgrid.coords[opbd,0]
        blat = self.hgrid.coords[opbd,1]
        #logger.info(f'blon min {np.min(blon)}, max {np.max(blon)}')
        NOP = len(blon)

        #create netcdf
        ntimes=self.rnday+1
        nComp1=1
        nComp2=2
        one=1
        #ndt=np.zeros([ntimes])

        if elev2D:
            #timeseries_el=np.zeros([ntimes,NOP,nComp1])
            #create netcdf 
            dst_elev = Dataset(outdir / 'elev2D.th.nc', 'w', format='NETCDF4')
            #dimensions
            dst_elev.createDimension('nOpenBndNodes', NOP)
            dst_elev.createDimension('one', one)
            dst_elev.createDimension('time', None)
            dst_elev.createDimension('nLevels', one)
            dst_elev.createDimension('nComponents', nComp1)

            #variables
            dst_elev.createVariable('time_step', 'f', ('one',))
            dst_elev['time_step'][:] = 86400

            dst_elev.createVariable('time', 'f', ('time',))
            #dst_elev['time'][:] = ndt

            dst_elev.createVariable('time_series', 'f', ('time', 'nOpenBndNodes', 'nLevels', 'nComponents'))
            #dst_elev['time_series'][:,:,:,:] = timeseries_el

        fname = 'aviso.nc' 
        ds=Dataset(fname)
        lon2 = ds['longitude'][:]
        lat2 = ds['latitude'][:]
        x2, y2 = transform_ll_to_cpp(lon2, lat2)

        time1 = ds['time']
        times=nc.num2date(time1,units=time1.units,only_use_cftime_datetimes=False)

        t0=time()
        for it, date in enumerate(self.timevector):

            print(f'Processing time {date}')
            #get time index
            time_idx = np.where( date == times)[0].item()
            print(f'time_idx is {time_idx}')

            #loop over each open boundary
            ind1 = 0
            ind2 = 0
            for ibnd in ocean_bnd_ids:

                #opbd = list(boundary.indexes)
                opbd = list(gdf.iloc[ibnd].indexes)
                ind1 = ind2
                ind2 = ind1 + len(opbd)
                #logger.info(f'ind1 = {ind1}, ind2 = {ind2}')
                blon = self.hgrid.coords[opbd,0]
                blat = self.hgrid.coords[opbd,1]
                xi,yi = transform_ll_to_cpp(blon, blat)
                bxy = np.c_[yi, xi]
          

                print('****Interpolation starts****')

                #ndt[it]=it*24*3600.

                if elev2D:
                    #ssh
                    ssh=np.squeeze(ds['adt'][time_idx, :,:])

                    dst_elev['time'][it] = it*24*3600.
                    try:
                        ssh_int = interp_to_points_2d(y2, x2, bxy, ssh)
                        dst_elev['time_series'][it,ind1:ind2,0,0] = ssh_int 
                    except:
                        print(f'Warning: no valid data for boundary {ibnd}, use previous day data!')
                        dst_elev['time_series'][it,ind1:ind2,0,0] = dst_elev['time_series'][it-1,ind1:ind2,0,0]

        print(f'Writing *th.nc takes {time()-t0} seconds')
