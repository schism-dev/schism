#!/usr/bin/env python3
import argparse
from datetime import datetime, timedelta
import pytz

import numpy as np
import pandas as pd
from pandas.tseries.offsets import Second
import netCDF4 as nc
from netCDF4 import Dataset


if __name__ == '__main__':

    #input paramters 
    argparser = argparse.ArgumentParser()
    argparser.add_argument('date',
        type=datetime.fromisoformat,
        help="The date format 'YYYY-MM-DD HH:MM:SS'",
    )

    args=argparser.parse_args()
    startdate=args.date
    print(f'startdate is {startdate}')

    #date = datetime.now() - timedelta(days=1)
    #startDT = datetime(date.year, date.month, date.day, 12)
    enddate = startdate + timedelta(days=1)
    datevectors = pd.date_range(start=startdate.strftime('%Y-%m-%d %H:00:00'), end=enddate.strftime('%Y-%m-%d %H:00:00'), tz='UTC')
    enddate2 = startdate + timedelta(days=5)
    datevectors2 = pd.date_range(start=startdate.strftime('%Y-%m-%d %H:00:00'), end=enddate2.strftime('%Y-%m-%d %H:00:00'), tz='UTC')

    rivers = ['St_lawrence']

    temp = {}

    #st Lawrence river - linear regression with airT (y=0.83x+2.817)
    point = (45.415, -73.623056)
    #ds = Dataset('./sflux/sflux_air_1.0001.nc')
    #ds = Dataset('./sflux/D.nc')
    ds = Dataset('stofs_3d_atl.t12z.gfs.rad.nc')

    lon = ds['lon'][0, :]
    lat = ds['lat'][:, 0]
    idxs = ((lat - point[0]) > 0) & ((lat - point[0]) < 0.2)
    lat_idx = np.where(idxs)[0]
    idxs = ((lon - point[1]) > 0) & ((lon - point[1]) < 0.2)
    lon_idx = np.where(idxs)[0]

    times = ds['time'][:]
    #sflux_startdate = datetime.strptime(ds['time'].units.split('since ')[-1], '%Y-%m-%d %H:%M:%S')
    sflux_startdate = pd.Timestamp(ds['time'].units.split('since ')[-1], tz='UTC')

    #timestamps = [sflux_startdate + timedelta(seconds=round(dt*86400)) for dt in times]
    timestamps = [sflux_startdate + round(dt*86400)*Second() for dt in times]

    #airT = np.squeeze(ds['stmp'][11:24*7:24, lat_idx, lon_idx] - 273.15)
    airT = np.squeeze(ds['stmp'][:, lat_idx, lon_idx] - 273.15)
    waterT = 0.83 * airT + 2.817
    #set waterT below zero to zero
    idxs = waterT < 0
    waterT[idxs] = 0
    df = pd.DataFrame(waterT, index=timestamps)

    #resample to hourly and fill nan if exists
    df2 = df.resample('h').mean()
    df2.fillna(method='bfill', inplace=True)

    #get time indexes for daily at T12
    timestamps = df2.index
    indices = np.where(np.isin(timestamps, datevectors2))[0]
    temp['St_lawrence'] = df2.iloc[np.array(indices)][0].values

    #write file
    data = []
    for i, date in enumerate(datevectors2):
        line = []
        dt = (date - datevectors[0]).total_seconds()
        print(f'time = {dt}')
        line.append(dt)
        for riv in rivers:
            line.append(temp[riv][i])

        data.append(line)

    newset = np.array(data)
    np.savetxt('TEM_1.th', newset, fmt='%.3f')
