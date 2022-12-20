'''
Sample Usage: python generate_station_timeseries.py --date YYYY-MM-DD-HH --input_dir ./outputs/ --output_dir ./
'''

from audioop import mul
from datetime import datetime, timedelta
from time import time
import argparse

import numpy as np
import pandas as pd
from netCDF4 import Dataset, stringtochar

from scipy import interpolate
import json


argparser = argparse.ArgumentParser()
argparser.add_argument('--date', type=datetime.fromisoformat, help='input file date')
argparser.add_argument('--input_dir', type=str, required=True)
argparser.add_argument('--output_dir', type=str, required=True)
args=argparser.parse_args()

# -------------------------- input paramters  ----------------------------------
# Input (1), command line input: date (e.g., 2000-01-01-12, yyyy-mm-dd-HH)
date=args.date

# Input (2), command line input: input file, fullpath to "staout_1"
input_dir=args.input_dir

# Input (3), command line input: output_dir, where *.nc will be saved.
output_dir=args.output_dir
# -------------------------- end input paramters  ----------------------------------

df=pd.read_csv("stations_noaa-coops_164.csv", index_col=[0])
#name=df['Name']
station_info=df['station_info']
lon=df['lon']
lat=df['lat']
nstation=len(station_info)
namelen=50

with open('./schism_staout.json') as d:
    var_dict = json.load(d)

#write to netcdf file
with Dataset(f"{output_dir}/schout_timeseries_at_obs_locations_{date.strftime('%Y%m%d')}.nc", "w", format="NETCDF4") as fout:
    for ivar, var in enumerate(var_dict):
        #read model output
        staout_fname = var_dict[var]['staout_fname']
        data=np.loadtxt(f"{input_dir}/{staout_fname}")
        time=data[:,0]
        nt=len(time)
        #print(nt)
        #print(nstation)
        model=np.ndarray(shape=(nt,nstation), dtype=float)
        model[:,:]=data[:,1:]

        out_dt = 360
        t_interp = np.arange(out_dt,time[-1]+out_dt/2, out_dt)
        f_interp = interpolate.interp1d(time, model, axis=0, fill_value='extrapolate')
        model = f_interp(t_interp)

        startdate=date  # -timedelta(days=1)

        #variables
        if ivar==0:
            #dimensions
            fout.createDimension('station', nstation)
            fout.createDimension('namelen', namelen)
            fout.createDimension('time', None)

            fout.createVariable('time', 'f8', ('time',))
            fout['time'].long_name="Time"
            fout['time'].units = f'seconds since {startdate.year}-{startdate.month}-{startdate.day} {startdate.hour}:00:00 UTC'
            fout['time'].base_date=f'{startdate.year}-{startdate.month}-{startdate.day} {startdate.hour}:00:00 UTC'
            fout['time'].standard_name="time"
            fout['time'][:] = t_interp

            fout.createVariable('station_name', 'c', ('station','namelen',))
            fout['station_name'].long_name="station name"
            names=[]
            names=np.empty((nstation,), 'S'+repr(namelen))
            for i in np.arange(nstation):
                names[i]=str(station_info[i])
            namesc=stringtochar(names)
            fout['station_name'][:]=namesc

            fout.createVariable('x', 'f8', ('station',))
            fout['x'].long_name="longitude"
            fout['x'].standard_name="longitude"
            fout['x'].units="degrees_east"
            fout['x'].positive="east"
            fout['x'][:]=lon

            fout.createVariable('y', 'f8', ('station',))
            fout['y'].long_name="latitude"
            fout['y'].standard_name="latitude"
            fout['y'].units="degrees_north"
            fout['y'].positive="north"
            fout['y'][:]=lat

            fout.title = 'SCHISM Model output'
            fout.source = 'SCHISM model output version v10'
            fout.references = 'http://ccrm.vims.edu/schismweb/'

        out_var = var_dict[var]['name']
        fout.createVariable(out_var, 'f8', ('time', 'station',), fill_value=-99999.)
        fout[out_var].long_name=var_dict[var]['long_name']
        fout[out_var].standard_name=var_dict[var]['stardard_name']
        fout[out_var].units=var_dict[var]['units']
        fout[out_var][:,:]=model
