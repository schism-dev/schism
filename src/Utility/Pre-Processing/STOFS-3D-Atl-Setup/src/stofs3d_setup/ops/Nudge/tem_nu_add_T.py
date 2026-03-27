#! /usr/bin/env python
import xarray as xr

ds = xr.open_dataset('TEM_nu.nc')
data = ds['tracer_concentration']
data += 1.0
ds['tracer_concentration'] = data
ds.to_netcdf('TEM_nu_Tplus1.nc')
