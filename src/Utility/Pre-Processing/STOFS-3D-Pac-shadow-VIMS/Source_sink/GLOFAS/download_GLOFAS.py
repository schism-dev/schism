import cdsapi
import datetime
'''cdsapi from https://cds.climate.copernicus.eu/how-to-api
requires CDS API personal token in
$HOME/.cdsapirc (or modify to use personal token as an argument to cdsapi.Client)
'''

dataset = "cems-glofas-forecast"
lon= [-180, 180]
lat= [-60, 90]
bnds = [lat[1], lon[0],lat[0],lon[1] ]
rnday = 7
fname = 'GLOFAS-global-today.nc'
client = cdsapi.Client(url='https://ewds.climate.copernicus.eu/api')
res = False
nday = 1
date = datetime.datetime.now() - datetime.timedelta(days=nday)
request = {
    "system_version": ["operational"],
    "hydrological_model": ["lisflood"],
    "product_type": ["control_forecast"],
    "variable": "river_discharge_in_the_last_24_hours",
    "year": ['%d' % date.year],
    "month": ['%02d' % date.month],
    "day": ['%02d' % date.day],
    "leadtime_hour": [ '%s' % (24*d) for d in range(1,rnday+1)
    ],
    "data_format": "netcdf",
    "download_format": "unarchived",
    "area": bnds
}
print(request)
req = client.retrieve(dataset, request, fname)

