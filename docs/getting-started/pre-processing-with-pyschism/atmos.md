#**Atmospheric forcing**
The sflux/ directory is required if nws=2 in *param.nml*. 

There are four types of files needed:

sflux_input.txt (<font color="red">required</font>): a namelist file

sflux_air_1.[XXXX].nc (<font color="red">required</font>): NetCDF files that have time (in days), wind speed at 10m above MSL (u, v), air temperatue and spedific humidity at 2m above MSL, sea level pressure;

sflux_prc_1.[XXXX].nc (<font color="blue">needed</font> if isconsv=1): NetCDF files that have time and precipitation rate;

sflux_rad_1.[XXXX].nc (<font color="blue">needed</font> if ihconsv=1): NetCDF files that have time, downward longwave and shortwave radiation fluxes.

PySCHISM supports three types of atmpsheric datasets.

Notes: startdate should be one day earlier than the actual run startdate, because there is no data at t00z. Acoordingly, add two more extra days for rnday in the following scripts.
##**ECMWF ERA5**
ERA5 provides hourly estimates of a large number of atmoshperic, land and oceanic climate variables. The data cover the Earth on a 30km grid and resolve the atmoshpere using 137 levels from the surface up to a height of 80km. The dataset covers from 1950 to present. PySCHISM downloads ERA5 data through CDS API service. To do so, you need to install CDS API package. [Here](https://cds.climate.copernicus.eu/api-how-to) is the instruction about how to install the package.

The python script to generate sflux file from ERA5 is as follows:
```python
from datetime import datetime
import pathlib

from pyschism.mesh.hgrid import Hgrid
from pyschism.forcing.nws.nws2.era5 import ERA5

if __name__ == '__main__'
    startdate=datetime(2022, 4, 1)
    rnday = 10
    hgrid=Hgrid.open('./hgrid.gr3',crs='EPSG:4326')
    bbox = hgrid.get_bbox('EPSG:4326', output_type='bbox')
    outdir = pathlib.Path('./')

    er=ERA5()
    er.write(outdir=outdir, start_date=startdate, rnday=rnday, air=True, rad=True, prc=True, bbox=bbox, output_interval=interval, overwrite=True)
```
##**GFS**
The Global Forecast System (GFS) is weather forecast model produced by the National Centers for Environmental Prediction (NCEP). PySCHISM uses data hosted on AWS S3 bucket.

The python script to generate sflux from GFS is as follows:
```python
from datetime import datetime

from pyschism.mesh.hgrid import Hgrid
from pyschism.forcing.nws.nws2.gfs2 import GFS

if __name__ == '__main__':
    startdate = datetime(2022, 3, 31)
    rnday = 10
    record = 1
    hgrid = Hgrid.open('./hgrid.gr3', crs='epsg:4326')
    pscr = '/sciclone/pscr/lcui01/GFS/'
    gfs = GFS(start_date=startdate, rnday=rnday, pscr=pscr, record=record, bbox=hgrid.bbox)
```
Parameter *startdate* should be one day earlier than the actual startdate in the *param.nml* becasue GFS have nodata at t00, *pscr* is the pre-generated directory to save dowloaded raw data (grib2), *record* is to specify how many days in each file. For hindcast, recommend *record=1*. For forecast, the maximum is 5, because GFS has 5-day forecast.

##**HRRR**
The HRRR is a NOAA real-time 3-km resolution, hourly updated, cloud-resolving, convection-allowing atmospheric model. The AWS archived data starts from August 2014 to present.

The python script to generate sflux from HRRR is as follows:
```python
from datetime import datetime

from pyschism.mesh.hgrid import Hgrid
from pyschism.forcing.nws.nws2.hrrr3 import HRRR

if __name__ == '__main__':
    startdate = datetime(2022, 3, 31)
    rnday = 10
    record = 1
    hgrid = Hgrid.open('../../../data/hgrid.gr3', crs='epsg:4326')
    pscr='/sciclone/pscr/lcui01/HRRR/'

    hrrr = HRRR(start_date=startdate, rnday=rnday, pscr=pscr, record=record, bbox=hgrid.bbox)
```
Parameter *startdate*, *pscr*, and *record* are defined the same as GFS', except for forecast, the maximum of record for HRRR is 2, because HRRR only has 2-day forecast.
