#**Boundary Condition**
Please refer to [this page](https://schism-dev.github.io/schism/master/input-output/bctides.html) for detailed horizontal B.C. and nudging options supported by SCHISM.

##**TPXO**
PySCHISM uses TPXO Global Tidal Model database. You may need to register an account on the website to download the data.

Generating *bctides.in*:   
```python
from datetime import datetime
from pyschism.mesh import Hgrid
from pyschism.forcing.bctides import Bctides, iettype, ifltype, isatype, itetype

if __name__ == '__main__':
    startdate = datetime(2022, 4, 1)
    rnday = 10
    hgrid = Hgrid.open("./hgrid.gr3", crs="epsg:4326")
    iet3 = iettype.Iettype3(constituents='major', database='tpxo')
    iet4 = iettype.Iettype4()
    iet5 = iettype.Iettype5(iettype3=iet3, iettype4=iet4)
    ifl3 = ifltype.Ifltype3(constituents='major', database='tpxo')
    ifl4 = ifltype.Ifltype4()
    ifl5 = ifltype.Ifltype5(ifltype3=ifl3, ifltype4=ifl4)
    isa3 = isatype.Isatype4()
    ite3 = itetype.Itetype4()
    bctides = Bctides(hgrid, iettype=iet5, ifltype=ifl5, isatype=isa3, itetype=ite3)
    bctides.write('./', startdate, rnday, bctides=True, elev2D=False, uv3D=False, tem3D=False, sal3D=False, overwrite=True)
```

##**HYCOM**
Generating *elev.2D.th.nc*, *SAL_3D.th.nc*, *TEM_3D.th.nc*, and *uv3D.th.nc*:
```python
from datetime import datetime

from pyschism.mesh.hgrid import Hgrid
from pyschism.forcing.hycom.hycom2schism import OpenBoundaryInventory

if __name__ == '__main__':
    start_date = datetime(2022, 4, 1)
    rnday = 10
    hgrid = Hgrid.open('./hgrid.gr3', crs='epsg:4326')
    vgrid = './vgrid.in'
    outdir = './'
    bnd = OpenBoundaryInventory(hgrid, vgrid)
    bnd.fetch_data(outdir, start_date, rnday, elev2D=True, TS=True, UV=True)
```
Generating *SAL_nu.nc* and *TEM_nu.nc* for nudging:
```python
from datetime import datetime

from pyschism.mesh import Hgrid
from pyschism.forcing.hycom.hycom2schism import Nudge

if __name__ == '__main__':
    start_date = datetime(2022, 4, 1)
    rnday = 10
    hgrid = Hgrid.open('./hgrid.gr3', crs='epsg:4326')
    vgrid = './vgrid.in'
    outdir = './'

    nudge=Nudge()
    nudge.fetch_data(outdir, hgrid, vgrid, start_date, rnday)
```
This script also generates nudging coefficient files *SAL_nu.gr3* and *TEM_nu.gr3*.