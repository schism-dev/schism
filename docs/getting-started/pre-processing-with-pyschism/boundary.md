#**Boundary Condition from HYCOM**
Please refer to [this page](https://schism-dev.github.io/schism/master/input-output/bctides.html) for detailed horizontal B.C. and nudging options supported by SCHISM.

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
