#**Bctides**
Please refer to [this page](https://schism-dev.github.io/schism/master/input-output/bctides.html) for detailed horizontal B.C. and nudging options supported by SCHISM.

PySCHISM supports both TPXO and FES2014 tidal database. Please download the data (you may need to register) and save it as:

For [TPXO](https://www.tpxo.net/tpxo-products-and-registration):   
    ~/.local/share/tpxo/h_tpxo9.v1.nc   
    ~/.local/share/tpxo/u_tpxo9.v1.nc


For [FES2014](https://www.aviso.altimetry.fr/en/data/products/auxiliary-products/global-tide-fes/description-fes2014.html):  
    ~/.local/share/fes2014/eastward_velocity/  
    ~/.local/share/fes2014/northward_velocity/  
    ~/.local/share/fes2014/ocean_tide_extrapolated/  

*Bctides class*:   
```python

arguments:
    hgrid: required, hgrid.ll (lon/lat)
    flags: requried, bctypes for each boundary [[iettype, ifltype, itetype, isatype], [...], [...], ...]
    constituents: optional, default is "major", which is eight major tidal constituents
    dabase: optional, default is "tpxo"
    add_earth_tidal: optional, default is True
    cutoff_depth: optional, default is 50.0
    ethconst: optional, needed when using constant elevation along the boundary (iettype=2), which is given as [v1, v2, v3, ...]
    vthconst: optional, needed when using constant discharge (ifltype=2)
    tthconst: optional, needed when using constant temperatur (itetype=2)
    sthconst: optional, needed when using constant temperature (isatype=2)
    tobc: optional, nudging factor, needed when itetype is not 0, [v1, v2, v3, ...]
    sobc: optional, nudging factor for salinity, needed when isatype is not 0, [v1, v2, v3, ...]
    relax: optional, needed when ifttype=-4, [rel1, rel2]
```
Below shows an example script to generate bctides.in for hgrid with three open boundaries (2 ocean + 1 river). Bctypes for ocean boundaries are [5, 5, 4, 4], and for river is [0, 1, 1, 2]:

```python

from datetime import datetime
import numpy as np

from pyschism.mesh import Hgrid
from pyschism.forcing.bctides import Bctides

if __name__ == '__main__':
    start_date = datetime(2017, 12, 1)
    rnday = 396
    bctypes = [[5, 5, 4, 4], [5, 5, 4, 4], [0, 1, 1, 2]]
    constituents = 'major'
    database = 'tpxo'
    earth_tidal_potential = True
    sthconst = [np.nan, np.nan, 0]
    tobc = [0.5, 0.5, 1]
    sobc = [0.5, 0.5, 1]
    outdir = './'
    hgrid = Hgrid.open("./hgrid.gr3", crs="epsg:4326")

    bctides = Bctides(
        hgrid = hgrid,
        flags = bctypes,
        constituents = constituents,
        database = database,
        add_earth_tidal = earth_tidal_potential,
        sthconst = sthconst,
        tobc = tobc,
        sobc = sobc,
    )

    bctides.write(
        outdir, 
        start_date=start_date, 
        rnday=rnday, 
        overwrite=True,
    )
```
