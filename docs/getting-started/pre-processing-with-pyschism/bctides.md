#**Bctides**
Please refer to [this page](https://schism-dev.github.io/schism/master/input-output/bctides.html) for detailed horizontal B.C. and nudging options supported by SCHISM.

PySCHISM supports both TPXO and FES2014 tidal database. Please download the data (you may need to register) and save it as:

For TPXO:   
    ~/.local/share/tpxo/

For FES2014:  
    ~/.local/share/fes2014/eastward_velocity/  
    ~/.local/share/fes2014/northward_velocity/  
    ~/.local/share/fes2014/ocean_tide_extrapolated/  

Generating *bctides.in*:   
```python
from datetime import datetime
from pyschism.mesh import Hgrid
from pyschism.forcing.bctides import Bctides, iettype, ifltype, isatype, itetype

if __name__ == '__main__':
    start_date = datetime(2017, 12, 1)
    rnday = 10
    outdir = './'

    hgrid = Hgrid.open("./hgrid.gr3", crs="epsg:4326")

    #elevation
    #constituents: major [M2, S2, N2, K2, K1, O1, P1, Q1]
    #database: 'tpxo' or 'fes2014'
    iet3 = iettype.Iettype3(constituents='major', database='tpxo')
    iet4 = iettype.Iettype4()
    iet5 = iettype.Iettype5(iettype3=iet3, iettype4=iet4)

    #velocity
    ifl3 = ifltype.Ifltype3(constituents='major', database='tpxo')
    ifl4 = ifltype.Ifltype4()
    ifl5 = ifltype.Ifltype5(ifltype3=ifl3, ifltype4=ifl4)

    #salinity
    isa4 = isatype.Isatype4() #time-varying

    #temperature
    ite4 = itetype.Itetype4() #time-varying


    #example1 - only one ocean boundary and type is [5, 5, 4, 4]
    #bctides = Bctides(hgrid, iettype=iet5, ifltype=ifl5, isatype=isa4, itetype=ite4)

    #example2 - two ocean open boundaries, one river boundary and type is:  
    #[5, 5, 4, 4] - first ocean boundary
    #[5, 5, 4, 4] - second ocean boudnary
    #[0, 0, 3, 3] - river boundary
    bctides = Bctides(
        hgrid,
        iettype = {'1': iet5, '2': iet5},
        ifltype = {'1': ifl5, '2': ifl5},
        isatype = {'1': isa4, '2':isa4, '3': isa2},
        itetype = {'1': ite4, '2':ite4, '3': ite2},
    )

    bctides.write(outdir, 
        start_date=start_date, 
        end_date=rnday, 
        bctides=True, 
        overwrite=True,
    )
```
