#**National Water Model**
NOAA NWM CONUS Retrospective Dataset is to provide streamflow when if_source=1 is defined in *param.nml*. This dataset covers from February 1979 to present, combination of different versions of NWM dataset. These four files are needed:   
source_sink.in    
vsource.th   
vsink.th    
msource.th

The python script to generate these files is as follows:
```python
from datetime import datetime

from pyschism.mesh import Hgrid
from pyschism.forcing.source_sink.nwm import NationalWaterModel, NWMElementPairings

if __name__ == '__main__':
    startdate = datetime(2022, 4, 4)
    rnday = 10
    hgrid = Hgrid.open("./hgrid.gr3", crs="epsg:4326")
    sources_pairings = pathlib.Path('./sources.json')
    sinks_pairings = pathlib.Path('./sinks.json')
    output_directory = pathlib.Path('./')

    cache = pathlib.Path(f'./{startdate.strftime("%Y%m%d")}')
    cache.mkdir(exist_ok=True, parents=True)

    if all([sources_pairings.is_file(), sinks_pairings.is_file()]) is False:
        pairings = NWMElementPairings(hgrid)
        sources_pairings.parent.mkdir(exist_ok=True, parents=True)
        pairings.save_json(sources=sources_pairings, sinks=sinks_pairings)
    else:
        pairings = NWMElementPairings.load_json(hgrid, sources_pairings, sinks_pairings)

    nwm=NationalWaterModel(pairings=pairings, cache=cache)
    nwm.write(output_directory, hgrid, startdate, rnday, overwrite=True)
```