from datetime import datetime
from time import time
import pathlib
import logging

#from pyschism.forcing import source_sink
from pyschism.forcing.source_sink.nwm import NationalWaterModel, NWMElementPairings
from pyschism.mesh import Hgrid


logging.basicConfig(
    format="[%(asctime)s] %(name)s %(levelname)s: %(message)s",
    force=True,
)
logging.captureWarnings(True)

log_level = logging.DEBUG
logging.getLogger('pyschism').setLevel(log_level)

if __name__ == '__main__':

    startdate = datetime(2022, 1, 1)
    rnday = 10
    hgrid = Hgrid.open("./hgrid.gr3", crs="epsg:4326")

    t0 = time()

    #source/sink json files, if not exist, it will call NWMElementPairings to generate.
    sources_pairings = pathlib.Path('./sources.json')
    sinks_pairings = pathlib.Path('./sinks.json')
    output_directory = pathlib.Path('./')

    #input directory which saves nc files
    cache = pathlib.Path(f'./{startdate.strftime("%Y%m%d")}')
    cache.mkdir(exist_ok=True, parents=True)

    # check if source/sink json file exists
    if all([sources_pairings.is_file(), sinks_pairings.is_file()]) is False:
        pairings = NWMElementPairings(hgrid)
        sources_pairings.parent.mkdir(exist_ok=True, parents=True)
        pairings.save_json(sources=sources_pairings, sinks=sinks_pairings)
    else:
        pairings = NWMElementPairings.load_json(
            hgrid, 
            sources_pairings, 
            sinks_pairings)

    #check nc files, if not exist will download
    nwm=NationalWaterModel(pairings=pairings, cache=cache)

    nwm.write(output_directory, hgrid, startdate, rnday, overwrite=True)
    print(f'It took {time()-t0} seconds to generate source/sink')
