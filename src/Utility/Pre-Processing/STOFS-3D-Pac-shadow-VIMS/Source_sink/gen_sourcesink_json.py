from datetime import datetime
from time import time
import pathlib
import logging

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

    hgrid = Hgrid.open("./hgrid_sub.gr3", crs="epsg:4326")
    lon = hgrid.nodes.coords[:,0]
    hgrid.nodes.coords[:,0] = (lon + 180) % 360 - 180
    t0 = time()

    #source/sink json files, if not exist, it will call NWMElementPairings to generate.
    sources_pairings = pathlib.Path('./sources_conus_v3.json')
    sinks_pairings = pathlib.Path('./sinks_conus_v3.json')
    output_directory = pathlib.Path('./')

    #input directory which saves downloaded nc files
    #cache = pathlib.Path(f'./{startdate.strftime("%Y%m%d")}')
    #cache.mkdir(exist_ok=True, parents=True)
    #cache = pathlib.Path('./NWM_v2.1')

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

    pairings.sources_gdf.to_file('sources_conus_v3.shp')
    pairings.sinks_gdf.to_file('sinks_conus_v3.shp')
    #check nc files, if not exist will download
    #nwm=NationalWaterModel(pairings=pairings, cache=cache)

    #nwm.write(output_directory, hgrid, startdate, rnday, overwrite=True)
    print(f'It took {time()-t0} seconds to generate source/sink json files')
