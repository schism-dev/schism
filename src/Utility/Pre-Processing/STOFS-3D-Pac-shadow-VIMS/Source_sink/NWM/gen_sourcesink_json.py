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

    startdate = datetime(2017, 12, 1)
    rnday = 10
    layers = ['conus', 'hawaii', 'alaska']
    #layers = ['alaska']
    nwm_v3 = '/sciclone/schism10/lcui01/schism20/ICOGS/ICOGS3D/Forecast/NWM/oper_3D/NWM_v3_final/NWM_v3_hydrofabric.gdb'

    t0 = time()

    for layer in layers:
        reach_layer = f'reaches_{layer}'
        hgrid = Hgrid.open(f"./hgrid_sub_convert_lon_{layer}.gr3", crs="epsg:4326")

        #source/sink json files, if not exist, it will call NWMElementPairings to generate.
        sources_pairings = pathlib.Path(f'./sources_{layer}.json')
        sinks_pairings = pathlib.Path(f'./sinks_{layer}.json')


        # check if source/sink json file exists
        if all([sources_pairings.is_file(), sinks_pairings.is_file()]) is False:
            pairings = NWMElementPairings(hgrid, nwm_file=nwm_v3, reach_layer=reach_layer)
            sources_pairings.parent.mkdir(exist_ok=True, parents=True)
            pairings.save_json(sources=sources_pairings, sinks=sinks_pairings)
        else:
            pairings = NWMElementPairings.load_json(
                hgrid, 
                sources_pairings, 
                sinks_pairings)

        pairings.sources_gdf.to_file(f'sources_{layer}_v3.shp')
        pairings.sinks_gdf.to_file(f'sinks_{layer}_v3.shp')

    print(f'It took {time()-t0} seconds to generate source/sink')
