'''
Generate source/sink files from NWM,
source/sink parings between elements and NWM streams (fid)
are saved in sources.json and sinks.json,
which are reusable as long as the horizontal coordinates
of hgrid does not change.
'''

import os
from datetime import datetime
from time import time
from pathlib import Path
import logging

from pyschism.forcing.source_sink.nwm import NationalWaterModel, NWMElementPairings
from pyschism.mesh import Hgrid


# Configure logging
logger = logging.getLogger(__name__)
logger.setLevel(logging.DEBUG)


def gen_sourcesink_nwm(hgrid_fname: str, startdate: datetime, rnday: float, cache_folder: Path = None):
    '''
    Generate source/sink files from NWM
    hgrid must be in epsg:4326
    '''

    t0 = time()
    logger.info('Reading hgrid')
    hgrid = Hgrid.open(hgrid_fname, crs="epsg:4326")
    logger.info("Reading hgrid took %.2f seconds", time()-t0)

    t0 = time()
    # source/sink json files, if not exist, it will call NWMElementPairings to generate.
    sources_pairings_file = Path.cwd() / 'sources.json'
    sinks_parings_file = Path.cwd() / 'sinks.json'
    output_directory = Path.cwd()

    # input directory which saves nc files
    if cache_folder is not None and os.path.exists(cache_folder):  # if cache folder exists, use it
        cache = Path(cache_folder)
    else:  # if cache folder does not exist, create it
        cache = Path(f'./{startdate.strftime("%Y%m%d")}')
        cache.mkdir(exist_ok=True, parents=True)

    # check if source/sink json file exists
    if all([sources_pairings_file.is_file(), sinks_parings_file.is_file()]) is False:
        logger.info('Generating source/sink pairings')
        pairings = NWMElementPairings(hgrid)
        sources_pairings_file.parent.mkdir(exist_ok=True, parents=True)
        logger.info('Saving source/sink pairings to json files')
        pairings.save_json(sources=sources_pairings_file, sinks=sinks_parings_file)
    else:
        logger.info('Loading source/sink pairings from existing json files')
        pairings = NWMElementPairings.load_json(hgrid, sources_pairings_file, sinks_parings_file)

    # check nc files, if not exist will download
    nwm = NationalWaterModel(pairings=pairings, cache=cache)

    nwm.write(output_directory, hgrid, startdate, rnday, overwrite=True)
    logger.info('Generating source/sink files took %.2f seconds', time()-t0)

    return cache.resolve()  # absolute path of the cache folder


def main():
    '''Sample usage'''
    working_dir = Path('/sciclone/schism10/feiye/STOFS3D-v7/Inputs/I12z/Source_sink/relocated_source_sink/')
    cache_folder = Path('/sciclone/schism10/feiye/STOFS3D-v7/Inputs/I12z/Source_sink/original_source_sink/20240305/')
    os.chdir(working_dir)
    gen_sourcesink_nwm(f'{working_dir}/hgrid.gr3', datetime(2024, 3, 5), 5, cache_folder=cache_folder)


if __name__ == '__main__':
    main()
    print('Done')
