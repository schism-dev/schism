from datetime import datetime, timedelta
from time import time
import multiprocessing as mp
import logging

import numpy as np
import pandas as pd

from pyschism.mesh.hgrid import Hgrid
from pyschism.forcing.nws.nws2.gfs2 import GFS
from pyschism.dates import nearest_cycle

logging.basicConfig(
    format="[%(asctime)s] %(name)s %(levelname)s: %(message)s",
    force=True,
)
logging.captureWarnings(True)

log_level = logging.DEBUG
logging.getLogger('pyschism').setLevel(log_level)

if __name__ == "__main__":
    """
    Needs pooch and cfgrib
    """

    t0 = time()
    startdate = datetime(2021, 9, 25)

    pscr = '/sciclone/pscr/feiye/GFS/'
    rnday = 1
    record = 1

    hgrid = Hgrid.open('../hgrid.gr3', crs='epsg:4326')

    gfs = GFS(start_date=startdate, rnday=rnday, pscr=pscr, record=record, bbox=hgrid.bbox)

    print(f'It took {(time()-t0)/60} mins to process {rnday} days, {record*24} records/day')
