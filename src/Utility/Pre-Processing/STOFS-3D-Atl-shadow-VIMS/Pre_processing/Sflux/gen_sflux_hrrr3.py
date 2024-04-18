from datetime import datetime, timedelta
from time import time
import logging

from pyschism.mesh.hgrid import Hgrid
from pyschism.forcing.nws.nws2.hrrr3 import HRRR
from pyschism.dates import nearest_cycle

logging.basicConfig(
    format="[%(asctime)s] %(name)s %(levelname)s: %(message)s",
    force=True,
)
logging.captureWarnings(True)

log_level = logging.DEBUG
logging.getLogger('pyschism').setLevel(log_level)

if __name__ == "__main__":
    '''
    Need to install cfgrib first
    '''
    t0 = time()
    rnday = 1
    record = 2
    startdate = datetime(2021, 9, 25)
    hgrid = Hgrid.open('../hgrid.gr3', crs='epsg:4326')

    pscr='/sciclone/pscr/feiye/HRRR/'
    hrrr = HRRR(start_date=startdate, rnday=rnday, pscr=pscr, record=record, bbox=hgrid.bbox)
    print(f'It took {(time()-t0)/60} mins to process {rnday} days, {record*24} files/day')
