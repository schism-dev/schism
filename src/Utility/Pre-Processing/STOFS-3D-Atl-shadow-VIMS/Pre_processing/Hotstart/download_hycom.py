from datetime import datetime
import logging

from pyschism.forcing.hycom.hycom2schism import DownloadHycom
from pyschism.mesh.hgrid import Hgrid

logging.basicConfig(
    format="[%(asctime)s] %(name)s %(levelname)s: %(message)s",
    force=True,
)
logger = logging.getLogger('pyschism')
logger.setLevel(logging.INFO)

if __name__ == '__main__':
    hgrid = Hgrid.open('./hgrid.gr3', crs='epsg:4326')
    date=datetime(2022, 1, 1)

    hycom = DownloadHycom(hgrid)
    hycom.fetch_data(date)
    
