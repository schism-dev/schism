#! /usr/bin/env python
from datetime import datetime, timedelta
import logging

from pyschism.mesh import Hgrid, Vgrid
from pyschism.forcing.hycom.hycom2schism import Nudge


logging.basicConfig(level=logging.INFO, force=True)


if __name__ == '__main__':

    hgrid=Hgrid.open('./hgrid.gr3', crs='epsg:4326')
    vgrid='./vgrid.in'
    outdir='./'
    start_date = datetime(2022, 1, 1)
    rnday=171

    nudge=Nudge()
    nudge.fetch_data(outdir, hgrid, vgrid, start_date, rnday)
