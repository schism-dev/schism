from datetime import datetime
from time import time
import pathlib

from pyschism.forcing.nws.nws2.era5 import ERA5

from pyschism.mesh.hgrid import Hgrid

if __name__ == '__main__':

    startdate=datetime(2017, 4, 26)
    rnday=10

    hgrid=Hgrid.open('./hgrid.gr3',crs='EPSG:4326')
    bbox = hgrid.get_bbox('EPSG:4326', output_type='bbox')

    outdir = pathlib.Path('./')

    #output interval
    interval = 1

    t0=time()

    er=ERA5()

    er.write(outdir=outdir, start_date=startdate, rnday=rnday, air=True, \
        rad=True, prc=True, bbox=bbox, output_interval=interval, overwrite=True)

    #write sflux_inputs.txt
    with open("./sflux_inputs.txt", "w") as f:
        f.write("&sflux_inputs\n/\n")

    print(f'It took {(time()-t0)/60} minutes to generate {rnday} days')
