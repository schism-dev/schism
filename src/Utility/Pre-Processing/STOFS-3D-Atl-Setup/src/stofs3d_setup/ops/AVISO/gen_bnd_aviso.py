from datetime import datetime

from pyschism.mesh.hgrid import Hgrid
from aviso2schism import OpenBoundaryInventory


if __name__ == "__main__":
    hgrid=Hgrid.open('../hgrid.gr3', crs='epsg:4326')
    print(hgrid.bbox)
    outdir='./'
    start_date=datetime(2014, 12, 1)
    rnday=396

    #boundary
    bnd=OpenBoundaryInventory(hgrid)
    bnd.fetch_data(outdir, start_date, rnday, ocean_bnd_ids=[0, 1], elev2D=True)
