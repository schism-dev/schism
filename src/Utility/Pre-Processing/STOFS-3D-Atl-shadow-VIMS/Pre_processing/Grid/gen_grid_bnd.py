from pylib_experimental.schism_file import cread_schism_hgrid
from pathlib import Path

hgrid_fname = '/sciclone/schism10/feiye/STOFS3D-v8/I02/Bathy_edit/DEM_loading/hgrid_max_dp.gr3'

hg = cread_schism_hgrid(hgrid_fname)

open_boundaries = {
    'ocean': [[-94.14734925, 29.59769601], [-88.96698021, 30.35921140]]
}

bxy = []
for key, bnd_xy in open_boundaries.items():
    bxy.append([bnd_xy[0][0], bnd_xy[1][0], bnd_xy[0][1], bnd_xy[1][1]])

hg.compute_bnd(bxy=bxy, method=1) 

hg.save(str(Path(hgrid_fname).with_suffix('.bnd.gr3')), fmt=1)
