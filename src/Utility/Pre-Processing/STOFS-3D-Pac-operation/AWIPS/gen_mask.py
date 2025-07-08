import os
import numpy as np
from pathlib import Path
from pylib import schism_grid  # from ZG's pylib: pip install pylibs-ocean


def gen_mask(gd:schism_grid, area_file='puertori'):
    xyz = np.loadtxt(area_file)
    out_mask_name = f'{Path(area_file).stem}_mask.txt'

    xyz[xyz[:, 0]>180, 0] -= 360

    ie, ip, acor = gd.compute_acor(np.c_[xyz[:, 0], xyz[:, 1]], fmt=0, out=0)
    ip[:] += 1  # node index start from 1
    ip[ie==-1, :] *= -1

    fmt = '%12d '*4 + '%20.10f '*3
    np.savetxt(out_mask_name, np.c_[np.array(range(len(ie)))+1, ip, acor], fmt=fmt)
    return out_mask_name, len(ie)

if __name__ == '__main__':
    # run script in the same folder as hgrid.gr3
    gd = schism_grid('./hgrid.gr3')
    gd.split_quads(angle_max=-1, angle_min=360)

    area_headers = {  # key is file name, value is header
        'alaska.xy': '',
        'conus_west.xy': '',
        'guam.xy': '',
        'hawaii.xy': '',
        'northpacific.xy': '',
    }
    for area_file, area_header in area_headers.items():
        output_name, npoints = gen_mask(gd, area_file=area_file)

        # insert header
        # if npoints != int(area_header.split()[0]):
        #     raise ValueError(f"Number of points in {area_file} does not match header: {area_header}")
        # os.system(f"sed -i '1i{area_header}' {output_name}")
