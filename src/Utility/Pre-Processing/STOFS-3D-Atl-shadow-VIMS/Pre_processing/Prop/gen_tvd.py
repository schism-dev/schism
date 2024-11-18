'''
Generate tvd.prop file for the high-order transport schemes in SCHISM
'''
import numpy as np

import socket
import pylib
from pylib import read_schism_bpfile, inside_polygon
if 'sciclone' in socket.gethostname():
    from pylib_experimental.schism_file import cread_schism_hgrid as schism_read
else:
    from pylib import schism_grid as schism_read


def gen_tvd_prop(hg: pylib.schism_grid, regions: list):
    '''
    Generate tvd.prop file for the high-order transport schemes in SCHISM

    Parameters
    ----------
    hg : schism_grid
        SCHISM grid object
    regions : list
        list of region files,
        in the final output tvd.prop: 1 for inside region, 0 for outside region

    Returns
    -------
    None

    '''
    hg.compute_ctr()
    tvd = np.zeros(hg.ne)

    for region in regions:
        bp = read_schism_bpfile(region, fmt=1)
        idx = inside_polygon(np.c_[hg.xctr, hg.yctr], bp.x, bp.y)
        tvd[idx == 1] = 1

    hg.write_prop('tvd.prop', value=tvd, fmt='{:d}')


def sample_usage():
    '''
    Sample usage of gen_tvd_prop
    '''
    gd = schism_read('hgrid.gr3')

    regions = [
        'tvd0_1.reg', 'tvd0_2.reg', 'tvd0_3.reg', 'tvd0_4.reg',
        'tvd0_5.reg', 'tvd0_6.reg', 'tvd0_7.reg',
        'upwind_east_Carribbean.rgn', 'upwind_west_Carribbean.rgn'
    ]

    gen_tvd_prop(gd, regions)


if __name__ == '__main__':
    sample_usage()
