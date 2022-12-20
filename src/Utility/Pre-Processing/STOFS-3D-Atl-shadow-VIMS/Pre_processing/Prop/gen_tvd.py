import os
import numpy as np

from pylib import schism_grid, read_schism_bpfile, inside_polygon

if __name__ == '__main__':

    #Read hgrid
    hgrid_name='hgrid.gr3'
    gr3_pkl = 'hgrid.pkl'
    if os.path.exists(gr3_pkl):
        print('Reading hgrid.pkl')
        gd = schism_grid(gr3_pkl)
    else:
        print('Reading hgrid.gr3')
        gd = schism_grid(hgrid_name)
        gd.save(gr3_pkl)


    #gd = read_schism_hgrid('hgrid.gr3')
    gd.compute_ctr()

    regions = ['tvd0_1.reg', 'tvd0_2.reg', 'tvd0_3.reg', 'tvd0_4.reg', \
        'tvd0_5.reg', 'tvd0_6.reg', 'tvd0_7.reg', \
        'upwind_east_Carribbean.rgn', 'upwind_west_Carribbean.rgn']

    #write constant (=1) *prop
    ab = np.ones(gd.ne)

    #reset values (=0) within region
    for region in regions:
        bp = read_schism_bpfile(region, fmt=1)
        sind = inside_polygon(np.c_[gd.xctr,gd.yctr], bp.x, bp.y)
        ab[sind==1] = 0

    gd.write_prop('tvd.prop',value=ab, fmt='{:d}')

