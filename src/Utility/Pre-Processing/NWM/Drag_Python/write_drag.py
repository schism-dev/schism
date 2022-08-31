#!/usr/bin/env python3
from pylib import *

def write_drag(grd,depths,mvalues,regions=None,rvalues=None,i_set_add_s=None,fname='drag.gr3'):
    '''
    write drag coefficient based on depth and specfied regions

    Input:
        grd:     grid name (*.gr3 or *.npz, where *.npz is python format)
        depths:  two hgrid depth(m) to distinguish river, land and the transition zone
        mvalues: lower and upper limits of drag values

        regions: list of region file names. e.g. regions=('GoME_1.reg','GoME_2.reg')
        rvalues: list of values for each region:  e.g. (0.1,0.5)
        i_set_add_s: identifier for setting or adding value,  0: set value; 1: add value
    '''
    #read hgrid
    if grd.endswith('.npz'):
        gd=loadz(grd).hgrid
    else:
        gd=read_schism_hgrid(grd)

    #compute drag coefficients
    mval=mvalues[0]+(gd.dp-depths[0])*(mvalues[1]-mvalues[0])/(depths[1]-depths[0])
    fpm=mval<mvalues[0]; mval[fpm]=mvalues[0]
    fpm=mval>mvalues[1]; mval[fpm]=mvalues[1]

    #set or add values in regions
    if regions is not None:
        for i_set_add, rvalue, region in zip(i_set_add_s, rvalues, regions):
            bp=read_schism_bpfile(region,fmt=1)
            sind=inside_polygon(c_[gd.x,gd.y], bp.x,bp.y)

            if i_set_add==0:
                print(f'setting {rvalue} drag in {region}')
                fp=sind==1
                mval[fp]=rvalue
            else:
                print(f'adding {rvalue} drag in {region}')
                sind2=gd.dp>depths[0]  # additional condition: deeper water, dp > -1 m
                fp=(sind & sind2)==1
                mval[fp]=mval[fp]+rvalue

    #save drag.gr3
    gd.dp=mval
    gd.write_hgrid(fname)

if __name__=="__main__":
    #------------------------------------------------------------------------------
    #input
    #------------------------------------------------------------------------------
    #general parameters
    grd='../hgrid.gr3'      #grid name (*.gr3 or *.npz, where *.npz is python format)
    depths=[-1,-3]       #two hgrid depth(m) to distinguish river, land and the transition zone
    mvalues=[0.0025,0.025]  #lower and upper limits of drag values

    #regions for certain values
    regions=['GoME+0.001.reg', 'Lake_Charles_0.reg']
    rvalues=[0.001, 0.0]
    i_set_add_s=[1, 0]  # 0: set value; 1: add value

    write_drag(grd,depths,mvalues,regions,rvalues,i_set_add_s,fname='../drag.gr3')
