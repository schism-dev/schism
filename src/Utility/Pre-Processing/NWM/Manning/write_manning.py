#!/usr/bin/env python3
from pylib import *

def write_manning(grd,depths,mvalues,regions=None,rvalues=None,fname='manning.gr3'):
    '''
    write manning coefficient based on depth and specfied regions

    Input:
        grd:     grid name (*.gr3 or *.npz, where *.npz is python format)
        depths:  two hgrid depth(m) to distinguish river, land and the transition zone
        mvalues: lower and upper limits of manning values

        regions: list of region file names. e.g. regions=('GoME_1.reg','GoME_2.reg')
        rvalues: list of values for each region:  e.g. (0.1,0.5)
    '''
    #read hgrid
    if grd.endswith('.npz'):
        gd=loadz(grd).hgrid
    else:
        gd=read_schism_hgrid(grd)

    #compute manning coefficients
    mval=mvalues[0]+(gd.dp-depths[0])*(mvalues[1]-mvalues[0])/(depths[1]-depths[0])
    fpm=mval<mvalues[0]; mval[fpm]=mvalues[0]
    fpm=mval>mvalues[1]; mval[fpm]=mvalues[1]

    #set values in regions
    if regions is not None:
        for m,region in enumerate(regions):
            print('modifying manning in {}'.format(region))
            bp=read_schism_bpfile(region,fmt=1)
            sind=inside_polygon(c_[gd.x,gd.y], bp.x,bp.y)
            fp=sind==1; mval[fp]=rvalues[m]

    #save manning.gr3
    gd.dp=mval
    gd.write_hgrid(fname)

if __name__=="__main__":
    #------------------------------------------------------------------------------
    #input
    #------------------------------------------------------------------------------
    #general parameters
    grd='hgrid.gr3'      #grid name (*.gr3 or *.npz, where *.npz is python format)
    depths=[-1,-3]       #two hgrid depth(m) to distinguish river, land and the transition zone
    mvalues=[0.02,0.05]  #lower and upper limits of manning values

    #regions for certain values
    regions=('GoME_1.reg','GoME_2.reg','Berwick.reg','BulkTerminal.reg','Pilottown.reg')
    rvalues=(0.2,0.2,0.005,0.005,0.005)

    write_manning(grd,depths,mvalues,regions,rvalues,fname='../manning.gr3')
