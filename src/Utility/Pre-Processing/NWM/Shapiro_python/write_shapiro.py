#!/usr/bin/env python3
from pylib import *


def write_shapiro(
    grd,
    shapiro_max=0.5, threshold_slope=0.5,
    depths=None, shapiro_vals1=None,
    regions=None, shapiro_vals2=None, i_set_add_s=None,
    fname='shapiro.gr3'
 ):
    '''
    write shapiro fileter strength value based on depth and specfied regions

    Input:
        grd:     grid name (*.gr3 or *.npz, where *.npz is python format)
        depths:  two hgrid depth(m) to distinguish river, land and the transition zone
        mvalues: lower and upper limits of shapiro values

        regions: list of region file names. e.g. regions=('GoME_1.reg','GoME_2.reg')
        rvalues: list of values for each region:  e.g. (0.1,0.5)
        i_set_add_s: identifier for setting or adding value,  0: set value; 1: add value
    '''
    #read hgrid
    if grd.endswith('.npz'):
        gd=loadz(grd).hgrid
    else:
        gd=read_schism_hgrid(grd)

    # compute bathymetry gradient on each node
    _, _, slope = gd.compute_gradient(fmt=2)

    # compute shapiro coefficients
    shapiro=shapiro_max*tanh(2*slope/threshold_slope)

    # further tweaks on shallow waters
    if len(depths) != len(shapiro_vals1):
        raise Exception(f'lengths of depths {len(depths)} and shapiro_vals1 {len(shapiro_vals1)} inconsistent')
    fp = gd.dp < depths[-1]
    shapiro[fp] = maximum(shapiro[fp], interp(gd.dp[fp], depths, shapiro_vals1))

    #set or add values in regions
    if regions is not None:
        for i_set_add, rvalue, region in zip(i_set_add_s, shapiro_vals2, regions):
            bp=read_schism_bpfile(region,fmt=1)
            sind=inside_polygon(c_[gd.x,gd.y], bp.x,bp.y).astype('bool')

            if i_set_add==0:
                print(f'setting {rvalue} shapiro in {region}')
                fp=sind
                shapiro[fp]=rvalue
            else:
                print(f'adding {rvalue} shapiro in {region}')
                sind2=(gd.dp>depths[0])  # additional condition: deeper water, dp > -1 m
                fp=(sind & sind2)
                shapiro[fp]=shapiro[fp]+rvalue

    #save shapiro.gr3
    gd.dp=shapiro
    gd.write_hgrid(fname)

if __name__=="__main__":
    rundir = '../'
    regdir = './'

    outfilename = f'{rundir}/shapiro.gr3'

    if os.path.exists(outfilename):
        os.remove(outfilename)

    write_shapiro(
        grd=f'{rundir}/hgrid.cpp',  # grid name (*.gr3 or *.npz, where *.npz is python format)
        shapiro_max=0.5,
        threshold_slope=0.5,
        depths=[-99999, 20, 50],  # tweaks in shallow waters
        shapiro_vals1=[0.2, 0.2, 0.05],  # tweaks in shallow waters
        regions=[f'{regdir}/coastal_0.2.cpp.reg',
                 f'{regdir}/coastal_0.5_1.cpp.reg',
                 f'{regdir}/coastal_0.5_2.cpp.reg'],  # tweaks in regions, the order matters
        shapiro_vals2=[0.2, 0.5, 0.5],  # tweaks in regions, the order matters
        i_set_add_s=[0, 0, 0],  # 0: set; 1: add
        fname=outfilename
    )
