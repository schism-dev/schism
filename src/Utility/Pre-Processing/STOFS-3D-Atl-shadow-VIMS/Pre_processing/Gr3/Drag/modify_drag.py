import numpy as np

from pylib import read_schism_hgrid, inside_polygon,read_schism_bpfile

if __name__ == '__main__':
    gd = read_schism_hgrid('hgrid.gr3')
    drag = read_schism_hgrid('drag.gr3')

    #set value to zero for Lake Charles
    bp = read_schism_bpfile('Lake_Charles_0.reg', fmt=1)
    sind = inside_polygon(np.c_[gd.x, gd.y], bp.x, bp.y)
    fp = sind==1
    drag.dp[fp] = 0

    #modify drag for GoME
    bfric1=0.01
    bfric2=0.02
    depths = [20, 5]
    bp = read_schism_bpfile('GoME2.reg', fmt=1)
    sind = inside_polygon(np.c_[gd.x, gd.y], bp.x, bp.y)

    #-1<dp<5
    sind2 = gd.dp>-1 
    sind3 = gd.dp<depths[-1]
    fp = (sind & sind2 & sind3) == 1
    drag.dp[fp] = bfric2

    #dp>20
    sind2 = gd.dp>depths[0]
    fp = (sind & sind2) == 1
    drag.dp[fp] = bfric1
 
    #5<dp<20
    mvalues = [bfric1, bfric2]
    sind2 = gd.dp>depths[-1]
    sind3 = gd.dp<depths[0]
    fp = (sind & sind2 & sind3) == 1
    mval = mvalues[0] + (gd.dp - depths[0])*(mvalues[1]-mvalues[0])/(depths[1]-depths[0])
    drag.dp[fp] = mval[fp]
    print(gd.dp[fp].min())
    print(gd.dp[fp].max())

    #limit Chatham area not excess than 0.0025
    bp = read_schism_bpfile('Chatham_max0.0025.reg', fmt=1)
    sind = inside_polygon(np.c_[gd.x, gd.y], bp.x, bp.y)
    fp = sind==1
    drag.dp[fp] = np.minimum(drag.dp[fp], 0.0025)
    
    drag.write_hgrid('drag_modified2.gr3')
