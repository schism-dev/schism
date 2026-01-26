import os

import numpy as np

from pylib import read_schism_hgrid, loadz, read_schism_bpfile, inside_polygon

if __name__ == '__main__':

    gd = read_schism_hgrid('hgrid.ll')

    ##get nodes inside ches_del_bay polgyon
    bp = read_schism_bpfile('chea_del_bay.reg', fmt=1)
    sind = inside_polygon(np.c_[gd.x, gd.y], bp.x, bp.y)

    idxs = np.where(sind == 1)[0]
    lon = gd.x[idxs]
    lat = gd.y[idxs]
    depth = -gd.dp[idxs]

    with open("hgrid_stofs3d_inland_ches_del.txt", "w") as f:
        f.write("\n".join(" ".join(map(str, line)) for line in zip(idxs+1, lon, lat, depth)))

    # get nodes inside stofs3d polygon but not in ches_del_bay polygon
    bp2 = read_schism_bpfile('stofs3d_inland2.reg', fmt=1)
    sind2 = inside_polygon(np.c_[gd.x, gd.y], bp2.x, bp2.y)
    #idxs2 = (sind2 == 1)
    idxs = np.where(~sind & sind2 == 1)[0]

    lon = gd.x[idxs]
    lat = gd.y[idxs]
    depth = -gd.dp[idxs]
  
    #seperate into multiple files to speed up conversion
    nfile = int(len(lon)/500000) + 1
    print(nfile)
 
    for i in np.arange(nfile):

        x1 = i*500000
        if i < 5:
            x2 = (i+1)*500000
        else:
            x2 =  gd.x.shape[0]
        print(f'file {i+1}, x1: {x1}, x2: {x2}')
        with open(f"hgrid_stofs3d_inland_{i+1}.txt", "w") as f:
            f.write("\n".join(" ".join(map(str, line)) for line in zip(idxs[x1:x2]+1, lon[x1:x2], lat[x1:x2], depth[x1:x2])))


