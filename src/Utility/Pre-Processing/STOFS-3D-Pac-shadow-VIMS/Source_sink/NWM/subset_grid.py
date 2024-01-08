import json
import numpy as np

from pylib import inside_polygon, read_schism_hgrid

def subset_grid(hgrid, name, bbox):

    #read hgrid
    gd = read_schism_hgrid(hgrid)
    gd.compute_ctr()

    #build polygon
    px = np.array([bbox[0], bbox[1], bbox[1], bbox[0]])
    py = np.array([bbox[3], bbox[3], bbox[2], bbox[2]])

    #indexes of elements inside the box
    eidxs = np.nonzero(inside_polygon(np.c_[gd.xctr, gd.yctr], px, py) == 1)[0]
    gd.elnode = gd.elnode[eidxs]

    #map sub_grid eid to global eid
    local_to_global = {}
    for i, ele in enumerate(eidxs):
        local_to_global[str(i+1)] = str(ele + 1)
    with open(f'local_to_global_{name}.json', 'w') as f: 
        json.dump(local_to_global, f)

    #construct subgrid
    nidxs = np.unique(gd.elnode)
    fpn = nidxs>=0
    nidxs = nidxs[fpn]
    gd.x = gd.x[nidxs]
    gd.y = gd.y[nidxs]
    gd.dp = gd.dp[nidxs]
    gd.ne = len(eidxs)
    gd.np = len(nidxs)
    gd.i34 = gd.i34[eidxs]

    #construct new element connectivity
    node2node = dict(zip(nidxs, np.arange(gd.np)))
    for i in np.arange(gd.elnode.size):
        if gd.elnode.ravel()[i]<0: continue
        gd.elnode.ravel()[i] = node2node[gd.elnode.ravel()[i]]

    #compute side
    gd.ns, isidenode, isdel = gd.compute_side(fmt=1)
    gd.save(f'hgrid_sub_{name}.gr3')


    gd.x = (gd.x + 180)%360 - 180
    gd.save(f'hgrid_sub_convert_lon_{name}.gr3')

if __name__ == '__main__':

    ##west coast
    #lon_min = 230.
    #lon_max = 360.
    #lat_min = 0.
    #lat_max = 70.

    ##Hawaii
    #lon_min = 198.
    #lon_max = 208.
    #lat_min = 0.
    #lat_max = 30.

    #Alaska
    #lon_min = 200.
    #lon_max = 223.
    #lat_min = 56.
    #lat_max = 65.

    hgrid = 'hgrid.gr3'

    layers = {
        'conus': [230, 360, 0, 70], 
        'hawaii': [198, 208, 0, 30],
        'alaska': [200, 223, 56, 65]
    }


    for name, bbox in layers.items():
        subset_grid(hgrid, name, bbox)
