import glob

import numpy as np

from pylib import read_schism_hgrid, inside_polygon

if __name__ == '__main__':

    gd = read_schism_hgrid('hgrid_dem_levee_loaded.gr3')

    #bnd files
    path = './VDatum_polygons/'
    files = glob.glob(f'{path}/*.bnd')

    print(files[0]) 


    for fname in files:
        print(fname)

        ##get nodes inside polgyon
        bp = np.loadtxt(fname)
        sind = inside_polygon(np.c_[gd.x, gd.y], bp[:, 0], bp[:, 1])

        idxs = np.where(sind == 1)[0]
        print(idxs.shape[0])
        if idxs.shape[0] == 0: continue
        lon = gd.x[idxs]
        lat = gd.y[idxs]
        #depth = -gd.dp[idxs] #for height
        depth = gd.dp[idxs] #for sounding

        name = fname.split('/')[-1].split('_')[0].lower()
        with open(f"hgrid_secofs_{name}.txt", "w") as f:
            f.write("\n".join(" ".join(map(str, line)) for line in zip(idxs+1, lon, lat, depth)))

