import os
import glob

import numpy as np
import pandas as pd

from pylib import read_schism_hgrid, loadz

if __name__ == '__main__':

    gd = read_schism_hgrid('hgrid_dem_levee_loaded.gr3')
    depth_navd = gd.dp.copy()

    datum = 'xGEOID20b'
    #files = glob.glob('./result_pacific/*.txt')
    files = glob.glob('./vdatum/result/*.txt')
    files.sort()

    depth_diff = np.empty(len(gd.dp))
    depth_diff[:] = np.NaN

    for fname in files:
        print(fname)

        with open(fname) as f:
            df = pd.read_csv(f, header=None, delim_whitespace=True, names=['id', 'lon', 'lat', 'depth'], na_values=-999999.0)

        idxs = df.index[~df['depth'].isnull()]
        depth_geoid = df['depth'][idxs]
        #df.head()

        #dp0 = gd.dp[df['id'][idxs.values]-1]
        #dp1 = -depth_geoid.values #reverse the depth to be consistent with schism

        gd.dp[df['id'][idxs.values]-1] = depth_geoid

        #depth_diff[df['id'][idxs.values]-1] = x - y

    gd.write_hgrid(f'hgrid_{datum}.gr3')

    depth_diff = gd.dp - depth_navd
    gd.dp = depth_diff
    gd.save(f'depth_NAVD-{datum}.gr3')

