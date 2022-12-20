import os

import numpy as np
import pandas as pd

from pylib import read_schism_hgrid, loadz

if __name__ == '__main__':

    if os.path.exists('grid.npz'):
        gd = loadz('grid.npz').hgrid
    else:
        gd = read_schism_hgrid('hgrid.gr3')
        gd.save('grid.npz') 

    #depth_navd88 = gd.dp

    #datums = ['LMSL', 'xGEOID20b']
    datums = ['xGEOID20b',]
    fids = ['nochesdel_navd', 'chesdel_navd']

    for  datum in datums:

        depth_diff = np.empty(len(gd.dp))
        depth_diff[:] = np.NaN

        for fid in fids:

            with open(f'hgrid_{fid}_positiveup_{datum}.txt', 'r') as f:
                df = pd.read_csv(f, header=None, delim_whitespace=True, names=['id', 'lon', 'lat', 'depth'], na_values=-999999.0)

            idxs = df.index[~df['depth'].isnull()]
            depth_lmsl = df['depth'][idxs]

            x = gd.dp[df['id'][idxs.values]-1]
            y = -depth_lmsl.values #reverse the depth to be consistent with schism

            gd.dp[df['id'][idxs.values]-1] = y

            depth_diff[df['id'][idxs.values]-1] = x - y


        gd.save(f'hgrid_{datum}.gr3')

        gd.dp = depth_diff
        gd.save(f'depth_NAVD-{datum}.gr3')



