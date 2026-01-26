import os

import numpy as np
import pandas as pd

from pylib import read_schism_hgrid, loadz

if __name__ == '__main__':

    gd = read_schism_hgrid('hgrid.ll')
    depth_navd88 = gd.dp

    #datums = ['LMSL', 'xGEOID20b']
    datum = 'xGEOID20b'
    fids = ['1', '2', '3', '4', 'ches_del']

    depth_diff = np.empty(len(gd.dp))
    depth_diff[:] = np.NaN

    for fid in fids:

        with open(f'hgrid_stofs3d_inland_{fid}_New.txt', 'r') as f:
            df = pd.read_csv(f, header=None, delim_whitespace=True, names=['id', 'lon', 'lat', 'depth'], na_values=-999999.0)

        idxs = df.index[~df['depth'].isnull()]
        depth_lmsl = df['depth'][idxs]
        df.head()


        x = gd.dp[df['id'][idxs.values]-1]
        y = -depth_lmsl.values #reverse the depth to be consistent with schism

        gd.dp[df['id'][idxs.values]-1] = y

        depth_diff[df['id'][idxs.values]-1] = x - y

    gd.write_hgrid(f'hgrid_{datum}.gr3')

    gd.dp = depth_diff
    gd.save(f'depth_NAVD-{datum}.gr3')

