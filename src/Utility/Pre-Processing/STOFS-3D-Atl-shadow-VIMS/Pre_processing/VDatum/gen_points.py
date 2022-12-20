import os

import numpy as np
import pandas as pd

from pylib import read_schism_hgrid, loadz, read_schism_bpfile, inside_polygon

if __name__ == '__main__':

    #dem files with navd datum
    navd = ["North_Carolina_USGS_3m", "al_ll","nc_ll","fl_ll","gulf_1_dem_usgs",
        "gulf_3_demcombined_ll","ge_ll","sc_ll", "cb_ll","db_ll",
        "new_england_topobathy_dem_3m_dd","Tile3_R3_DEMv2","cb_bay_dem_v3.1_ll",
        "ncei19_2019v1","tx_ncei19","ncei19_MS_LA", "ncei19_MA_NH_NE", 
        "ncei19_AL_nwFL", "ncei19_FL"]

    #build a dictionary between dem_id and {'navd', 'msl'}
    ids = {}
    names = ['cb', 'db', 'gulf_3_demcombined', 'new_england'] #files are NAVD but cannot be distinguished by regex with '_' 
    #names = ['gulf_3_demcombined', 'new_england']
    
    f = open('hgrid.ll.new_dem_name', 'r')
    lines = f.readlines()
    for line in lines:
        dem_id = (line.split(':')[0])
        dem_name = '_'.join(line.split(':')[1].strip().split('_')[:-1])
        if dem_name.strip() in navd or dem_name.strip().startswith(names[0], 0, len(names[0]))  or \
            dem_name.strip().startswith(names[1], 0, len(names[1])):
            #ids[dem_id] = ['navd', dem_name]
            ids[dem_id] = 'navd'
        else:
            #ids[dem_id] = ['msl', dem_name]
            ids[dem_id] = 'msl'

    f.close()

    #read hgrid.ll.new_dem_id
    f = open('hgrid.ll.new_dem_id', 'r') 
    lines = f.readlines()
    NP = len(lines)
    nodes = np.arange(0, NP)

    # get node ids which receive depth from NAVD dem files
    navd_ids = []
    for i, line in enumerate(lines):
        dem_id = line.split('\n')[0]
        #print(ids[dem_id.strip()])
        if ids[dem_id.strip()] == 'navd': 
            navd_ids.append(1)
        else:
            navd_ids.append(0)
    print(len(navd_ids))

    #read hgrid
    if os.path.exists('grid.npz'):
        gd = loadz('grid.npz').hgrid 
    else:
        gd = read_schism_hgrid('hgrid.ll.new')
        gd.save('grid.npz')

    #exclude Ches/Del BAY
    bp1 = read_schism_bpfile('chea_del_bay.reg', fmt=1)
    sind1 = inside_polygon(np.c_[gd.x, gd.y], bp1.x, bp1.y)
    idxs1 = (sind1 == 0)

    #nodes receiving NAVD dem
    idxs = np.where(np.array(navd_ids) & idxs1)

    #save nodes with navd to txt (without chesapeake/delware bay)              
    lon = gd.x[idxs]
    lat = gd.y[idxs]
    depth = gd.dp[idxs]

    with open("hgrid_nochesdel_navd_positiveup.txt", "w") as f:
        f.write("\n".join(" ".join(map(str, line)) for line in zip(idxs[0] + 1, lon, lat, -depth)))

    #chesapeake/delware bay
    idxs1 = (sind1 == 1)
    idxs = np.where(np.array(navd_ids) & idxs1)
    lon = gd.x[idxs]
    lat = gd.y[idxs]
    depth = gd.dp[idxs]

    with open("hgrid_chesdel_navd_positiveup.txt", "w") as f:
        f.write("\n".join(" ".join(map(str, line)) for line in zip(idxs[0]+1, lon, lat, -depth)))
