#!/usr/bin/env python3

'''
This script is used to edit the DEM for hgrid.
This contains several steps:
'''
from pathlib import Path
import os

import numpy as np

from pylib import schism_grid, grd2sms, sms2grd
from pylib_essentials.schism_file import read_schism_hgrid_cached
from Levee.set_levees import set_levees
from NCF.load_NCF import load_NCF
from xGEOID.convert2xgeoid import convert2xgeoid
from Chart.load_chart import load_chart
from Dredge.dredge_auto_channels import dredge_auto_channels
from SetFeederDp.set_feeder_dp import set_feeder_dp

wdir = '/sciclone/schism10/feiye/STOFS3D-v7/Inputs/I11/Bathy_edit/'

# ----------------- load the hgrid -------------------
hgrid_obj = schism_grid(f'{wdir}/DEM_loading/hgrid_max_dp.gr3')  # after pload_depth.py
dp_dem = hgrid_obj.dp.copy()
print("finished reading the DEM-loaded hgrid")

# ----------------- load NCF (National Channel Framework) -------------------
hgrid_obj = load_NCF(hgrid_obj=hgrid_obj, NCF_shpfile=Path('/sciclone/schism10/Hgrid_projects/Charts/Mississippi/channel_quarter_NCF.shp'), expansion=20)
hgrid_obj.dp = np.maximum(dp_dem, hgrid_obj.dp)  # change the depth only if it is deeper than the original depth
grd2sms(hgrid_obj, f'{wdir}/NCF/hgrid.ll.dem_NCF.2dm')
print("finished loading NCF")

# ----------------- set levees -------------------
os.chdir(f'{wdir}/Levee') #to set the directory
hgrid_obj = set_levees(hgrid_obj=hgrid_obj, wdir=f'{wdir}/Levee/')
grd2sms(hgrid_obj, f'{wdir}/Levee/hgrid_dem_NCF_levee.2dm')
print("finished setting levees")


# ----------------- convert to xGEOID -------------------
hgrid_obj = convert2xgeoid(wdir=f'{wdir}/xGEOID/', hgrid_obj=hgrid_obj)  # the output is also saved at hgrid_xGEOID20b.gr3
temp_dp = hgrid_obj.dp.copy()
print("finihsed converting the vdatum to xGEOID")

# # ----------------- load Chart -------------------
# hgrid_obj = load_chart(
#     hgrid_obj=hgrid_obj,
#     sounding_shpfile=Path('/sciclone/schism10/Hgrid_projects/Charts/Savanna_Cooper/savannah_cooper_sounding_3_xyz_edited_xgeoid.shp'),
#     region_shpfile=Path('/sciclone/schism10/Hgrid_projects/Charts/Savanna_Cooper/secofs_chart_loading_zones.shp'), crs_region='esri:102008'
# )
# hgrid_obj.dp = np.maximum(temp_dp, hgrid_obj.dp)  # change the depth only if it is deeper than the original depth
# hgrid_obj.write_hgrid(f'{wdir}/Chart/hgrid_dem_levee_loaded_NCF_loaded_xGEOID20b_chart_loaded.gr3')
# print("finished loading Chart")
# 
# # ----------------- load dredging -------------------
# dredge_depth = 2 # set the dredge depth
# hgrid_obj = dredge_auto_channels(hgrid_obj=hgrid_obj, dredge_polygon_file=Path('/sciclone/schism10/Hgrid_projects/SECOFS/new20_JZ/total_river_polys_clipped_test.shp'), dredge_depth=dredge_depth)
# hgrid_obj.write_hgrid(f'{wdir}/Dredge/hgrid_dredged_{dredge_depth}m.gr3')
# print("finished loading dredging depth")

# ----------------- Set feeder channel depth -----------------
# A grid without feeder is needed to identify which feeder points are outside and should be deepened
# Only the boundary matters, the interior of the grid doesn't matter,
# so if you don't have a grid without feeders, you can just generate a simplified grid with the lbnd_ocean map
gd_no_feeder = sms2grd('/sciclone/schism10/Hgrid_projects/STOFS3D-v7/v20.0/no_feeder.2dm')
gd_no_feeder.proj(prj0='esri:102008', prj1='epsg:4326')
temp_dp = hgrid_obj.dp.copy()
hgrid_obj = set_feeder_dp(
    feeder_info_dir='/sciclone/schism10/Hgrid_projects/STOFS3D-v7/v20.0/Feeder/',
    hgrid_obj=hgrid_obj, hgrid_obj_no_feeder=gd_no_feeder
)
dp_diff = temp_dp - hgrid_obj.dp
print(f'min deepening: {min(dp_diff)}; max deepening: {max(dp_diff)} \n')

# ------------------- save final product -------------------
hgrid_obj.save(f'{wdir}/hgrid.ll', fmt=0)