#!/usr/bin/env python3

'''
This script is used to edit the DEM for hgrid.
This contains several steps:
'''
from pathlib import Path
import os

import numpy as np

from pylib import schism_grid
from Levee.set_levees import set_levees
from NCF.load_NCF import load_NCF
from xGEOID.convert2xgeoid import convert2xgeoid
from Chart.load_chart import load_chart
from Dredge.dredge_auto_channels import dredge_auto_channels

#wdir = '/sciclone/schism10/feiye/SECOFS/Inputs/I04a/DEM_edit/'
wdir = '/sciclone/schism10/Hgrid_projects/TMP/DEM_edit' #test

# load the hgrid
hgrid_obj = schism_grid(f'{wdir}/DEM_loading/hgrid.ll.dem_loaded.gr3')  # after pload_depth.py
hgrid_obj_dem = hgrid_obj.dp.copy()
print("finished the loading the hgrid process")

# set levees
os.chdir(f'{wdir}/Levee') #to set the directory 
hgrid_obj = set_levees(hgrid_name=hgrid_obj.source_file)
hgrid_obj.write_hgrid(f'{wdir}/Levee/hgrid.ll.dem_loaded_levee.gr3')
hgrid_obj_dem_levee = hgrid_obj.dp.copy()
print("finished the setting levees process")

# load NCF (National Channel Framework)
hgrid_obj = load_NCF(hgrid_obj=hgrid_obj, NCF_shpfile=Path('/sciclone/schism10/Hgrid_projects/Charts/Mississippi/channel_quarter_NCF.shp'))
hgrid_obj.dp = np.maximum(hgrid_obj_dem_levee, hgrid_obj.dp)  # change the depth only if it is deeper than the original depth
hgrid_obj.write_hgrid(f'{wdir}/NCF/hgrid.ll.dem_loaded_levee_NCF.gr3')
print("finished the loading the NCF process")

# convert to xGEOID
hgrid_obj = convert2xgeoid(wdir=f'{wdir}/xGEOID/', hgrid_obj=hgrid_obj)
hgrid_obj_dem_levee_xgeoid = hgrid_obj.dp.copy()
print("finihsed the converting the vdatum to xGEOID process")

# load Chart
hgrid_obj = load_chart(
    hgrid_obj=hgrid_obj,
    sounding_shpfile=Path('/sciclone/schism10/Hgrid_projects/Charts/Savanna_Cooper/savannah_cooper_sounding_3_xyz_edited_xgeoid.shp'),
    region_shpfile=Path('/sciclone/schism10/Hgrid_projects/Charts/Savanna_Cooper/secofs_chart_loading_zones.shp'), crs_region='esri:102008'
)
hgrid_obj.dp = np.maximum(hgrid_obj_dem_levee_xgeoid, hgrid_obj.dp)  # change the depth only if it is deeper than the original depth
hgrid_obj.write_hgrid(f'{wdir}/Chart/hgrid_dem_levee_loaded_NCF_loaded_xGEOID20b_chart_loaded.gr3')
print("finished the loading Chart process")

# load dredging
dredge_depth = 2 # set the dredge depth
hgrid_obj = dredge_auto_channels(hgrid_obj=hgrid_obj, dredge_polygon_file=Path('/sciclone/schism10/Hgrid_projects/SECOFS/new20_JZ/total_river_polys_clipped_test.shp'), dredge_depth=dredge_depth)
hgrid_obj.write_hgrid(f'{wdir}/Dredge/hgrid_dredged_{dredge_depth}m.gr3')
print("finished the loading dredging depth process")
