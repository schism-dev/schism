#!/usr/bin/env python3

'''
This script is used to edit the DEM for hgrid.
This contains several steps, which can be turned on/off by setting the tasks dictionary.
'''

from pathlib import Path
import os

import numpy as np

from pylib import grd2sms, sms2grd  # , schism_grid
from pylib_experimental.schism_file import cread_schism_hgrid
from Levee.set_levees import set_levees
from NCF.load_NCF import load_NCF
from xGEOID.convert2xgeoid import convert2xgeoid
from Chart.load_chart import load_chart
from Dredge.dredge_auto_channels import dredge_auto_channels
from SetFeederDp.set_feeder_dp import set_feeder_dp

# ----------------- inputs -------------------
WDIR = '/sciclone/schism10/feiye/STOFS3D-v8/I04b/Bathy_edit/'
hgrid_obj = cread_schism_hgrid(f'{WDIR}/DEM_loading/hgrid_max_dp.gr3')  # after pload_depth.py

tasks = {
    'NCF': True,
    'Levee': True,
    'xGEOID': True,
    'Chart': False,
    'Dredge': False,
    'Feeder': False,
}
# ----------------- end inputs -------------------

print("Finished reading the DEM-loaded hgrid.\n")
hgrid_base_name = 'hgrid_ll_dem'

if tasks['NCF']:  # load NCF (National Channel Framework)
    hgrid_base_name += '_NCF'
    initial_dp = hgrid_obj.dp.copy()  # save dp before processing
    hgrid_obj = load_NCF(
        hgrid_obj=hgrid_obj, buf=20,
        NCF_shpfile=Path('/sciclone/schism10/Hgrid_projects/Charts/Mississippi/channel_quarter_NCF.shp'))
    hgrid_obj.dp = np.maximum(initial_dp, hgrid_obj.dp)  # change the depth only if it is deeper than the original depth
    print("finished loading NCF.\n")

if tasks['Levee']:  # set levees
    hgrid_base_name += '_levee'
    os.chdir(f'{WDIR}/Levee')  # to set the directory
    hgrid_obj = set_levees(hgrid_obj=hgrid_obj, wdir=f'{WDIR}/Levee/')
    grd2sms(hgrid_obj, f'{WDIR}/Levee/{hgrid_base_name}.2dm')
    print("Finished setting levees.\n")

if tasks['xGEOID']:  # convert to xGEOID
    hgrid_base_name += '_xGEOID'
    hgrid_obj, _ = convert2xgeoid(
        wdir=f'{WDIR}/xGEOID/', hgrid_obj=hgrid_obj, diag_output=f'{WDIR}/{hgrid_base_name}.2dm')
    print("Finihsed converting the vdatum to xGEOID.\n")

if tasks['Chart']:  # load Chart, the Chart has been converted in xGEOID
    # if the chart is in NAVD88, this step should be done before the xGEOID conversion
    initial_dp = hgrid_obj.dp.copy()  # save dp before processing
    hgrid_obj = load_chart(
        hgrid_obj=hgrid_obj,
        sounding_shpfile=Path('/sciclone/schism10/Hgrid_projects/Charts/'
                              'Savanna_Cooper/savannah_cooper_sounding_3_xyz_edited_xgeoid.shp'),
        region_shpfile=Path('/sciclone/schism10/Hgrid_projects/Charts/'
                            'Savanna_Cooper/secofs_chart_loading_zones.shp'),
        crs_region='esri:102008'
    )
    hgrid_obj.dp = np.maximum(initial_dp, hgrid_obj.dp)  # change the depth only if it is deeper than the original depth
    hgrid_base_name += '_chart'
    hgrid_obj.write_hgrid(f'{WDIR}/Chart/{hgrid_base_name}.gr3')
    print("Finished loading Chart.\n")

if tasks['Dredge']:  # dredge the channels made by RiverMapper
    DREDGE_DEPTH = 2  # set the dredge depth
    hgrid_obj = dredge_auto_channels(
        hgrid_obj=hgrid_obj,
        dredge_polygon_file=Path('/sciclone/schism10/Hgrid_projects/SECOFS/'
                                 'new20_JZ/total_river_polys_clipped_test.shp'),
        dredge_depth=DREDGE_DEPTH)
    hgrid_base_name += f'_dredged_{DREDGE_DEPTH}m'
    hgrid_obj.write_hgrid(f'{WDIR}/Dredge/{hgrid_base_name}.gr3')
    print("Finished loading dredging depth.\n")

if tasks['Feeder']:  # set feeder dp
    # A grid without feeder is needed to identify which feeder points are outside and should be deepened
    # Only the boundary matters, the interior of the grid doesn't matter,
    # so if you don't have a grid without feeders, you can just generate a simplified grid with the lbnd_ocean map
    gd_no_feeder = sms2grd('/sciclone/schism10/Hgrid_projects/STOFS3D-v7/v20.0/no_feeder.2dm')
    gd_no_feeder.proj(prj0='esri:102008', prj1='epsg:4326')
    initial_dp = hgrid_obj.dp.copy()
    hgrid_obj = set_feeder_dp(
        feeder_info_dir='/sciclone/schism10/Hgrid_projects/STOFS3D-v7/v20.0/Feeder/',
        hgrid_obj=hgrid_obj, hgrid_obj_no_feeder=gd_no_feeder
    )
    dp_diff = initial_dp - hgrid_obj.dp
    print(f'min deepening: {min(dp_diff)}; max deepening: {max(dp_diff)} \n')
    hgrid_base_name += '_feeder_deepened'
    hgrid_obj.save(f'{WDIR}/Feeder/{hgrid_base_name}.gr3')
    print("Finished setting feeder dp.\n")

# ------------------- save final product -------------------
hgrid_obj.save(f'{WDIR}/hgrid.ll', fmt=1)
print(f"Finished saving the final product {WDIR}/hgrid.ll\n")
