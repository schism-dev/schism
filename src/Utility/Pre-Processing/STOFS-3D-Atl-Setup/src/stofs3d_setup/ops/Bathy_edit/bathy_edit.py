#!/usr/bin/env python3

'''
This script is used to edit the depth of a DEM-loaded hgrid.

Multiple steps may be needed, which can be turned on/off by setting the "TASKS",
see "def sample_usage".

Each step imports scripts under the subdirectories.

You can run this script inside the SCHISM git folder after
setting paths in the sample_usage function.
The script will copy itself and the subdirectories being used
to the working directory to keep a record.

Alternatively, you can copy the entire script folder (Bathy_edit/)
to your working directory, then set paths in the sample_usage function.

Some larger files are not included in the Git repository,
you need to specify the paths in the "LARGE_FILES" dictionary below.

Note that xGEOID is just a wrapper around vdatum.jar,
which takes a lot of memory and quite slow.
Running xGEOID on a large grid may take hours.
Use viz.sciclone.wm.edu; other SciClone subclusters often run out of memory.
This is deprecated, use Felicio's workflow after testing.
'''

import os
from pathlib import Path
import shutil
from typing import Optional, List

import numpy as np

from pylib import grd2sms
from pylib import schism_grid as schism_read


IMPLEMENTED_TASKS = [  # order matters
    'Regional_tweaks',  # set minimum depth in regions specified in regional_tweaks
    'NCF',  # load NCF (National Channel Framework) maintained depth
    'Levee',  # set levees height based on National Levee Database
    'Levee_dev',  # a development version of levee setting, not fully tested,
                  # using 2025 MTG levee data and not forcing min height
    'xGEOID',  # convert from NAVD88 to xGEOID, use viz (gulf has memory issues); deprecated, use Felicio's workflow
    'xGEOID_cmvd',  # convert from NAVD88 to xGEOID using Felicio's alternative tool cmvd, much faster than vdatum
    'xGEOID_from_diff',  # convert from NAVD88 to xGEOID based on dp difference, must be the same (x, y)
    'Chart',  # load chart depth, the chart has been converted to xGEOID
    'Dredge',  # dredge the channels made by RiverMapper, relative, vertical datum doesn't matter
    'Ensure_channel_connectivity',  # ensure channel connectivity by ensure a minimum elevation drop
                                    # from bank to thalweg for river transects defined by RiverMapper,
                                    # relative, vertical datum doesn't matter
    'Feeder',  # set feeder channel depth, relative, vertical datum doesn't matter
    'Temporary_Fix_v7p2',  # tweak depths around Bayou Lafourche
    'Temporary_Fix_v7.2.1',  # tweak depth around Philadelphia International Airport and Minas Basin
    'Temporary_Fix_v7.4',  # tweak depth around Philadelphia International Airport
    'Temporary_Fix_secofs',  # tweak depths around Bayou Lafourche
]

TASKS = {
    'Temporary_Fix_secofs',
    'Regional_tweaks',
    'NCF',
    'Levee_dev',
    'xGEOID_cmvd',
    'Ensure_channel_connectivity',
    'Temporary_Fix_v7p2',
}

# Test refactor note:
# Keep production behavior identical to bathy_edit.py, but gather hardwired
# paths and frequently tuned parameters here so future cleanup can be reviewed
# without changing task logic.
WORKFLOW_CONSTANTS = {
    'Dredge': {
        'dredge_depth': 2,
        'dredge_polygon_file': Path(
            '/sciclone/schism10/Hgrid_projects/SECOFS/'
            'new20_JZ/total_river_polys_clipped_test.shp'
        ),
    },
    'Ensure_channel_connectivity': {
        'min_channel_depth': 1.0,
        'measured_from_high_bank': True,
        'river_extra_info_map_file': (
            '/sciclone/schism10/Hgrid_projects/STOFS3D-v7/v19_RiverMapper/Outputs/'
            'bora_v19.1.v19_ie_v18_3_nwm_clipped_in_cudem_missing_tiles_20-core/'
            'total_river_arcs_extra.map'
        ),
        'region_gdf_file': (
            '/sciclone/schism10/Hgrid_projects/STOFS3D-v7.4/v32e/Clip/outputs/'
            'watershed.shp'
        ),
        'exclude_region_gdf_file_list': [
            (
                '/sciclone/schism10/feiye/STOFS3D-v8/I15_v7/Bathy_edit/'
                'RiverArc_Dredge/watershed_ME.shp'
            )
        ],
    },
    'Feeder': {
        'hgrid_without_feeders': Path('/sciclone/schism10/feiye/STOFS3D-v8/R13p_v7/hgrid.gr3'),
        'feeder_info_dir': Path('/sciclone/schism10/Hgrid_projects/STOFS3D-v8/v31/Feeder/'),
    },
    'Temporary_Fix_v7p2': {
        'reference_hgrid_file': (
            '/sciclone/schism10/Hgrid_projects/STOFS3D-v7.4/v32f_test/Bathy_edit/'
            'DEM_loading_v7.2_original/hgrid.ll.dem_loaded.mpi.gr3'
        ),
    },
    'Temporary_Fix_secofs': {
        'reference_hgrid_file': (
            '/sciclone/schism10/hjyoo/task/task10_Atlantic/RUN100g/src/Bathy_edit_org/'
            'DEM_loading/DEM_loading_for_SECOFS_v1/hgrid.ll.dem_loaded.mpi.gr3'
        ),
    },
    'sample_usage': {
        'wdir': Path('/sciclone/schism10/hjyoo/task/task10_Atlantic/RUN100g/src/Bathy_edit_org/'),
        'hgrid_fname': Path(
            '/sciclone/schism10/hjyoo/task/task10_Atlantic/RUN100g/src/Bathy_edit_org/'
            'DEM_loading/DEM_loading_for_STOFS_v7p4/hgrid.ll.dem_loaded.mpi.upper_Hudson.gr3'
        ),
    },
}

# larger files not included in the Git repository, need to be copied to the working directory
LARGE_FILES = {
    'NCF': [
        Path('/sciclone/schism10/Hgrid_projects/NCF_patched_James/'),
        Path('/work/noaa/nosofs/feiye/STOFS-3D-Atl-Example-Setup/DATA/NCF_patched_James/')
    ],
    'Levee': [
        Path('/sciclone/schism10/Hgrid_projects/Levees/Levee_v3/FEMA_regions/FEMA_region_levees/'),
        Path('/work/noaa/nosofs/feiye/STOFS-3D-Atl-Example-Setup/DATA/Levee/')
    ],
    'Levee_dev': [
        Path('/sciclone/schism10/Hgrid_projects/Levees/Levee_2025/FEMA_regions/FEMA_region_levees/'),
        Path('/work/noaa/nosofs/feiye/STOFS-3D-Atl-Example-Setup/DATA/Levee_2025/')  # not uploaded yet
    ],
    'Chart': [
        Path('/sciclone/schism10/Hgrid_projects/Charts/Savanna_Cooper/SECOFS/'),
        Path('/work/noaa/nosofs/feiye/STOFS-3D-Atl-Example-Setup/DATA/Charts/Savanna_Cooper/SECOFS/')
    ]
}


def handle_tasks(user_tasks: Optional[List[str]]) -> List[str]:
    '''
    Handle task validation and ensure tasks follow the allowed order.

    Args:
    - user_tasks: A list of tasks provided by the user.
    Returns:
    - A list of valid tasks in the order of IMPLEMENTED_TASKS.
    Raises:
    - ValueError: If an undefined task is provided.
    '''
    if user_tasks is None:
        user_tasks = []

    # Convert user_tasks to a set for easy comparison
    user_tasks_set = set(user_tasks)

    # Check for undefined tasks
    invalid_tasks = user_tasks_set - set(IMPLEMENTED_TASKS)
    if invalid_tasks:
        raise ValueError(f"Undefined tasks: {', '.join(invalid_tasks)}")

    # Return tasks in the same order as IMPLEMENTED_TASKS
    ordered_tasks = [task for task in IMPLEMENTED_TASKS if task in user_tasks]

    return ordered_tasks


def copy_large_files(src_paths: List[Path], dest_dir: Path):
    '''
    Copy larger files not in the Git repository to the working directory.

    Args:
    - wdir: The working directory where the files will be copied.
    - src_paths: A list of paths to the source files.
    - dest_dir: The destination directory where the files will be copied.
    '''
    copy_success = False
    for src_path in src_paths:
        if src_path.is_dir():
            shutil.copytree(src_path, dest_dir, dirs_exist_ok=True)
            print(f'Copied {src_path} to {dest_dir}')
            copy_success = True
            break

    if not copy_success:
        raise FileNotFoundError(f'None of the source paths exist: {src_paths}')


def prepare_dir(wdir: Path, tasks: str):
    '''
    Make a copy of the scripts in the working directory.
    This includes the scripts for each requested task and
    the larger files not in the Git repository.
    '''
    script_dir = Path(__file__).parent
    if script_dir == wdir:
        print('The script is already in the working directory; no need to copy.')
    else:
        existing_task_dirs = [wdir / task for task in tasks if (wdir / task).exists()]
        if existing_task_dirs:
            existing_dirs = '\n'.join(f'  {path}' for path in existing_task_dirs)
            raise FileExistsError(
                'Refusing to prepare the working directory because these requested '
                f'task directories already exist:\n{existing_dirs}\n'
                'Move or remove them, or choose a new working directory. No files were copied.'
            )

        print(f'Copying the script and the subdirectories to {wdir}')
        os.makedirs(wdir, exist_ok=True)
        os.system(f'cp {__file__} {wdir}')  # copy the script itself to the working directory
        for task in tasks:
            if Path(f'{script_dir}/{task}/').exists():
                shutil.copytree(
                    f'{script_dir}/{task}/', f'{wdir}/{task}/',
                    symlinks=False, dirs_exist_ok=True)
                print(f'Copied {task} to {wdir}')
            else:
                raise FileNotFoundError(f'Task directory not found: {script_dir}/{task}/')

    # copy larger files not in the Git repository
    for task, file_path_list in LARGE_FILES.items():
        if task in tasks:
            copy_large_files(file_path_list, f'{wdir}/{task}/')


def bathy_edit(wdir: Path, hgrid_fname: Path, tasks: list = None):
    '''
    Edit the DEM for hgrid.

    - wdir: The working directory where the final hgrid.ll file will be saved.
    - hgrid_fname: The path to the DEM-loaded hgrid, must be in lon/lat
    - tasks: A list of tasks to be performed, the order matters.
    '''
    tasks = handle_tasks(tasks)  # check and order tasks
    prepare_dir(wdir, tasks)  # save a copy of the scripts to the working directory

    hgrid_obj = schism_read(str(hgrid_fname))  # after pload_depth.py
    print("Finished reading the DEM-loaded hgrid.\n")

    hgrid_base_name = 'hgrid_ll_dem'
    initial_dp = hgrid_obj.dp.copy()  # save dp before processing

    if 'Temporary_Fix_secofs' in tasks:  # replace depth in SECOFS domain (SJR, Savannah, Charleston)
        from Temporary_Fix_secofs.temp_fix_secofs import temp_fix_secofs
        task_cfg = WORKFLOW_CONSTANTS['Temporary_Fix_secofs']
        hgrid_base_name += '_temp_fix_secofs'
        hgrid_obj = temp_fix_secofs(
            hgrid_obj,
            wdir=f'{wdir}/Temporary_Fix_secofs/',
            reference_hgrid_file=task_cfg['reference_hgrid_file']
        )
        hgrid_obj.grd2sms(
            f'{wdir}/Temporary_Fix_secofs/{hgrid_base_name}.2dm'
        )
        print("Finished setting temporary fix for SECOFS.\n")

    if 'Regional_tweaks' in tasks:  # set minimum depth in regions
        from Regional_tweaks.regional_tweaks import tweak_hgrid_depth, tweak_shp_hgrid_depth
        hgrid_obj = tweak_hgrid_depth(
            hgrid=hgrid_obj, regions_dir=f'{wdir}/Regional_tweaks/regions/')
        hgrid_obj = tweak_shp_hgrid_depth(
            hgrid=hgrid_obj, regions_dir=f'{wdir}/Regional_tweaks/shp_regions/')        
        grd2sms(hgrid_obj, f'{wdir}/Regional_tweaks/{hgrid_base_name}_tweaks.2dm')
        initial_dp = hgrid_obj.dp.copy()  # treat the regional tweaks as the initial dp
        print("Finished setting regional tweaks and updating initial dp.\n")

    if 'NCF' in tasks:  # load NCF (National Channel Framework)
        from NCF.load_NCF import load_NCF
        hgrid_base_name += '_NCF'
        hgrid_obj = load_NCF(
            hgrid_obj=hgrid_obj, buf=20, NCF_shpfile=f'{wdir}/NCF/channel_quarter_NCF.shp')
        hgrid_obj.dp = np.maximum(initial_dp, hgrid_obj.dp)
        grd2sms(hgrid_obj, f'{wdir}/NCF/{hgrid_base_name}.2dm')
        print("finished loading NCF.\n")

    if 'Levee' in tasks:  # set levees
        err_msg = (
            "The 'Levee' task is deprecated. Use 'Levee_dev' instead. "
            "Use 'min_levee_height_meters' to set the minimum levee height in meters (positive upward). "
            "2 m is a reasonable value, but it needs more testing in the shadow forecast. "
            "The operation uses 7 m to be safe."
        )
        raise NotImplementedError(err_msg)
        # from Levee.set_levees_operation import set_levees
        # hgrid_base_name += '_levee'
        # os.chdir(f'{wdir}/Levee')  # to set the directory
        # hgrid_obj = set_levees(min_levee_height_meters=7, hgrid_obj=hgrid_obj, wdir=f'{wdir}/Levee/')
        # grd2sms(hgrid_obj, f'{wdir}/Levee/{hgrid_base_name}.2dm')
        # print("Finished setting levees.\n")

    if 'Levee_dev' in tasks:  # set levees
        from Levee.set_levees_dev import set_levees
        hgrid_base_name += '_levee'
        os.chdir(f'{wdir}/Levee_dev')  # to set the directory
        hgrid_obj = set_levees(hgrid_obj=hgrid_obj, wdir=f'{wdir}/Levee_dev/', min_levee_height_meters=7)
        grd2sms(hgrid_obj, f'{wdir}/Levee_dev/{hgrid_base_name}.2dm')
        print("Finished setting levees (dev).\n")

    if 'xGEOID' in tasks:  # convert from NAVD88 to xGEOID
        from xGEOID.convert2xgeoid import convert2xgeoid
        hgrid_base_name += '_xGEOID'
        hgrid_obj, _ = convert2xgeoid(
            wdir=f'{wdir}/xGEOID/', hgrid_obj=hgrid_obj,
            diag_output=f'{wdir}/{hgrid_base_name}.2dm')
        print("Finihsed converting the vdatum to xGEOID.\n")

    if 'xGEOID_cmvd' in tasks:  # convert from NAVD88 to xGEOID using Felicio's alternative tool cmvd
        from xGEOID_cmvd.xGEOID_cmvd import convert2xgeoid_cmvd
        hgrid_base_name += '_xGEOID_cmvd'
        hgrid_obj = convert2xgeoid_cmvd(hgrid_obj)
        hgrid_obj.grd2sms(f'{wdir}/xGEOID_cmvd/{hgrid_base_name}.2dm')
        print("Finished converting the vdatum to xGEOID using cmvd.\n")

    if 'xGEOID_from_diff' in tasks:  # convert from NAVD88 to xGEOID based on dp difference
        hgrid_diff = schism_read(f'{wdir}/xGEOID_from_diff/hgrid_xGEOID20b_dp_minus_NAVD_dp.gr3')
        hgrid_obj.dp += hgrid_diff.dp
        hgrid_base_name += '_xGEOID_from_diff'
        hgrid_obj.write_hgrid(f'{wdir}/xGEOID_from_diff/{hgrid_base_name}.gr3')
        print("Finished converting the vdatum to xGEOID based on dp difference.\n")

    if 'Chart' in tasks:  # load Chart, the Chart has been converted in xGEOID
        from Chart.load_chart import load_chart
        # if the chart is in NAVD88, this step should be done before the xGEOID conversion
        initial_dp = hgrid_obj.dp.copy()  # save dp before processing
        hgrid_obj = load_chart(
            hgrid_obj=hgrid_obj,
            sounding_shpfile=Path(f'{wdir}/Chart/savannah_cooper_sounding_3_xyz_edited_xgeoid.shp'),
            region_shpfile=Path(f'{wdir}/Chart/secofs_chart_loading_zones.shp'),
            crs_region='esri:102008'
        )
        # change the depth only if it is deeper than the original depth
        hgrid_obj.dp = np.maximum(initial_dp, hgrid_obj.dp)
        hgrid_base_name += '_chart'
        hgrid_obj.write_hgrid(f'{wdir}/Chart/{hgrid_base_name}.gr3')
        print("Finished loading Chart.\n")

    if 'Dredge' in tasks:  # dredge the channels made by RiverMapper
        from Dredge.dredge_auto_channels import dredge_auto_channels
        task_cfg = WORKFLOW_CONSTANTS['Dredge']
        dredge_depth = task_cfg['dredge_depth']
        hgrid_obj = dredge_auto_channels(
            hgrid_obj=hgrid_obj,
            dredge_polygon_file=task_cfg['dredge_polygon_file'],
            dredge_depth=dredge_depth)
        hgrid_base_name += f'_dredged_{dredge_depth}m'
        hgrid_obj.write_hgrid(f'{wdir}/Dredge/{hgrid_base_name}.gr3')
        print("Finished loading dredging depth.\n")

    if 'Ensure_channel_connectivity' in tasks:  # dredging river transects defined by RiverMapper
        from Ensure_channel_connectivity.ensure_channel_connectivity import ensure_channel_connectivity
        task_cfg = WORKFLOW_CONSTANTS['Ensure_channel_connectivity']
        hgrid_obj = ensure_channel_connectivity(
            hgrid_obj,
            min_channel_depth=task_cfg['min_channel_depth'],
            measured_from_high_bank=task_cfg['measured_from_high_bank'],
            river_extra_info_map_file=task_cfg['river_extra_info_map_file'],
            region_gdf_file=task_cfg['region_gdf_file'],
            exclude_region_gdf_file_list=task_cfg['exclude_region_gdf_file_list'],
            output_dir=f'{wdir}/Ensure_channel_connectivity/'
        )
        hgrid_base_name += '_channel_connectivity_ensured'
        hgrid_obj.write_hgrid(f'{wdir}/Ensure_channel_connectivity/{hgrid_base_name}.gr3')
        print("Finished ensuring channel connectivity.\n")

    if 'Feeder' in tasks:  # set feeder channel depth
        from SetFeederDp.set_feeder_dp import set_feeder_dp
        task_cfg = WORKFLOW_CONSTANTS['Feeder']
        # A grid without feeder is needed to identify which feeder points are outside and should be deepened
        # Only the boundary matters, the interior of the grid doesn't matter,
        # so if you don't have a grid without feeders, you can just generate a simplified grid with the lbnd_ocean map
        gd_no_feeder = schism_read(str(task_cfg['hgrid_without_feeders']))
        # gd_no_feeder.proj(prj0='esri:102008', prj1='epsg:4326')
        initial_dp = hgrid_obj.dp.copy()
        hgrid_obj = set_feeder_dp(
            feeder_info_dir=str(task_cfg['feeder_info_dir']),
            hgrid_obj=hgrid_obj, hgrid_obj_no_feeder=gd_no_feeder
        )
        dp_diff = initial_dp - hgrid_obj.dp
        print(f'min deepening: {min(dp_diff)}; max deepening: {max(dp_diff)} \n')
        hgrid_base_name += '_feeder_deepened'
        hgrid_obj.save(f'{wdir}/Feeder/{hgrid_base_name}.gr3')
        print("Finished setting feeder dp.\n")

    if 'Temporary_Fix_v7p2' in tasks:  # tweak depths around Bayou Lafourche
        from Temporary_Fix_v7p2.temp_fix_v7p2 import temp_fix_v7p2
        task_cfg = WORKFLOW_CONSTANTS['Temporary_Fix_v7p2']
        hgrid_base_name += '_temp_fix_v7.2'
        hgrid_obj = temp_fix_v7p2(
            hgrid_obj,
            wdir=f'{wdir}/Temporary_Fix_v7p2/',
            reference_hgrid_file=task_cfg['reference_hgrid_file']
        )
        hgrid_obj.grd2sms(f'{wdir}/Temporary_Fix_v7p2/{hgrid_base_name}.2dm')
        print("Finished setting temporary fix for v7p2.\n")

    if 'Temporary_Fix_v7.2.1' in tasks:  # tweaks around Philadelphia International Airport and Bay of Fundy
        from Regional_tweaks.regional_tweaks import shape_tweak
        hgrid_base_name += '_temp_fix_v7.2.1'
        gpkg_file = f'{wdir}/Temporary_Fix_v7.2.1/v7.2.1_fix.gpkg'
        hgrid_obj, hgrid_tweaked_idx = shape_tweak(hgrid_obj, gpkg_file)
        grd2sms(hgrid_obj, f'{wdir}/Temporary_Fix_v7.2.1/{hgrid_base_name}.2dm')
        print("Finished setting shape tweaks.\n")

    if 'Temporary_Fix_v7.4' in tasks:
        from Regional_tweaks.regional_tweaks import shape_tweak
        hgrid_base_name += '_temp_fix_v7.4'
        gpkg_file = f'{wdir}/Temporary_Fix_v7.4/v7.4_fix.gpkg'
        hgrid_obj, hgrid_tweaked_idx = shape_tweak(hgrid_obj, gpkg_file)
        grd2sms(hgrid_obj, f'{wdir}/Temporary_Fix_v7.4/{hgrid_base_name}.2dm')
        print("Finished setting shape tweaks for v7.4.\n")

    # ------------------- save final product -------------------
    hgrid_obj.save(f'{wdir}/hgrid_dem_edit.ll', fmt=1)
    print(f"Finished saving the final product {wdir}/hgrid_dem_edit.ll\n")


def sample_usage():
    '''
    Sample usage of the bathy_edit function.

    Recommended to copy he entire Bathy_edit/ folder to your working directory to keep
    a record of the scripts and the larger files used.
    Then set the paths below and run this function.
    '''
    WDIR = WORKFLOW_CONSTANTS['sample_usage']['wdir'].resolve(strict=True)
    HGRID_FNAME = WORKFLOW_CONSTANTS['sample_usage']['hgrid_fname']
    bathy_edit(wdir=WDIR, hgrid_fname=HGRID_FNAME, tasks=TASKS)


if __name__ == '__main__':
    sample_usage()
    print('Done')
