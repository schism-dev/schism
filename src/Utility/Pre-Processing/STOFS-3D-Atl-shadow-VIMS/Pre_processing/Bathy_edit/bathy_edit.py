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

from pylib import grd2sms, sms2grd
try:  # c++ function to speed up the grid reading
    from pylib_experimental.schism_file import xread_schism_hgrid as schism_read
except ImportError:
    from pylib import schism_grid as schism_read


IMPLEMENTED_TASKS = [  # order matters
    'Regional_tweaks',  # set minimum depth in regions specified in regional_tweaks
    'NCF',  # load NCF (National Channel Framework) maintained depth
    'Levee',  # set levees height based on National Levee Database
    'xGEOID',  # convert from NAVD88 to xGEOID
    'Chart',  # load chart depth, the chart has been converted to xGEOID
    'Dredge',  # dredge the channels made by RiverMapper, relative, datum doesn't matter
    'Feeder',  # set feeder channel depth, relative, datum doesn't matter
]

DEFAULT_TASKS = {'Regional_tweaks', 'NCF', 'Levee'}

# larger files not included in the Git repository, need to be copied to the working directory
LARGE_FILES = {
    'NCF': [
        Path('/sciclone/schism10/Hgrid_projects/NCF/'),
        Path('/work/noaa/nosofs/feiye/STOFS-3D-Atl-Example-Setup/DATA/NCF/')
    ],
    'Levee': [
        Path('/sciclone/schism10/Hgrid_projects/Levees/Levee_v3/FEMA_regions/'
             'FEMA_region_levees/'),
        Path('/work/noaa/nosofs/feiye/STOFS-3D-Atl-Example-Setup/DATA/Levee/')
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
        print(f'Copying the script and the subdirectories to {wdir}')
        for task in tasks:
            shutil.copytree(
                f'{script_dir}/{task}/', f'{wdir}/{task}/',
                symlinks=False, dirs_exist_ok=True)

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

    if 'Regional_tweaks' in tasks:  # set minimum depth in regions
        from Regional_tweaks.regional_tweaks import tweak_hgrid_depth
        hgrid_obj = tweak_hgrid_depth(
            hgrid=hgrid_obj, regions_dir=f'{wdir}/Regional_tweaks/regions/')
        initial_dp = hgrid_obj.dp.copy()  # treat the regional tweaks as the initial dp
        print("Finished setting regional tweaks and updating initial dp.\n")

    if 'NCF' in tasks:  # load NCF (National Channel Framework)
        from NCF.load_NCF import load_NCF
        hgrid_base_name += '_NCF'
        hgrid_obj = load_NCF(
            hgrid_obj=hgrid_obj, buf=20, NCF_shpfile=f'{wdir}/NCF/channel_quarter_NCF.shp')
        hgrid_obj.dp = np.maximum(initial_dp, hgrid_obj.dp)
        print("finished loading NCF.\n")

    if 'Levee' in tasks:  # set levees
        from Levee.set_levees import set_levees
        hgrid_base_name += '_levee'
        os.chdir(f'{wdir}/Levee')  # to set the directory
        hgrid_obj = set_levees(hgrid_obj=hgrid_obj, wdir=f'{wdir}/Levee/')
        grd2sms(hgrid_obj, f'{wdir}/Levee/{hgrid_base_name}.2dm')
        print("Finished setting levees.\n")

    if 'xGEOID' in tasks:  # convert from NAVD88 to xGEOID
        from xGEOID.convert2xgeoid import convert2xgeoid
        hgrid_base_name += '_xGEOID'
        hgrid_obj, _ = convert2xgeoid(
            wdir=f'{wdir}/xGEOID/', hgrid_obj=hgrid_obj,
            diag_output=f'{wdir}/{hgrid_base_name}.2dm')
        print("Finihsed converting the vdatum to xGEOID.\n")

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
        DREDGE_DEPTH = 2  # set the dredge depth
        hgrid_obj = dredge_auto_channels(
            hgrid_obj=hgrid_obj,
            dredge_polygon_file=Path(
                '/sciclone/schism10/Hgrid_projects/SECOFS/'
                'new20_JZ/total_river_polys_clipped_test.shp'),
            dredge_depth=DREDGE_DEPTH)
        hgrid_base_name += f'_dredged_{DREDGE_DEPTH}m'
        hgrid_obj.write_hgrid(f'{wdir}/Dredge/{hgrid_base_name}.gr3')
        print("Finished loading dredging depth.\n")

    if 'Feeder' in tasks:  # set feeder channel depth
        from SetFeederDp.set_feeder_dp import set_feeder_dp
        # A grid without feeder is needed to identify which feeder points are outside and should be deepened
        # Only the boundary matters, the interior of the grid doesn't matter,
        # so if you don't have a grid without feeders, you can just generate a simplified grid with the lbnd_ocean map
        gd_no_feeder = sms2grd('/sciclone/schism10/Hgrid_projects/STOFS3D-v8/v23.1/hgrid.2dm')
        gd_no_feeder.proj(prj0='esri:102008', prj1='epsg:4326')
        initial_dp = hgrid_obj.dp.copy()
        hgrid_obj = set_feeder_dp(
            feeder_info_dir='/sciclone/schism10/Hgrid_projects/STOFS3D-v7/v23.3/Feeder/',
            hgrid_obj=hgrid_obj, hgrid_obj_no_feeder=gd_no_feeder
        )
        dp_diff = initial_dp - hgrid_obj.dp
        print(f'min deepening: {min(dp_diff)}; max deepening: {max(dp_diff)} \n')
        hgrid_base_name += '_feeder_deepened'
        hgrid_obj.save(f'{wdir}/Feeder/{hgrid_base_name}.gr3')
        print("Finished setting feeder dp.\n")

    # ------------------- save final product -------------------
    hgrid_obj.save(f'{wdir}/hgrid.ll', fmt=1)
    print(f"Finished saving the final product {wdir}/hgrid.ll\n")


def sample_usage():
    '''
    Sample usage of the bathy_edit function.
    '''
    WDIR = Path('/sciclone/schism10/feiye/STOFS3D-v8/I10/Bathy_edit2/')
    HGRID_FNAME = Path(  # Typically, this is the DEM-loaded hgrid
        '/sciclone/schism10/feiye/STOFS3D-v8/I10/Bathy_edit/'
        'DEM_loading/hgrid.ll.dem_loaded.mpi.gr3'
    )
    TASKS = ['xGEOID']

    bathy_edit(wdir=WDIR, hgrid_fname=HGRID_FNAME, tasks=TASKS)


if __name__ == '__main__':
    sample_usage()
    print('Done')
