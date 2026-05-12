#!/usr/bin/env python3
'''
Sample usage of the driver

Recommended usage:
 - Make a copy of this file under the same directory as stofs3d_atl_driver.py.
 - Based on the sample_setup() function,
   add your customized setup function (e.g., my_first_run, my_second_run, etc.)


Please read the comments in the sample_setup() function
'''

from datetime import datetime
from stofs3d_atl_driver import stofs3d_atl_driver
from stofs3d_atl_config import ConfigStofs3dAtlantic


def La_v8_2024_reforecast():
    '''
    STOFS3D Atl v8 setup for a local test on Louisianna

    Pre-requisites:
        hgrid.gr3 (in lon/lat, with boundaries) must be prepared before running this script
    '''

    # define the project_dir, the current run folders will be created under it
    project_dir = '/sciclone/schism10/feiye/STOFS3D-v8/'

    # run ID. e.g, 00, 01a, 10b, etc.; the run folders are named using this ID as follows:
    # I{runid}: input directory for holding model input files and scripts.
    #           This script and all packages called by this script will be copied to
    #           I{runid}/Pre_processing_scripts_backup/ for future reference.
    # R{runid}: run directory, where the run will be submitted to queue,
    # O{runid}: output directory for holding raw outputs and post-processing,
    # These three folders will be put under project_dir
    runid = '19'

    hgrid_path = '/sciclone/schism10/feiye/STOFS3D-v8/I19/hgrid.gr3'

    # get a configuration preset and make changes if necessary
    config = ConfigStofs3dAtlantic.v8_louisianna()
    config.rnday = 35
    config.startdate = datetime(2024, 3, 5)
    config.nwm_cache_folder = (
        '/sciclone/schism10/feiye/STOFS3D-v8/I09/Source_sink/original_source_sink/20240305/'
    )  # set to None if you don't have cache and a cache will be generated.
    # You can find it under I{runid}/Source_sink/original_source_sink/{starting_date}/
    # which can be reused for future runs with a same simulation period

    # alternatively, you can set the parameters from scratch on a
    # new instance of ConfigStofs3dAtlantic:
    # config = ConfigStofs3dAtlantic()
    # You can use any of the preset as a template.

    # Define the scr_dir if you want to save schism outputs to a scratch directory as
    # scr_dir/R{runid}/outputs/
    # a symlink will be created in the run directory to link to outputs/ on scratch as
    # project_dir/R{runid}/outputs/ -> scr_dir/R{runid}/outputs/
    # Set to None if you don't want to use a scratch directory.
    scr_dir = '/sciclone/scr10/feiye/'

    # swithes to generate different input files
    input_files = {
        'bctides': False,
        'vgrid': False,  # already generated
        'gr3': False,
        'nudge_gr3': False,
        'shapiro': False,
        'drag': False,
        'elev_ic': False,
        'source_sink': False,  # already generated
        'hotstart.nc': False,
        '3D.th.nc': False,
        'elev2D.th.nc': False,
        '*nu.nc': False,
        '*.prop': True,
    }

    # call the driver
    stofs3d_atl_driver(
        hgrid_path=hgrid_path, config=config,
        project_dir=project_dir, runid=runid,
        scr_dir=scr_dir, input_files=input_files
    )


def La_Ida_v8():
    '''
    STOFS3D Atl v8 setup for a local test on Louisianna

    Pre-requisites:
        hgrid.gr3 (in lon/lat, with boundaries) must be prepared before running this script
    '''

    # define the project_dir, the current run folders will be created under it
    project_dir = '/sciclone/schism10/feiye/STOFS3D-v8/'

    # run ID. e.g, 00, 01a, 10b, etc.; the run folders are named using this ID as follows:
    # I{runid}: input directory for holding model input files and scripts.
    #           This script and all packages called by this script will be copied to
    #           I{runid}/Pre_processing_scripts_backup/ for future reference.
    # R{runid}: run directory, where the run will be submitted to queue,
    # O{runid}: output directory for holding raw outputs and post-processing,
    # These three folders will be put under project_dir
    runid = '09z'

    hgrid_path = '/sciclone/schism10/feiye/STOFS3D-v8/I09d/09c.gr3'

    # get a configuration preset and make changes if necessary
    config = ConfigStofs3dAtlantic.v8_louisianna()
    config.rnday = 55
    config.startdate = datetime(2021, 8, 1)
    config.nwm_cache_folder = None
    # config.nwm_cache_folder = (
    #     '/sciclone/schism10/feiye/STOFS3D-v8/I09g/'
    #     'Source_sink/original_source_sink/20210801'
    # )  # set to None if you don't have cache and a cache will be generated.
    # You can find it under I{runid}/Source_sink/original_source_sink/{starting_date}/
    # which can be reused for future runs with a same simulation period

    # alternatively, you can set the parameters from scratch on a
    # new instance of ConfigStofs3dAtlantic:
    # config = ConfigStofs3dAtlantic()
    # You can use any of the preset as a template.

    # Define the scr_dir if you want to save schism outputs to a scratch directory as
    # scr_dir/R{runid}/outputs/
    # a symlink will be created in the run directory to link to outputs/ on scratch as
    # project_dir/R{runid}/outputs/ -> scr_dir/R{runid}/outputs/
    # Set to None if you don't want to use a scratch directory.
    scr_dir = '/sciclone/scr10/feiye/'

    # swithes to generate different input files
    input_files = {
        'bctides': False,
        'vgrid': False,  # already generated
        'gr3': False,
        'nudge_gr3': False,
        'shapiro': False,
        'drag': False,
        'elev_ic': False,
        'source_sink': True,  # already generated
        'hotstart.nc': False,
        '3D.th.nc': False,
        'elev2D.th.nc': False,
        '*nu.nc': False,
        '*.prop': False,
    }

    # call the driver
    stofs3d_atl_driver(
        hgrid_path=hgrid_path, config=config,
        project_dir=project_dir, runid=runid,
        scr_dir=scr_dir, input_files=input_files
    )


def v7p2_Isabel():
    '''
    STOFS3D Atl v7.2 setup for Isabel case
    '''

    project_dir = '/sciclone/schism10/feiye/STOFS3D-v7.2/'
    run_id = '16x'
    hgrid_path = '/sciclone/schism10/feiye/STOFS3D-v8/R15c_v7/hgrid.gr3'
    scr_dir = '/sciclone/scr10/feiye/'

    config = ConfigStofs3dAtlantic.v7p2()
    config.rnday = 23
    config.startdate = datetime(2003, 9, 8)

    input_files = {
        'bctides': True,
        'vgrid': True,
        'gr3': True,
        'nudge_gr3': True,
        'shapiro': True,
        'drag': True,
        'elev_ic': True,
        'source_sink': True,
        'hotstart.nc': False,
        '3D.th.nc': True,
        'elev2D.th.nc': True,
        '*nu.nc': True,
        '*.prop': True,
    }

    stofs3d_atl_driver(
        hgrid_path=hgrid_path, config=config,
        project_dir=project_dir, runid=run_id,
        scr_dir=scr_dir, input_files=input_files
    )


def v7p2_2018_hindcast():
    '''
    STOFS3D Atl v7.2
    '''

    project_dir = '/sciclone/schism10/feiye/STOFS3D-v8/'
    run_id = '15n_v7'
    hgrid_path = '/sciclone/schism10/feiye/STOFS3D-v8/I15a_v7/FeederDp/hgrid.feeder_dp.ll'
    scr_dir = None

    config = ConfigStofs3dAtlantic.v7p2()
    config.rnday = 396
    config.startdate = datetime(2017, 12, 1)
    config.existing_source_json_path = (
        '/sciclone/schism10/feiye/STOFS3D-v8/I15k_v7/Source_sink/original_source_sink/'
    )
    config.nwm_cache_folder = (
        # '/sciclone/schism10/feiye/STOFS3D-v8/I15/Source_sink/original_source_sink/20171201/'
        '/sciclone/schism10/feiye/STOFS3D-v8/NWM/CONUS/netcdf/CHRTOUT/for_2018_hindcast/'
    )

    input_files = {
        'bctides': False,
        'vgrid': False,
        'gr3': False,
        'nudge_gr3': False,
        'shapiro': False,
        'drag': False,
        'elev_ic': False,
        'source_sink': True,
        'hotstart.nc': False,
        '3D.th.nc': False,
        'elev2D.th.nc': False,
        '*nu.nc': False,
        '*.prop': False,
    }

    stofs3d_atl_driver(
        hgrid_path=hgrid_path, config=config,
        project_dir=project_dir, runid=run_id,
        scr_dir=scr_dir, input_files=input_files
    )


def v7p2_2018_hindcast_atlas():
    '''
    STOFS3D Atl v7.2
    '''

    project_dir = '/sciclone/schism10/feiye/STOFS3D-v8/'
    run_id = '15o_v7'
    hgrid_path = '/sciclone/schism10/feiye/STOFS3D-v8/I15a_v7/FeederDp/hgrid.feeder_dp.ll'
    scr_dir = None

    config = ConfigStofs3dAtlantic.v7p2()
    config.rnday = 396
    config.startdate = datetime(2017, 12, 1)
    config.existing_source_json_path = (
        '/sciclone/schism10/feiye/STOFS3D-v8/I15k_v7/Source_sink/original_source_sink/'
    )
    config.nwm_cache_folder = (
        # '/sciclone/schism10/feiye/STOFS3D-v8/I15/Source_sink/original_source_sink/20171201/'
        '/sciclone/schism10/feiye/STOFS3D-v8/NWM/CONUS/netcdf/CHRTOUT/for_2018_hindcast/'
    )
    config.replace_nwm_with_usgs = True

    input_files = {
        'bctides': False,
        'vgrid': False,
        'gr3': False,
        'nudge_gr3': False,
        'shapiro': False,
        'drag': False,
        'elev_ic': False,
        'source_sink': True,
        'hotstart.nc': False,
        '3D.th.nc': False,
        'elev2D.th.nc': False,
        '*nu.nc': False,
        '*.prop': False,
    }

    stofs3d_atl_driver(
        hgrid_path=hgrid_path, config=config,
        project_dir=project_dir, runid=run_id,
        scr_dir=scr_dir, input_files=input_files
    )



def fariborz_case():
    '''
    Fariborz's case setup
    '''
    project_dir = '/sciclone/schism10/feiye/TEMP/clip_by_polygon/'
    runid = '01'
    hgrid_path = '/sciclone/schism10/feiye/TEMP/clip_by_polygon/hgrid.gr3'
    config = ConfigStofs3dAtlantic.v7p2_subset()
    config.rnday = 1
    config.startdate = datetime(2018, 9, 1)

    scr_dir = None

    input_files = {
        'bctides': False,
        'vgrid': False,
        'gr3': False,
        'nudge_gr3': False,
        'shapiro': False,
        'drag': False,
        'elev_ic': False,
        'source_sink': True,
        'hotstart.nc': False,
        '3D.th.nc': False,
        'elev2D.th.nc': False,
        '*nu.nc': False,
        '*.prop': False,
    }

    stofs3d_atl_driver(
        hgrid_path=hgrid_path, config=config,
        project_dir=project_dir, runid=runid,
        scr_dir=scr_dir, input_files=input_files
    )


def my_second_run():
    '''
    Placeholder for another setup.
    Change the name of this function as needed.
    '''
    return None


if __name__ == '__main__':
    La_v8_2024_reforecast()
