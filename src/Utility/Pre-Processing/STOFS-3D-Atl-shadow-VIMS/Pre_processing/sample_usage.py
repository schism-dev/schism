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


def sample_setup():
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
    config.nwm_cache_folder = (
        '/sciclone/schism10/feiye/STOFS3D-v8/I09g/'
        'Source_sink/original_source_sink/20210801'
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

    # call the driver
    stofs3d_atl_driver(
        hgrid_path=hgrid_path, config=config,
        project_dir=project_dir, runid=runid,
        scr_dir=scr_dir, input_files=input_files
    )


def my_first_run():
    '''
    Placeholder for another setup.
    Change the name of this function as needed.
    '''


def my_second_run():
    '''
    Placeholder for another setup.
    Change the name of this function as needed.
    '''


if __name__ == '__main__':
    sample_setup()
