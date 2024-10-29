#!/usr/bin/env python3
'''
This is the driver script to set up the STOFS-3D-ATL model.
Simpler tasks such as generating uniform *.gr3 are included in this script.
More complex tasks such as generating source_sink are imported as modules from subfolders.

Usage:
    Under the current schism git directory,

    python stofs3d-atl-driver.py

    You can discard any changes when you git pull again,
    because a copy of the script and imported modules will be automatically saved in the model
    input folder to keep a record of the model setup.

    Alternatively, you can copy the current directory (Pre_processing/) to your preferred
    location and run the script there.

Default settings for the STOFS-3D-ATL model are set in the ConfigStofs3dAtlantic class.
Read the docstring of the main function for how to use a configuration preset.
Also, set paths at the beginning of the main function for each new setup.

The script will prepare the model folders and keep a record of itself in the
model input folder.
'''


# ------------------------- Import modules ---------------------------
import os
import subprocess
import socket
from pathlib import Path
import shutil
from datetime import datetime
import copy

import numpy as np
from scipy.spatial import cKDTree
import shapefile  # from pyshp

# self-defined modules
import pylib
from pylib import inside_polygon, read_schism_reg
from pylib_experimental.schism_file import source_sink
if 'sciclone' in socket.gethostname():
    from pylib_experimental.schism_file import cread_schism_hgrid as schism_read
    print('Using c++ function to accelerate hgrid reading')
else:
    from pylib import schism_grid as schism_read
    print('Using python function to read hgrid')
from pyschism.mesh import Hgrid

# Import from the sub folders. These are not from installed packages.
from Source_sink.NWM.gen_sourcesink_nwm import gen_sourcesink_nwm
from Source_sink.Constant_sinks.set_constant_sink import set_constant_sink
from Source_sink.Relocate.relocate_source_feeder import relocate_sources2
from Source_sink.Relocate.relocate_source_feeder import v19p2_mandatory_sources_coor
from Prop.gen_tvd import gen_tvd_prop
from Bctides.bctides.bctides import Bctides  # temporary, bctides.py will be merged into pyschism
# from pyschism.forcing.bctides import Bctides

# Global variables:
# use the full path of the Pre-Processing dir inside your schism repo
# e.g, script_path = '/my_dir/schism/src/Utility/Pre-Processing/'
# If you are running this script in the schism repo, you can also use the following line:
script_path = os.path.dirname(os.path.realpath(__file__))
print(f"script_path: {script_path}")

DRIVER_PRINT_PREFIX = '\n-----------------STOFS3D-ATL driver:---------------------\n'


# ---------------------------------------------------------------------
#                               Classes
# ---------------------------------------------------------------------
class ConfigStofs3dAtlantic():
    '''A class to handle the configuration of STOFS-3D-ATL model,
    i.e., processing the parameters and storing the factory settings.
    '''
    def __init__(
        self,
        startdate=datetime(2017, 12, 1),  # start date of the model
        rnday=60,  # number of days to run the model
        nudging_zone_width=1.5,  # in degrees
        nudging_day=1.0,  # in days
        shapiro_zone_width=2.5,  # in degrees
        shapiro_tilt=2.0,  # more abrupt transition in the shapiro zone
        bctides_flags=None,  # a list of lists, each sublist is a set of flags for an open boundary
        nwm_cache_folder=None,
        relocate_source=True,
        feeder_info_file=None,  # file containing feeder info,
                                # made by make_feeder_channel.py in RiverMapper
        gr3_values=None,
        tvd_regions=None
    ):

        self.startdate = startdate
        self.rnday = rnday
        self.nudging_zone_width = nudging_zone_width
        self.nudging_day = nudging_day
        self.shapiro_zone_width = shapiro_zone_width
        self.shapiro_tilt = shapiro_tilt
        self.relocate_source = relocate_source
        self.nwm_cache_folder = nwm_cache_folder
        self.feeder_info_file = feeder_info_file
        if bctides_flags is None:
            self.bctides_flags = [
                [3, 3, 0, 0],  # tides for elev and vel
                [3, 3, 0, 0],  # tides for elev and vel
                [3, 3, 0, 0],  # tides for elev and vel
                [3, 3, 0, 0],  # tides for elev and vel
                [3, 3, 0, 0],  # tides for elev and vel
            ]
        else:
            self.bctides_flags = bctides_flags

        if gr3_values is None:
            self.gr3_values = {  # uniform gr3 values
                'albedo': 0.1,
                'diffmax': 1.0,
                'diffmin': 1e-6,
                'watertype': 1.0,
                'windrot_geo2proj': 0.0
            }
        else:
            self.gr3_values = gr3_values

        if tvd_regions is None:
            self.tvd_regions = []
        else:
            self.tvd_regions = tvd_regions

    @classmethod
    def v6(cls):
        '''Factory method to create a configuration for STOFS3D-v6'''
        return cls(
            nudging_zone_width=.3,  # very wide nudging zone
            shapiro_zone_width=11.5,  # very wide shapiro zone
            shapiro_tilt=3.5,  # very abrupt transition in the shapiro zone
            feeder_info_file=(
                '/sciclone/schism10/feiye/STOFS3D-v5/Inputs/v14/Parallel/'
                'SMS_proj/feeder/feeder.pkl'),
            bctides_flags=[[5, 5, 4, 4]],
            tvd_regions=[
                'tvd0_1.reg', 'tvd0_2.reg', 'tvd0_3.reg', 'tvd0_4.reg',
                'tvd0_5.reg', 'tvd0_6.reg', 'tvd0_7.reg',
                'upwind_east_Carribbean.rgn', 'upwind_west_Carribbean.rgn'
            ]
        )

    @classmethod
    def v7(cls):
        '''Factory method to create a configuration for STOFS3D-v7'''
        return cls(
            nudging_zone_width=7.3,  # default nudging zone
            shapiro_zone_width=11.5,  # default shapiro zone
            shapiro_tilt=3.5,  # default abrupt transition in the shapiro zone
            feeder_info_file=(
                '/sciclone/schism10/Hgrid_projects/STOFS3D-v7/v20.0/Feeder/'
                'feeder_heads_bases.xy'),
            nwm_cache_folder=Path('/sciclone/schism10/whuang07/schism20/NWM_v2.1/'),
            bctides_flags=[
                [3, 3, 0, 0],  # Atlantic Ocean
                [3, 3, 0, 0],  # Gulf of St. Lawrence
                [0, 1, 0, 0],  # St. Lawrence River
            ],
            tvd_regions=[
                'tvd0_1.reg', 'tvd0_2.reg', 'tvd0_3.reg', 'tvd0_4.reg',
                'tvd0_5.reg', 'tvd0_6.reg', 'tvd0_7.reg',
                'upwind_east_Carribbean.rgn', 'upwind_west_Carribbean.rgn'
            ]
        )

    @classmethod
    def v7_hercules_test(cls):
        '''Factory method to create a configuration for STOFS3D-v7'''
        return cls(
            nudging_zone_width=7.3,  # default nudging zone
            shapiro_zone_width=11.5,  # default shapiro zone
            shapiro_tilt=3.5,  # default abrupt transition in the shapiro zone
            relocate_source=True,  # need the feeder info file
            feeder_info_file=(
                '/work/noaa/nosofs/feiye/STOFS-3D-Atl-Example-Setup/DATA/'
                'Feeder_channels/feeder_heads_bases.xy'),
            nwm_cache_folder=None,
            bctides_flags=[
                [3, 3, 0, 0],  # Atlantic Ocean
                [3, 3, 0, 0],  # Gulf of St. Lawrence
                [0, 1, 0, 0],  # St. Lawrence River
            ],
            tvd_regions=[
                'tvd0_1.reg', 'tvd0_2.reg', 'tvd0_3.reg', 'tvd0_4.reg',
                'tvd0_5.reg', 'tvd0_6.reg', 'tvd0_7.reg',
                'upwind_east_Carribbean.rgn', 'upwind_west_Carribbean.rgn'
            ]
        )

    @classmethod
    def v8_louisianna(cls):
        '''Factory method to create a configuration for STOFS3D-v8's local test in Louisianna'''
        return cls(
            nudging_zone_width=0,  # default nudging zone
            shapiro_zone_width=0,  # default shapiro zone
            shapiro_tilt=0,  # default abrupt transition in the shapiro zone
            feeder_info_file='',
            relocate_source=False,
            nwm_cache_folder=Path('/sciclone/schism10/whuang07/schism20/NWM_v2.1/'),
            bctides_flags=[[5, 5, 4, 4]]
        )


# ---------------------------------------------------------------------
#                         Utility functions
# ---------------------------------------------------------------------
def mkcd_new_dir(path):
    '''Make a new directory and change to it'''
    remove_folder(path)
    os.makedirs(path)
    os.chdir(path)


def remove_folder(path):
    '''Remove a folder and all its contents'''
    try:
        shutil.rmtree(path)
        print(f"{path} and all contents removed successfully")
    except FileNotFoundError:
        print(f"{path} does not exist, no need to remove.")
    except OSError as error:
        print(f"Error removing {path}: {error}")
        raise


def try_remove(file):
    '''Try to remove a file, ignore if it does not exist'''
    try:
        os.remove(file)
    except FileNotFoundError:
        print(f"{file} does not exist, no need to remove.")
    except OSError:
        print(f"Error removing: {file}")
        raise


def find_points_in_polyshp(pt_xy, shapefile_names):
    '''Find points inside polygons defined in shapefiles'''
    ind = np.zeros(pt_xy[:, 0].shape)
    for shapefile_name in shapefile_names:
        # shapefile # records = sf.records() # sf.shapeType # len(sf) # s = sf.shape(idx)
        sf = shapefile.Reader(shapefile_name)
        shapes = sf.shapes()

        for i, shp in enumerate(shapes):
            poly_xy = np.array(shp.points).T
            print(f'shape {i+1} of {len(shapes)}, {poly_xy[:, 0]}')
            ind += inside_polygon(pt_xy, poly_xy[0], poly_xy[1])  # 1: true; 0: false

    ind = ind.astype('bool')
    return ind


def prep_run_dir(parent_dir, runid, scr_dir=None):
    '''
    Prepare the run directory and return the path to the model input folder
    - parent_dir: the parent directory where the run directory is created
    - runid: the run directory name
    - scr_dir: the parent directory where outputs are stored, e.g.:
        scr_dir/R{runid}/outputs/,
        if None, use the run directory
        parent_dir/R{runid}/outputs/

    Model inputs and the scripts are saved in the Input dir, e.g., I01;
    The run directory is saved in the Run dir, e.g., R01;
    The outputs are saved on the scratch disk and linked under parent dir. e.g.,
    parent_dir/O01/outputs/.
    Other post-processed outputs can be put under parent_dir/O01/
    '''
    if scr_dir is None:
        scr_dir = parent_dir

    model_input_path = f'{parent_dir}/I{runid}'
    rundir = f'{parent_dir}/R{runid}'
    output_dir = f'{scr_dir}/R{runid}/outputs'

    os.makedirs(rundir, exist_ok=True)
    os.makedirs(model_input_path, exist_ok=True)
    os.makedirs(output_dir, exist_ok=True)
    os.makedirs(f'{parent_dir}/O{runid}', exist_ok=True)

    os.chdir(rundir)
    subprocess.run(
        f'ln -sf {output_dir} .',
        shell=True, stdout=subprocess.DEVNULL,
        stderr=subprocess.DEVNULL, check=False)  # ignore errors, existing folder okay

    os.chdir(f'{parent_dir}/O{runid}')
    os.system(f'ln -sf {output_dir} .')
    os.system(f'ln -sf {model_input_path}/hgrid.gr3 .')

    return model_input_path


# ---------------------------------------------------------------------
# Simple pre-processing functions that do not need subfolders
# ---------------------------------------------------------------------

def gen_nudge_coef(hgrid: pylib.schism_grid, rlmax=1.5, rnu_day=0.25, open_bnd_list=None):
    """
    Set up nudge zone within rlmax distance from the ocean boundary:
    - hgrid must be in lon/lat coordinates
    - rlmax can be a uniform value, e.g., rl_max = 1.5,
        or a 2D array of the same size as the hgrid's number of nodes,
        e.g., rl_max = np.zeros(NP) with further tuning of the nudging zone width.
    - If there are more than one ocean boundary,
        specify the boundary index in open_bnd_list as a list
    """
    if open_bnd_list is None:
        open_bnd_list = [0, 1]

    # maximum nudging strength
    rnu_max = 1.0 / rnu_day / 86400.0

    nudge_coeff = np.zeros((hgrid.np, ), dtype=float)
    for open_bnd_idx in open_bnd_list:
        print(
            f'boundary {open_bnd_idx}: {len(hgrid.iobn[open_bnd_idx])} open boundary nodes')
        bnd_node_idx = hgrid.iobn[open_bnd_idx]

        # distance between each node to its nearest open boundary node
        distance, _ = cKDTree(np.c_[
            hgrid.x[bnd_node_idx], hgrid.y[bnd_node_idx]
        ]).query(np.c_[hgrid.x, hgrid.y])

        # nudge coefficient based on the current open boundary
        this_nudge_coeff = np.clip((1.0-distance/rlmax)*rnu_max, 0.0, rnu_max)
        nudge_coeff = np.maximum(this_nudge_coeff, nudge_coeff)

    return nudge_coeff


def gen_drag(hgrid: pylib.schism_grid):
    '''generate drag coefficient based on the depth and regions'''

    # 1) overall: depth based
    grid_depths = [-3, -1]
    drag_coef = [0.025, 0.0025]
    # linear interpolation with constant extrapolation of nearest end values
    drag = np.interp(hgrid.dp, grid_depths, drag_coef, left=drag_coef[0], right=drag_coef[-1])

    # 2) tweak: regions with constant drag
    region_files = [
        f'{script_path}/Gr3/Drag/Lake_Charles_0.reg',
    ]
    region_tweaks = [0.0]
    for region_tweak, region_file in zip(region_tweaks, region_files):
        reg = read_schism_reg(region_file)
        idx = inside_polygon(np.c_[hgrid.x, hgrid.y], reg.x, reg.y).astype(bool)
        drag[idx] = region_tweak

    # 4) tweak: for nodes inside GoME and dp > -1 m, set drag between 0.01 and 0.02 based on depth
    grid_depths = [5, 20]
    drag_coef = [0.02, 0.01]
    reg = read_schism_reg(f'{script_path}/Gr3/Drag/GoME2.reg')
    idx = inside_polygon(np.c_[hgrid.x, hgrid.y], reg.x, reg.y).astype(bool) & (hgrid.dp > -1)
    drag[idx] = np.interp(
        hgrid.dp[idx], grid_depths, drag_coef,
        left=drag_coef[0], right=drag_coef[-1])

    # 5) tweak: limit some areas to be no more than 0.0025
    for region in [
        f'{script_path}/Gr3/Drag/Portland_max0.0025.reg',
        f'{script_path}/Gr3/Drag/Chatham_max0.0025.reg'
    ]:
        bp = read_schism_reg(region)
        idx = inside_polygon(np.c_[hgrid.x, hgrid.y], bp.x, bp.y).astype(bool)
        drag[idx] = np.minimum(drag[idx], 0.0025)

    # 3) tweak: regions with constant drag but only for river nodes
    region_tweaks = {
        'Mississippi_0': {  # all segments of the Mississippi River
            'drag': 1e-8,
            'region_file': f'{script_path}/Gr3/Drag/Mississippi_0.reg'
        },
        'Mississippi_downstream_0': {  # only the downstream part of the Mississippi River
            'drag': 0,
            'region_file': f'{script_path}/Gr3/Drag/Mississippi_downstream_0.reg'
        },
        'Eastport_0': {
            'drag': 0.0,
            'region_file': f'{script_path}/Gr3/Drag/Eastport_0.reg'
        }
    }
    for _, tweak in region_tweaks.items():
        print(f"Applying drag {tweak['drag']} in {tweak['region_file']}")
        reg = read_schism_reg(tweak['region_file'])
        idx = inside_polygon(np.c_[hgrid.x, hgrid.y], reg.x, reg.y).astype(bool)
        river = hgrid.dp > 0
        drag[idx & river] = tweak['drag']

    return drag


def gen_shapiro_strength(hgrid: pylib.schism_grid, init_shapiro_dist: np.ndarray = None, tilt=2.5):
    '''generate shapiro strength based on the depth and regions'''

    shapiro_min = 100
    shapiro_max = 1000

    # ---------- dependency on the ocean nudging zone ------------
    if init_shapiro_dist is None:
        # no depency on nudging zone, will be set to shapiro_min eventually
        shapiro_ocean_bnd = np.zeros((hgrid.np, ), dtype=float)
    else:
        # scale nudging strength (actually shapiro strength) so that the max equals shapiro_max
        shapiro_ocean_bnd = shapiro_max * (init_shapiro_dist / max(init_shapiro_dist))
    print(f'initial shapiro_nudge: {max(shapiro_ocean_bnd)}, {min(shapiro_ocean_bnd)}')
    # scale nudging transition by a factor, cutoff shapiro at shapiro_max
    # > 1 for more abrupt transition
    shapiro_ocean_bnd = np.minimum(shapiro_max, shapiro_ocean_bnd * tilt)

    # ---------- dependency on regions ------------
    shapiro_region = np.zeros(hgrid.np)
    region_tweaks = {
        'coastal_buffer_1:': {
            'filename': f'{script_path}/Gr3/Shapiro/coastal_0.2.lonlat.reg',
            'strength': 200
        },
        'coastal_buffer_2:': {
            'filename': f'{script_path}/Gr3/Shapiro/coastal_0.5_1.lonlat.reg',
            'strength': 300
        },
        'coastal_buffer_3:': {
            'filename': f'{script_path}/Gr3/Shapiro/coastal_0.5_2.lonlat.reg',
            'strength': 300
        }
    }
    for _, tweak in region_tweaks.items():
        reg = read_schism_reg(tweak['filename'])
        idx = inside_polygon(np.c_[hgrid.x, hgrid.y], reg.x, reg.y).astype(bool)
        shapiro_region[idx] = tweak['strength']

    # -------- combine all shapiro tweaks ---------
    shapiro = shapiro_region + shapiro_ocean_bnd
    shapiro = np.clip(shapiro, shapiro_min, shapiro_max)

    return shapiro


def gen_elev_ic(hgrid=None, h0=0.1, city_shape_fnames=None, base_elev_ic=None):
    '''
    set initial elevation: 0 in the ocean, h0 below ground on higher grounds and in cities

    Inputs:
    - city_shape_fname: the shapefile must be in the same projection as the hgrid
    - base_elev_ic: the base elevation ic file, if not None,
                    the initial ocean elevation will be set to this value
    '''
    if city_shape_fnames is None:
        in_city = np.ones(hgrid.dp.shape, dtype=bool)
    else:
        in_city = find_points_in_polyshp(
            pt_xy=np.c_[hgrid.x, hgrid.y], shapefile_names=city_shape_fnames)

    high_ground = hgrid.dp < 0

    if base_elev_ic is not None:
        elev_ic = schism_read(base_elev_ic).dp
    else:
        elev_ic = np.zeros(hgrid.dp.shape, dtype=float)

    ind = np.logical_or(high_ground, in_city)
    hgrid.save('city_high_ground.gr3', value=ind.astype(int))  # diagnostic

    # set initial elevation to h0 below ground on higher grounds and in cities
    elev_ic[ind] = - hgrid.dp[ind] - h0

    return elev_ic


# ---------------------------------------------------------------------
#       Main function to generate inputs for STOFS-3D-ATL
# ---------------------------------------------------------------------
def main():
    '''
    Main function to generate inputs for STOFS3D-ATL.

    Set the configuration parameters at the beginning of this function,
    i.e., between "---------input---------" and "---------end input---------".

    hgrid.gr3 (in lon/lat, with boundaries) must be prepared before running this script
    '''

    # -----------------input---------------------
    # hgrid generated by SMS, pre-processed, and converted to *.gr3
    hgrid_path = ('/work/noaa/nosofs/feiye/STOFS-3D-Atl-Example-Setup/'
                  '2D_Setup/hgrid_levee_xGEOID20b.gr3')

    # get a configuration preset and make changes if necessary
    # alternatively, you can set the parameters directly on an
    # new instance of ConfigStofs3dAtlantic
    config = ConfigStofs3dAtlantic.v7_hercules_test()
    config.rnday = 3
    config.startdate = datetime(2024, 3, 5)

    # define the project dir, where the run folders are located
    project_dir = '/work/noaa/nosofs/feiye/STOFS-3D-Atl-Example-Setup/2D_Setup/'

    # run ID. e.g, 00, 01a, 10b, etc.; the run folders are named using this ID as follows:
    # I{runid}: input directory for holding model input files and scripts;
    # R{runid}: run directory, where the run will be submitted to queue;
    # O{runid}: output directory for holding raw outputs and post-processing.
    # under project_dir
    runid = '00'

    # swithes to generate differ
    input_files = {
        'bctides': True,
        'vgrid': False,
        'gr3': True,
        'nudge_gr3': True,
        'shapiro': True,
        'drag': True,
        'elev_ic': True,
        'source_sink': True,
        'hotstart.nc': False,
        '*D.th.nc': False,
        '*nu.nc': False,
        '*.prop': True,
    }
    # -----------------end input---------------------

    print(f'{DRIVER_PRINT_PREFIX}reading hgrid from {hgrid_path} ...')
    hgrid = schism_read(hgrid_path)

    # -----------------begin generating model inputs---------------------

    # define and make the model_input_path, the run_dir and the output dir
    model_input_path = prep_run_dir(project_dir, runid)

    # make a copy of the script itself to the model_input_path
    os.system(f'cp -r {script_path} {model_input_path}/')
    # make a copy of the hgrid to the model_input_path
    os.system(f'cp {hgrid_path} {model_input_path}/hgrid.gr3')

    # -----------------bctides---------------------
    if input_files['bctides']:
        hgrid_pyschism = Hgrid.open(hgrid_path, crs='epsg:4326')
        sub_dir = 'Bctides'
        print(f'{DRIVER_PRINT_PREFIX}Generating bctides.in ...')
        mkcd_new_dir(f'{model_input_path}/{sub_dir}')
        Bctides(
            hgrid=hgrid_pyschism,  # needs hgrid in pyschism's Hgrid class
            flags=config.bctides_flags,
            database='fes2014'
        ).write(
            output_directory=f'{model_input_path}/{sub_dir}',
            start_date=config.startdate,
            rnday=config.rnday,
        )

    # -----------------Vgrid---------------------
    if input_files['vgrid']:
        sub_dir = 'Vgrid'
        print(f'{DRIVER_PRINT_PREFIX}Generating vgrid.in ...')
        mkcd_new_dir(f'{model_input_path}/{sub_dir}')
        os.system('ln -sf ../hgrid.gr3 .')

        # call a fortran program to generate vgrid.in
        print(f'compile the fortran program {script_path}/Vgrid/gen_vqs if necessary')
        # fortran_script_path = '/path/to/your/fortran_script.f90'
        # fortran_command = ['gfortran', fortran_script_path, '-o', 'fortran_script']
        fortran_process = subprocess.Popen(
            f'{script_path}/Vgrid/gen_vqs', stdin=subprocess.PIPE
        )  # , stdout=subprocess.PIPE, stderr=subprocess.PIPE)

        # the command line argument 0 means no outputs of a sample vgrid along a given transect
        fortran_process.communicate(input=str(0).encode())

        print(f'{DRIVER_PRINT_PREFIX}converting the format of vgrid.in ...')
        # rename vgrid.in to vgrid.in.old
        try_remove(f'{model_input_path}/{sub_dir}/vgrid.in.old')
        os.rename('vgrid.in', 'vgrid.in.old')

        # convert the format of the vgrid.in file using a fortran script
        subprocess.run(
            [f'{script_path}/Vgrid/change_vgrid'], check=True
        )  # , stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)

        # rename vgrid.in.new to vgrid.in
        try_remove(f'{model_input_path}/{sub_dir}/vgrid.in')
        os.rename('vgrid.in.new', 'vgrid.in')

        os.chdir(model_input_path)
        os.system(f'ln -sf {sub_dir}/vgrid.in .')

    # -----------------spatially uniform Gr3---------------------
    if input_files['gr3']:
        sub_dir = 'Gr3'
        print(f'{DRIVER_PRINT_PREFIX}Generating spatially uniform gr3 files ...')
        mkcd_new_dir(f'{model_input_path}/{sub_dir}')

        for name, value in config.gr3_values.items():
            print(f'{DRIVER_PRINT_PREFIX}Generating {name}.gr3 ...')
            try_remove(f'{model_input_path}/{sub_dir}/{name}.gr3')
            hgrid.write(f'{model_input_path}/{sub_dir}/{name}.gr3', value=value)

        os.chdir(model_input_path)

    # -------------------------------------------------
    # ---------begin spatially varying Gr3 ------------

    # -----------------nudge.gr3---------------------
    if input_files['nudge_gr3']:
        sub_dir = 'Nudge_gr3'
        print(f'{DRIVER_PRINT_PREFIX}Generating nudge.gr3 ...')
        mkcd_new_dir(f'{model_input_path}/{sub_dir}')

        # generate nudging coefficient based on the proximity to the open boundaries
        nudge_coef = gen_nudge_coef(
            hgrid, rlmax=config.nudging_zone_width, rnu_day=config.nudging_day)

        # write nudge.gr3
        hgrid.save(f'{model_input_path}/{sub_dir}/nudge.gr3', value=nudge_coef)

        # link nudge.gr3 to SAL_nudge.gr3 and TEM_nudge.gr3
        try_remove(f'{model_input_path}/{sub_dir}/SAL_nudge.gr3')
        os.system('ln -s nudge.gr3 SAL_nudge.gr3')
        try_remove(f'{model_input_path}/{sub_dir}/TEM_nudge.gr3')
        os.system('ln -s nudge.gr3 TEM_nudge.gr3')

        os.chdir(model_input_path)

    # -----------------shapiro.gr3---------------------
    if input_files['shapiro']:
        sub_dir = 'Shapiro'
        print(f'{DRIVER_PRINT_PREFIX}Generating shapiro.gr3 ...')
        mkcd_new_dir(f'{model_input_path}/{sub_dir}')

        if config.shapiro_zone_width > 0:
            # use the same method for nudge.gr3 to generate a buffer zone along the open boundaries
            distribute_coef = gen_nudge_coef(hgrid, rlmax=config.shapiro_zone_width)
        else:
            distribute_coef = None
        # using a distribution coefficient simliar to the nudging coefficient
        shapiro = gen_shapiro_strength(
            hgrid, init_shapiro_dist=distribute_coef, tilt=config.shapiro_tilt)

        hgrid.save(f'{model_input_path}/{sub_dir}/shapiro.gr3', value=shapiro)

        os.chdir(model_input_path)

    # -----------------drag.gr3---------------------
    if input_files['drag']:
        sub_dir = 'Drag'
        print(f'{DRIVER_PRINT_PREFIX}Generating drag.gr3 ...')
        mkcd_new_dir(f'{model_input_path}/{sub_dir}')
        drag = gen_drag(hgrid)

        hgrid.save(f'{model_input_path}/{sub_dir}/drag.gr3', value=drag)

        os.chdir(model_input_path)

    # -----------------elev_ic---------------------
    if input_files['elev_ic']:
        sub_dir = 'Elev_ic'
        print(f'{DRIVER_PRINT_PREFIX}Generating elev_ic.gr3 ...')
        mkcd_new_dir(f'{model_input_path}/{sub_dir}')
        elev_ic = gen_elev_ic(
            hgrid, h0=0.1,
            city_shape_fnames=[f'{script_path}/Hotstart/LA_urban_polys_lonlat_v2.shp'],
            base_elev_ic=None  # '/sciclone/schism10/feiye/STOFS3D-v7/Inputs/
                               # Iv7_Francine_update/Reinit_hot/elev0.gr3'
        )

        hgrid.save(f'{model_input_path}/{sub_dir}/elev_ic.gr3', value=elev_ic)
        os.symlink(
            f'{model_input_path}/{sub_dir}/elev_ic.gr3',
            f'{model_input_path}/{sub_dir}/elev.ic')

        os.chdir(model_input_path)

    # ----- end spatially varying Gr3 -----------------
    # -------------------------------------------------

    # -----------------source_sink---------------------
    if input_files['source_sink']:
        sub_dir = 'Source_sink'
        print(f'{DRIVER_PRINT_PREFIX}Generating source_sink.in ...')

        # mkcd_new_dir(f'{model_input_path}/{sub_dir}')

        # generate source_sink files by intersecting NWM river segments
        # with the model land boundary
        mkcd_new_dir(f'{model_input_path}/{sub_dir}/original_source_sink/')
        os.symlink(f'{model_input_path}/hgrid.gr3', 'hgrid.gr3')
        gen_sourcesink_nwm(
            startdate=config.startdate, rnday=config.rnday,
            cache_folder=config.nwm_cache_folder)

        # relocate source locations to resolved river channels, the result is the "base" source/sink
        if config.relocate_source:
            mkcd_new_dir(f'{model_input_path}/{sub_dir}/relocated_source_sink/')
            os.symlink(f'{model_input_path}/hgrid.gr3', 'hgrid.gr3')
            # this will generate relocated sources.json and sinks.json
            relocate_sources2(
                old_ss_dir=f'{model_input_path}/{sub_dir}/original_source_sink/',
                feeder_info_file=config.feeder_info_file,
                hgrid_fname=f'{model_input_path}/hgrid.gr3',
                outdir=f'{model_input_path}/{sub_dir}/relocated_source_sink/',
                max_search_radius=2100, mandatory_sources_coor=v19p2_mandatory_sources_coor,
                allow_neglection=False
            )
            # regenerate source/sink based on relocated sources.json and sinks.json
            os.symlink('../original_source_sink/sinks.json', 'sinks.json')  # dummy
            gen_sourcesink_nwm(
                startdate=config.startdate, rnday=config.rnday,
                cache_folder=config.nwm_cache_folder)
            relocated_ss = source_sink.from_files(
                f'{model_input_path}/{sub_dir}/relocated_source_sink/')
            # remove sinks
            os.unlink('vsink.th')
            base_ss = source_sink(
                vsource=relocated_ss.vsource, vsink=None, msource=relocated_ss.msource)
            base_ss.writer(f'{model_input_path}/{sub_dir}/relocated_source_sink/')
        else:
            base_ss = source_sink.from_files(f'{model_input_path}/{sub_dir}/original_source_sink/')

        # set constant sinks (pumps and background sinks)
        mkcd_new_dir(f'{model_input_path}/{sub_dir}/constant_sink/')
        # copy *.shp to the current directory
        os.system(f'cp {script_path}/Source_sink/Constant_sinks/levee_4_pump_polys.* .')
        hgrid_utm = copy.deepcopy(hgrid)
        hgrid_utm.proj(prj0='epsg:4326', prj1='epsg:26918')
        background_ss = set_constant_sink(
            wdir=f'{model_input_path}/{sub_dir}/constant_sink/', hgrid_utm=hgrid_utm)

        # assemble source/sink files and write to model_input_path
        total_ss = base_ss + background_ss
        total_ss.writer(f'{model_input_path}/{sub_dir}/')

        # write diagnostic outputs
        hgrid.compute_ctr()
        np.savetxt(
            f'{model_input_path}/{sub_dir}/vsource.xyz',
            np.c_[hgrid.xctr[total_ss.source_eles-1],
                  hgrid.yctr[total_ss.source_eles-1],
                  total_ss.vsource.df.mean().values]
        )
        np.savetxt(
            f'{model_input_path}/{sub_dir}/vsink.xyz',
            np.c_[hgrid.xctr[total_ss.sink_eles-1],
                  hgrid.yctr[total_ss.sink_eles-1],
                  total_ss.vsink.df.mean().values]
        )

        os.chdir(model_input_path)

    # -----------------*prop---------------------
    if input_files['*.prop']:
        sub_dir = 'Prop'
        print(f'{DRIVER_PRINT_PREFIX}Generating *prop files ...')
        mkcd_new_dir(f'{model_input_path}/{sub_dir}')

        os.system(f'cp {script_path}/Prop/* .')
        gen_tvd_prop(hgrid, config.tvd_regions)

    # -----------------hotstart.nc---------------------
    if input_files['hotstart.nc']:
        sub_dir = 'Hotstart'
        print(f'{DRIVER_PRINT_PREFIX}Generating hotstart.nc ...')
        mkcd_new_dir(f'{model_input_path}/{sub_dir}')

        # generate hotstart.nc
        raise NotImplementedError('hotstart.nc is not implemented yet')

    # -----------------*D.th.nc---------------------
    if input_files['*D.th.nc']:
        sub_dir = 'D.th'
        print(f'{DRIVER_PRINT_PREFIX}Generating *D.th.nc files ...')
        mkcd_new_dir(f'{model_input_path}/{sub_dir}')

        # generate *D.th.nc
        raise NotImplementedError('*D.th.nc is not implemented yet')

    # -----------------*nu.nc---------------------
    if input_files['*nu.nc']:
        sub_dir = 'Nu'
        print(f'{DRIVER_PRINT_PREFIX}Generating *nu.nc files ...')
        mkcd_new_dir(f'{model_input_path}/{sub_dir}')

        # generate *nu.nc
        raise NotImplementedError('*nu.nc is not implemented yet')


def test():
    '''Temporary tests'''
    print("--------------------- Temporary tests ---------------------")

    # inputs
    model_input_path = '/sciclone/schism10/feiye/STOFS3D-v7/Inputs/I12z/'
    config = ConfigStofs3dAtlantic.v7()
    config.startdate = datetime(2024, 3, 5)
    config.rnday = 5
    config.nwm_cache_folder = Path(
        '/sciclone/schism10/feiye/STOFS3D-v7/Inputs/I12z/'
        'Source_sink/original_source_sink/20240305/')
    hgrid = schism_read('/sciclone/schism10/feiye/STOFS3D-v7/Inputs/I12y/hgrid.gr3')
    # end inputs

    sub_dir = 'Source_sink'

    # generate source_sink files by intersecting NWM river segments
    # with the model land boundary
    os.chdir(f'{model_input_path}/{sub_dir}/relocated_source_sink/')
    gen_sourcesink_nwm(
        startdate=config.startdate, rnday=config.rnday,
        cache_folder=config.nwm_cache_folder)
    relocated_ss = source_sink.from_files(
        f'{model_input_path}/{sub_dir}/relocated_source_sink/')

    # set constant sinks (pumps and background sinks)
    mkcd_new_dir(f'{model_input_path}/{sub_dir}/constant_sink/')
    # copy *.shp to the current directory
    os.system(f'cp {script_path}/Source_sink/Constant_sinks/levee_4_pump_polys.* .')
    hgrid_utm = copy.deepcopy(hgrid)
    hgrid_utm.proj(prj0='epsg:4326', prj1='epsg:26918')
    background_ss = set_constant_sink(
        wdir=f'{model_input_path}/{sub_dir}/constant_sink/', hgrid_utm=hgrid_utm)

    # assemble source/sink files and write to model_input_path
    total_ss = relocated_ss + background_ss
    total_ss.writer(f'{model_input_path}/{sub_dir}/')

    # write diagnostic outputs
    hgrid.compute_ctr()
    np.savetxt(
        f'{model_input_path}/{sub_dir}/vsource.xyz',
        np.c_[
            hgrid.xctr[total_ss.source_eles-1],
            hgrid.yctr[total_ss.source_eles-1],
            total_ss.vsource.df.mean().values
        ]
    )
    np.savetxt(
        f'{model_input_path}/{sub_dir}/vsink.xyz',
        np.c_[
            hgrid.xctr[total_ss.sink_eles-1],
            hgrid.yctr[total_ss.sink_eles-1],
            total_ss.vsink.df.mean().values
        ]
    )

    os.chdir(model_input_path)


if __name__ == '__main__':
    # test()
    main()

    print('Done')
