#!/usr/bin/env python3
'''
This is the driver script to set up the STOFS-3D-ATL model.
Simpler tasks such as generating uniform *.gr3 are included in this script.
More complex tasks such as generating source_sink are imported as modules from subfolders.

Some scripts are written in Fortran and C++ for speed.
You may need to compile them before running this driver script.
See compile instructions at the beginning of each script.

- Vgrid/: gen_vqs.f90, change_vgrid.f90

Usage:
    Under the current schism git directory,

    python stofs3d-atl-driver.py

    You can discard any changes when you git pull again,
    because a copy of the script and imported modules will be automatically saved in the model
    input folder to keep a record of the model setup.

    Alternatively, you can copy the current directory (Pre_processing/) to your preferred
    location and run the script there.

Default settings for the STOFS-3D-ATL model are set in the ConfigStofs3dAtlantic class
in the stofs3d_atl_config.py file.

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
import logging

import numpy as np
from scipy.spatial import cKDTree
import shapefile  # from pyshp
import netCDF4

# self-defined modules
import pylib
from pylib import inside_polygon, read_schism_reg, read_schism_vgrid
from pylib_experimental.schism_file import source_sink
if 'sciclone' in socket.gethostname():
    from pylib_experimental.schism_file import cread_schism_hgrid as schism_read
    print('Using c++ function to accelerate hgrid reading')
else:
    from pylib import schism_grid as schism_read
    print('Using python function to read hgrid')
from pyschism.mesh import Hgrid, Vgrid
from pyschism.forcing.hycom.hycom2schism import OpenBoundaryInventory as OpenBoundaryInventoryHYCOM
from AVISO.aviso2schism import OpenBoundaryInventory as OpenBoundaryInventoryAVISO
from pyschism.forcing.hycom.hycom2schism import Nudge

# Import from the sub folders. These are not from installed packages.
from Source_sink.NWM.gen_sourcesink_nwm import gen_sourcesink_nwm
from Source_sink.Constant_sinks.set_constant_sink import set_constant_sink
from Source_sink.Relocate.relocate_source_feeder import relocate_sources2
from Source_sink.Relocate.relocate_source_feeder import v19p2_mandatory_sources_coor
from Prop.gen_tvd import gen_tvd_prop
from Bctides.bctides.bctides import Bctides  # temporary, bctides.py will be merged into pyschism
# from pyschism.forcing.bctides import Bctides

# Import configuration
from stofs3d_atl_config import ConfigStofs3dAtlantic

# Global variables:
# use the full path of the Pre-Processing dir inside your schism repo
# e.g, script_path = '/my_dir/schism/src/Utility/Pre-Processing/'
# If you are running this script in the schism repo, you can also use the following line:
script_path = os.path.dirname(os.path.realpath(__file__))
print(f"script_path: {script_path}")

DRIVER_PRINT_PREFIX = '\n-----------------STOFS3D-ATL driver:---------------------\n'


# ---------------------------------------------------------------------
#                         Utility functions
# ---------------------------------------------------------------------
def mkcd_new_dir(path, remove=True):
    '''Make a new directory and change to it'''
    if remove:
        remove_folder(path)
    os.makedirs(path, exist_ok=True)
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

    return model_input_path, rundir, output_dir


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


def gen_3dbc(hgrid_fname, vgrid_fname, outdir, start_date, rnday):
    """
    Generate 3D boundary conditions based on the HYCOM data

    Sometimes the server can be busy, and the connection may be lost.
    """
    # check input files
    if not os.path.exists(hgrid_fname):
        raise FileNotFoundError(f'{hgrid_fname} not found')
    if not os.path.exists(vgrid_fname):
        raise FileNotFoundError(f'{vgrid_fname} not found')
    # enable pyschism's logging
    logging.basicConfig(
        format="[%(asctime)s] %(name)s %(levelname)s: %(message)s",
        force=True,
    )
    logger = logging.getLogger('pyschism')
    logger.setLevel(logging.INFO)
    
    max_attempts = 10
    restart = False  # first attempt doesn't need "restart"

    hgrid = Hgrid.open(hgrid_fname, crs='epsg:4326')
    bnd = OpenBoundaryInventoryHYCOM(hgrid, vgrid_fname)
 
    nattempts = max_attempts
    while nattempts > 0:
        try:
            if nattempts < max_attempts:
                restart = True  # continue from the last attempt
            bnd.fetch_data(
                outdir, start_date, rnday,
                elev2D=True, TS=True, UV=True,
                ocean_bnd_ids=[0, 1], restart=restart
            )
            break  # Exit the loop if fetch_data is successful
        except Exception as e:
            nattempts -= 1
            print(f"Attempt failed: {e}. Retrying... ({nattempts} attempts left)")

    if nattempts == 0:
        raise Exception(f"Failed to fetch 3D boundary conditions after {max_attempts} attempts")


def gen_elev2d(hgrid_fname, outdir, start_date, rnday, uniform_shift=0.0):
    '''
    Generate 2D elevation boundary conditions based on AVISO (CMEMS) data

    Needs to pre-download data from AVISO website, because
    some clusters like Hercules and SciClone do not allow direct downloading
    via copernicusmarine.

    See instructions in AVISO/README and sample download script in AVISO/download_aviso*.py

    Save the downloaded data in the current diretory as aviso.nc

    A uniform shift can be applied to the elevation data, e.g., -0.42 m in STOFS-3D v7
    '''
    if not os.path.exists(hgrid_fname):
        raise FileNotFoundError(f'{hgrid_fname} not found')
    # enable pyschism's logging
    logging.basicConfig(
        format="[%(asctime)s] %(name)s %(levelname)s: %(message)s",
        force=True,
    )
    logger = logging.getLogger('pyschism')
    logger.setLevel(logging.INFO)

    while True:  # search for aviso.nc until it is provided
        if not os.path.exists('aviso.nc'):
            print('\n' + '-'*50)
            print(f'aviso.nc not found, please download the data to {outdir}')
            print(f'See instructions in {script_path}/AVISO/README and '
                  f'{script_path}/AVISO/download_aviso*.py')
            input('Press Enter to continue after downloading aviso.nc')
            print('\n' + '-'*50)
        else:
            break
        
    hgrid = Hgrid.open(hgrid_fname, crs='epsg:4326')
    bnd = OpenBoundaryInventoryAVISO(hgrid)
    bnd.fetch_data(outdir, start_date, rnday, ocean_bnd_ids=[0, 1], elev2D=True)

    # modify the elevation data based on the uniform shift
    if np.abs(uniform_shift) > 1e-6:
        print(f'shifting {uniform_shift}')

        with netCDF4.Dataset(f'{outdir}/elev2D.th.nc', mode="r+") as nc_file:
            time_series = nc_file.variables["time_series"]
            print("Original values:", time_series[:10])
            time_series[:] = time_series[:] + uniform_shift
            print("Updated values:", time_series[:10])

        # os.system(f'mv {outdir}/elev2D.th.nc {outdir}/elev2D.th.nc.aviso')
        # ds = xr.open_dataset(f'{outdir}/elev2D.th.nc.aviso', engine="netcdf4")
        # ds['time_series'] += uniform_shift
        # ds.to_netcdf(f'{outdir}/elev2D.th.nc', engine="netcdf4")


def gen_nudge_stofs(hgrid_fname, vgrid_fname, outdir, start_date, rnday):
    '''
    Generate nudge coefficient,
    adapted from pyschism's sample script
    '''
    hgrid = Hgrid.open(hgrid_fname, crs='epsg:4326')

    # ocean_bnd_ids - segment indices, starting from zero.
    nudge = Nudge(hgrid=hgrid, ocean_bnd_ids=[0, 1])

    # rlmax - max relax distance in m or degree
    # rnu_day - max relax strength in days
    # restart = True will append to the existing nc file, works when first try doesn't break.
    nudge.fetch_data(outdir, vgrid_fname, start_date, rnday, restart=False, rnu_day=1, rlmax=7.3)


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
    hgrid_path = ('/sciclone/schism10/feiye/STOFS3D-v8/I11/hgrid.gr3')
    # hgrid_path = ('/work/noaa/nosofs/feiye/Runs/R16x/hgrid.gr3')

    # get a configuration preset and make changes if necessary
    # alternatively, you can set the parameters directly on an
    # new instance of ConfigStofs3dAtlantic
    config = ConfigStofs3dAtlantic.v8()
    config.rnday = 3
    config.startdate = datetime(2024, 3, 5)
    config.nwm_cache_folder = Path(
        '/sciclone/schism10/feiye/STOFS3D-v8/I09/Source_sink/original_source_sink/20240305/')

    # define the project dir, where the run folders are located
    project_dir = '/sciclone/schism10/feiye/STOFS3D-v8/'
    # project_dir = '/work/noaa/nosofs/feiye/Runs/'

    # run ID. e.g, 00, 01a, 10b, etc.; the run folders are named using this ID as follows:
    # I{runid}: input directory for holding model input files and scripts;
    # R{runid}: run directory, where the run will be submitted to queue;
    # O{runid}: output directory for holding raw outputs and post-processing.
    # under project_dir
    runid = '11y'

    # swithes to generate different input files
    input_files = {
        'bctides': False,
        'vgrid': False,
        'gr3': False,
        'nudge_gr3': False,
        'shapiro': False,
        'drag': False,
        'elev_ic': False,
        'source_sink': False,
        'hotstart.nc': False,
        '3D.th.nc': False,
        'elev2D.th.nc': True,
        '*nu.nc': False,
        '*.prop': False,
    }
    # -----------------end input---------------------

    print(f'{DRIVER_PRINT_PREFIX}reading hgrid from {hgrid_path} ...')
    hgrid = schism_read(hgrid_path)

    # -----------------begin generating model inputs---------------------

    # define and make the model_input_path, the run_dir and the output dir
    model_input_path, run_dir, _ = prep_run_dir(project_dir, runid)

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
            bc_flags=config.bc_flags,
            bc_const=config.bc_const,
            bc_relax=config.bc_relax,
            database='fes2014',
        ).write(
            output_directory=f'{model_input_path}/{sub_dir}',
            start_date=config.startdate,
            rnday=config.rnday,
        )

        os.chdir(run_dir)
        os.system(f'ln -sf ../I{runid}/{sub_dir}/bctides.in .')
        os.chdir(model_input_path)

    # -----------------Vgrid---------------------
    if input_files['vgrid']:
        sub_dir = 'Vgrid'
        print(f'{DRIVER_PRINT_PREFIX}Generating vgrid.in ...')
        mkcd_new_dir(f'{model_input_path}/{sub_dir}')
        os.system('ln -sf ../hgrid.gr3 .')

        # call a fortran program to generate vgrid.in
        print(f'compile the fortran program {script_path}/Vgrid/gen_vqs if necessary')
        fortran_process = subprocess.Popen(
            f'{script_path}/Vgrid/gen_vqs', stdin=subprocess.PIPE
        )

        # the command line argument 0 means no outputs of a sample vgrid along a given transect
        fortran_process.communicate(input=str(0).encode())

        print(f'{DRIVER_PRINT_PREFIX}converting the format of vgrid.in ...')
        vg = read_schism_vgrid(f'{model_input_path}/{sub_dir}/vgrid.in')
        os.rename('vgrid.in', 'vgrid.in.old')
        vg.save(f'{model_input_path}/{sub_dir}/vgrid.in')

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

        os.chdir(run_dir)
        os.system(f'ln -sf ../I{runid}/{sub_dir}/*.gr3 .')
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

        os.chdir(run_dir)
        os.system(f'ln -sf ../I{runid}/{sub_dir}/*_nudge.gr3 .')
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

        os.chdir(run_dir)
        os.system(f'ln -sf ../I{runid}/{sub_dir}/shapiro.gr3 .')
        os.chdir(model_input_path)

    # -----------------drag.gr3---------------------
    if input_files['drag']:
        sub_dir = 'Drag'
        print(f'{DRIVER_PRINT_PREFIX}Generating drag.gr3 ...')
        mkcd_new_dir(f'{model_input_path}/{sub_dir}')
        drag = gen_drag(hgrid)

        hgrid.save(f'{model_input_path}/{sub_dir}/drag.gr3', value=drag)

        os.chdir(run_dir)
        os.system(f'ln -sf ../I{runid}/{sub_dir}/drag.gr3 .')
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

        os.chdir(run_dir)
        os.system(f'ln -sf ../I{runid}/{sub_dir}/elev.ic .')
        os.chdir(model_input_path)

    # ----- end spatially varying Gr3 -----------------
    # -------------------------------------------------

    # -----------------source_sink---------------------
    if input_files['source_sink']:
        sub_dir = 'Source_sink'
        print(f'{DRIVER_PRINT_PREFIX}Generating source_sink.in ...')

        # mkcd_new_dir(f'{model_input_path}/{sub_dir}')

        # # generate source_sink files by intersecting NWM river segments
        # # with the model land boundary
        # mkcd_new_dir(f'{model_input_path}/{sub_dir}/original_source_sink/')
        # os.symlink(f'{model_input_path}/hgrid.gr3', 'hgrid.gr3')
        # gen_sourcesink_nwm(
        #     startdate=config.startdate, rnday=config.rnday,
        #     cache_folder=config.nwm_cache_folder)

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
                max_search_radius=2100, mandatory_sources_coor=config.mandatory_sources_coor,
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

        os.chdir(run_dir)
        os.system(f'ln -sf ../I{runid}/{sub_dir}/source.nc .')
        os.chdir(model_input_path)

    # -----------------*prop---------------------
    if input_files['*.prop']:
        sub_dir = 'Prop'
        print(f'{DRIVER_PRINT_PREFIX}Generating *prop files ...')
        mkcd_new_dir(f'{model_input_path}/{sub_dir}')

        os.system(f'cp {script_path}/Prop/* .')
        gen_tvd_prop(hgrid, config.tvd_regions)

        os.chdir(run_dir)
        os.system(f'ln -sf ../I{runid}/{sub_dir}/tvd.prop .')

    # -----------------hotstart.nc---------------------
    if input_files['hotstart.nc']:
        sub_dir = 'Hotstart'
        print(f'{DRIVER_PRINT_PREFIX}Generating hotstart.nc ...')
        mkcd_new_dir(f'{model_input_path}/{sub_dir}')

        # generate hotstart.nc
        raise NotImplementedError('hotstart.nc is not implemented yet')

    # -----------------3D.th.nc---------------------
    if input_files['3D.th.nc']:
        sub_dir = '3Dth'
        print(f'{DRIVER_PRINT_PREFIX}Generating *3D.th.nc files ...')
        mkcd_new_dir(f'{model_input_path}/{sub_dir}')

        # generate *D.th.nc
        gen_3dbc(
            hgrid_fname=f'{model_input_path}/hgrid.gr3',
            vgrid_fname=f'{model_input_path}/vgrid.in',
            outdir=f'{model_input_path}/{sub_dir}',
            start_date=config.startdate, rnday=config.rnday
        )

        os.chdir(run_dir)
        os.system(f'ln -sf ../I{runid}/{sub_dir}/*3D.th.nc .')
    
    # -----------------elev2D.th.nc---------------------
    if input_files['elev2D.th.nc']:
        sub_dir = 'Elev2Dth'
        print(f'{DRIVER_PRINT_PREFIX}Generating elev2D.th.nc ...')
        mkcd_new_dir(f'{model_input_path}/{sub_dir}', remove=False)

        gen_elev2d(
            hgrid_fname=f'{model_input_path}/hgrid.gr3',
            outdir=f'{model_input_path}/{sub_dir}',
            start_date=config.startdate, rnday=config.rnday,
            uniform_shift=config.elev2d_uniform_shift
        )

        os.chdir(run_dir)
        os.system(f'ln -sf ../I{runid}/{sub_dir}/elev2D.th.nc .')

    # -----------------*nu.nc---------------------
    if input_files['*nu.nc']:
        sub_dir = 'Nudge'
        print(f'{DRIVER_PRINT_PREFIX}Generating *nu.nc files ...')
        mkcd_new_dir(f'{model_input_path}/{sub_dir}')

        os.system(f'cp {script_path}/Nudge/* .')
        gen_nudge_stofs(
            hgrid_fname=f'{model_input_path}/hgrid.gr3',
            vgrid_fname=f'{model_input_path}/vgrid.in',
            outdir=f'{model_input_path}/{sub_dir}',
            rnday=config.rnday, start_date=config.startdate
        )

        os.chdir(run_dir)
        os.system(f'ln -sf ../I{runid}/{sub_dir}/*nu.nc .')


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
