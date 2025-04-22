'''
This is the driver script to set up the STOFS-3D-ATL model.
Simpler tasks such as generating uniform *.gr3 are included in this script.
More complex tasks such as generating source_sink are imported as modules from subfolders.

Some scripts are written in Fortran and C++ for speed.
You may need to compile them before running this driver script.
See compile instructions at the beginning of each script:

- Vgrid/: gen_vqs.f90, change_vgrid.f90

Usage:
    see sample_usage.py in the same folder.

The script will prepare the model folders and keep a record of itself in the
model input folder.

For the author, there are a few "todo" items in the script.
1) The second call to gen_sourcesink_nwm is not optimal, because it downloads
    the NWM data again if config.nwm_cache_folder is not set.
    This cost time for longer runs.
    It is better to directly use generated sources/sinks.
    As a temporary solution, run the script with config.relocate_source = False,
    locate the downloaded data in {model_input_path}/Source_sink/original_source_sink/,
    symlink all downloaded NWM data in a single folder and set config.nwm_cache_folder.
    The symlink is necessary because NWM data may be stored in subfolders of different years.

2) hotstart.nc

Temporary changes to be tested:
1) admit more BlueTopo in DEM_loading, see bluetopo_region.shp
2) adopted regional depth tweaks from SECOFS, added DEFAULT_REGIONAL_TWEAKS2 in regional_tweaks.py
3) adopted drag tweaks from SECOFS, added a few regions in gen_drag()
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
import json

# self-defined modules
import pylib
from pylib import inside_polygon, read_schism_reg, read_schism_vgrid
from pylib_experimental.schism_file import source_sink, TimeHistory
if 'gulf' in socket.gethostname():
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
from Source_sink.Relocate.relocate_source_feeder import relocate_sources2, remove_duplicate_dict_values
from Source_sink.Relocate.relocate_source_feeder import v19p2_mandatory_sources_coor, find_duplicate_dict_values
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


def gen_3dbc(
    hgrid_fname, vgrid_fname, outdir, start_date, rnday,
    ocean_bnd_ids=None, max_attempts=10
):
    """
    Generate 3D boundary conditions based on the HYCOM data

    Sometimes the server can be busy, and the connection may be lost.
    """
    if ocean_bnd_ids is None:
        raise ValueError("ocean_bnd_ids must be specified, e.g., [0, 1]")

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

    hgrid = Hgrid.open(hgrid_fname, crs='epsg:4326')
    bnd = OpenBoundaryInventoryHYCOM(hgrid, vgrid_fname)

    restart = False  # First attempt doesn't need "restart"
    last_error = None  # To store the last exception
    for attempt in range(max_attempts):
        try:
            bnd.fetch_data(
                outdir, start_date, rnday,
                elev2D=True, TS=True, UV=True,
                ocean_bnd_ids=ocean_bnd_ids,
                restart=restart
            )
            return  # Successful fetch
        except Exception as e:  # Replace with specific exceptions
            last_error = e
            restart = True
            print(f"Attempt failed: {e}. Retrying..."
                  f" ({max_attempts-attempt-1} attempts left)")

    # Raise an exception after exhausting all retries
    raise Exception(
        f"Failed to fetch 3D boundary conditions after {max_attempts} attempts. "
        f"Last error: {last_error}"
    )


def gen_elev2d(
    hgrid_fname, outdir, start_date, rnday,
    ocean_bnd_ids=None, uniform_shift=0.0
):
    '''
    Generate 2D elevation boundary conditions based on AVISO (CMEMS) data

    Needs to pre-download data from AVISO website, because
    some clusters like Hercules and SciClone do not allow direct downloading
    via copernicusmarine.

    See instructions in AVISO/README and sample download script in AVISO/download_aviso*.py

    Save the downloaded data in the current diretory as aviso.nc

    A uniform shift can be applied to the elevation data, e.g., -0.42 m in STOFS-3D v7
    '''
    if ocean_bnd_ids is None:
        raise ValueError("ocean_bnd_ids must be specified, e.g., [0, 1]")
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
    bnd.fetch_data(
        outdir, start_date, rnday, ocean_bnd_ids=ocean_bnd_ids, elev2D=True)

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


def gen_nudge_stofs(
    hgrid_fname, vgrid_fname, outdir, start_date, rnday,
    ocean_bnd_ids=None, max_retry=10, hycom_download_dir=None
):
    '''
    Generate nudge coefficient,
    adapted from pyschism's sample script

    Note: fixed rnu_day=1, rlmax=7.3
    '''
    if ocean_bnd_ids is None:
        raise ValueError("ocean_bnd_ids must be specified, e.g., [0, 1]")

    hgrid = Hgrid.open(hgrid_fname, crs='epsg:4326')

    # ocean_bnd_ids - segment indices, starting from zero.
    nudge = Nudge(hgrid=hgrid, ocean_bnd_ids=ocean_bnd_ids)

    # for nudge.fetch_data:
    # rlmax - max relax distance in m or degree
    # rnu_day - max relax strength in days
    # restart = True will append to the existing nc file, works when first try doesn't break.
    restart = False
    retries_left = max_retry  # Keep track of retries separately
    while retries_left > 0:
        try:
            nudge.fetch_data(
                outdir, vgrid_fname, start_date, rnday, restart=restart, rnu_day=1, rlmax=7.3,
                hycom_download_dir=hycom_download_dir
            )
            return  # Successful fetch
        except Exception as e:  # Replace with specific exceptions as needed
            retries_left -= 1
            print(f"Attempt failed: {e}. Retrying... ({retries_left} attempts left)")
            restart = True

    # Include the last exception in the error message
    raise Exception(f"Failed to fetch nudge coefficient after {max_retry} attempts. Last error: {e}")


def gen_drag(hgrid: pylib.schism_grid):
    '''generate drag coefficient based on the depth and regions'''

    # - overall: depth based
    grid_depths = [-3, -1]
    # default [0.025, 0.0025]; [0.005, 0.0025] for STOFS-3D v8 R20e/f; [0.02, 0.001] for R13r_v7
    drag_coef = [0.025, 0.0025]
    # linear interpolation with constant extrapolation of nearest end values
    drag = np.interp(hgrid.dp, grid_depths, drag_coef, left=drag_coef[0], right=drag_coef[-1])

    # - tweak: regions with constant drag
    region_files = [
        f'{script_path}/Gr3/Drag/Lake_Charles_0.reg',
        # f'{script_path}/Gr3/Drag/Mayport.reg',
    ]
    region_tweaks = [0.0]
    for region_tweak, region_file in zip(region_tweaks, region_files):
        reg = read_schism_reg(region_file)
        idx = inside_polygon(np.c_[hgrid.x, hgrid.y], reg.x, reg.y).astype(bool)
        drag[idx] = region_tweak

    # - tweak: for nodes inside GoME and dp > -1 m, set drag between 0.01 and 0.02 based on depth
    grid_depths = [5, 20]
    drag_coef = [0.02, 0.01]
    reg = read_schism_reg(f'{script_path}/Gr3/Drag/GoME2.reg')
    idx = inside_polygon(np.c_[hgrid.x, hgrid.y], reg.x, reg.y).astype(bool) & (hgrid.dp > -1)
    drag[idx] = np.interp(
        hgrid.dp[idx], grid_depths, drag_coef,
        left=drag_coef[0], right=drag_coef[-1])

    # - tweak: limit some areas to be no more than 0.0025
    for region in [
        f'{script_path}/Gr3/Drag/Portland_max0.0025.reg',
        f'{script_path}/Gr3/Drag/Chatham_max0.0025.reg'
    ]:
        bp = read_schism_reg(region)
        idx = inside_polygon(np.c_[hgrid.x, hgrid.y], bp.x, bp.y).astype(bool)
        drag[idx] = np.minimum(drag[idx], 0.0025)

    # - tweak: regions with constant drag but only for river nodes
    region_tweaks = {
        # 'MayPort_to_Wacha_0.001': {
        #     'drag': 0.001,
        #     'region_file': f'{script_path}/Gr3/Drag/MayPort_to_Wacha_0.001.reg'
        # },
        'NY.reg': {
            'drag': 0.001,
            'region_file': f'{script_path}/Gr3/Drag/NY.reg'
        },
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
        },
        # 'Upper_Delaware_Bay': {  # starting from STOFS-3D v8 R20e/f
        #     'drag': 0.001,
        #     'region_file': f'{script_path}/Gr3/Drag/Upper_Delaware_Bay.reg'
        # },
        'Hudson': {  # starting from STOFS-3D v8 R20e/f
            'drag': 0.001,
            'region_file': f'{script_path}/Gr3/Drag/Hudson.reg'
        },
        'Dalhgren_0': {  # from SECOFS
            'drag': 0.0005,
            'region_file': f'{script_path}/Gr3/Drag/Dalhgren_0.reg'
        },
        'Wilmington_0': {  # from SECOFS
            'drag': 0.000,
            'region_file': f'{script_path}/Gr3/Drag/Wilmington_0.reg'
        },
        'Fernandina_0': {  # from SECOFS, increased drag
            'drag': 0.00,
            'region_file': f'{script_path}/Gr3/Drag/Fernandina_0.reg'
        },
        'Fort_Myers_0': {  # from SECOFS, slightly increased drag
            'drag': 0.000,
            'region_file': f'{script_path}/Gr3/Drag/Fort_Myers_0.reg'
        },
        'Tampa_Bay_0': {  # from SECOFS, increased drag
            'drag': 0.002,
            'region_file': f'{script_path}/Gr3/Drag/Tampa_Bay_0.reg'
        },
        'StJohns_0': {  # from SECOFS
            'drag': 0.0,
            'region_file': f'{script_path}/Gr3/Drag/StJohns_0.reg'
        },
        # 'Virginia_Key_0.0035': {
        #     'drag': 0.01,
        #     'region_file': f'{script_path}/Gr3/Drag/Virginia_Key_0.0035.reg'
        # },
        # 'GoMX_east_0.001': {
        #     'drag': 0.001,
        #     'region_file': f'{script_path}/Gr3/Drag/GoMX_east_0.001.reg'
        # },
        'Oyster_Landing_0.0': {
            'drag': 0.0,
            'region_file': f'{script_path}/Gr3/Drag/Oyster_Landing_0.reg'
        },
        'Wachapreague_0.0': {
            'drag': 0.0,
            'region_file': f'{script_path}/Gr3/Drag/Wachapreague_0.reg'
        },
        'Manchester_0.001': {
            'drag': 0.00,
            'region_file': f'{script_path}/Gr3/Drag/Manchester_0.001.reg'
        },
        'Sabine_Lake_0.001': {
            'drag': 0.001,
            'region_file': f'{script_path}/Gr3/Drag/Sabine_Lake_0.001.reg'
        },
        'High_Island_0': {
            'drag': 0.0,
            'region_file': f'{script_path}/Gr3/Drag/High_Island_0.reg'
        },
        'West_Fowl_River_0': {
            'drag': 0.0,
            'region_file': f'{script_path}/Gr3/Drag/West_Fowl_River_0.reg'
        },
        'Chickasaw_0': {
            'drag': 0.0,
            'region_file': f'{script_path}/Gr3/Drag/Chickasaw_0.reg'
        },
        'Panama_City_0.001.reg': {
            'drag': 0.001,
            'region_file': f'{script_path}/Gr3/Drag/Panama_City_0.001.reg'
        },
        'Money_Pt_0.reg': {
            'drag': 0.0,
            'region_file': f'{script_path}/Gr3/Drag/Money_Pt_0.reg'
        },
        'Nantucket_Island_0.reg': {
            'drag': 0.0,
            'region_file': f'{script_path}/Gr3/Drag/Nantucket_Island_0.reg'
        },
        'Atchafalaya_0.reg': {
            'drag': 0.002,
            'region_file': f'{script_path}/Gr3/Drag/Atchafalaya_0.reg'
        },
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


def gen_relocated_source(original_source_sink_dir, relocated_source_sink_dir):
    """
    Generate relocated vsource.th for STOFS-3D-ATL
    based on sources.json and sinks.json generated by relocate_sources2
    """

    original_ss = source_sink.from_files(source_dir=original_source_sink_dir)

    original_ele_nwm_mapping = json.load(
        open(f'{original_source_sink_dir}/sources.json', 'r')
    )
    # reverse the mapping
    original_nwm_ele_mappping = {}
    for ele, nwm_fids in original_ele_nwm_mapping.items():
        for nwm_fid in nwm_fids:
            if nwm_fid in original_nwm_ele_mappping:
                print(
                    f"Warning: NWM feature ID {nwm_fid} is already mapped to "
                    f"{original_nwm_ele_mappping[nwm_fid]}. Overwriting with {ele}"
                )  # if an nwm_fid is mapped to multiple elements, keep the last one.

                # This can happen if the original source_sink is generated by pyschism,
                # with an NWM segment weaving in and out of multiple elements,
                # resulting in multiple sources and sinks maped to the same nwm_fid.
                # However, the relocation process will ignore all sinks, so only one source
                # should be mapped to the same nwm_fid. And the exact location is somewhat arbitrary.
                # In a rare case, a river can go in as one nwm_fid and out as another nwm_fid,
                # so make sure this doesn't happen by properly aligning the schism land boundary
                # and avoid NWM segments weaving in and out of multiple elements.
            original_nwm_ele_mappping[nwm_fid] = ele

    relocated_ele_nwm_mapping = json.load(
        open(f'{relocated_source_sink_dir}/sources.json', 'r')
    )

    relocated_ele_mapping = {}
    for ele, nwm_fids in relocated_ele_nwm_mapping.items():
        relocated_ele_mapping[ele] = [original_nwm_ele_mappping[nwm_fid] for nwm_fid in nwm_fids]
        # remove duplicates
        relocated_ele_mapping[ele] = list(set(relocated_ele_mapping[ele]))

    vsource_data = np.zeros((original_ss.vsource.n_time, len(relocated_ele_mapping)), dtype=float)
    for i, [relocated_ele, original_eles] in enumerate(relocated_ele_mapping.items()):
        for original_ele in original_eles:
            vsource_data[:, i] += original_ss.vsource.df[original_ele].values

    # assemble relocated vsource.th and msource.th
    relocated_vsource = TimeHistory(
        data_array=np.c_[original_ss.vsource.time, vsource_data], columns=list(relocated_ele_mapping.keys()))
    
    relocated_msource_data = np.ones((original_ss.vsource.n_time, len(relocated_ele_mapping)), dtype=float)
    relocated_msource_list = [
        TimeHistory(
            data_array=np.c_[original_ss.vsource.time, relocated_msource_data*-9999],
            columns=list(relocated_ele_mapping.keys())
        ),  # Temperature
        TimeHistory(
            data_array=np.c_[original_ss.vsource.time, relocated_msource_data*0],
            columns=list(relocated_ele_mapping.keys())
        )  # Salinity
    ]

    return relocated_vsource, relocated_msource_list


# ---------------------------------------------------------------------
#       Main function to generate inputs for STOFS-3D-ATL
# ---------------------------------------------------------------------
def stofs3d_atl_driver(
    hgrid_path: str,
    config: ConfigStofs3dAtlantic,
    project_dir: str, runid: str, scr_dir: str,
    input_files: dict = None,
):
    '''
    Main function to generate inputs for STOFS3D-ATL.
    '''

    if input_files is None:
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
            'tvd.prop': True,
        }

    print(f'{DRIVER_PRINT_PREFIX}reading hgrid from {hgrid_path} ...')
    hgrid = schism_read(hgrid_path)

    # -----------------begin generating model inputs---------------------

    # define and make the model_input_path, the run_dir and the output dir
    model_input_path, run_dir, _ = prep_run_dir(project_dir, runid, scr_dir=scr_dir)

    # make a copy of the script itself to the model_input_path
    os.system(f'cp -r {script_path} {model_input_path}/Pre_processing_scripts_backup')
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
    # Note:

    # Normal case: main hgrid has pseudo channels as feeders.
    # An "hgrid_without_feeders" needs to be provided in the stofs3d_atl_config.py
    # to generate the original sources/sinks.
    # This is necessary because the feeder channels are not real channels
    # and "hgrid_without_feeders" gives the correct source locations.
    # Moreover, the land boundary of "hgrid_without_feeders" is set to avoid
    # NWM segments weaving in and out of the schism domain.
    # The script will then relocate the original sources to the head of each feeder channel.

    # Special case: main hgrid has no feeders or has real channels as feeders.
    # In this case, the original sources/sinks are generated based on the main hgrid directly,
    # and the option "hgrid_without_feeders" should be set to None.

    if input_files['source_sink']:
        sub_dir = 'Source_sink'
        print(f'{DRIVER_PRINT_PREFIX}Generating source_sink.in ...')

        if config.hgrid_without_feeders is not None:
            print(
                'Normal case: beside the main hgrid, '
                f'an hgrid without feeders is provided: {config.hgrid_without_feeders}, '
                'assuming the main hgrid has feeders!'
            )
            main_hgrid_has_feeder = True
        else:
            print(
                'Special case: only the main hgrid is provided,'
                'assuming the main hgrid has NO feeders!'
                'Caution: you should not search for NWM source/sink on a grid with feeders!'
            )
            main_hgrid_has_feeder = False
            
        # '''  comment out the following code to skip generating original source_sink files
        # Generate original source_sink files
        mkcd_new_dir(f'{model_input_path}/{sub_dir}')

        # generate source_sink files by intersecting NWM river segments
        # with the model land boundary
        mkcd_new_dir(f'{model_input_path}/{sub_dir}/original_source_sink/')
        if config.hgrid_without_feeders is not None:
            os.symlink(f'{config.hgrid_without_feeders}', 'hgrid.gr3')
        else:
            os.symlink(f'{model_input_path}/hgrid.gr3', 'hgrid.gr3')

        gen_sourcesink_nwm(
            startdate=config.startdate, rnday=config.rnday,
            cache_folder=config.nwm_cache_folder)
        # '''
        
        # A single NWM segment weaving in and out will create duplicate sources/sinks
        # , so it is not necessary to remove duplicates here.
        # 
        # # find any duplicate sources
        # with open(
        #     f'{model_input_path}/{sub_dir}/original_source_sink/sources.json',
        #     'r', encoding='utf-8'
        # ) as f:
        #     old_sources2fids = json.load(f)
        # fid_list = [fid for fids in old_sources2fids.values() for fid in fids]
        # if len(fid_list) != len(set(fid_list)):
        #     print(f'Number of duplicated fids: {len(fid_list) - len(set(fid_list))}')
        #     # raise ValueError('Duplicated fids in new2fid')

        #     for fid in set(fid_list):
        #         if fid_list.count(fid) > 1:
        #             print(f'Duplicated fid: {fid}')

        #     # backup the original source_sink files
        #     os.system(f'cp -r {model_input_path}/{sub_dir}/original_source_sink/ '
        #               f'{model_input_path}/{sub_dir}/original_source_sink_0/')

        #     # remove duplicated sources in the original source_sink files
        #     old_sources2fids = remove_duplicate_dict_values(old_sources2fids)
        #     # remove keys with empty values
        #     old_sources2fids = {k: v for k, v in old_sources2fids.items() if v}

        #     # regenerate old sources based on updated old_sources2fids
        #     with open(
        #         f'{model_input_path}/{sub_dir}/original_source_sink/sources.json',
        #         'w', encoding='utf-8'
        #     ) as f:
        #         json.dump(old_sources2fids, f, indent=4)
        #     gen_sourcesink_nwm(
        #         startdate=config.startdate, rnday=config.rnday,
        #         cache_folder=config.nwm_cache_folder)


        # Relocate source locations to resolved river channels, the result is the "base" source/sink.
        # Set appropriate no_feeder option, mandatory_sources_coor, and
        # feeder_info_file in stofs3d_atl_config.py
        if config.relocate_source:
            # relocate
            mkcd_new_dir(f'{model_input_path}/{sub_dir}/relocated_source_sink/')
            os.symlink(f'{model_input_path}/hgrid.gr3', 'hgrid.gr3')

            # this will generate relocated sources.json and sinks.json
            relocate_sources2(
                old_ss_dir=f'{model_input_path}/{sub_dir}/original_source_sink/',
                feeder_info_file=config.feeder_info_file,
                hgrid_fname=f'{model_input_path}/hgrid.gr3',
                outdir=f'{model_input_path}/{sub_dir}/relocated_source_sink/',
                max_search_radius=2100, mandatory_sources_coor=config.mandatory_sources_coor,
                allow_neglection=False, main_hgrid_has_feeder=main_hgrid_has_feeder,
            )

            # regenerate vsource.th based on relocated sources.json
            regenerate_using_pyschism = True

            if regenerate_using_pyschism:
                # gen_sourcesin_nwm requires sinks.json
                os.system(f'ln -sf {model_input_path}/{sub_dir}/original_source_sink/sinks.json .')
                gen_sourcesink_nwm(  # with existing sources.json and sinks.json
                    startdate=config.startdate, rnday=config.rnday,
                    cache_folder=config.nwm_cache_folder
                )
                relocated_ss = source_sink.from_files(
                    source_dir=f'{model_input_path}/{sub_dir}/relocated_source_sink/',
                )  # sinks will be discarded later, only sources will be used
                base_ss = source_sink(
                    vsource=relocated_ss.vsource, vsink=None, msource=relocated_ss.msource
                )
            else:
                # In the case the original sources have been adjusted by USGS obs,
                # Don't call gen_sourcesink_nwm again, use gen_relocated_source instead.
                # This may create minor duplicates in the sources, check the print message to
                # see if the duplicated sources are small.
                # This doesn't matter for operation, since only sources.json is used.
                # Todo: remove the duplicated sources;
                relocated_vsource, relocated_msource_list = gen_relocated_source(
                    original_source_sink_dir=f'{model_input_path}/{sub_dir}/original_source_sink/',
                    relocated_source_sink_dir=f'{model_input_path}/{sub_dir}/relocated_source_sink/',
                )
                base_ss = source_sink(
                    vsource=relocated_vsource, vsink=None, msource=relocated_msource_list)
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

        # the final source/sink uses the relocated sources and constant sinks (no json needed)
        os.chdir(f'{model_input_path}/{sub_dir}')
        os.system('ln -sf ./relocated_source_sink/sources.json .')

        # write diagnostic outputs
        hgrid.compute_ctr()
        with open(f'{model_input_path}/{sub_dir}/vsource.xyz', 'w') as f:
            f.write('lon lat vsource\n')
            for i, ele in enumerate(total_ss.source_eles):
                f.write(f'{hgrid.xctr[ele-1]} {hgrid.yctr[ele-1]} {total_ss.vsource.df.iloc[:, i].mean()}\n')
        with open(f'{model_input_path}/{sub_dir}/vsink.xyz', 'w') as f:
            f.write('lon lat vsink\n')
            for i, ele in enumerate(total_ss.sink_eles):
                f.write(f'{hgrid.xctr[ele-1]} {hgrid.yctr[ele-1]} {total_ss.vsink.df.iloc[:, i].mean()}\n')

        os.chdir(run_dir)
        os.system(f'ln -sf ../I{runid}/{sub_dir}/source.nc .')
        os.chdir(model_input_path)

        print('Done generating source/sink files.')

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
        # todo: implement this function
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
            start_date=config.startdate, rnday=config.rnday,
            ocean_bnd_ids=config.ocean_bnd_ids,
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
            ocean_bnd_ids=config.ocean_bnd_ids,
            uniform_shift=config.elev2d_uniform_shift,
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
            rnday=config.rnday, start_date=config.startdate,
            ocean_bnd_ids=config.ocean_bnd_ids,
            # use pre-downloaded files to speed up
            hycom_download_dir='/sciclone/schism10/feiye/STOFS3D-v8/I15_v7/HYCOM_files_2018_hindcast/',
            # hycom_download_dir=None,
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


def sample2():
    '''Subsetting example workflow for STOFS3D Atlantic v7.2,
    assuming the original source_sink files have been generated by gen_sourcesink_nwm'''

    import geopandas as gpd
    from shapely import get_coordinates
    import json
    from Source_sink.Relocate.relocate_source_feeder import v19p2_for_sms_v27_mandatory_sources_coor
    
    region_gdf = gpd.read_file(f'/sciclone/schism10/feiye/TEMP/clip_by_polygon/hires/hires.shp').to_crs('EPSG:4326')
    region_list = [get_coordinates(p) for p in region_gdf.explode(index_parts=True).exterior]

    wdir = '/sciclone/schism10/feiye/TEMP/clip_by_polygon/I01/Source_sink/'

    # read original source_sink files
    original_ss = source_sink.from_files(source_dir=f'{wdir}/original_source_sink/')

    # subset
    inside_ss, outside_ss = original_ss.clip_by_polygons(
        hgrid=schism_read(f'{wdir}/original_source_sink/hgrid.gr3'),  # can also use pylib's read()
        polygons_xy=region_list
    )

    # relocate sources in high-res regions (inside_ss)
    # reassign sources.json and sinks.json
    for i, json_file in enumerate(['sources.json', 'sinks.json']):
        original_ele2nwm_fid = json.load(open(f'{wdir}/original_source_sink/{json_file}', 'r'))
        ele2nwm_fid = {
            str(int(k)): v for k, v in original_ele2nwm_fid.items()
            if int(k) in inside_ss.source_sink_in.ip_group[i]
        }
        with open(f'{wdir}/inside_source_sink/{json_file}', 'w') as f:
            json.dump(ele2nwm_fid, f, indent=4)
    # generate sources.json for relocated sources
    relocate_sources2(
        old_ss_dir=f'{wdir}/inside_source_sink/',
        feeder_info_file='/sciclone/schism10/Hgrid_projects/STOFS3D-v8/v31/Feeder/feeder_heads_bases.xy',
        hgrid_fname=f'{wdir}/hgrid.gr3', outdir=f'{wdir}/relocated_source_sink/',
        max_search_radius=2100, mandatory_sources_coor=v19p2_for_sms_v27_mandatory_sources_coor,
        allow_neglection=True, region_list=region_list
    )
    # generate relocated source/sink
    os.chdir(f'{wdir}/relocated_source_sink/')
    os.system(f'ln -sf {wdir}/inside_source_sink/sinks.json .')
    os.system(f'ln -sf {wdir}/original_source_sink/hgrid.gr3 .')
    gen_sourcesink_nwm(  # with existing sources.json and sinks.json
        startdate=datetime(2018, 9, 1), rnday=1,
        cache_folder=f'{wdir}/original_source_sink/20180901/'
    )
    relocated_ss = source_sink.from_files(
        source_dir=f'{wdir}/relocated_source_sink/',
    )  # sinks will be discarded later, only sources will be used
    # discard sinks from relocated_ss
    relocated_ss = source_sink(
        vsource=relocated_ss.vsource, vsink=None, msource=relocated_ss.msource
    )

    # combine source/sinks inside and outside the high-res regions
    combined_ss = outside_ss + relocated_ss
    combined_ss.writer(f'{wdir}/combined_source_sink/')

    print('Done!')


if __name__ == '__main__':
    sample2()
    print('Done')
