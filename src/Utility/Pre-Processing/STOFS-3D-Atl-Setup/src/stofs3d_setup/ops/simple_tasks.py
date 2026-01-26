
"""
Simple tasks with few steps and dependencies.
When a tasks becomes more complex, it should be moved to a separate module.
"""


# ------------------------- Import modules ---------------------------
import os
import time
from pathlib import Path
import logging

import numpy as np
from scipy.spatial import cKDTree
import netCDF4


# self-defined modules
import pylib
from pylib import inside_polygon, read_schism_reg
from pyschism.mesh import Hgrid as pyschism_Hgrid

from pyschism.forcing.hycom.hycom2schism import OpenBoundaryInventory as OpenBoundaryInventoryHYCOM
from ..ops.AVISO.aviso2schism import OpenBoundaryInventory as OpenBoundaryInventoryAVISO
from pyschism.forcing.hycom.hycom2schism import Nudge

# Import from the sub folders. These are not from installed packages.
from ..utils.utils import find_points_in_polyshp

# Global variables:
# use the full path of the Pre-Processing dir inside your schism repo
# e.g, script_path = '/my_dir/schism/src/Utility/Pre-Processing/'
# If you are running this script in the schism repo, you can also use the following line:
script_path = Path(__file__).resolve().parent.parent / "ops"
print(f"script_path: {script_path}")


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
    ocean_bnd_ids=None, max_attempts=50,
    hycom_download_dir=None
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

    hgrid = pyschism_Hgrid.open(hgrid_fname, crs='epsg:4326')
    bnd = OpenBoundaryInventoryHYCOM(hgrid, vgrid_fname)

    restart = False  # First attempt doesn't need "restart"
    last_error = None  # To store the last exception
    for attempt in range(10000000):  # effectively infinite attempts
        if attempt > 10:
            time.sleep(1200)  # wait for 20 minutes before retrying
        elif attempt > max_attempts:  # let user decide to continue or not
            input("Press Enter to continue retrying...")

        try:
            bnd.fetch_data(
                outdir, start_date, rnday,
                elev2D=True, TS=True, UV=True,
                ocean_bnd_ids=ocean_bnd_ids,
                restart=restart,
                hycom_download_dir=hycom_download_dir
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

    hgrid = pyschism_Hgrid.open(hgrid_fname, crs='epsg:4326')
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

    hgrid = pyschism_Hgrid.open(hgrid_fname, crs='epsg:4326')

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
    raise Exception(f"Failed to fetch nudge coefficient after {max_retry} attempts due to repeated errors.")

    print('Done generating nudging netcdf files.')


def gen_soil(hgrid: pylib.schism_grid):
    """
    Prepare soil_conductivity.gr3, which specifies the soil thermal conductivity in W/m^2/K,
    and soil_thick.gr3 (in meters).
    """

    soil_conductivity = np.zeros_like(hgrid.dp)
    soil_thick = np.zeros_like(hgrid.dp)

    shallow = hgrid.dp < 20
    soil_conductivity[shallow] = 5  # W/m^2/K for shallow areas
    soil_thick[:] = 1.0

    return [soil_conductivity, soil_thick]


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
        elev_ic = pylib.schism_grid(base_elev_ic).dp
    else:
        elev_ic = np.zeros(hgrid.dp.shape, dtype=float)

    ind = np.logical_or(high_ground, in_city)
    hgrid.save('city_high_ground.gr3', value=ind.astype(int))  # diagnostic

    # set initial elevation to h0 below ground on higher grounds and in cities
    elev_ic[ind] = - hgrid.dp[ind] - h0

    return elev_ic
