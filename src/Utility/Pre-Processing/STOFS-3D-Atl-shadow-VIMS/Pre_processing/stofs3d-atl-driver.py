'''
This is the driver script to run the STOFS3D-ATL model.
'''


# Import modules
import os
import numpy as np
from glob import glob
import subprocess
from scipy.spatial import cKDTree
import geopandas as gpd
import copy
from datetime import datetime
from pyproj import Transformer
from pathlib import Path
import shutil
import shapefile
import hgrid_pybind

# self-defined modules
from pylib_essentials.schism_file import read_schism_reg, source_sink
from pylib_essentials.schism_file import schism_grid, read_schism_hgrid_cached, cread_schism_hgrid
from pylib_essentials.utility_functions import inside_polygon
from pyschism.mesh import Hgrid

# import from the sub folders, not from installed packages
from Source_sink.NWM.gen_sourcesink import gen_sourcesink
from Source_sink.Constant_sinks.set_constant_sink import set_constant_sink
from Source_sink.Relocate.relocate_source_feeder import relocate_sources, v19p2_mandatory_sources_coor
from Bctides.bctides.bctides import Bctides  # temporary, bctides.py will be merged into pyschism
# from pyschism.forcing.bctides import Bctides

# Global variables:
# use the full path of the Pre-Processing dir inside your schism repo
# e.g, script_path = '/my_dir/schism/src/Utility/Pre-Processing/'
# If you are running this script in the schism repo, you can also use the following line:
script_path = Path('./').absolute()

# Classes and functions:
class Config_stofs3d_atlantic():
    '''A class to handle the configuration of STOFS-3D-ATL model,
    i.e., processing the parameters and storing the factory settings.
    '''
    def __init__(self,
        startdate = datetime(2017, 12, 1),  # start date of the model
        rnday = 60,  # number of days to run the model
        nudging_zone_width = 1.5,  # in degrees
        nudging_day = 1.0,  # in days
        shapiro_zone_width = 2.5,  # in degrees
        shapiro_tilt = 2.0,  # more abrupt transition in the shapiro zone
        nwm_cache_folder = None,
        feeder_info_file = None,  # the file that contains the feeder info, made by make_feeder_channel.py in RiverMapper
        gr3_values = {  # uniform gr3 values
            'albedo': 0.1,
            'diffmax': 1.0,
            'diffmin': 1e-6,
            'watertype': 1.0,
            'windrot_geo2proj': 0.0
        }
    ):
        self.startdate = startdate
        self.rnday = rnday
        self.nudging_zone_width = nudging_zone_width
        self.nudging_day = nudging_day
        self.shapiro_zone_width = shapiro_zone_width
        self.shapiro_tilt = shapiro_tilt
        self.nwm_cache_folder = nwm_cache_folder
        self.feeder_info_file = feeder_info_file
        self.gr3_values = gr3_values

    @classmethod
    def v6(cls):
        return cls(
            nudging_zone_width = 7.3,  # very wide nudging zone
            shapiro_zone_width = 11.5,  # very wide shapiro zone
            shapiro_tilt = 3.5,  # very abrupt transition in the shapiro zone
            feeder_info_file = f'/sciclone/schism10/feiye/STOFS3D-v5/Inputs/v14/Parallel/SMS_proj/feeder/feeder.pkl'
        )

    @classmethod
    def v7(cls):
        return cls(
            nudging_zone_width = 7.3,  # default nudging zone
            shapiro_zone_width = 11.5,  # default shapiro zone
            shapiro_tilt = 3.5,  # default abrupt transition in the shapiro zone
            feeder_info_file = f'/sciclone/schism10/Hgrid_projects/STOFS3D-v7/v20.0/Feeder/feeder_heads_bases.xy',
            nwm_cache_folder = Path('/sciclone/schism10/whuang07/schism20/NWM_v2.1/')
        )

def mkcd_new_dir(path):
    remove_folder(path)
    os.makedirs(path)
    os.chdir(path)

def remove_folder(path):
    import shutil
    try:
        shutil.rmtree(path)
        print(f"{path} and all contents removed successfully")
    except FileNotFoundError:
        print(f"{path} does not exist, and that's okay.")
    except Exception as error:
        print(f"Error removing {path}: {error}")

def try_remove(file):
    try:
        os.remove(file)
    except FileNotFoundError:
        pass
    except OSError as e:
        print("Error removing path:", file)
        print("Error message:", str(e))

def find_points_in_polyshp(pt_xy, shapefile_names):
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

def gen_nudge_coef(hgrid:schism_grid, rlmax = 1.5, rnu_day=0.25, open_bnd_list=[0, 1]):
        """
        Set up nudge zone within rlmax distance from the ocean boundary:
        - schism_grid (from pylib) must be in lon/lat coordinates
        - rlmax can be a uniform value, e.g., rl_max = 1.5,
          or a 2D array of the same size as the hgrid's number of nodes,
          e.g., rl_max = np.zeros(NP) with further tuning of the nudging zone width.
        - If there are more than one ocean boundary, specify the boundary index in open_bnd_list as a list
        """

        # maximum nudging strength
        rnu_max = 1.0 / rnu_day / 86400.0

        nudge_coeff = np.zeros((hgrid.np, ), dtype=float)
        for open_bnd_idx in open_bnd_list:
            print(f'boundary {open_bnd_idx}: {len(hgrid.iobn[open_bnd_idx])} open boundary nodes')
            bnd_node_idx = hgrid.iobn[open_bnd_idx]

            # distance between each node to its nearest open boundary node
            distance, _ = cKDTree(np.c_[hgrid.x[bnd_node_idx], hgrid.y[bnd_node_idx]]).query(np.c_[hgrid.x, hgrid.y])

            # nudge coefficient based on the current open boundary
            this_nudge_coeff = np.clip((1.0-distance/rlmax)*rnu_max, 0.0, rnu_max)
            nudge_coeff = np.maximum(this_nudge_coeff, nudge_coeff)

        return nudge_coeff

def gen_drag(hgrid:schism_grid):
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
    drag[idx] = np.interp(hgrid.dp[idx], grid_depths, drag_coef, left=drag_coef[0], right=drag_coef[-1])

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
        'Mississippi_0': {
            'drag': 1e-8,
            'region_file': f'{script_path}/Gr3/Drag/Mississippi_0.reg'  # all segments of the Mississippi River
        },
        'Mississippi_downstream_0': {
            'drag': 0,
            'region_file': f'{script_path}/Gr3/Drag/Mississippi_downstream_0.reg'  # only the downstream part of the Mississippi River
        },
        'Eastport_0': {
            'drag': 0.0,
            'region_file': f'{script_path}/Gr3/Drag/Eastport_0.reg'
        }
    }
    for keys, tweak in region_tweaks.items():
        print(f"Applying drag {tweak['drag']} in {tweak['region_file']}")
        reg = read_schism_reg(tweak['region_file'])
        idx = inside_polygon(np.c_[hgrid.x, hgrid.y], reg.x, reg.y).astype(bool)
        river = hgrid.dp > 0
        drag[idx & river] = tweak['drag']

    return drag

def gen_shapiro_strength(hgrid:schism_grid, init_shapiro_dist:np.ndarray=None, tilt=2.5):

    shapiro_min=100; shapiro_max=1000

    # ---------- dependency on the ocean nudging zone ------------
    if init_shapiro_dist is None:
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
    region_files = [
        f'{script_path}/Gr3/Shapiro/coastal_0.2.lonlat.reg',
        f'{script_path}/Gr3/Shapiro/coastal_0.5_1.lonlat.reg',
        f'{script_path}/Gr3/Shapiro/coastal_0.5_2.lonlat.reg'
    ]  # the order matters
    region_tweaks = [200, 300, 300]
    for region_tweak, region_file in zip(region_tweaks, region_files):
        reg = read_schism_reg(region_file)
        idx = inside_polygon(np.c_[hgrid.x, hgrid.y], reg.x, reg.y).astype(bool)
        shapiro_region[idx] = region_tweak

    # -------- combine all shapiro tweaks ---------
    shapiro = shapiro_region + shapiro_ocean_bnd
    shapiro = np.clip(shapiro, shapiro_min, shapiro_max)

    return shapiro

def gen_elev_ic(hgrid=None, h0=0.1, city_shape_fnames=None):
    '''
    set initial elevation: 0 in the ocean, h0 below ground on higher grounds and in cities
    city_shape_fname: the shapefile must be in the same projection as the hgrid
    '''
    if city_shape_fnames is None:
        in_city = np.ones(hgrid.dp.shape, dtype=bool)
    else:
        in_city = find_points_in_polyshp(pt_xy=np.c_[hgrid.x, hgrid.y], shapefile_names=city_shape_fnames)

    high_ground = hgrid.dp < 0

    elev_ic = np.zeros(hgrid.dp.shape, dtype=float)

    ind = np.logical_or(high_ground, in_city)

    # set initial elevation to h0 below ground on higher grounds and in cities
    elev_ic[ind] = - hgrid.dp[ind] - h0

    return elev_ic

def prep_run_dir(parent_dir, runid):
    model_input_path = f'{parent_dir}/Inputs/I{runid}'
    os.makedirs(model_input_path, exist_ok=True)
    os.makedirs(f'/sciclone/scr10/feiye/R{runid}/outputs', exist_ok=True)

    rundir = f'{parent_dir}/Runs/R{runid}'
    output_dir = f'/sciclone/scr10/feiye/R{runid}/outputs'

    os.makedirs(rundir, exist_ok=True)
    os.chdir(rundir)
    os.system(f'ln -sf {output_dir} .')

    os.makedirs(f'{parent_dir}/Outputs/O{runid}', exist_ok=True)
    os.chdir(f'{parent_dir}/Outputs/O{runid}')
    os.system(f'ln -sf {output_dir} .')
    os.system(f'ln -sf {model_input_path}/hgrid.gr3 .')
    os.system(f'ln -sf {model_input_path}/vgrid.in .')

    return model_input_path

def main(is_test=False):
    # prepare a configuration
    config = Config_stofs3d_atlantic.v7()
    config.nwm_cache_folder = None
    config.rnday = 37
    config.startdate = datetime(2024, 3, 5)

    driver_print_prefix = '-----------------STOFS3D-ATL driver:---------------------\n'
    # define the path where the model inputs are generated
    project_dir = '/sciclone/schism10/feiye/STOFS3D-v7/'
    runid = '15f'

    # define and make the model_input_path, the run_dir and the output dir
    model_input_path = prep_run_dir(project_dir, runid)
    
    # make a copy of the script itself to the model_input_path
    os.system(f'cp {script_path}/stofs3d-atl-driver.py {model_input_path}')

    # hgrid must be prepared before running this script
    # temporary tests
    if is_test:
        hgrid = cread_schism_hgrid(hgrid_pybind.HGrid(f'{model_input_path}/hgrid.gr3'))
        test(config, model_input_path, hgrid)
        return

    hgrid = cread_schism_hgrid(f'{model_input_path}/hgrid.gr3')
    hgrid_pyschism = Hgrid.open(f'{model_input_path}/hgrid.gr3', crs='epsg:4326')

    # set a dict to indicate input files to be generated
    input_files = {
        'bctides': False,
        'vgrid': False,
        'gr3': False,
        'nudge_gr3': False,
        'shapiro': False,
        'drag': True,
        'elev_ic': False,
        'source_sink': False,
        # 'hotstart.nc``
        # '*D.th.nc'
        # '*nu.nc'
        # '*.prop'
    }

    # -----------------bctides---------------------
    if input_files['bctides']:
        sub_dir = 'Bctides'
        print(f'{driver_print_prefix}Generating bctides.in ...')
        mkcd_new_dir(f'{model_input_path}/{sub_dir}')
        Bctides(
            hgrid=hgrid_pyschism,  # needs hgrid in pyschism's Hgrid class
            flags=[[5,3,0,0],[3,3,0,0],[0,1,0,0]],
            database='fes2014'
        ).write(
            output_directory=f'{model_input_path}/{sub_dir}',
            start_date=config.startdate,
            rnday=config.rnday,
        )

    # -----------------Vgrid---------------------
    if input_files['vgrid']:
        sub_dir = 'Vgrid'
        print(f'{driver_print_prefix}Generating vgrid.in ...')
        mkcd_new_dir(f'{model_input_path}/{sub_dir}')
        os.system(f'ln -sf ../hgrid.gr3 .')

        # call a fortran program to generate vgrid.in
        # compile the fortran program if necessary
        # fortran_script_path = '/path/to/your/fortran_script.f90'
        # fortran_command = ['gfortran', fortran_script_path, '-o', 'fortran_script']
        fortran_process = subprocess.Popen(f'{script_path}/Vgrid/gen_vqs', stdin=subprocess.PIPE)  # , stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        # the command line argument 0 means no outputs of a sample vgrid along a given transect
        fortran_process.communicate(input=str(0).encode())

        print(f'{driver_print_prefix}converting the format of vgrid.in ...')
        # rename vgrid.in to vgrid.in.old
        try_remove(f'{model_input_path}/{sub_dir}/vgrid.in.old')
        os.rename('vgrid.in', 'vgrid.in.old')

        # convert the format of the vgrid.in file using a fortran script
        subprocess.run([f'{script_path}/Vgrid/change_vgrid'], check=True)  # , stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
        # rename vgrid.in.new to vgrid.in
        try_remove(f'{model_input_path}/{sub_dir}/vgrid.in')
        os.rename('vgrid.in.new', 'vgrid.in')

        os.chdir(model_input_path)
        os.system(f'ln -sf {sub_dir}/vgrid.in .')

    # -----------------spatially uniform Gr3---------------------
    if input_files['gr3']:
        sub_dir = 'Gr3'
        print(f'{driver_print_prefix}Generating spatially uniform gr3 files ...')
        mkcd_new_dir(f'{model_input_path}/{sub_dir}')

        for name, value in config.gr3_values.items():
            print(f'{driver_print_prefix}Generating {name}.gr3 ...')
            try_remove(f'{model_input_path}/{sub_dir}/{name}.gr3')
            hgrid.write(f'{model_input_path}/{sub_dir}/{name}.gr3', value=value)

        os.chdir(model_input_path)

    # >>> spatially varying Gr3 >>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    # -----------------nudge.gr3---------------------
    if input_files['nudge_gr3']:
        sub_dir = 'Nudge_gr3'
        print(f'{driver_print_prefix}Generating nudge.gr3 ...')
        mkcd_new_dir(f'{model_input_path}/{sub_dir}')

        # generate nudging coefficient based on the proximity to the open boundaries
        nudge_coef = gen_nudge_coef(hgrid, rlmax=config.nudging_zone_width, rnu_day=config.nudging_day)

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
        print(f'{driver_print_prefix}Generating shapiro.gr3 ...')
        mkcd_new_dir(f'{model_input_path}/{sub_dir}')

        # use the same method for nudge.gr3 to generate a buffer zone along the open boundaries
        nudge_coef = gen_nudge_coef(hgrid, rlmax=config.shapiro_zone_width)
        # using nudging coefficient as a proxy for the distribution of shapiro strength
        shapiro = gen_shapiro_strength(hgrid, nudge_coef, tilt=config.shapiro_tilt)

        hgrid.save(f'{model_input_path}/{sub_dir}/shapiro.gr3', value=shapiro)

        os.chdir(model_input_path)

    # -----------------drag.gr3---------------------
    if input_files['drag']:
        sub_dir = 'Drag'
        print(f'{driver_print_prefix}Generating drag.gr3 ...')
        mkcd_new_dir(f'{model_input_path}/{sub_dir}')
        drag = gen_drag(hgrid)

        hgrid.save(f'{model_input_path}/{sub_dir}/drag.gr3', value=drag)

        os.chdir(model_input_path)

    # -----------------elev_ic---------------------
    if input_files['elev_ic']:
        sub_dir = 'Elev_ic'
        print(f'{driver_print_prefix}Generating elev_ic.gr3 ...')
        mkcd_new_dir(f'{model_input_path}/{sub_dir}')
        elev_ic = gen_elev_ic(
            hgrid, h0=0.1,
            city_shape_fnames=[f'{script_path}/Hotstart/LA_urban_polys_lonlat.shp']
        )

        hgrid.save(f'{model_input_path}/{sub_dir}/elev_ic.gr3', value=elev_ic)
        os.symlink(f'{model_input_path}/{sub_dir}/elev_ic.gr3', f'{model_input_path}/{sub_dir}/elev.ic')

        os.chdir(model_input_path)

    # <<< end spatially varying Gr3 <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

    # -----------------source_sink---------------------
    if input_files['source_sink']:
        sub_dir = 'Source_sink'
        # file_list = ['source_sink.in', 'vsource.th', 'msource.th', 'vsink.th']
        print(f'{driver_print_prefix}Generating source_sink.in ...')
        mkcd_new_dir(f'{model_input_path}/{sub_dir}')

        # generate source_sink files by intersecting NWM river segments with the model land boundary
        mkcd_new_dir(f'{model_input_path}/{sub_dir}/original_source_sink/')
        os.symlink(f'{model_input_path}/hgrid.gr3', 'hgrid.gr3')
        gen_sourcesink(startdate=config.startdate, rnday=config.rnday, cache_folder=config.nwm_cache_folder)

        # relocate source locations to resolved river channels
        mkcd_new_dir(f'{model_input_path}/{sub_dir}/relocated_source_sink/')
        os.symlink(f'{model_input_path}/hgrid.gr3', 'hgrid.gr3')
        relocated_ss = relocate_sources(
            old_ss_dir=f'{model_input_path}/{sub_dir}/original_source_sink/',
            feeder_info_file=config.feeder_info_file,
            hgrid_fname=f'{model_input_path}/hgrid.gr3',
            outdir=f'{model_input_path}/{sub_dir}/relocated_source_sink/',
            max_search_radius=2100, mandatory_sources_coor=v19p2_mandatory_sources_coor, relocate_map=None,
            allow_neglection=False
        )

        # set constant sinks (pumps and background sinks)
        mkcd_new_dir(f'{model_input_path}/{sub_dir}/constant_sink/')
        # copy *.shp to the current directory
        os.system(f'cp {script_path}/Source_sink/Constant_sinks/levee_4_pump_polys.* .')
        hgrid_utm = copy.deepcopy(hgrid)
        hgrid_utm.proj(prj0='epsg:4326', prj1='epsg:26918')
        background_ss = set_constant_sink(wdir=f'{model_input_path}/{sub_dir}/constant_sink/', hgrid_utm=hgrid_utm)

        # assemble source/sink files and write to model_input_path
        total_ss = relocated_ss + background_ss
        total_ss.writer(f'{model_input_path}/{sub_dir}/')

        # write diagnostic outputs
        hgrid.compute_ctr()
        np.savetxt(
            f'{model_input_path}/{sub_dir}/vsource.xyz',
            np.c_[hgrid.xctr[total_ss.source_eles-1], hgrid.yctr[total_ss.source_eles-1], total_ss.vsource.df.mean().values]
        )
        np.savetxt(
            f'{model_input_path}/{sub_dir}/vsink.xyz',
            np.c_[hgrid.xctr[total_ss.sink_eles-1], hgrid.yctr[total_ss.sink_eles-1], total_ss.vsink.df.mean().values]
        )

        os.chdir(model_input_path)

def test(config, model_input_path, hgrid):
    sub_dir = 'Source_sink'
    relocated_ss = relocate_sources(
        old_ss_dir=f'{model_input_path}/{sub_dir}/original_source_sink/',
        feeder_info_file=config.feeder_info_file,
        hgrid_fname=f'{model_input_path}/hgrid.gr3',
        outdir=f'{model_input_path}/{sub_dir}/relocated_source_sink/',
        max_search_radius=2100, mandatory_sources_coor=v19p2_mandatory_sources_coor, relocate_map=None,
        allow_neglection=False
    )
    # set constant sinks (pumps and background sinks)
    mkcd_new_dir(f'{model_input_path}/{sub_dir}/constant_sink/')
    # copy *.shp to the current directory
    os.system(f'cp {script_path}/Source_sink/Constant_sinks/levee_4_pump_polys.* .')
    hgrid_utm = copy.deepcopy(hgrid)
    hgrid_utm.proj(prj0='epsg:4326', prj1='epsg:26918')
    background_ss = set_constant_sink(wdir=f'{model_input_path}/{sub_dir}/constant_sink/', hgrid_utm=hgrid_utm)

    # assemble source/sink files and write to model_input_path
    total_ss = relocated_ss + background_ss
    total_ss.writer(f'{model_input_path}/{sub_dir}/')

    # write diagnostic outputs
    hgrid.compute_ctr()
    np.savetxt(
        f'{model_input_path}/{sub_dir}/vsource.xyz',
        np.c_[hgrid.xctr[total_ss.source_eles-1], hgrid.yctr[total_ss.source_eles-1], total_ss.vsource.df.mean().values]
    )
    np.savetxt(
        f'{model_input_path}/{sub_dir}/vsink.xyz',
        np.c_[hgrid.xctr[total_ss.sink_eles-1], hgrid.yctr[total_ss.sink_eles-1], total_ss.vsink.df.mean().values]
    )

    os.chdir(model_input_path)

if __name__ == '__main__':
    main(is_test=False)
