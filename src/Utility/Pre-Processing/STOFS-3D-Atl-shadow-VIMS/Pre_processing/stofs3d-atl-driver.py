'''
This is the driver script to run the STOFS3D-ATL model.
'''

# Import modules
import os
import numpy as np
from glob import glob
import subprocess
from scipy.spatial import cKDTree
import copy

# self-defined modules
from spp_essentials.Hgrid_extended import read_schism_hgrid_cached as schism_grid  # only for testing purposes
# support both the experimental and original versions of pylib; the former is imported by default
try:
    # from pylib_essentials.schism_file import schism_grid
    from pylib_essentials.schism_file import read_schism_reg
    from pylib_utils.utility_functions import inside_polygon
except:
    from pylib import inside_polygon, schism_grid, read_schism_reg
# modules to be moved to pylib
from schism_py_pre_post.Grid.SourceSinkIn import source_sink

# import from the sub folders, not the installed package
from Source_sink.Constant_sinks.set_constant_sink import set_constant_sink

# Global variables:
# define script path; scripts are available from schism's git repo
script_path = '/sciclone/data10/feiye/SCHISM_REPOSITORY/schism/src/Utility/Pre-Processing/STOFS-3D-Atl-shadow-VIMS/Pre_processing/'

# Classes and functions:
class Config_stofs3d_atlantic():
    '''A class to handle the configuration of STOFS-3D-ATL model,
    i.e., processing the parameters and storing the factory settings.
    '''
    def __init__(self,
        nudging_zone_width = 1.5,  # in degrees
        nudging_day = 1.0,  # in days
        shapiro_zone_width = 2.5,  # in degrees
        shapiro_tilt = 2.0,  # more abrupt transition in the shapiro zone
    ):
        # see a description of the parameters in the function make_river_map()
        self.nudging_zone_width = nudging_zone_width
        self.nudging_day = nudging_day
        self.shapiro_zone_width = shapiro_zone_width
        self.shapiro_tilt = shapiro_tilt

    @classmethod
    def v6(cls):
        return cls(
            nudging_zone_width = 7.3,  # very wide nudging zone
            shapiro_zone_width = 11.5,  # very wide shapiro zone
            shapiro_tilt = 3.5,  # very abrupt transition in the shapiro zone
        )

def prep_and_change_to_subdir(path, target_file):
    try_mkdir(path)
    os.chdir(path)

    if type(target_file) is str:
        target_file = [target_file]
    elif type(target_file) is list:
        if len(target_file) == 0:  # remove all if empty list
            target_file = glob(f'{path}/*')
    else:
        raise TypeError("target_file must be a string or a list of strings")

    for file in target_file:
        try_remove(file)

def try_mkdir(folder_path):
    try:
        os.makedirs(folder_path)
        print("Folder created successfully:", folder_path)
    except FileExistsError:
        print("Folder already exists:", folder_path)
    except OSError as e:
        print("Error creating folder:", folder_path)
        print("Error message:", str(e))
    
def try_remove(file):
    try:
        os.remove(file)
            # print("Folder removed successfully:", path)
    except FileNotFoundError:
        pass
    except OSError as e:
        print("Error removing path:", file)
        print("Error message:", str(e))

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

def main():
    # prepare a configuration
    config = Config_stofs3d_atlantic.v6()

    driver_print_prefix = '-----------------STOFS3D-ATL driver:---------------------\n'
    # define the path where the model inputs are generated
    model_input_path = '/sciclone/schism10/feiye/STOFS3D-v6/Inputs/I16/'
    # make the directory if it does not exist
    try_mkdir(model_input_path)

    # hgrid must be prepared before running this script
    hgrid = schism_grid(f'{model_input_path}hgrid.gr3')

    # set a dict to indicate input files to be generated
    input_files = {
        'vgrid': False,
        'gr3': False,
        'nudge_gr3': False,
        'shapiro': False,
        'drag': False,
        'source_sink': True,
    }

    # -----------------Vgrid---------------------
    if input_files['vgrid']:
        sub_dir = 'Vgrid'
        print(f'{driver_print_prefix}Generating vgrid.in ...')
        prep_and_change_to_subdir(f'{model_input_path}/{sub_dir}', [])
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
        try_mkdir(f'{model_input_path}/{sub_dir}')
        os.chdir(f'{model_input_path}/{sub_dir}')
        gr3_values = {
            'albedo': 0.1,
            'diffmax': 1.0,
            'diffmin': 1e-6,
            'watertype': 7.0,
            'windrot_geo2proj': 0.0
        }
        for name, value in gr3_values.items():
            print(f'{driver_print_prefix}Generating {name}.gr3 ...')
            try_remove(f'{model_input_path}/{sub_dir}/{name}.gr3')
            hgrid.write(f'{model_input_path}/{sub_dir}/{name}.gr3', value=value)

        os.chdir(model_input_path)

    # >>> spatially varying Gr3 >>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    # -----------------nudge.gr3---------------------
    if input_files['nudge_gr3']:
        sub_dir = 'Nudge_gr3'
        print(f'{driver_print_prefix}Generating nudge.gr3 ...')
        prep_and_change_to_subdir(f'{model_input_path}/{sub_dir}', 'nudge.gr3')

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
        prep_and_change_to_subdir(f'{model_input_path}/{sub_dir}', 'shapiro.gr3')

        # use the same method for nudge.gr3 to generate a buffer zone along the open boundaries
        nudge_coef = gen_nudge_coef(hgrid, rlmax=config.shapiro_zone_width)
        # using nudging coefficient as a proxy for the distribution of shapiro strength
        shapiro = gen_shapiro_strength(hgrid, nudge_coef, tilt=config.shapiro_tilt)

        hgrid.save(f'{model_input_path}/{sub_dir}/shapiro.gr3', value=shapiro)

        os.chdir(model_input_path)

    # -----------------drag.gr3---------------------
    if input_files['drag']:
        pass
    # <<< end spatially varying Gr3 <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    
    # -----------------source_sink---------------------
    if input_files['source_sink']:
        sub_dir = 'Source_sink'
        file_list = ['source_sink.in', 'vsource.th', 'msource.th', 'vsink.th']
        print(f'{driver_print_prefix}Generating source_sink.in ...')
        prep_and_change_to_subdir(f'{model_input_path}/{sub_dir}', file_list)

        # generate source_sink files by intersecting NWM river segments with the model land boundary
        # prep_and_change_to_subdir(f'{model_input_path}/{sub_dir}/original_source_sink/', [])
        pass
        
        # relocate source locations to resolved river channels
        # prep_and_change_to_subdir(f'{model_input_path}/{sub_dir}/relocated_source_sink/', [])
        # os.symlink(f'{model_input_path}/hgrid.gr3', 'hgrid.gr3')
        # subprocess.run(
        #     [f'{script_path}/Source_sink/Relocate/relocate_source_feeder.py'],
        #     cwd=f'{model_input_path}/{sub_dir}/relocated_source_sink/'
        # )
        relocated_ss = source_sink(source_dir=f'{model_input_path}/{sub_dir}/relocated_source_sink/')

        # set constant sinks (pumps and background sinks)
        prep_and_change_to_subdir(f'{model_input_path}/{sub_dir}/constant_sink/', [])
        # copy *.shp to the current directory
        os.system(f'cp {script_path}/Source_sink/Constant_sinks/levee_4_pump_polys.* .')

        hgrid_utm = copy.deepcopy(hgrid)
        hgrid_utm.proj(prj0='epsg:4326', prj1='epsg:26918')  # needed because the levee shapefile is in 26918 (lon/lat would lead to problems due to truncation errors)
        background_ss = set_constant_sink(wdir=f'{model_input_path}/{sub_dir}/constant_sink/', hgrid_utm=hgrid_utm)

        total_ss = relocated_ss + background_ss
        total_ss.writer(f'{model_input_path}/{sub_dir}/')
        
        os.chdir(model_input_path)

if __name__ == '__main__':
    main()