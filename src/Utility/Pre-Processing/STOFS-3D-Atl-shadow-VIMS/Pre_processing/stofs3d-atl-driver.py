'''
This is the driver script to run the STOFS3D-ATL model.
'''

# Import modules
import os
import sys
import shutil
import subprocess
# from pylib_essentials.schism_file import schism_grid
from spp_essentials.Hgrid_extended import read_schism_hgrid_cached as schism_grid
import os

def prep_and_change_to_subdir(path, target_file):
    try_mkdir(path)
    os.chdir(path)
    try_remove(target_file)

def try_mkdir(folder_path):
    try:
        os.makedirs(folder_path)
        print("Folder created successfully:", folder_path)
    except FileExistsError:
        print("Folder already exists:", folder_path)
    except OSError as e:
        print("Error creating folder:", folder_path)
        print("Error message:", str(e))
    
def try_remove(path):
    try:
        if os.path.isfile(path):
            os.remove(path)
            # print("File removed successfully:", path)
        elif os.path.isdir(path):
            shutil.rmtree(path)
            # print("Folder removed successfully:", path)
        else:
            pass
            # print("Path does not exist:", path)
    except OSError as e:
        print("Error removing path:", path)
        print("Error message:", str(e))

driver_print_prefix = '-----------------STOFS3D-ATL driver:---------------------\n'


# define script path; scripts are available from schism's git repo
script_path = '/sciclone/data10/feiye/SCHISM_REPOSITORY/schism/src/Utility/Pre-Processing/STOFS-3D-Atl-shadow-VIMS/Pre_processing/'

# define the path where the model inputs are generated
model_input_path = '/sciclone/schism10/feiye/STOFS3D-v6/Inputs/I99/'
# make the directory if it does not exist
try_mkdir(model_input_path)

# hgrid must be prepared before running this script
hgrid = schism_grid(f'{model_input_path}hgrid.gr3')

# set a dict to indicate input files to be generated
input_files = {
    'vgrid': False,
    'gr3': True,
}

# -----------------Vgrid---------------------
if input_files['vgrid']:
    print(f'{driver_print_prefix}Generating vgrid.in ...')
    prep_and_change_to_subdir(f'{model_input_path}/Vgrid', 'vgrid.in')
    try_remove('hgrid.gr3')
    os.symlink(f'{model_input_path}hgrid.gr3', 'hgrid.gr3')
    # call a fortran program to generate vgrid.in
    # compile the fortran program if necessary
    # fortran_script_path = '/path/to/your/fortran_script.f90'
    # fortran_command = ['gfortran', fortran_script_path, '-o', 'fortran_script']
    fortran_process = subprocess.Popen(f'{script_path}/Vgrid/gen_vqs', stdin=subprocess.PIPE)  # , stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    # the command line argument 0 means no outputs of a sample vgrid along a given transect
    fortran_process.communicate(input=str(0).encode())

    print(f'{driver_print_prefix}converting the format of vgrid.in ...')
    # rename vgrid.in to vgrid.in.old
    try_remove('vgrid.in.old')
    os.rename('vgrid.in', 'vgrid.in.old')
    # convert the format of the vgrid.in file using a fortran script
    subprocess.run([f'{script_path}/Vgrid/change_vgrid'], check=True)  # , stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
    # rename vgrid.in.new to vgrid.in
    try_remove('vgrid.in')
    os.rename('vgrid.in.new', 'vgrid.in')
    os.chdir(model_input_path)

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

# -----------------spatially varying Gr3---------------------
# -----------------shapiro.gr3---------------------
if input_files['shapiro']:
    sub_dir = 'Shapiro'
    print(f'{driver_print_prefix}Generating shapiro.gr3 ...')
    prep_and_change_to_subdir(f'{model_input_path}/{sub_dir}', 'shapiro.gr3')

    os.chdir(model_input_path)




