"""
Utility functions for various tasks
"""
import os
import shutil
import shapefile
import numpy as np
import subprocess
from pylib import inside_polygon


# List of states in the Atlantic STOFS3D domain
STOFS3D_ATL_STATES = [
    'TX', 'LA', 'MS', 'AL', 'FL', 'GA', 'SC', 'NC', 'VA',
    'DC', 'MD', 'DE', 'PA', 'NJ', 'NY', 'CT', 'RI', 'MA', 'NH', 'ME'
]


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


