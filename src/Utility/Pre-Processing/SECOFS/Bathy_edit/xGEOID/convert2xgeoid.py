#!/usr/bin/env python
import os
import subprocess
from pathlib import Path
from glob import glob

import re
import numpy as np
import pandas as pd

from pylib import read_schism_hgrid, inside_polygon

def extract_trailing_numbers(s):
    match = re.search(r'\d+$', s)
    return match.group() if match else None

def prep_folder(wdir):
    os.chdir(wdir)

    # remove any previous folders
    os.system(f'rm -rf {wdir}/vdatum/region*')

    # make a clean copy of vdatum
    os.makedirs(f'{wdir}/vdatum', exist_ok=True)
    os.makedirs(f'{wdir}/vdatum/result', exist_ok=True)
    source_path = '/sciclone/schism10/Hgrid_projects/DEMs/vdatum/vdatum/'  # this is a clean copy of vdatum, all files under this folder are symlinks
    os.chdir(f'{wdir}/vdatum')
    os.system(f'ln -sf {source_path}/* .')
    os.chdir(wdir)

    # copy over the polygons
    source_path = '/sciclone/schism10/hjyoo/task/task6_SECOFS/simulation/Whole_Domain/Grid/Script/xGEOID/VDatum_polygons/'
    os.system(f'cp -rL {source_path} .')

def generate_input_txt(hgrid_obj, wdir):
    grouping = {'al': 4, 'fl': 4, 'ga': 4, 'nc': 4, 'de': 5, 'md': 5, 'nj': 5, 'va': 5}
    polygon_files = glob(f'{wdir}/VDatum_polygons/*.bnd')

    unique_groups = np.unique(list(grouping.values()))
    for group in unique_groups:
        os.makedirs(f'{wdir}/vdatum/region{group}', exist_ok=True)
        os.system(f'rm -rf {wdir}/vdatum/region{group}/*')

    for fname in polygon_files:
        prefix = Path(fname).name[:2].lower()
        group_id = grouping[prefix]

        ##get nodes inside polgyon
        bp = np.loadtxt(fname)
        sind = inside_polygon(np.c_[hgrid_obj.x, hgrid_obj.y], bp[:, 0], bp[:, 1])

        idxs = np.where(sind == 1)[0]
        print(idxs.shape[0])
        if idxs.shape[0] == 0: continue
        lon = hgrid_obj.x[idxs]
        lat = hgrid_obj.y[idxs]
        depth = hgrid_obj.dp[idxs] #for sounding

        name = Path(fname).name.split('_')[0].lower()
        with open(f"{wdir}/vdatum/region{group_id}/hgrid_secofs_{name}.txt", "w") as f:
            f.write("\n".join(" ".join(map(str, line)) for line in zip(idxs+1, lon, lat, depth)))
        
def file_check(input_fname, result_fname):
    '''
    Check if the result file is complete such that
    it has all lines processed from the input file.
    If not, put the failed lines into a file.
    '''
    input_fname = Path(input_fname)
    result_fname = Path(result_fname)

    arr1 = np.loadtxt(input_fname)
    arr2 = np.loadtxt(result_fname)
    failed_indices = np.where(~np.isin(arr1[:, 0], arr2[:, 0]))[0]

    # put the failed lines into a file and try again
    if len(failed_indices) > 0:
        # save the failed portion of the original input file to a failed folder for later inspection
        failed_folder = f'{input_fname.parent}_failed'
        if not os.path.exists(failed_folder):
            os.mkdir(failed_folder)
        failed_file = Path(f'{failed_folder}/{input_fname.name}')

        # write the failed lines to a file
        with open(input_fname, 'r') as file:
            lines = file.readlines()
        with open(failed_file, 'w') as file:
            for index in failed_indices:
                file.write(lines[index])

        return failed_file
    else:
        return None

def replace_depth(wdir, hgrid_obj):
    depth_navd = hgrid_obj.dp.copy()
    destination_datum = 'xGEOID20b'

    files = glob(f'{wdir}/vdatum/result/*.txt')
    files.sort()

    depth_diff = np.zeros_like(hgrid_obj.dp)
    depth_diff[:] = np.NaN

    for fname in files:
        print(fname)

        with open(fname) as f:
            df = pd.read_csv(f, header=None, delim_whitespace=True, names=['id', 'lon', 'lat', 'depth'], na_values=-999999.0)

        idxs = df.index[~df['depth'].isnull()]
        depth_geoid = df['depth'][idxs]

        hgrid_obj.dp[df['id'][idxs.values]-1] = depth_geoid  # id is node id, which starts from 1, so need to subtract 1 for dp

    hgrid_obj.write_hgrid(f'{wdir}/hgrid_{destination_datum}.gr3')

    depth_diff = hgrid_obj.dp - depth_navd

    return hgrid_obj, depth_diff


def convert2xgeoid(wdir, hgrid_obj):
    prep_folder(wdir=wdir)
    generate_input_txt(hgrid_obj=hgrid_obj, wdir=wdir)

    vdatum_folder = f"{wdir}/vdatum/" 
    region_folders = glob(f"{vdatum_folder}/region*")  # outputs from generate_input_txt(), e.g., region4, region5, etc.
    regions = [int(extract_trailing_numbers(folder)) for folder in region_folders]
    regions = sorted(regions)  # must be sorted,
    # e.g., if region4 is being processed and an input file (*.txt) has the first half of the points in region4,
    # then the outputs will only contain the first half (valid points).
    # the rest (invalid points) will be used as input for the next region (region5).
    # However, if region5 is processed before region4, the outputs will be empty, i.e., neglecting the valid points in region4.

    # clear the result folder
    os.makedirs(f'{vdatum_folder}/result', exist_ok=True)
    os.system(f"rm -rf {vdatum_folder}/result/*")

    # vdatum.jar needs to be in the same folder, and it treates all upper case letters as lower case in the input file name
    os.chdir(vdatum_folder)
    for i, region_id in enumerate(regions): # calling vdatum.jar for each region
        input_fnames = glob(f"region{region_id}/*.txt")

        # Starting the processes
        processes = [subprocess.Popen(f"java -jar vdatum.jar  ihorz:NAD83_2011 ivert:navd88:m:sounding ohorz:igs14:geo:deg overt:xgeoid20b:m:sounding -file:txt:space,1,2,3,skip0:{fname}:result region:{region_id}", shell=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL) for fname in input_fnames]

        # Wait for all processes to complete
        for i, process in enumerate(processes):
            exit_code = process.wait()
            if exit_code == 0:
                print(f"Process {i+1} ({input_fnames[i]}) completed successfully")
            else:
                print(f"Process {i+1} ({input_fnames[i]}) failed with exit code {exit_code}")
        
        # Check if the output files are correct
        for input_fname in input_fnames:
            failed_file = file_check(input_fname, f'{vdatum_folder}/result/{Path(input_fname).name}')  # check if the output file is complete
            if failed_file is not None:
                if region_id == regions[-1]:
                    raise Exception(f"No other groups to try: Failed to convert {failed_file}")
                else:
                    print(f"Failed to convert {failed_file} in region{region_id}, sending the failed portion to the next region")
                    # check if the next region also has the failed file, and remove it if it does because it probably has
                    # an overlap with the failed file from the previous region. The overlap with previous region leads to no outputs.
                    fname = f"region{region_id+1}/{failed_file.name}"
                    if os.path.exists(fname):
                        os.system(f"rm -rf {fname}")
                    # copy the failed file (a portion of the original input file) to the next region's input folder,
                    # rename the file so that the previous results are not overwritten
                    os.system(f"cp {failed_file} region{region_id+1}/{failed_file.stem}_failed_in_region{region_id}.txt")  # copy the failed file (a portion of the original input file) to the next region's input folder


    hgrid_obj,depth_diff = replace_depth(wdir=wdir, hgrid_obj=hgrid_obj)
    #hgrid_obj_diff = hgrid_obj.copy() # for save diff between NAVD and xGEOID
    #hgrid_obj_diff.dp = depth_diff
    #hgrid_obj_diff.save(f'{wdir}/depth_NAVD-{destination_datum}.gr3')
    
    return hgrid_obj

if __name__ == "__main__":
    wdir = '/sciclone/schism10/Hgrid_projects/TMP/DEM_edit/xGEOID/'
    hgrid_obj = read_schism_hgrid(f'{wdir}/hgrid.gr3')
    hgrid_obj = convert2xgeoid(wdir=wdir, hgrid_obj=hgrid_obj)
    
    pass
