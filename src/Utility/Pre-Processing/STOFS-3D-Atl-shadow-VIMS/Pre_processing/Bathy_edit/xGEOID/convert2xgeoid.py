#!/usr/bin/env python
import os
import subprocess
from pathlib import Path
from glob import glob

import re
import numpy as np
import pandas as pd

from pylib import read_schism_hgrid, inside_polygon, read_schism_bpfile

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
    ##get nodes inside ches_del_bay polgyon
    bp = read_schism_bpfile(f'{wdir}/chea_del_bay.reg', fmt=1)
    idx = inside_polygon(np.c_[hgrid_obj.x, hgrid_obj.y], bp.x, bp.y).astype(bool)
    id = np.where(idx)[0] + 1  # add 1 to the index to make it 1-based

    lon = hgrid_obj.x[idx]
    lat = hgrid_obj.y[idx]
    depth = hgrid_obj.dp[idx]

    print(f'file 0, chesbay and delaware bay')
    with open(f"{wdir}/hgrid_stofs3d_inland_ches_del.txt", "w") as f:
        f.write("\n".join(" ".join(map(str, line)) for line in zip(id, lon, lat, -depth)))  # reverse the depth sign for "height" in vdatum

    # get nodes inside stofs3d polygon but not in ches_del_bay polygon
    bp2 = read_schism_bpfile(f'{wdir}/stofs3d_inland2.reg', fmt=1)
    idx2 = inside_polygon(np.c_[hgrid_obj.x, hgrid_obj.y], bp2.x, bp2.y).astype(bool)
    idx = np.logical_and(~idx, idx2)
    id = np.where(idx)[0] + 1  # add 1 to the index to make it 1-based

    lon = hgrid_obj.x[idx]
    lat = hgrid_obj.y[idx]
    depth = hgrid_obj.dp[idx]
  
    #seperate into multiple files to speed up conversion
    n_sub = 500000
    nfile = int(len(lon)/n_sub) + 1
    nfile = np.ceil(len(lon)/n_sub).astype(int)
    print(nfile)
    for i in np.arange(nfile):
        x1 = i * n_sub
        x2 = min(x1 + n_sub, len(lon))
        print(f'file {i+1}, x1: {x1}, x2: {x2}')
        with open(f"{wdir}/hgrid_stofs3d_inland_{i+1}.txt", "w") as f:
            f.write("\n".join(" ".join(map(str, line)) for line in zip(id[x1:x2], lon[x1:x2], lat[x1:x2], -depth[x1:x2])))   # reverse the depth sign for "height" in vdatum 
        
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

        hgrid_obj.dp[df['id'][idxs.values]-1] = -depth_geoid  # id is node id, which starts from 1, so need to subtract 1 for dp's index

    hgrid_obj.write_hgrid(f'{wdir}/hgrid_{destination_datum}.gr3')

    depth_diff = hgrid_obj.dp - depth_navd

    return hgrid_obj, depth_diff


def point_conversion(x, y, z, print_info=''):
    from tqdm import tqdm
    # sample:
    # lib\java_home\openjdk-11.0.2\bin\java -jar vdatum.jar ihorz:NAD27:geo:deg ivert:DTL:US_ft:height ohorz:NAD83_2011:geo:deg overt:NAVD88:m:height -deg2dms -pt:-97.30965,26.3897528,3.545 region:4
    vdatum_folder = '/sciclone/schism10/Hgrid_projects/DEMs/vdatum/vdatum/'
    z_convention = 'height'  # "sounding": positive downwards; "height": positive upwards
    z_converted = z.copy()
    regions = [4, 5]

    for i in tqdm(range(len(x)), desc=f"{print_info} Processing"):
        success = False
        for region in regions:
            result = subprocess.run(f"java -jar {vdatum_folder}/vdatum.jar ihorz:NAD83_2011 ivert:navd88:m:{z_convention} ohorz:igs14:geo:deg overt:xgeoid20b:m:{z_convention} -pt:{x[i]},{y[i]},{z[i]} region:{region}", shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
            if result.returncode == 0:
                z_converted[i] = float(result.stdout.decode().split()[417])
                print(f"{print_info}Point {i+1} ({x[i]}, {y[i]}, {z[i]}) successfully converted to {z_converted[i]}")
                success = True
                continue  # found the correct region, skip the rest
        
        if success == False:
            print(f"{print_info}Point {i+1} ({x[i]}, {y[i]}, {z[i]}) failed to convert")
            # print("swapping default regions")
            # regions = [regions[1], regions[0]]

        # if not np.isclose(z[i], float(result.stdout.decode().split()[416]), atol=1e-4):  # set toloreance to 1e-4, i.e., 0.1 mm
            # raise Exception(f"Input z and output z do not match for point {i+1} ({x[i]}, {y[i]}, {z[i]})")
    
    return z_converted


def convert2xgeoid(wdir, hgrid_obj):
    print("Converting to xGEOID20b ... \n")
    z_convention = 'height'  # "sounding": positive downwards; "height": positive upwards

    prep_folder(wdir=wdir)
    generate_input_txt(hgrid_obj=hgrid_obj, wdir=wdir)

    vdatum_folder = f"{wdir}/vdatum/" 
    os.system(f"mv {wdir}/*.txt {vdatum_folder}")  # move the input files to the vdatum folder, since vdautm.jar only reads files in the current folder

    # clear the result folder
    os.makedirs(f'{vdatum_folder}/result', exist_ok=True)
    os.system(f"rm -rf {vdatum_folder}/result/*")

    # vdatum.jar needs to be in the same folder, and it treates all upper case letters as lower case in the input file name
    os.chdir(vdatum_folder)

    # the first group should have no failed files, since they are strictly in region 4
    # this may not be true for other domains, so manually go over the workflow first before using the script
    input_fnames = glob("*_[0-9].txt")
    # Starting the processes
    processes = [subprocess.Popen(f"java -jar vdatum.jar  ihorz:NAD83_2011 ivert:navd88:m:{z_convention} ohorz:igs14:geo:deg overt:xgeoid20b:m:{z_convention} -file:txt:space,1,2,3,skip0:{fname}:result region:4", shell=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL) for fname in input_fnames]

    # This group should have failed files, since they are not strictly in region 4 or 5
    fname = f'{vdatum_folder}/hgrid_stofs3d_inland_ches_del.txt'
    input_fnames.append(fname)
    processes.append(subprocess.Popen(f"java -jar vdatum.jar  ihorz:NAD83_2011 ivert:navd88:m:{z_convention} ohorz:igs14:geo:deg overt:xgeoid20b:m:{z_convention} -file:txt:space,1,2,3,skip0:{Path(fname).name}:result region:5", shell=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL))

    # Wait for all processes to complete
    for i, process in enumerate(processes):
        exit_code = process.wait()
        if exit_code == 0:
            print(f"Process {i+1} ({input_fnames[i]}) completed successfully")
        else:
            print(f"Process {i+1} ({input_fnames[i]}) failed with exit code {exit_code}")

    # Check if the output files are complete
    retry_inputs = []
    for input_fname in input_fnames:
        failed_file = file_check(input_fname, f'{vdatum_folder}/result/{Path(input_fname).name}')  # check if the output file is complete
        if failed_file is not None:
            print(f"Failed to convert {failed_file}, sending the failed portion to another region")
            retry_input = f"{Path(failed_file).stem}_retry.txt"
            os.system(f"cp {failed_file} {retry_input}")
            retry_inputs.append(retry_input)
    
    # Retry the failed files in other regions
    for retry_input in retry_inputs:
        process = subprocess.Popen(f"java -jar vdatum.jar  ihorz:NAD83_2011 ivert:navd88:m:{z_convention} ohorz:igs14:geo:deg overt:xgeoid20b:m:{z_convention} -file:txt:space,1,2,3,skip0:{retry_input}:result region:4", shell=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
        exit_code = process.wait()
        if exit_code == 0:
            print(f"{retry_input} completed successfully")
        else:
            raise Exception(f"{retry_input} failed with exit code {exit_code}")

    hgrid_obj,depth_diff = replace_depth(wdir=wdir, hgrid_obj=hgrid_obj)
    
    return hgrid_obj

if __name__ == "__main__":
    wdir = '/sciclone/schism10/Hgrid_projects/TMP/DEM_edit/xGEOID/'
    xyz = np.loadtxt('/sciclone/schism10/Hgrid_projects/TMP/DEM_edit/xGEOID/vdatum/region4_failed/hgrid_secofs_nccoast11.txt')
    xyz = xyz[:, 1:]
    point_conversion(xyz[:, 0], xyz[:, 1], xyz[:, 2])
    
    hgrid_obj = read_schism_hgrid(f'{wdir}/hgrid.gr3')
    hgrid_obj = convert2xgeoid(wdir=wdir, hgrid_obj=hgrid_obj)
    
    pass
