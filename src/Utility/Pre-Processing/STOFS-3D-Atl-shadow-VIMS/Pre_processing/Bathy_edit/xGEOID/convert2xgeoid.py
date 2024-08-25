#!/usr/bin/env python
import os
import subprocess
from pathlib import Path
from glob import glob

from tqdm import tqdm
import numpy as np
import pandas as pd

from pylib import read_schism_hgrid, inside_polygon, read_schism_bpfile, grd2sms


def prep_folder(wdir):
    """
    Prepare the working directory for vdatum conversion
    """
    os.chdir(wdir)

    # remove any previous folders
    os.system(f'rm -rf {wdir}/vdatum/region*')

    # make a clean copy of vdatum
    os.makedirs(f'{wdir}/vdatum', exist_ok=True)
    os.makedirs(f'{wdir}/vdatum/result', exist_ok=True)

    # this is a clean copy of vdatum, all files under this folder are symlinks
    source_path = '/sciclone/schism10/Hgrid_projects/DEMs/vdatum/vdatum/'

    os.chdir(f'{wdir}/vdatum')
    os.system(f'ln -sf {source_path}/* .')
    os.chdir(wdir)

    # copy over the polygons
    source_path = ('/sciclone/schism10/hjyoo/task/task6_SECOFS/'
                   'simulation/Whole_Domain/Grid/Script/xGEOID/VDatum_polygons/')
    os.system(f'cp -rL {source_path} .')


def generate_input_txt(hgrid_obj, wdir, n_sub=500000):
    '''
    This function reproduces the manual steps for generating input files for vdatum.jar

    Generate input files for vdatum.jar by
    breaking the nodes into multiple files/regions to speed up conversion
    The nodes are divided into multiple files, each with n_sub nodes.
    If there are any nodes in the Chesapeake Bay and Delaware Bay,
    they are put into a separate file, because a different region is needed
    when calling vdatum.
    '''
    generated_input_files = []

    # ------- get nodes inside ches_del_bay polgyon -------------------------------
    bp = read_schism_bpfile(f'{wdir}/chea_del_bay.reg', fmt=1)
    idx = inside_polygon(np.c_[hgrid_obj.x, hgrid_obj.y], bp.x, bp.y).astype(bool)
    node_id = np.where(idx)[0] + 1  # add 1 to the index to make it 1-based

    lon = hgrid_obj.x[idx]
    lat = hgrid_obj.y[idx]
    depth = hgrid_obj.dp[idx]
    if len(node_id) > 0:
        generated_input_files.append('hgrid_stofs3d_inland_ches_del.txt')
        print('file 0, chesbay and delaware bay')
        with open(f"{wdir}/hgrid_stofs3d_inland_ches_del.txt", "w", encoding='utf-8') as f:
            # reverse the depth sign for "height" in vdatum
            f.write("\n".join(" ".join(map(str, line)) for line in zip(node_id, lon, lat, -depth)))
    else:
        print('No nodes inside chesbay and delaware bay')

    # ------------- get nodes inside stofs3d polygon but not in ches_del_bay polygon --------
    bp2 = read_schism_bpfile(f'{wdir}/stofs3d_inland2.reg', fmt=1)
    idx2 = inside_polygon(np.c_[hgrid_obj.x, hgrid_obj.y], bp2.x, bp2.y).astype(bool)
    idx = np.logical_and(~idx, idx2)
    node_id = np.where(idx)[0] + 1  # add 1 to the index to make it 1-based

    lon = hgrid_obj.x[idx]
    lat = hgrid_obj.y[idx]
    depth = hgrid_obj.dp[idx]

    # seperate into multiple files to speed up conversion
    nfile = int(len(lon)/n_sub) + 1
    nfile = np.ceil(len(lon)/n_sub).astype(int)
    print(f"nfile: {nfile}")
    for i in np.arange(nfile):
        x1 = i * n_sub
        x2 = min(x1 + n_sub, len(lon))
        print(f'file {i+1}, x1: {x1}, x2: {x2}')
        generated_input_files.append(f'hgrid_stofs3d_inland_{i+1}.txt')
        with open(f"{wdir}/hgrid_stofs3d_inland_{i+1}.txt", "w", encoding='utf-8') as f:
            # reverse the depth sign for "height" in vdatum
            f.write("\n".join(" ".join(map(str, line))
                    for line in zip(node_id[x1:x2], lon[x1:x2], lat[x1:x2], -depth[x1:x2])))

    if len(generated_input_files) == 0:
        raise ValueError('No nodes found in the domain')

    return generated_input_files


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
        with open(input_fname, 'r', encoding='utf-8') as file:
            lines = file.readlines()
        with open(failed_file, 'w', encoding='utf-8') as file:
            for index in failed_indices:
                file.write(lines[index])

        return failed_file
    else:
        return None


def replace_depth(wdir, hgrid_obj):
    '''
    Replace the original grid depths with the converted depths,
    which are saved in multiple result files (as outputs of vdatum.jar)

    '''
    depth_navd = hgrid_obj.dp.copy()
    destination_datum = 'xGEOID20b'

    files = glob(f'{wdir}/vdatum/result/*.txt')
    files.sort()

    depth_diff = np.zeros_like(hgrid_obj.dp)
    depth_diff[:] = np.NaN

    for fname in files:
        print(fname)

        with open(fname, encoding='utf-8') as f:
            df = pd.read_csv(
                f, header=None, delim_whitespace=True, names=['id', 'lon', 'lat', 'depth'], na_values=-999999.0)

        idxs = df.index[~df['depth'].isnull()]
        depth_geoid = df['depth'][idxs]
        # node id starts from 1, subtract 1 for dp's index
        hgrid_obj.dp[df['id'][idxs.values]-1] = -depth_geoid

    hgrid_obj.write_hgrid(f'{wdir}/hgrid_{destination_datum}.gr3')

    depth_diff = hgrid_obj.dp - depth_navd
    hgrid_obj.write_hgrid(f'{wdir}/hgrid_{destination_datum}_dp_minus_NAVD_dp.gr3')

    return hgrid_obj, depth_diff


def point_conversion(x, y, z, print_info=''):
    """
    Wrapper function for vdatum.jar to convert a single point from NAVD88 to xGEOID20b
    """

    # sample_command = "lib\java_home\openjdk-11.0.2\bin\java -jar vdatum.jar ihorz:NAD27:geo:deg "
    #                  "ivert:DTL:US_ft:height ohorz:NAD83_2011:geo:deg overt:NAVD88:m:height "
    #                  "-deg2dms -pt:-97.30965,26.3897528,3.545 region:4"

    vdatum_folder = '/sciclone/schism10/Hgrid_projects/DEMs/vdatum/vdatum/'
    z_convention = 'height'  # "sounding": positive downwards; "height": positive upwards
    z_converted = z.copy()
    regions = [4, 5]

    for i in tqdm(range(len(x)), desc=f"{print_info} Processing"):
        success = False
        for region in regions:
            result = subprocess.run(
                f"java -jar {vdatum_folder}/vdatum.jar ihorz:NAD83_2011 "
                f"ivert:navd88:m:{z_convention} ohorz:igs14:geo:deg overt:xgeoid20b:m:{z_convention} "
                f"-pt:{x[i]},{y[i]},{z[i]} region:{region}",
                shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE, check=True)
            if result.returncode == 0:
                z_converted[i] = float(result.stdout.decode().split()[417])
                print(f"{print_info}Point {i+1} ({x[i]}, {y[i]}, {z[i]}) successfully converted to {z_converted[i]}")
                success = True
                continue  # found the correct region, skip the rest

        if not success:
            print(f"{print_info}Point {i+1} ({x[i]}, {y[i]}, {z[i]}) failed to convert")
            # print("swapping default regions")
            # regions = [regions[1], regions[0]]

        # set toloreance to 1e-4, i.e., 0.1 mm
        # if not np.isclose(z[i], float(result.stdout.decode().split()[416]), atol=1e-4):
            # raise Exception(f"Input z and output z do not match for point {i+1} ({x[i]}, {y[i]}, {z[i]})")

    return z_converted


def convert2xgeoid(wdir, hgrid_obj, diag_output=None):
    '''
    Wrapper function for vdatum.jar to convert the hgrid depth to xGEOID20b
    Review the workflow and test it manually before using the script,
    because the script may not work for all cases.

    return:
    hgrid_obj: updated hgrid object with the converted depths
    depth_diff: converted depth - original depth
    '''

    print("Converting to xGEOID20b ... \n")
    z_convention = 'height'  # "sounding": positive downwards; "height": positive upwards

    prep_folder(wdir=wdir)

    # this generates the input files for vdatum.jar, including
    # hgrid_stofs3d_inland_?.txt and hgrid_stofs3d_inland_ches_del.txt
    generated_input_files = generate_input_txt(hgrid_obj=hgrid_obj, wdir=wdir, n_sub=100000)

    # see if the input files are complete

    vdatum_folder = f"{wdir}/vdatum/"
    # move the input files to the vdatum folder, since vdautm.jar only reads files in the current folder
    os.system(f"mv {wdir}/*.txt {vdatum_folder}")

    # clear the result folder
    os.makedirs(f'{vdatum_folder}/result', exist_ok=True)
    os.system(f"rm -rf {vdatum_folder}/result/*")

    # vdatum.jar needs to be in the same folder,
    # and it treates all upper case letters as lower case in the input file name
    os.chdir(vdatum_folder)

    # the first group should have no failed files, since they are strictly in region 4
    # this may not be true for other domains, so manually go over the workflow first before using the script
    input_fnames = glob("*_[0-9].txt")
    # Starting the processes
    processes = [subprocess.Popen(
            f"java -jar vdatum.jar  ihorz:NAD83_2011 ivert:navd88:m:{z_convention} "
            f"ohorz:igs14:geo:deg overt:xgeoid20b:m:{z_convention} "
            f"-file:txt:space,1,2,3,skip0:{fname}:result region:4",
            shell=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL
    ) for fname in input_fnames]

    if 'hgrid_stofs3d_inland_ches_del.txt' in generated_input_files:
        # This group should have failed files, since they are not strictly in region 4 or 5
        # This may be due to slight misalignment of the polygons provided by vdatum
        fname = f'{vdatum_folder}/hgrid_stofs3d_inland_ches_del.txt'
        input_fnames.append(fname)
        processes.append(subprocess.Popen(
            f"java -jar vdatum.jar ihorz:NAD83_2011 ivert:navd88:m:{z_convention} "
            f"ohorz:igs14:geo:deg overt:xgeoid20b:m:{z_convention} "
            f"-file:txt:space,1,2,3,skip0:{Path(fname).name}:result region:5",
            shell=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL))

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
        failed_file = file_check(input_fname, f'{vdatum_folder}/result/{Path(input_fname).name}')
        if failed_file is not None:
            print(f"Failed to convert {failed_file}, sending the failed portion to another region")
            retry_input = f"{Path(failed_file).stem}_retry.txt"
            os.system(f"cp {failed_file} {retry_input}")
            retry_inputs.append(retry_input)

    # Retry the failed files in Region 4
    for retry_input in retry_inputs:
        process = subprocess.Popen(
            f"java -jar vdatum.jar  ihorz:NAD83_2011 ivert:navd88:m:{z_convention} "
            f"ohorz:igs14:geo:deg overt:xgeoid20b:m:{z_convention} "
            f"-file:txt:space,1,2,3,skip0:{retry_input}:result region:4",
            shell=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
        exit_code = process.wait()
        if exit_code == 0:
            print(f"{retry_input} completed successfully")
        else:
            raise ValueError(f"{retry_input} failed with exit code {exit_code}")

    hgrid_obj, depth_diff = replace_depth(wdir=wdir, hgrid_obj=hgrid_obj)

    # diagnostic output
    if diag_output is not None:
        if Path(diag_output).suffix == '.2dm':
            grd2sms(hgrid_obj, diag_output)
        elif Path(diag_output).suffix == '.gr3':
            hgrid_obj.write_hgrid(diag_output, fmt=1)
        else:
            print(f"Unsupported output format: {Path(diag_output).suffix}\n"
                  "Assuming .gr3 format")
            hgrid_obj.write_hgrid(diag_output, fmt=1)

    return hgrid_obj, depth_diff


if __name__ == "__main__":
    # sample usage
    WORKING_DIR = '/sciclone/schism10/Hgrid_projects/TMP/DEM_edit/xGEOID/'
    xyz = np.loadtxt('/sciclone/schism10/Hgrid_projects/TMP/'
                     'DEM_edit/xGEOID/vdatum/region4_failed/hgrid_secofs_nccoast11.txt')
    xyz = xyz[:, 1:]
    point_conversion(xyz[:, 0], xyz[:, 1], xyz[:, 2])

    # sample usage
    hg = read_schism_hgrid(f'{WORKING_DIR}/hgrid.gr3')
    hg = convert2xgeoid(wdir=WORKING_DIR, hgrid_obj=hg)

    print('Done!')
