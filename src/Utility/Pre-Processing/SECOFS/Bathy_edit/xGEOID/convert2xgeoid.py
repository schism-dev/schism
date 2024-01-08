import numpy as np
import subprocess
import os
from glob import glob
from pathlib import Path
import time

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

if __name__ == "__main__":
    # ---------- inputs ----------
    vdatum_folder = "/sciclone/schism10/feiye/SECOFS/Inputs/I03a/Bathy_edit/xGEOID/vdatum/"
    regions = [4, 5]  # must be in ascending order
    # --------------------

    # clear the result folder
    os.system(f"rm -rf {vdatum_folder}/result/*")

    # vdatum.jar needs to be in the same folder, and it treates all upper case letters as lower case in the input file name
    os.chdir(vdatum_folder)
    for i, region_id in enumerate(regions):
        input_fnames = glob(f"region{region_id}/*.txt")

        # Starting the processes
        processes = [subprocess.Popen(f"java -jar vdatum.jar  ihorz:NAD83_2011 ivert:navd88:m:height ohorz:igs14:geo:deg overt:xgeoid20b:m:height -file:txt:space,1,2,3,skip0:{fname}:result region:{region_id}", shell=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL) for fname in input_fnames]

        # Wait for all processes to complete
        for i, process in enumerate(processes):
            exit_code = process.wait()
            if exit_code == 0:
                print(f"Process {i+1} completed successfully")
            else:
                print(f"Process {i+1} failed with exit code {exit_code}")
        
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


        pass

