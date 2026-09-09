#!/usr/bin/env python3

'''
Automatically submit and monitor the run status.
Restart the run from the latest hotstart point if the run stops.

(1) Copy this script, run_test and run_comb into rundir

(2) Inside rundir: prep inputs as before
    The "hotstart.nc" (if it exists under ihot=1 or 2) will be overwritten by this script by symlinks ("ln -sf");
    if you want to keep it, "mv hotstart.nc hotstart.nc.0" then "ln -s hotstart.nc.0 hotstart.nc"

(3) Under the run dir,
      "python auto_hotstart.py >& scrn.out"
    or
      "./auto_hotstart.py > scrn.out &"
      "tail -f scrn.out"
    The latter prints jobs status and run speeds to the screen continuously.

    This can be done at any stage of the run, for example:
        Before launching a run:
            - The script can submit the initial job and monitor its progress.
        After an interruption:
            - This is often used to continue a run in the same folder, for example, after the run stops due to a time
            limit or a manual scancel.
        After manually submitting a run (i.e., while the run is pending or ongoing):
            - If resuming a run copied or moved from another location, it is preferable to launch it manually first
            to ensure it starts properly, unless the full output files are transferred.
            Sometimes, only staout* and flux.out are available for resuming a run.

    Note: the script will use the last part of the current dir as runid, e.g., the runid of "RUN13a" or "R13a" will be "13a"
          , make sure you don't have duplicate run ids.

(4) Sometimes a run can hang, e.g., on Hercules, the script will automatically detect this and interrupt the run then relaunch from
    the nearest hotstarting point

(5) Start by specifying inputs below, between the lines "inputs" and "end inputs"
'''

import os
import time
import subprocess
import glob
from copy import copy
import re
from datetime import datetime
from pathlib import Path

#----------input---------------------------------
#(1) Copy this script, run_test and run_comb into rundir
#(2) Inside rundir: prep inputs as before
#    The "hotstart.nc" (if it exists under ihot=1 or 2) will be overwritten by this script by symlinks;
#    if you want to keep it, "mv hotstart.nc hotstart.nc.0" then "ln -s hotstart.nc.0 hotstart.nc"
#(3) Under the run dir, "python auto_hotstart.py >& scrn.out"
#Script will use the last part of current dir as runid
#This script can also be launched inside an on-going run


# ----------------inputs---------------------
JOB_SCHEDULER = 'slurm'  # 'pbs' or 'slurm'

rundir = os.getcwd()  # use os.getcwd if launch from rundir; otherwise specify a string yourself
last_stack = None  # if None, the script will try to find the last stack number in the file "param.nml"
                   # make sure the run can finish the specified rnday in param.nml (i.e., the forcing covers the whole period);
                   # otherwise, change the "rnday" in param.nml or specify another number here
# ----------------end inputs---------------------


# ---------------------  embedded functions -----------------------
def decor_print(msg, prefix='', suffix=''):
    '''
    Decorative print with prefix and suffix
    '''
    print(f'{prefix}{msg}{suffix}', flush=True)


def Replace_string_in_file(fname, str_orig_pattern, str_replace):
    '''
    Match and replace a string in a file
    '''
    if '~' in fname:
        fname = fname.replace('~', os.path.expanduser('~'))
    fname = os.path.abspath(fname)
    with open(fname, "rt", encoding='utf-8') as fin:
        with open("tmp.txt", "wt", encoding='utf-8') as fout:
            for line in fin:
                fout.write(re.sub(str_orig_pattern, str_replace, line))
    os.system(f"mv tmp.txt {fname}")


def ReplaceJobName(fname, job_name, job_scheduler='slurm'):
    import fileinput

    if job_scheduler == 'slurm':
        pattern = r"(SBATCH\s+-J\s+)(\S+)"
    else:
        raise ValueError('job_scheduler must be either "slurm"; pbs not implemented yet')

    replacement = rf'SBATCH -J {job_name}'

    # Use fileinput to edit the file in place
    match_found = False
    with fileinput.FileInput(fname, inplace=True) as file:
        for line in file:
            match = re.search(pattern, line)
            if match:
                match_found = True
            modified_line = re.sub(pattern, replacement, line)
            print(modified_line, end='')

    if not match_found:
        raise ValueError(f'Job name specification not found in {fname}')


def Get_hotstart_step(run_out_dir):
    if '~' in run_out_dir:
        run_out_dir = run_out_dir.replace('~', os.path.expanduser('~'))
    run_out_dir = os.path.abspath(run_out_dir)

    hot_files = glob.glob(f"{run_out_dir}/hotstart_000000_*.nc")

    hot_steps = []
    for hot_file in hot_files:
        sub_str = hot_file.split("_")
        sub_str = sub_str[-1].split(".")
        hot_steps.append(int(sub_str[0]))

    hot_steps.sort()
    return hot_steps


def Get_var_from_file(fname, var_name, reverse_search=False):
    if '~' in fname:
        fname = fname.replace('~', os.path.expanduser('~'))
    fname = os.path.abspath(fname)

    pattern = fr'{var_name}\s*=\s*(\d+)'
    with open(fname, "rt") as fin:
        lines = fin.readlines()  # Read all lines into a list

    if reverse_search:
        lines = [line for line in reversed(lines)]

    for line in lines:
        match = re.search(pattern, line)
        if match:
            var_value = match.group(1)  # Extract the matched group
            return var_value

    print(f'Variable {var_name} not found in {fname}')
    return None


def save_interruption_diagnostics(outputs_dir):
    '''
    Save diagnostic files from an interrupted run using the next unused index.
    '''
    outputs_dir = Path(outputs_dir).expanduser().resolve()
    diagnostic_files = [
        'fatal.error',
        'mirror.out',
    ]

    existing_indices = []
    for diagnostic_file in diagnostic_files:
        pattern = re.compile(rf'^{re.escape(diagnostic_file)}\.(\d+)$')
        for saved_file in outputs_dir.glob(f'{diagnostic_file}.*'):
            match = pattern.match(saved_file.name)
            if match:
                existing_indices.append(int(match.group(1)))

    diagnostic_index = max(existing_indices, default=0) + 1
    saved_files = []
    for diagnostic_file in diagnostic_files:
        source_file = outputs_dir / diagnostic_file
        if not source_file.exists():
            decor_print(f'{source_file} does not exist; skipping diagnostic save.')
            continue

        target_file = outputs_dir / f'{diagnostic_file}.{diagnostic_index}'
        if target_file.exists():
            raise FileExistsError(f'Refusing to overwrite existing diagnostic file: {target_file}')

        target_file.write_bytes(source_file.read_bytes())
        saved_files.append(str(target_file))

    if saved_files:
        decor_print(f'Saved interruption diagnostics: {", ".join(saved_files)}')
    else:
        decor_print('No interruption diagnostics were saved.')
# --------------------- end embedded functions -----------------------


def auto_hotstart(job_scheduler='slurm', rundir=None, last_stack=None):
    '''
    Main function to submit and monitor the run status
    '''
    my_print_prefix = ''
    my_print_suffix = ''

    # ---------------------  get background information -----------------------
    # options not exposed to users
    run_id_normal = ''
    rundir_normal = ''  # f'/scratch1/02786/{user_name}/RUN{run_id_normal}/'
    # end options not exposed to users

    user_name = os.getlogin()  # getlogin() or specify a string yourself

    if job_scheduler == 'slurm':
        batch_cmd = 'sbatch'
        queue_query_str = f"squeue -u {user_name}"
    elif job_scheduler == 'pbs':
        batch_cmd = 'qsub'
        queue_query_str = f"qstat -u {user_name}"
    else:
        raise ValueError('job_scheduler must be either "slurm" or "pbs"')

    run_id = os.path.basename(rundir)
    if len(run_id) > 8:
        raise ValueError('run_id must be 8 characters or less, otherwise it may be truncated by the job scheduler')
    decor_print(f'RUN job name : {run_id}')


    combine_job_name = copy(run_id)
    if combine_job_name.startswith("RUN"):
        combine_job_name = combine_job_name[3:]  # Remove "RUN"
    elif combine_job_name.startswith("R"):
        combine_job_name = combine_job_name[1:]  # Remove "R"
    combine_job_name = f'C{combine_job_name}'

    if len(combine_job_name) > 8:
        raise ValueError(f'combine_job_name ({combine_job_name}) must be 8 characters or less, otherwise it may be truncated by the job scheduler')
    decor_print(f'Combine job name : {combine_job_name}')

    # get the last stack number in the file "param.nml"
    # Define the regular expression pattern to match "rnday =" followed by a number
    if last_stack is None:
        rnday = float(Get_var_from_file(f'{rundir}/param.nml', 'rnday'))
        dt = float(Get_var_from_file(f'{rundir}/param.nml', 'dt'))
        ihfskip = int(Get_var_from_file(f'{rundir}/param.nml', 'ihfskip'))

        last_stack = int(rnday * 86400 / dt / ihfskip)

    # ----------------- prepare job names in batch script --------------------
    ReplaceJobName(f'{rundir}/run_test', run_id, job_scheduler)
    ReplaceJobName(f'{rundir}/run_comb', combine_job_name, job_scheduler)

    os.chdir(f'{rundir}')

    # --------------------- monitor the run -----------------------
    previous_time_step = -1  # initialize
    while (not os.path.exists(f'{rundir}/outputs/schout_000000_{last_stack+1}.nc')) and (not os.path.exists(f'{rundir}/outputs/out2d_{last_stack}.nc')):
        print(f'\n{"%" * 80}\n  Local time: {datetime.now()}\n{"%" * 80}\n')  # step header

        job_status = subprocess.getoutput(queue_query_str)
        print(job_status)

        # extract the line corresponding to the current run_id
        current_job_status = None
        pattern = rf'\b{re.escape(run_id)}\b'  # Ensures exact match with word boundaries
        for line in job_status.splitlines():
            if re.search(pattern, line):  # Search for an exact match
                current_job_status = line.strip()
                print(f'\n{"-"*100}\ncurrent job status: {current_job_status}\n{"-"*100}\n')
                break


        if rundir_normal != '':
            if re.search(rf"{run_id_normal}\s+{user_name}\s+R", current_job_status) is not None:
                decor_print(f'{run_id_normal} running')
            elif re.search(rf"{run_id_normal}\s+{user_name}\s+PD", current_job_status) is not None:
                decor_print(f'{run_id_normal} pending')
            else:
                decor_print(f'{run_id_normal} does not exist')

        if current_job_status is not None:
            job_id = re.search(r'\b\d+\b', current_job_status).group()
            decor_print(f'job_id: {job_id}')

            # check run status
            if re.search(rf"{run_id}\s+{user_name}\s+R", current_job_status) is not None:  # job id in squeue

                # make sure mirror.out is written after the job starts;
                # if it is not written after 60s, then the run probably hangs
                time.sleep(60)

                # see if the run hangs
                hanged = False
                # check if mirror.out is written
                if not os.path.exists(f'{rundir}/outputs/mirror.out'):
                    hanged = True
                else:
                    time_step = Get_var_from_file(f'{rundir}/outputs/mirror.out', 'TIME STEP', reverse_search=True)
                    if time_step is not None:
                        time_step = int(time_step)
                        if time_step == previous_time_step:
                            decor_print(f'{run_id} hangs at Time Step')
                            hanged = True
                        else:
                            decor_print(f'Run advancing, time step: {time_step}')
                            previous_time_step = time_step
                if hanged:
                    decor_print(f'{run_id}: killing job and wait for automatic relaunching')
                    print(f'scancel {job_id}')
                    os.system(f'scancel {job_id}')
                    # time.sleep(20)  # wait for scancel
                    # print(f'{batch_cmd} run_test')
                    # os.system(f'{batch_cmd} run_test')
                    # decor_print(f'{run_id} hanged and was resubmitted ...')
                else:
                    decor_print(f'{run_id} running, wait ...')
            else:
                previous_time_step = -1  # resetting for unforeseen cases
                decor_print(f'{run_id} queueing, wait ...')
            time.sleep(120)
        else:  # i.e., current_job_status is None (job id not in squeue)
            # check if the run finishes normally
            # open a file and check if the last line contains "Run completed successfully"
            if os.path.exists(f'{rundir}/outputs/mirror.out'):
                with open(f'{rundir}/outputs/mirror.out', 'r') as file:
                    lines = file.readlines()  # Read all lines in the file
                    last_line = lines[-1] if lines else ''  # Get the last line if the file is not empty
                    # Check if the last line contains the desired phrase
                    if "Run completed successfully" in last_line:
                        decor_print("The last line indicates that the run completed successfully.")
                        break
                    else:
                        decor_print("The last line does not indicate a successful completion, try combining the last hotstart.nc then restart the run.")
                        save_interruption_diagnostics(f'{rundir}/outputs')
                # combine hotstart
                hot_steps = Get_hotstart_step(f'{rundir}/outputs/')
                if len(hot_steps) == 0:
                    raise Exception('No hotstart files generated before run stopped.')
                hot_step = hot_steps[-1]

                hot_combined = f'{rundir}/outputs/hotstart_it={hot_step}.nc'
                decor_print(f'{run_id} stopped, last hotstart to combine: {hot_combined}')

                os.chdir(f'{rundir}/outputs/')
                Replace_string_in_file(
                    f'{rundir}/run_comb',
                    r'~/bin/combine_hotstart7 -i\s+\d+',
                    f'~/bin/combine_hotstart7 -i {hot_step}'
                )

                decor_print(f'{batch_cmd} run_comb')
                os.system(f'{batch_cmd} {rundir}/run_comb')  # run_comb under rundir, but submitted inside outputs

                # restore run_cmb template
                # os.system(f'cat {rundir}/run_comb')
                # Replace_string_in_file(f'{rundir}/run_comb', f'-i {hot_step}', '-i 0000')

                # wait for combine hotstart to finish
                while combine_job_name in subprocess.getoutput(queue_query_str):
                    time.sleep(20)
                    decor_print(f'waiting for job {combine_job_name} to finish')
                if not os.path.exists(hot_combined):
                    raise Exception(f'Failed generating {hot_combined}')

                time.sleep(20)

                # link combined hotstart.nc
                os.chdir(f'{rundir}')
                decor_print(f'linking {hot_combined}')
                os.system(f'rm {rundir}/hotstart.nc')
                os.symlink(hot_combined, f'{rundir}/hotstart.nc')

                # dealing with the paired normal run
                if rundir_normal != '':
                    if re.search(rf"{run_id_normal}\s+{user_name}\s+R", current_job_status) is not None:
                        decor_print(f'{run_id_normal} has started, skipping syncing files with it')
                    else:
                        decor_print(f'copying staout* and *.out to {rundir_normal}')
                        os.system(f'cp {rundir}/outputs/staout* {rundir_normal}/outputs/')
                        os.system(f'cp {rundir}/outputs/*.out {rundir_normal}/outputs/')

                Replace_string_in_file(f'{rundir}/param.nml', r'ihot\s*=\s*\d+', 'ihot = 2')
                time.sleep(20)
            # -- end if  mirror.out exists

            # submit new job
            os.chdir(f'{rundir}')
            decor_print(f'{batch_cmd} run_test')
            os.system(f'{batch_cmd} run_test')

        # restore template run_test
            # Replace_string_in_file('~/bin/run_test', f'RUN{run_id}', 'RUNxxx')

    decor_print('Autohotstart Task Completed.')


if __name__ == '__main__':
    auto_hotstart(job_scheduler=JOB_SCHEDULER, rundir=rundir, last_stack=last_stack)
