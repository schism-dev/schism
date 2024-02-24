import os
import time
import subprocess
import glob
import re

#----------input---------------------------------
#(1) Copy this script, run_test and run_comb into rundir
#(2) Inside rundir: prep inputs as before
#    The "hotstart.nc" (if it exists under ihot=1 or 2) will be overwritten by this script by symlinks; 
#    if you want to keep it, "mv hotstart.nc hotstart.nc.0" then "ln -s hotstart.nc.0 hotstart.nc"
#(3) Under the run dir, "python auto_hotstart.py >& scrn.out"
#Script will use the last part of current dir as runid
#This script can also be launched inside an on-going run

job_scheduler = 'slurm' # 'pbs' or 'slurm'

rundir = os.getcwd()  # use os.getcwd if launch from rundir; otherwise specify a string yourself
last_stack = None  # if None, the script will try to find the last stack number in the file "param.nml"
                   # make sure the run can finish the specified rnday in param.nml (i.e., the forcing covers the whole period);
                   # otherwise, change the "rnday" in param.nml or specify another number here

#----------end input---------------------------------


# ---------------------  embedded functions -----------------------

def Replace_string_in_file(fname, str_orig, str_replace):
    if '~' in fname:
        fname = fname.replace('~', os.path.expanduser('~'))
    fname = os.path.abspath(fname)
    with open(fname, "rt") as fin:
        with open("tmp.txt", "wt") as fout:
            for line in fin:
                fout.write(line.replace(str_orig, str_replace))
    os.system(f"mv tmp.txt {fname}")

def ReplaceJobName(fname, job_name, job_scheduler='slurm'):
    import fileinput

    if job_scheduler == 'slurm':
        pattern = r"(SBATCH\s+-J\s+)(\S+)"
    else:
        raise Exception('job_scheduler must be either "slurm"; pbs not implemented yet')

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
        raise Exception(f'Job name specification not found in {fname}')

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

def Get_var_from_file(fname, var_name):
    if '~' in fname:
        fname = fname.replace('~', os.path.expanduser('~'))
    fname = os.path.abspath(fname)

    pattern = fr'{var_name}\s*=\s*(\d+)'
    with open(fname, "rt") as fin:
        for line in fin:
            match = re.search(pattern, line)
            if match:
                var_value = match[1]
                return var_value

    raise Exception(f'Variable {var_name} not found in {fname}')

# --------------------- end embedded functions -----------------------

my_print_decor = '\n%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%\n'

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
    raise Exception('job_scheduler must be either "slurm" or "pbs"')

run_id = os.path.basename(rundir)
if len(run_id) > 8:
    raise Exception('run_id must be 8 characters or less, otherwise it may be truncated by the job scheduler')
print(f'{my_print_decor}RUN job name : {run_id}{my_print_decor}', flush=True)

combine_job_name = f'{run_id}_cmb'
if "RUN" in combine_job_name:
    combine_job_name = combine_job_name.replace("RUN", "")  # shorten the name if it contains "RUN"
if len(combine_job_name) > 8:
    raise Exception('combine_job_name must be 8 characters or less, otherwise it may be truncated by the job scheduler')
print(f'{my_print_decor}Combine job name : {combine_job_name}{my_print_decor}', flush=True)

# get the last stack number in the file "param.nml"
# Define the regular expression pattern to match "rnday =" followed by a number
if last_stack is None:
    rnday = float(Get_var_from_file(f'{rundir}/param.nml', 'rnday'))
    dt = float(Get_var_from_file(f'{rundir}/param.nml', 'dt'))
    ihfskip = int(Get_var_from_file(f'{rundir}/param.nml', 'ihfskip'))

    last_stack = int(rnday * 86400 / dt / ihfskip)

# -------------------- Launch the init run --------------------
ReplaceJobName(f'{rundir}/run_test', run_id, job_scheduler)
ReplaceJobName(f'{rundir}/run_comb', combine_job_name, job_scheduler)

os.chdir(f'{rundir}')
# print(f'{my_print_decor}Launch first run_test{my_print_decor}', flush=True)
# os.system(f'{batch_cmd} run_test')

# --------------------- monitor the run -----------------------
while (not os.path.exists(f'{rundir}/outputs/schout_000000_{last_stack+1}.nc')) and (not os.path.exists(f'{rundir}/outputs/out2d_{last_stack}.nc')):

    job_status = subprocess.getoutput(queue_query_str)
    print(job_status)

    if rundir_normal != '':
        if re.search(rf"{run_id_normal}\s+{user_name}\s+R", job_status) is not None:
            print(f'{my_print_decor}{run_id_normal} running{my_print_decor}')
        elif re.search(rf"{run_id_normal}\s+{user_name}\s+PD", job_status) is not None:
            print(f'{my_print_decor}{run_id_normal} pending{my_print_decor}')
        else:
            print(f'{my_print_decor}{run_id_normal} does not exist{my_print_decor}')

    if run_id in job_status:
        if re.search(rf"{run_id}\s+{user_name}\s+R", job_status) is not None:
            print(f'{my_print_decor}{run_id} running, wait ...{my_print_decor}', flush=True)
        else:
            print(f'{my_print_decor}{run_id} queueing, wait ...{my_print_decor}', flush=True)
        time.sleep(120)
    else:
        # check if the run finishes normally
        # open a file and check if the last line contains "Run completed successfully"
        with open(f'{rundir}/outputs/mirror.out', 'r') as file:
            lines = file.readlines()  # Read all lines in the file
            last_line = lines[-1] if lines else ''  # Get the last line if the file is not empty
            # Check if the last line contains the desired phrase
            if "Run completed successfully" in last_line:
                print("The last line indicates that the run completed successfully.")
                break
            else:
                print("The last line does not indicate a successful completion, try combining the last hotstart.nc the restart the run.")

        # combine hotstart
        hot_steps = Get_hotstart_step(f'{rundir}/outputs/')
        if len(hot_steps) == 0:
            raise Exception('No hotstart files generated before run stopped.')
        hot_step = hot_steps[-1]

        hot_combined = f'{rundir}/outputs/hotstart_it={hot_step}.nc'
        print(f'{my_print_decor}{run_id} stopped, last hotstart to combine: {hot_combined}{my_print_decor}', flush=True)

        os.chdir(f'{rundir}/outputs/')
        Replace_string_in_file(f'{rundir}/run_comb', '-i 0000', f'-i {hot_step}')

        print(f'{my_print_decor}{batch_cmd} run_comb{my_print_decor}', flush=True)
        os.system(f'{batch_cmd} {rundir}/run_comb')

	# restore run_cmb template
        os.system(f'cat {rundir}/run_comb')
        Replace_string_in_file(f'{rundir}/run_comb', f'-i {hot_step}', '-i 0000')

	# wait for combine hotstart to finish
        while combine_job_name in subprocess.getoutput(queue_query_str):
            time.sleep(20)
            print(f'{my_print_decor}waiting for job {combine_job_name} to finish{my_print_decor}', flush=True)
        if not os.path.exists(hot_combined):
            raise Exception(f'Failed generating {hot_combined}')
        time.sleep(20)

        # link combined hotstart.nc
        os.chdir(f'{rundir}')
        print(f'{my_print_decor}linking {hot_combined}{my_print_decor}', flush=True)
        os.system(f'rm {rundir}/hotstart.nc')
        os.symlink(hot_combined, f'{rundir}/hotstart.nc')
        
        # dealing with the paired normal run
        if rundir_normal != '':
            if re.search(rf"{run_id_normal}\s+{user_name}\s+R", job_status) is not None:
                print(f'{my_print_decor}{run_id_normal} has started, skipping syncing files with it{my_print_decor}', flush=True)
            else:
                print(f'{my_print_decor}copying staout* and *.out to {rundir_normal}{my_print_decor}', flush=True)
                os.system(f'cp {rundir}/outputs/staout* {rundir_normal}/outputs/')
                os.system(f'cp {rundir}/outputs/*.out {rundir_normal}/outputs/')

        Replace_string_in_file(f'{rundir}/param.nml', 'ihot = 1', 'ihot = 2')

        # Replace_string_in_file('~/bin/run_test', 'RUNxxx', f'RUN{run_id}')

        # submit new job
        os.chdir(f'{rundir}')
        print(f'{my_print_decor}{batch_cmd} run_test{my_print_decor}', flush=True)
        os.system(f'{batch_cmd} run_test')

	# restore template run_test
        # Replace_string_in_file('~/bin/run_test', f'RUN{run_id}', 'RUNxxx')

print(f'{my_print_decor}Task completed. {my_print_decor}', flush=True)
