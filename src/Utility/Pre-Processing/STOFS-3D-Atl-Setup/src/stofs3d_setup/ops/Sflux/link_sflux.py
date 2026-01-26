from glob import glob
import os

folder_pattern = "[0-9][0-9][0-9][0-9][0-9][0-9][0-9][0-9]"
products = ['gfs', 'hrrr' ]  # order matters
os.makedirs("sflux", exist_ok=True)
os.chdir("sflux")

for k, product in enumerate(products):
    files = sorted(glob(f"../{folder_pattern}/{product}*.nc"))
    for i, file in enumerate(files):
        for var in ['air', 'prc', 'rad']:
            os.symlink(file, f"sflux_{var}_{k+1}.{i+1:04d}.nc")

with open("sflux_inputs.txt", "w") as f:
    f.write("&sflux_inputs\n/\n")

pass