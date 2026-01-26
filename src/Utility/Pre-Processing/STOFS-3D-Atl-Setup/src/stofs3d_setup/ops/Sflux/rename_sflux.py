"""
Rename sflux file names for later versions of SCHISM.
For example, rename "sflux_air_1.0011.nc" to "sflux_air_1.11.nc".
"""

import os
from pathlib import Path
import re

def rename_sflux_files(input_dir:Path, output_dir:Path) -> None:
    """
    Rename sflux files in the input directory and link them to the output directory.
    
    Parameters:
    - input_dir: Directory containing the original sflux files.
    - output_dir: Directory where renamed files will be saved.
    """
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)
    
    # remove existing files in the output directory
    for file in output_dir.glob('sflux*.*'):
        print(f'Removing existing file: {file}')
        file.unlink()

    for old_file_path in input_dir.glob('sflux*.nc'):
        old_filename = Path(old_file_path).name
        new_filename = re.sub(r'\.(\d+)\.', lambda m: f".{int(m.group(1))}.", old_filename)
        new_file_path = output_dir / new_filename
        print(f'Linking {old_file_path} to {new_file_path}')
        os.symlink(old_file_path, new_file_path)
    
    os.symlink(input_dir / 'sflux_inputs.txt', output_dir / 'sflux_inputs.txt')


if __name__ == "__main__":
    input_directory = Path('../sflux.0/')
    output_directory = Path('/sciclone/schism10/feiye/STOFS3D-v8/R15b1_v7/sflux/')

    os.chdir(output_directory)
    rename_sflux_files(input_directory, output_directory)
    print('Renaming completed.')