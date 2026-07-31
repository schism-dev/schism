"""
Generate hotstart.nc file for STOFS3D model
"""

import os
import subprocess

from glob import glob
from pylib import read
from pyschism.mesh.hgrid import Hgrid
from pyschism.forcing.hycom.hycom2schism import DownloadHycom
from .tweak_stofs3d_hotstart import tweak_stofs3d_hotstart


LEVEE_POLYGON_SHAPEFILE = "levee_pump_polys_2026_with_poly_type_lonlat.shp"


def gen_hotstart_nc(
    wdir=None,
    start_date=None,
    fortran_exe=None,
    aviso_path=None,
    python_exe=None,
):
    """
    Generate hotstart.nc file for STOFS3D model in the specified working directory.
    This function assumes that necessary input files are already present in the working directory.
    """

    script_path = os.path.dirname(os.path.abspath(__file__))
    aviso_script_path = os.path.abspath(os.path.join(script_path, '..', 'AVISO'))
    if wdir is None:
        wdir = os.getcwd() + os.sep

    if fortran_exe is None and python_exe is not None:
        print("Warning: gen_hotstart_nc(python_exe=...) is deprecated; use fortran_exe=... instead")
        fortran_exe = python_exe

    while True:
        if fortran_exe is None:
            print(
                "If you have successfully compiled schism source code, the executable is at "
                "schism/build/bin/gen_hot_3Dth_from_hycom."
            )
            fortran_exe = input('Please provide the path to the Fortran executable to run the script: ').strip()
        if os.path.exists(fortran_exe):
            break
        print(f'Fortran executable not found: {fortran_exe}')
        fortran_exe = None

    if aviso_path is not None and not os.path.exists('aviso.nc'):
        os.symlink(str(aviso_path), 'aviso.nc')

    # check for aviso.nc
    while True:  # search for aviso.nc until it is provided
        if not os.path.exists('aviso.nc'):
            print('\n' + '-'*50)
            print(f'aviso.nc not found, please download the data to {wdir}')
            print(f'See instructions in {aviso_script_path}/README and '
                  f'{aviso_script_path}/download_adt.py')
            input('Press Enter to continue after downloading aviso.nc')
            print('\n' + '-'*50)
        else:
            break

    # link grid files
    os.system('ln -sf ../hgrid.gr3 .')
    os.system('ln -sf ../hgrid.gr3 hgrid.ll')
    os.system('ln -sf ../vgrid.in .')

    # download hycom files
    hgrid = Hgrid.open('./hgrid.gr3', crs='epsg:4326')
    hycom = DownloadHycom(hgrid)
    hycom.fetch_data(start_date)
    hycom_files = glob('hycom_*.nc')
    if len(hycom_files) != 1:
        raise RuntimeError('There should be exactly one HYCOM file downloaded')
    for f in ['TS_1.nc', 'UV_1.nc', 'SSH_1.nc']:
        os.system(f'ln -sf {hycom_files[0]} {f}')

    # make estuary.gr3 and fill with 0
    hg = read('hgrid.gr3')
    hg.write('estuary.gr3', value=0.0)

    # write the input file for the external Fortran script
    with open('gen_hot_3Dth_from_nc.in', 'w') as f:
        f.write('1 !1: include vel and elev. in hotstart.in (and *[23D].th start from non-0 values); 0: only T,S\n')
        f.write('20 33 !T,S values for estuary points defined in estuary.gr3, overwritten by obs if any\n')
        f.write('20 33 !T,S values for pts outside bg grid in nc\n')
        f.write('86400. !time step in .nc [sec]\n')
        f.write('1 1  !# of open bnds that need *3D.th; list of IDs\n')
        f.write('9999 !# of days needed\n')
        f.write('1 ! # of HYCOM stacks (e.g., 3 if there are SSH_1.nc, SSH_2.nc, SSH_3.nc)\n')

    # generate hotstart.nc using the external fortran script
    subprocess.run([fortran_exe], check=True)

    # rename hotstart.nc for further processing
    os.system('mv hotstart.nc hotstart.nc.hycom')

    # tweak hotstart.nc based on coastal observations if any
    levee_polygon_basename = os.path.splitext(LEVEE_POLYGON_SHAPEFILE)[0]
    os.system(f'cp {script_path}/{levee_polygon_basename}.* {wdir}/')  # copy levee polygon shapefiles to working directory
    completed = tweak_stofs3d_hotstart(
        wdir=wdir,
        hotstart_date_str=start_date.strftime('%Y-%m-%d'),
        city_shapefile_names=[LEVEE_POLYGON_SHAPEFILE],
        aviso_file='aviso.nc',
    )

    return completed
