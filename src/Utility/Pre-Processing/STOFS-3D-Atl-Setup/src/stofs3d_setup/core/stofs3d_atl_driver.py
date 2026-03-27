'''
This is the driver script to set up the STOFS-3D-ATL model.
Simpler tasks such as generating uniform *.gr3 are included in this script.
More complex tasks such as generating source_sink are imported as modules from subfolders.

Some scripts are written in Fortran and C++ for speed.
You may need to compile them before running this driver script.
See compile instructions at the beginning of each script:

- Vgrid/: gen_vqs.f90, change_vgrid.f90

Usage:
    see sample_usage.py in the same folder.

The script will prepare the model folders and keep a record of itself in the
model input folder.

For the author, there are a few "todo" items in the script.
1) The second call to gen_sourcesink_nwm is not optimal, because it downloads
    the NWM data again if config.nwm_cache_folder is not set.
    This cost time for longer runs.
    It is better to directly use generated sources/sinks.
    As a temporary solution, run the script with config.relocate_source = False,
    locate the downloaded data in {model_input_path}/Source_sink/original_source_sink/,
    symlink all downloaded NWM data in a single folder and set config.nwm_cache_folder.
    The symlink is necessary because NWM data may be stored in subfolders of different years.

2) hotstart.nc

Temporary changes to be tested:
1) admit more BlueTopo in DEM_loading, see bluetopo_region.shp
2) adopted regional depth tweaks from SECOFS, added DEFAULT_REGIONAL_TWEAKS2 in regional_tweaks.py
3) adopted drag tweaks from SECOFS, added a few regions in gen_drag()
'''


# ------------------------- Import modules ---------------------------
import os
from pathlib import Path
import socket

from datetime import timedelta

# self-defined modules
from pylib import read_schism_vgrid
if 'gulf' in socket.gethostname():
    from pylib_experimental.schism_file import cread_schism_hgrid as schism_read
    print('Using c++ function to accelerate hgrid reading')
else:
    from pylib import schism_grid as schism_read
    print('Using python function to read hgrid')
from pyschism.mesh import Hgrid as pyschism_Hgrid

# Import from the sub folders. These are not from installed packages.
from ..ops.simple_tasks import gen_nudge_coef, gen_shapiro_strength, gen_soil, gen_drag, gen_elev_ic
from ..ops.simple_tasks import gen_3dbc, gen_elev2d, gen_nudge_stofs
from ..ops.Vgrid.gen_vqs import gen_vqs
from ..ops.River.gen_Canada_river_flux_th import gen_Canada_river_flux_th
from ..utils.utils import mkcd_new_dir, try_remove, prep_run_dir
from ..ops.Prop.gen_tvd import gen_tvd_prop
from ..ops.Bctides.bctides.bctides import Bctides  # temporary, bctides.py will be merged into pyschism
# from pyschism.forcing.bctides import Bctides
from ..ops.Source_sink.assemble_source_sink import assemble_source_sink
from ..ops.Hotstart.gen_hotstart_nc import gen_hotstart_nc

# Import configuration
from ..config.stofs3d_atl_config import ConfigStofs3dAtlantic

# Global variables:
# use the full path of the Pre-Processing dir inside your schism repo
# e.g, script_path = '/my_dir/schism/src/Utility/Pre-Processing/'
# If you are running this script in the schism repo, you can also use the following line:
script_path = Path(__file__).resolve().parent.parent / "ops"
print(f"script_path: {script_path}")

DRIVER_PRINT_PREFIX = '\n-----------------STOFS3D-ATL driver:---------------------\n'


# ---------------------------------------------------------------------
#       Main function to generate inputs for STOFS-3D-ATL
#       Only house keeping here, the core functions are imported
# ---------------------------------------------------------------------
def stofs3d_atl_driver(
    hgrid_path: str,
    vgrid_path: str,
    config: ConfigStofs3dAtlantic,
    project_dir: str, runid: str, scr_dir: str,
    input_files: dict = None,
):
    '''
    Main function to generate inputs for STOFS3D-ATL.
    '''

    if input_files is None:
        input_files = {
            'vgrid': True,
            'bctides': True,
            'gr3': True, 'nudge_gr3': True, 'shapiro': True, 'drag': True, 'elev_ic': True,
            'tvd.prop': True,
            'flux_th': True, 'source_sink': True,
            'hotstart.nc': False,
            '3D.th.nc': True, 'elev2D.th.nc': True, '*nu.nc': True,
        }

    print(f'{DRIVER_PRINT_PREFIX}reading hgrid from {hgrid_path} ...')
    hgrid = schism_read(str(hgrid_path))

    # -----------------begin generating model inputs---------------------

    # define and make the model_input_path, the run_dir and the output dir
    model_input_path, run_dir, _ = prep_run_dir(project_dir, runid, scr_dir=scr_dir)

    # make a copy of the script itself to the model_input_path
    os.system(f'cp -rf {script_path} {model_input_path}/Pre_processing_scripts_backup')
    # make a copy of the hgrid to the model_input_path
    os.system(f'cp {hgrid_path} {model_input_path}/hgrid.gr3')

    # ------------------vgrid---------------------
    if vgrid_path is not None:
        if os.path.exists(vgrid_path):
            print(f'{DRIVER_PRINT_PREFIX}linking provided vgrid from {vgrid_path} ...')
            os.system(f'ln -sf {vgrid_path} {model_input_path}/vgrid.in')
        else:
            raise FileNotFoundError(f'Provided vgrid_path {vgrid_path} does not exist!')
    else:
        if input_files['vgrid']:
            sub_dir = 'Vgrid'
            print(f'{DRIVER_PRINT_PREFIX}Generating vgrid.in ...')
            mkcd_new_dir(f'{model_input_path}/{sub_dir}')
            os.system('ln -sf ../hgrid.gr3 .')

            gen_vqs(hgrid_file=hgrid_path, output_dir=f'{model_input_path}/{sub_dir}')

            # # call a fortran program to generate vgrid.in
            # print(f'compile the fortran program {script_path}/Vgrid/gen_vqs if necessary')
            # fortran_process = subprocess.Popen(
            #     f'{script_path}/Vgrid/gen_vqs', stdin=subprocess.PIPE
            # )
            # # the command line argument 0 means no outputs of a sample vgrid along a given transect
            # fortran_process.communicate(input=str(0).encode())

            print(f'{DRIVER_PRINT_PREFIX}converting the format of vgrid.in ...')
            vg = read_schism_vgrid(f'{model_input_path}/{sub_dir}/vgrid.in')
            os.rename('vgrid.in', 'vgrid.in.old')
            vg.save(f'{model_input_path}/{sub_dir}/vgrid.in')

            os.chdir(model_input_path)
            os.system(f'ln -sf {sub_dir}/vgrid.in .')

    # -----------------bctides---------------------
    if input_files['bctides']:
        hgrid_pyschism = pyschism_Hgrid.open(hgrid_path, crs='epsg:4326')
        sub_dir = 'Bctides'
        print(f'{DRIVER_PRINT_PREFIX}Generating bctides.in ...')
        mkcd_new_dir(f'{model_input_path}/{sub_dir}')
        Bctides(
            hgrid=hgrid_pyschism,  # needs hgrid in pyschism's Hgrid class
            bc_flags=config.bc_flags,
            bc_const=config.bc_const,
            bc_relax=config.bc_relax,
            database='fes2014',
        ).write(
            output_directory=f'{model_input_path}/{sub_dir}',
            start_date=config.startdate,
            rnday=config.rnday,
        )

        os.chdir(run_dir)
        os.system(f'ln -sf ../I{runid}/{sub_dir}/bctides.in .')
        os.chdir(model_input_path)

    # -----------------elev2D.th.nc---------------------
    if input_files['elev2D.th.nc']:
        sub_dir = 'Elev2Dth'
        print(f'{DRIVER_PRINT_PREFIX}Generating elev2D.th.nc ...')
        mkcd_new_dir(f'{model_input_path}/{sub_dir}', remove=False)

        gen_elev2d(
            hgrid_fname=f'{model_input_path}/hgrid.gr3',
            outdir=f'{model_input_path}/{sub_dir}',
            start_date=config.startdate, rnday=config.rnday,
            ocean_bnd_ids=config.ocean_bnd_ids,
            uniform_shift=config.elev2d_uniform_shift,
        )

        os.chdir(run_dir)
        os.system(f'ln -sf ../I{runid}/{sub_dir}/elev2D.th.nc .')

    # -----------------spatially uniform Gr3---------------------
    if input_files['gr3']:
        sub_dir = 'Gr3'
        print(f'{DRIVER_PRINT_PREFIX}Generating spatially uniform gr3 files ...')
        mkcd_new_dir(f'{model_input_path}/{sub_dir}')

        for name, value in config.gr3_values.items():
            print(f'{DRIVER_PRINT_PREFIX}Generating {name}.gr3 ...')
            try_remove(f'{model_input_path}/{sub_dir}/{name}.gr3')
            hgrid.write(f'{model_input_path}/{sub_dir}/{name}.gr3', value=value)

        os.chdir(run_dir)
        os.system(f'ln -sf ../I{runid}/{sub_dir}/*.gr3 .')
        os.chdir(model_input_path)

    # -------------------------------------------------
    # ---------begin spatially varying Gr3 ------------

    # -----------------nudge.gr3---------------------
    if input_files['nudge_gr3']:
        sub_dir = 'Nudge_gr3'
        print(f'{DRIVER_PRINT_PREFIX}Generating nudge.gr3 ...')
        mkcd_new_dir(f'{model_input_path}/{sub_dir}')

        # generate nudging coefficient based on the proximity to the open boundaries
        nudge_coef = gen_nudge_coef(
            hgrid, rlmax=config.nudging_zone_width, rnu_day=config.nudging_day)

        # write nudge.gr3
        hgrid.save(f'{model_input_path}/{sub_dir}/nudge.gr3', value=nudge_coef)

        # link nudge.gr3 to SAL_nudge.gr3 and TEM_nudge.gr3
        try_remove(f'{model_input_path}/{sub_dir}/SAL_nudge.gr3')
        os.system('ln -s nudge.gr3 SAL_nudge.gr3')
        try_remove(f'{model_input_path}/{sub_dir}/TEM_nudge.gr3')
        os.system('ln -s nudge.gr3 TEM_nudge.gr3')

        os.chdir(run_dir)
        os.system(f'ln -sf ../I{runid}/{sub_dir}/*_nudge.gr3 .')
        os.chdir(model_input_path)

    # -----------------shapiro.gr3---------------------
    if input_files['shapiro']:
        sub_dir = 'Shapiro'
        print(f'{DRIVER_PRINT_PREFIX}Generating shapiro.gr3 ...')
        mkcd_new_dir(f'{model_input_path}/{sub_dir}')

        if config.shapiro_zone_width > 0:
            # use the same method for nudge.gr3 to generate a buffer zone along the open boundaries
            distribute_coef = gen_nudge_coef(hgrid, rlmax=config.shapiro_zone_width)
        else:
            distribute_coef = None
        # using a distribution coefficient simliar to the nudging coefficient
        shapiro = gen_shapiro_strength(
            hgrid, init_shapiro_dist=distribute_coef, tilt=config.shapiro_tilt)

        hgrid.save(f'{model_input_path}/{sub_dir}/shapiro.gr3', value=shapiro)

        os.chdir(run_dir)
        os.system(f'ln -sf ../I{runid}/{sub_dir}/shapiro.gr3 .')
        os.chdir(model_input_path)

    # -----------------soil*.gr3---------------------
    if input_files['soil']:
        sub_dir = 'Soil'
        print(f'{DRIVER_PRINT_PREFIX}Generating soil_conductivity.gr3 and soil_thick.gr3 ...')
        mkcd_new_dir(f'{model_input_path}/{sub_dir}')

        soil_conductivity, soil_thick = gen_soil(hgrid)

        hgrid.save(f'{model_input_path}/{sub_dir}/soil_conductivity.gr3', value=soil_conductivity)
        hgrid.save(f'{model_input_path}/{sub_dir}/soil_thick.gr3', value=soil_thick)

        os.chdir(run_dir)
        os.system(f'ln -sf ../I{runid}/{sub_dir}/soil_conductivity.gr3 .')
        os.system(f'ln -sf ../I{runid}/{sub_dir}/soil_thick.gr3 .')
        os.chdir(model_input_path)

    # -----------------drag.gr3---------------------
    if input_files['drag']:
        sub_dir = 'Drag'
        print(f'{DRIVER_PRINT_PREFIX}Generating drag.gr3 ...')
        mkcd_new_dir(f'{model_input_path}/{sub_dir}')
        drag = gen_drag(hgrid)

        hgrid.save(f'{model_input_path}/{sub_dir}/drag.gr3', value=drag)

        os.chdir(run_dir)
        os.system(f'ln -sf ../I{runid}/{sub_dir}/drag.gr3 .')
        os.chdir(model_input_path)

    # -----------------elev_ic---------------------
    if input_files['elev_ic']:
        sub_dir = 'Elev_ic'
        print(f'{DRIVER_PRINT_PREFIX}Generating elev_ic.gr3 ...')
        mkcd_new_dir(f'{model_input_path}/{sub_dir}')
        elev_ic = gen_elev_ic(
            hgrid, h0=0.1,
            city_shape_fnames=[f'{script_path}/Hotstart/LA_urban_polys_lonlat_v2.shp'],
            base_elev_ic=None  # '/sciclone/schism10/feiye/STOFS3D-v7/Inputs/
                               # Iv7_Francine_update/Reinit_hot/elev0.gr3'
        )

        hgrid.save(f'{model_input_path}/{sub_dir}/elev_ic.gr3', value=elev_ic)
        os.symlink(
            f'{model_input_path}/{sub_dir}/elev_ic.gr3',
            f'{model_input_path}/{sub_dir}/elev.ic')

        os.chdir(run_dir)
        os.system(f'ln -sf ../I{runid}/{sub_dir}/elev.ic .')
        os.chdir(model_input_path)

    # ----- end spatially varying Gr3 -----------------
    # -------------------------------------------------

    # -------flux.th at St. Lawrence River boundary--------
    if input_files['flux_th']:
        sub_dir = 'Flux_th'
        print(f'{DRIVER_PRINT_PREFIX}Generating flux.th at St. Lawrence River boundary ...')
        mkcd_new_dir(f'{model_input_path}/{sub_dir}')

        # generate flux.th
        gen_Canada_river_flux_th(
            start_datetime=config.startdate,
            end_datetime=config.startdate + timedelta(days=config.rnday+1),
        )
        os.chdir(run_dir)
        os.system(f'ln -sf ../I{runid}/{sub_dir}/flux.th .')
        os.chdir(model_input_path)

    # -----------------source_sink---------------------
    # Note:

    # Normal case: main hgrid has pseudo channels as feeders.
    # An "hgrid_without_feeders" needs to be provided in the stofs3d_atl_config.py
    # to generate the original sources/sinks.
    # This is necessary because the feeder channels are not real channels
    # and "hgrid_without_feeders" gives the correct source locations.
    # Moreover, the land boundary of "hgrid_without_feeders" is set to avoid
    # NWM segments weaving in and out of the schism domain.
    # The script will then relocate the original sources to the head of each feeder channel.

    # Special case: main hgrid has no feeders or has real channels as feeders.
    # In this case, the original sources/sinks are generated based on the main hgrid directly,
    # and the option "hgrid_without_feeders" should be set to None.

    if input_files['source_sink']:
        sub_dir = 'Source_sink'
        print(f'{DRIVER_PRINT_PREFIX}Generating source_sink.in ...')

        mkcd_new_dir(f'{model_input_path}/{sub_dir}')
        os.chdir(f'{model_input_path}/{sub_dir}')

        assemble_source_sink(config, hgrid, model_input_path=model_input_path,
                             wdir=f'{model_input_path}/{sub_dir}')

        os.chdir(run_dir)
        for f in ['source_sink.in', 'source.nc', 'vsource.th', 'msource.th', 'vsink.th']:
            os.system(f'ln -sf ../I{runid}/{sub_dir}/{f} .')
        os.chdir(model_input_path)

        print('Done generating source/sink files.')

    # -----------------*prop---------------------
    if input_files['*.prop']:
        sub_dir = 'Prop'
        print(f'{DRIVER_PRINT_PREFIX}Generating *prop files ...')
        mkcd_new_dir(f'{model_input_path}/{sub_dir}')

        os.system(f'cp {script_path}/Prop/* .')
        gen_tvd_prop(hgrid, config.tvd_regions)

        os.chdir(run_dir)
        os.system(f'ln -sf ../I{runid}/{sub_dir}/tvd.prop .')

    # -----------------hotstart.nc---------------------
    if input_files['hotstart.nc']:
        sub_dir = 'Hotstart'
        print(f'{DRIVER_PRINT_PREFIX}Generating hotstart.nc ...')
        mkcd_new_dir(f'{model_input_path}/{sub_dir}')

        # generate hotstart.nc
        gen_hotstart_nc(
            start_date=config.startdate
        )

        os.chdir(run_dir)
        os.system(f'ln -sf ../I{runid}/{sub_dir}/hotstart.nc .')

    # -----------------*nu.nc---------------------
    if input_files['*nu.nc']:
        sub_dir = 'Nudge'
        print(f'{DRIVER_PRINT_PREFIX}Generating *nu.nc files ...')
        mkcd_new_dir(f'{model_input_path}/{sub_dir}')

        os.system(f'cp {script_path}/Nudge/* .')
        gen_nudge_stofs(
            hgrid_fname=f'{model_input_path}/hgrid.gr3',
            vgrid_fname=f'{model_input_path}/vgrid.in',
            outdir=f'{model_input_path}/{sub_dir}',
            rnday=config.rnday, start_date=config.startdate,
            ocean_bnd_ids=config.ocean_bnd_ids,
            # use pre-downloaded files to speed up
            hycom_download_dir='/sciclone/schism10/feiye/STOFS3D-v8/HYCOM_DOWNLOADS_Nudge/',
            # hycom_download_dir=None,
        )

        os.chdir(run_dir)
        os.system(f'ln -sf ../I{runid}/{sub_dir}/*nu.nc .')

    # -----------------3D.th.nc---------------------
    if input_files['3D.th.nc']:
        sub_dir = '3Dth'
        print(f'{DRIVER_PRINT_PREFIX}Generating *3D.th.nc files ...')
        mkcd_new_dir(f'{model_input_path}/{sub_dir}')

        # generate *D.th.nc
        gen_3dbc(
            hgrid_fname=f'{model_input_path}/hgrid.gr3',
            vgrid_fname=f'{model_input_path}/vgrid.in',
            outdir=f'{model_input_path}/{sub_dir}',
            start_date=config.startdate, rnday=config.rnday,
            ocean_bnd_ids=config.ocean_bnd_ids,
            hycom_download_dir='/sciclone/schism10/feiye/STOFS3D-v8//HYCOM_DOWNLOADS_3Dth/',
        )

        os.chdir(run_dir)
        for f in ['SAL_3D.th.nc', 'TEM_3D.th.nc', 'uv3D.th.nc']:
            os.system(f'rm {f}')
            os.system(f'ln -s ../I{runid}/{sub_dir}/{f} .')


if __name__ == '__main__':
    print('Done')
