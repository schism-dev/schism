Pre-processing scripts used for coupling SCHISM to National Water Model (NWM) in NOAA's Water Initiative project
(1) More details are in the online documentation: http://ccrm.vims.edu/yinglong/feiye/Workshop_20190701/TEMP/Doc/main.html

(2) hgrid and vgrid are needed for most scripts, so make them first

(3) When you use these scripts for the first time on a new cluster,
    some of them may need to be recomplied, including:

./SCRIPTS/*.cpp
./SCRIPTS/Elev_IC/gen_elev.f90
./SCRIPTS/Shapiro/gen_slope_filter.f90
./SCRIPTS/HYCOM_nudge/gen_nudge_from_hycom.f90  
./SCRIPTS/Vgrid/gen_vqs.f90
./SCRIPTS/Hot/gen_hot_3Dth_from_hycom.f90  
./SCRIPTS/Rough3/gen_source2.f90
./SCRIPTS/DEM/interpolate_depth_structured2.f90
./SCRIPTS/DEM_USGS/interpolate_depth_structured2.f90

A sample compiling cmd is provided in the beginning few lines in each source code.
The binaries provided here are compiled with ifort and gcc.


(4) When setting up a new run from scratch, it is recommended to copy the following sub-folders to your run directory:
	[$your_script_dir] cp -r Grid_manipulation/ Prop/ Shapiro/ Rough3/ Elev_IC/ Hot/ Vgrid/ NWM_Coupling/ Load_Bathy/ $your_run_dir 

  , and creat symbolic links for a few large data files:
	[$your_script_dir] cd $your_run_dir
	[$your_run_dir] ln -s $your_script_dir/HYCOM_IRENE_PERIOD/ .
	[$your_run_dir] ln -s $your_script_dir/DEM*/ .

  If you decide not to follow this, some paths in the "auto.pl" under each sub-folder (e.g. Rough3/auto.pl) may need to be manually changed.
  Look for the "dirs" section in each "auto.pl".


(5) Load_bathy does not work on Stampede2, because the module 'proj' (PROJ.4 - Cartographic Projections Library) is not available
