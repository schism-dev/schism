Pre-processing scripts used for coupling SCHISM to National Water Model (NWM) in NOAA's Water Initiative project


(1) Also refer to the online documentation: http://ccrm.vims.edu/yinglong/feiye/Workshop_20190701/TEMP/Doc/main.html

(2) hgrid and vgrid are needed for most scripts

(3) Some scripts may need to be recomplied, including:

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


(4) When setting up a new run from scratch, it is recommended to copy the following scripts to your run directory:
	[$your_script_dir] cp -r Prop/ Shapiro/ Rough3/ Elev_IC/ Hot/ Vgrid/ NWM_Coupling/ Load_Bathy/ $your_run_dir 

    and creat symbolic links for a few others:
	[$your_script_dir] cd $your_run_dir
	[$your_run_dir] ln -s $your_script_dir/SCRIPTS/auto_edit_prop .
	[$your_run_dir] ln -s $your_script_dir/SCRIPTS/auto_edit_region .
	[$your_run_dir] ln -s $your_script_dir/HYCOM_IRENE_PERIOD/ .
	[$your_run_dir] ln -s $your_script_dir/DEM*/ .

    ; otherwise some paths in the "auto.pl" under each sub-folder (e.g. Rough3/auto.pl) may need to be manually changed.
