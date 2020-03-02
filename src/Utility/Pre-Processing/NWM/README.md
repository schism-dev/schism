Pre-processing scripts used for coupling SCHISM to National Water Model (NWM) for NOAA's Water Initiative project


(1) More details are in the online documentation: http://ccrm.vims.edu/yinglong/feiye/Workshop_20190701/TEMP/Doc/main.html


(2) hgrid and vgrid are needed for most scripts, so make them first


(3) When you use these scripts for the first time on a new cluster,
    some of them may need to be recomplied (cmake will be added soon), including:
```
"./Grid_manipulation/*.cpp",
"./Grid_manipulation/*.f90",
"./Elev_IC/gen_elev.f90",
"./Shapiro/gen_slope_filter.f90",
"./HYCOM_nudge/gen_nudge_from_hycom.f90",
"./Vgrid/gen_vqs.f90",
"./Hot/gen_hot_3Dth_from_hycom.f90",
"./Hot/modify_hot.f90",
"./Manning/gen_source2.f90",
"./Drag/gen_source2.f90",
"./Rough2D/gen_source2.f90",
"./Rough3D/gen_source2.f90",
"./DEM/interpolate_depth_structured2.f90",
"./DEM_USGS/interpolate_depth_structured2.f90"
```
(Rough2D and Rough3D will be replaced by Roughness and Drag in the future)

A sample compiling cmd is provided in the beginning few lines in each source code.
The sample binaries provided here were compiled with ifort and gcc.

Also, make sure you have appropriate netcdf libraries; as a test, see if you can successfully compile:
./HYCOM_nudge/gen_nudge_from_hycom.f90


(4) When setting up a new run from scratch, it is recommended to copy the following sub-folders to your run directory:

	[$your_script_dir] cp -r Grid_manipulation/ Prop/ Shapiro/ HYCOM_nudge/ Rough3/ Elev_IC/ Hot/ Vgrid/ NWM_Coupling/ Load_Bathy/ $your_run_dir 

  , and creat symbolic links for a few large data files:

	[$your_script_dir] cd $your_run_dir

	[$your_run_dir] ln -s $your_script_dir/HYCOM_IRENE_PERIOD/ .

	[$your_run_dir] ln -s $your_script_dir/DEM*/ .

  If you decide not to follow this, some paths in the "auto.pl" under each sub-folder (e.g. Rough3/auto.pl) may need to be manually changed.
  Look for the "dirs" section in each "auto.pl".


(5) If you get these scripts from Github, then you have to mannually download some large data files.
    The download links are provided in:

    DEM_USGS/REMOVED_FILE_LOCATION (not necessay if you are not working on the hgrid)

    HYCOM_IRENE_PERIOD/REMOVED_FILE_LOCATION

    If wget doesn't work on some clusters, just paste the link in your web browser.


(6) Specifically for Stampede2:

    (6.1) Load_bathy does not work on Stampede2, because the module 'proj' (PROJ.4 - Cartographic Projections Library) is not available;

    (6.2) request an interative session for Hot/ and HYCOM_nudge/, since they can take a while (around 1 hour);

    (6.3) You can also submit any "auto.pl" to queue (2 hours on dev queue should be enough).



