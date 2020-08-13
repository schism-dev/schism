Pre-processing scripts used for coupling SCHISM to National Water Model (NWM) for NOAA's Water Initiative project

(1) More details are provided in the online documentation for the Delaware Bay setup:
     http://ccrm.vims.edu/yinglong/feiye/Workshop_20190701/TEMP/Doc/main.html
     The Delaware Bay setup (finalized in 2019) is being used as a tutorial, so it won't be updated regularly;
     see the "Appendix" at the end of this file for the updates since then.


(2) hgrid and vgrid are needed for most scripts, so make them first


(3) When you use these scripts for the first time on a new cluster,
    some of them may need to be recomplied, including:
```
"./Grid_manipulation/*.cpp",
"./Grid_manipulation/*.f90",
"./Elev_IC/gen_elev.f90",
"./Shapiro/gen_slope_filter.f90",
"./Vgrid/gen_vqs.f90",
"./Manning/gen_source2.f90",
"./Hot_and_3Dth/gen_hot_3Dth_from_hycom.f90",
"./Hot_and_3Dth/modify_hot.f90",
"./HYCOM_nudge/gen_nudge_from_hycom.f90",
```
For your convenience, the compiled binaries are provided in Git,
which should be usable on most systems with c++/fortran compilers.
However, the last 3 scripts requires netcdf libraries and 
probably won't work until you re-compile them on your local machine/cluster.
In the future, we may put all the scripts that need to be compiled under a central folder
and implement cmake.
This has been put on hold because there is a plan to convert all scripts to python.

A sample compiling cmd is provided in the beginning few lines in each source code.
The existing binaries provided in Git were compiled with ifort and g++.


(4) When setting up a new run from scratch, it is recommended to copy the following sub-folders under NWM/ to your run directory:

	[$your_script_dir] cp -rL Grid_manipulation/ Prop/ Shapiro/ HYCOM_nudge/ Manning/ Drag/  Elev_IC/ Hot/ Vgrid/ NWM_Coupling/ Load_Bathy_mpi/ $your_run_dir 
    Note that the -L option is recommended to follow the symbolic links
    Many scripts will use tools in Grid_manipulation/.
    In the future, symbolic links to Grid_manipulation/$sometools will be provided in subfolders (such as Prop/) so that
    users do not need to explicitly copy Grid_manipulation/ to the run dir. 
    This has been put on hold because there is a plan to convert all scripts to python.

  Then, creat symbolic links for the HYCOM files of the period you need, e.g.:
	[$your_script_dir] cd $your_run_dir
	[$your_run_dir] ln -s $your_git_dir/HYCOM_IRENE_PERIOD/ .


(5) If you get these scripts from Github, then you have to mannually download some large data files to your own machine/cluster.
    The large data files (such as HYCOM_IRENE_PERIOD/) is not saved on Git and you'll see a broken link at the end of Step (4)
    When this happens, look for the "REMOVED_FILE_LOCATION" file under the same folder as the large data file you're trying to link.
    The download links are provided in "REMOVED_FILE_LOCATION".
    Contact the VIMS group if any "REMOVED_FILE_LOCATION" file is missing, which means the data have not been shared publicly yet.


(6) Specifically for Stampede2:

    (6.1) request an interative session for Hot_and_3Dth/ and HYCOM_nudge/, since they can take a while (around 1 hour);
          do "./auto.pl" under the interative session.

    (6.2) You can also submit the "./auto.pl" to queue (2 hours on the dev/debug queue should be enough).



Appendix A. Updates since the 2019 Delaware Bay setup

  Online documentation needed.
  Compared to the previous documentation (http://ccrm.vims.edu/yinglong/feiye/Workshop_20190701/TEMP/Doc/main.html),
  the following items are changed:

    · param.nml:
         air-sea exchange
         ielm_transport

    · Coupling with NWM
         Fei's matlab script -> Wei's Fortran script
         , automation not available

    · *.prop:
         For the automation, default *.reg files are no longer provided 
           because the script is intended for different events;
           instead, user-specified regions (*.reg or *.rng) are needed in the subfolders
           "Prop/tvd/" and "Prop/fluxflag/".
         Non-tvd regions are needed in "Prop/tvd/", any basename will do, e.g., "DE_Bay.reg";
         Regions for recording flux transects are needed in "Prop/fluxflag/",
           the basenames of a transect should follow this convention:
           for "some_river", specify two regions at each transect, e.g.,
           "some_river+.reg" and "some_river-.reg"
           , where the flux is measured from the "+" region to the "-" region.
           The automation script will glob all "*+.reg" and "*-.reg" and assign values properly.
         Sample regions are provided in Prop/Samples*/

    · Bottom friction:
         Spatially varying values based on specified regions temporarily commented out,
         since the latest setup is purely dependant on depth.
         2D uses Manning/, 3D uses Drag/, the methods are the same:
           dp>-1 m: val1; dp<-3 m: val2; linear transition in between.
         The default values are no longer provided, since these can vary for different runs/events
         , but samples are provided in Manning/Sample*/ and Drag/Sample*/

    · elev.ic (same as the elev i.c. in hotstart.nc if ihot=1)
         2e-2 m above ground -> 1e-6 m under ground

    · shapiro.gr3
         previously:
         1.Based on slope, up to 0.5;
         2.Based on user-specified regions: 0.5 in ChesBay; 0.1 and 0.2 in DE Bay
         now:
         1.Based on slope, up to 0.5;
         2.Based on depth: 0.05~0.2 when dp=50m~20m; 0.2 when dp<20m
         3.Based on local regions: 0.5 in all coastal bays with a small transition zone of max(0.2, shapiro)

    · hotstart.nc:
         The procdure itself is not changed, but more local regions are added.
         Three Steps:
         1. HYCOM
         2. Modifications 
           Elev:
             same as elev.ic
           Sal/Tem in coastal bays:
             a) USGS obs are used for the regions influenced by Harvey and Florence;
             b) approximated values are used in DE Bay;
             c) CBP obs are used for ChesBay;
             For a) and b):
               · The obs/approximation is pre-defined in a small background grid (*.gr3)
                   for each region of interest (*.reg), then interpolated onto the full model grid.
               · Since Sal/Tem varies for different periods/regions,
                   it's the users' responsibility to ensure the integrity of
                   the background grid and the region of interest,
                   and specify them properly in the "inputs" section of the "auto.pl".
               · The convention of the "inputs" section is self-explanatory.
               · USGS obs is not handled automatically because the data points are scarce and
                   the spatial/temporal availability is not guaranteed at each station.
                   Standalone Python scripts are used for generating a draft of the background grid
                   , with nodal values based on all available USGS stations in the region of interest at the specified time.
                   However, due to the scarcity of the obs points, users need to mannually add points to properly delineate 
                   rivers and bays based on his/her best estimate.
             For c):
               · The script handles ChesBay obs automatically since the spatial/temporal availability is good.
           Velocity:
             Not specified
         3. Further limiting watershed salinity based on depth: 
           dp>-1 m: existing value; dp<-3 m: 0 PSU; linear transition in between.
