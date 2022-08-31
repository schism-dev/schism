## Visualization with Matlab
The directory Utility/Vis_Matlab/ has matlab scripts that can visualize outputs along a horizontal slab (at a fixed
 z level or at a sigma level) or vertical transects. In particular, `SCHISM_SLAB2.m`  and `SCHISM_TRANSECT2.m` for for 
 the new scribed outputs, while `SCHISM_SLAB.m` and `SCHISM_TRANSECT.m` are for the old outputs (schout*.nc).

## Visualization with Python
[more coming]

## Visualization with VisIT
The most comprehensive way to visualize SCHISM nc4 outputs is via VisIT.

Shortly after v5.9.0, we have successfully implemented a new mode of efficient I/O using dedicated 'scribes'.
 At runtime, the user needs to specify number of scribe cores (= # of 3D outputs variables 
(vectors counted as 2) plus 1), and the code, compiled without `OLDIO`, will output 
 combined netcdf outputs for each 3D variable and also all 2D variables in `out2d*.nc`. 
Sample 3D outputs are: `salinity_*.nc`, `horizontalVelX_*.nc` etc - note that vectors variable names end with `X,Y`.

You can download newer versions of VisIT plugins c/o Dr. Jon Shu, DWR by following these steps:

**On Windows 7 or 10**

1. First download VisIT from [LLNL](https://wci.llnl.gov/simulation/computer-codes/visit/downloads) site and install. Use default dir (and record it), e.g. `C:\Users\username\AppData\Local\Programs\LLNL\VisIt*`
2. Make sure MS visualc++ 2012 x64 is installed. If not, google it and install and restart (this is required for using multi-core VisIT)
3. Download pre-built plug-in, developed at California Dept of Water Resource
    * [For VisIT v2.13.3](https://cadwr.box.com/s/tiuaeb87ka1tb7i185l0cpf7g3da2h63)
    * [For VisIT v3.1.4](https://cadwr.box.com/s/jrj83yycyumhkf9kqod7teb7biqovmta)
    
    You need to put dlls on `C:\Users\username\Documents\VisIt\databases`. Alternatively, please email Jon Shu at Qiang.shu[at]water[dot]ca[dot]gov and let him know the version of visit you are using.

4. After these steps, you'd be able to read in SCHISM outputs in ViSIT; look for `SCHISM`, `gr3` format from the dropdown list. To load in vectors, select only the `X` file.

**On Linux Systems**

Newer versions can be found at the master branch of [github](https://github.com/schism-dev/schism_visit_plugin).

!!!note
    Note that the new plugins also work with the old I/O (combined `schout*.nc`).
