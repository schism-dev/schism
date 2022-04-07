You can find some useful post-processing tools in the src/Utility directory. 

## Combining scripts
`combine_output11.f90` is used to combine process-specific netcdf to global netcdf. 

`combine_hotstart7.f90` is used combine process-specific hotstart outputs (`outputs/hotstart_0*.nc`) into `hotstart.nc`. 

`combine_gr3.f90` is used to combine process-specific `maxelev_*` and `maxdahv_*` (ASCII) into `maxelev.gr3` or `maxdahv.gr3`.

`autocombine_MPI_elfe.pl` is a simple perl wrapper script that automatically combines all available outputs during or after the run.

## One-way nesting
`OneWayNestScripts/interpolate_variables7.f90`: The purpose of this script is to generate `elev2D.th.nc`, `SAL_3D.th.nc`, `TEM_3D.th.nc` and/or `uv3D.th.nc` from a large-domain run to be used in a small-domain run. This is of limited utility now because `uv3D.th.nc` for the sub-tidal component can be generated using tidal packages such as FES.

To prepare for the nesting, first do a 2D barotropic run for a larger or same grid, with only elevation b.c. Note that 2D model is inherently more stable than 3D model, and to further enhance stability, make sure you use `indvel=1 (ishapiro=ihorcon=0)`, `thetai=1`, and also use a large Manning’s $n$ (e.g., 0.025 or larger) near the boundary. Once this is done

1. use interpolate_variables7.f90 to generate *[23]D.th.nc for the small-domain run;
2. use the new *[23]D.th.nc as inputs for the small-domain run.

Note that `interpolate_variables.in` in the directory are sample inputs for the script. A common mistake is that the parent elements of some open boundary nodes in `fg.gr3` (i.e. ‘small-domain’ hgrid) become dry and the script then fails to find a wet element to interpolate from. So make sure all open boundary nodes in `fg.gr3` are located in wet region of `bg.gr3`; this is especially important for those nodes near the coast. You can use xmgredit5 to ascertain this. If necessary, modify `fg.gr3` by moving some nodes to deeper depths.

## Extraction
`read_output*.f90`: This group of scripts read multiple nc4 outputs and extract time series of a point, a slab, a transect etc. They share similar code structure and can be used to understand the nc4 output format as well as how to do your own processing. You may start from `read_output*_xyz.f90`. After you are familiar with these scripts, you can easily customize them for your own purpose.