All .gr3 and .prop inputs can be visualized/generated using xmgredit5. For other inputs, you can find some useful pre-proc tools in the src/Utility directory. In general, all FORTRAN scripts have a header that explains its purpose, how to use it (inputs/outputs) and sample compilation command.

## Grid conversion
`2dm2gr3_m2m.pl` and `grd2sms.pl` : These 2 perl scripts are used to convert between .2dm (SMS) and .gr3.

## Interpolation
`interpolate_unstructured.f90`: This efficient script can be used to interpolate depths from .gr3 onto another 
.gr3 quickly. It uses a bucket search algorithm along either x or y direction.

## $LSC^2$ scripts
`gen_vqs.f90` and `plot_VQS.m`: The FORTRAN script can be used as a start for creating a LSC2 grid and you need to use the matlab script to plot out transects to see if the vertical grid makes sense. You may lump multiple transects into 1 transect.bp. If you want to use this type of vgrid, make sure you go through the FORTRAN carefully and understand the details.

## Nudging scripts
`gen_nudge.f90`, `gen_nudge2.f90`: the two scripts generate either a 
simple elliptic nudging zone or a zone with fixed distance from boundary as *_nudge.gr3. 
`gen_nudge_from_hycom.f90`: This
script generates the actual nudge data for tracers `*_nu.nc` from HYCOM (you may modify this 
for other gridded data sources from other structured-grid models).

## Hotstart
`Gen_Hotstart/change_hotstart4.f90`: This simple script shows you the internal structure of hotstart.nc and how to manipulate it.

`Gen_Hotstart/gen_*_from_hycom.f90`: These scripts show you how to create hotstart.nc and *.th.nc from gridded outputs 
 like HYCOM. It essentially consists of 3D interpolations.

## Sflux_nc
The matlab scripts inside Sflux_nc dir show you the structure of sflux*.nc as well as how to generate your own files.
