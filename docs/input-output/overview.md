## Type of inputs for SCHISM
Most SCHISM inputs can be visualized with ACE/xmgredit5 or xmgr5 tools; see src/Utility/ACE for instructions on how to install these tools. Other useful tools will be mentioned in the later chapters.

SCHISM input files can be broadly categorized into following 7 groups:

1. `*.gr3`, `hgrid.ll`: node centered spatial data and mesh connectivity. These file can be visualized using `ACE/xmgredit5`;
2. `*.th`: time history files in ASCII format. The ASCII files can be visualized using `ACE/xmgr5`;
3. `*.ic`: initial condition files. Some of these files use `.gr3` format, while others have simple ASCII format;
4. `*.prop`: element-centered spatial data and properties; can be visualized using `ACE/xmgredit5`;
5. `*.nc`: netcdf4 inputs, including time history (`*.th.nc`), hotstart (`hotstart.nc`), and nudging inputs (`*_nu.nc`);
6. `*.nml`: main parameter input (`param.nml`);
7. `*.in`: role-specific input files with individual formats. ASCII inputs include vertical grid (`vgrid.in`), B.C. input (`bctides.in`), and `hydraulics.in` (for hydraulics module) etc;
8. `sflux/`: atmospheric and heat flux files in netcdf format (CF convention v1.0). These files can be visualized using standard tools like ncview, ferret etc;
9. Inputs from modules: `.nml`, `.inp` etc.

## Mandatory inputs
These inputs are required for all SCHISM simulations:

1. Horizontal grid (`hgrid.gr3`)
2. Vertical grid (`vgrid.in`)
3. Parameter input (`param.nml`)
4. B.C. input (`bctides.in`)
5. Bottom friction input (`drag.gr3`, or `rough.gr3` or `manning.gr3`)

We’ll explain these inputs in detail below. Comments/explanations are usually preceded by '!'.
