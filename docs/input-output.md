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

We’ll explain these inputs in detail below. Comments/explanations are in **bold**.

### hgrid.gr3
The format of this file is shown below. It has 4 part.

First part of the file is the information - 
```
hgrid.gr3 ! alphanumeric description; ignored by code
60356 31082 ! # of elements and nodes in the horizontal grid
```

Second part is the node info.
```
1 402672.000000 282928.000000 2.0000000e+01 ! node #, x,y, depth
2 402416.000000 283385.000000 2.0000000e+01
3 402289.443000 282708.750000 2.0000000e+01
4 402014.597000 283185.897000 2.0000000e+01
.............................................
31082 331118.598253 112401.547031 2.3000000e-01 !last node
```

Third part is the connectivity table.
```
1 4 1 2 3 101 ! element #, element type (triangle or quad), nodes 1-4
2 3 2 4 3
3 3 4 5 3
...........................................
60356 3 26914 30943 26804 !last element
```

The last part is the list of open and land boundary segments. This part is needed for hgrid.gr3 only; not needed
for other .gr3 files.

```
3 = Number of open boundaries
95 = Total number of open boundary nodes
3 = Number of nodes for open boundary 1
29835 ! first node on this segment
29834 ! 2nd node on this segment
.
.
.
30001 !last node on this segment
90 = Number of nodes for open boundary 2
.
.
.
16 = number of land boundaries (including islands)
1743 = Total number of land boundary nodes
753 0 = Number of nodes for land boundary 1 ('0' means the exterior land boundary)
30381 ! first node on this segment
.......................................
1 !last node on this segment
741 0 ! Number of nodes for land boundary 2 ('0' means the exterior boundary)
.
.
.
10 1 = Number of nodes for island boundary 1 ('1' means island)
29448 ! first node on this island
.
.
.
29449 !last node on this island (note this is different from the first node ‘29448’ above)
```

Note: 
1. B.C. part can be generated with xmgredit5 $\rightarrow$ GridDEM $\rightarrow$ Create open/land boundaries; it can also be generated with SMS;
2. If you have no open boundary, you can create two land boundary segments that are linked to each other. Likewise, if you have no land boundary, you should create two open boundary segments that are connected to each other;
3. Although not required, we recommend you follow the following convention when generating the boundary segments. For the exterior boundary (open+land), go in counter-clockwise direction. With xmgredit5, the island boundaries are automatically created once you have finished designating all open and land segments on the exterior boundary.
4. Note that this format is the same as fort.14 of ADCIRC; Keep an eye on the Land boundary sagment, where instead of `741 0` in the above example SMS will produce `741 10`. 
5. If WWM is used, the land boundary flags (cf. red texts above) are required, and also there must not be any open boundary segments on any island. Since WWM can only handle triangles, the mixed grid needs to be converted to a pure triangular grid for WWM using a pre-processing script.

### vgrid.in
See [Numetical Formulation](/schism/numerical-formulation) chapter for details of different types of vgrid supported in SCHISM. Following are a few example `vgrid.in`.

#### An example of SZ grid 

```
2 !ivcor (1: LSC2; 2: SZ)
54 18 100. !nvrt(=Nz); kz (# of Z-levels); hs (transition depth between S and Z)
Z levels !Z-levels in the lower portion
1 -5000. !level index, z-coordinates
2 -2300.
3 -1800.
4 -1400.
5 -1000.
6 -770.
7 -570.
8 -470.
9 -390.
10 -340.
11 -290.
12 -240.
13 -190.
14 -140.
15 -120.
16 -110.
17 -105.
18 -100. !z-coordinate of the last Z-level must match -h_s
S levels !S-levels below
30. 0.7 10. ! constants used in S-transformation: hc, theta_b, theta_f
18 -1.
 !first S-level (sigma-coordinate must be -1)
19 -0.972222
 !levels index, sigma-coordinate
20 -0.944444
.
.
.
54 0.
 !last sigma-coordinate must be 0
```

Notes:
 - Global outputs are from the bottom (`kbp`, variable in space) to surface (level `nvrt`) at each node;
 - The code will crash if the surface elevation falls below $–h_c$ so make sure $h_c$ is sufficiently large (there is a hardwired lower bound for this around 5m in the code).

#### An example of pure S grid
If a "pure S" model is desired, use only 1 Z-level and set hs to a very large number (e.g., 1.e6)
above. For example, using vgrid.in below leads to a 2D model.

```
2 !ivcor
2 1 1.e6
Z levels
1 -1.e6
S levels
40. 1. 1.e-4
1 -1.
2 0.
```

#### An example of $LSC^2$ grid
This type of grid requires some user experience and can be generated using scripts (e.g., `Utility/Pre-Processing/gen_vqs.f90`).

```
1 !ivcor (1: LSC2; 2: SZ)
39 !nvrt(=Nz)
1 29 -1. -0.893491 -0.793885 -0.694399 -0.595016 -0.495718 -0.396491 -0.297320 -0.198191 -0.099089 0. !node ID, bottom level index, sigma coordinate from bottom (-1) to surface (0)
....
19520 12 -1.0 -0.985926 -0.918190 -0.854562 -0.794749 -0.738477 -0.685488 -0.635540 -0.588405 -0.543870 -0.501733 -0.461804 -0.423903 -0.387860 -0.353515 -0.320714 -0.289311 -0.259169 -0.230153 -0.202137 -0.174997 -0.148616 -0.122878 -0.097671 -0.072887 -0.048417 -0.024156 0.0 !@last node
```

### bctides.in
### param.nml

## Optional inputs
### .gr3, hgrid.ll
### .th (ASCII)
### .th.nc (netcdf4)
### .prop
### .ic
### .nc
### .in
### slufx

## Outputs
### Run infor output (mirror.out)
### Global output
### Station outputs
### Warning and fatal messages