See [Numetical Formulation](/schism/numerical-formulation) chapter for details of different types of vgrid supported in SCHISM. Following are a few example `vgrid.in`.

### An example of SZ grid 

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

### An example of pure S grid
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

### An example of $LSC^2$ grid
This type of grid requires some user experience and can be generated using scripts (e.g., `Utility/Pre-Processing/gen_vqs.f90`).

```
1 !ivcor (1: LSC2; 2: SZ)
39 !nvrt(=Nz)
1 29 -1. -0.893491 -0.793885 -0.694399 -0.595016 -0.495718 -0.396491 -0.297320 -0.198191 -0.099089 0. !node ID, bottom level index, sigma coordinate from bottom (-1) to surface (0)
....
19520 12 -1.0 -0.985926 -0.918190 -0.854562 -0.794749 -0.738477 -0.685488 -0.635540 -0.588405 -0.543870 -0.501733 -0.461804 -0.423903 -0.387860 -0.353515 -0.320714 -0.289311 -0.259169 -0.230153 -0.202137 -0.174997 -0.148616 -0.122878 -0.097671 -0.072887 -0.048417 -0.024156 0.0 !@last node
```