All SCHISM outputs (except system outputs) can be found in `outputs/` directory.

## Run info output (mirror.out)
This is a mirror image of now-defunct screen output. Below is a sample:

```
Barotropic model without ST calculation
 # of tracers in each module:            1            1            0
            0            0            0            0            0            0      0
 Total # of tracers=            2
 Index ranges of each module:            1            1            2
            2            3            2            3            2            3
            2            3            2            3            2            3
            2            3            2            3            2
 # of global outputs=           27
 done reading param.in; s2_mxnbt in param.in =    3.000000000000000
 lhas_quad=  T
 mnei, mnei_p =             4            9
 lhas_quad=  T
Global Grid Size (ne,np,ns,nvrt):        108       130       237         2

**********Augmented Subdomain Sizes**********
 rank     nea      ne     neg     nea2     neg2     npa      np     npg     npa2     npg2     nsa      ns     nsg     nsa2     nsg2
    0      30      14      16      40      10      43      24      19      43       0      72      37      35      72       0
    1      28      14      14      38      10      40      23      17      40       0      67      36      31      67       0
    2      26      13      13      32       6      38      23      15      38       0      63      35      28      63       0
    3      23      13      10      30       7      34      22      12      34       0      56      34      22      56       0
    4      23      13      10      30       7      34      22      12      34       0      56      34      22      56       0
    5      28      14      14      38      10      40      23      17      40       0      67      36      31      67       0
    6      26      13      13      32       6      38      23      15      38       0      63      35      28      63       0
    7      30      14      16      40      10      43      24      19      43       0      72      37      35      72       0

**********Global Boundary Sizes**********
    nope    neta   nland    nvel
       1      13       1      31

**********Augmented Subdomain Boundary Sizes**********
    rank    nope    neta   nland    nvel
       0       0       0       1       7
       1       0       0       1       7
       2       0       0       1      10
       3       0       0       1      10
       4       1       4       1       7
       5       1       7       0       0
       6       1       5       1       6
       7       1       7       0       0
Max. & min. sidelength=     19934.91537849107         7973.938973963872
 done init (1)...
 done init. tracers..
 done initializing cold start
Done initializing variables
Done initializing outputs
 done computing initial vgrid...
 done computing initial nodal vel...
 done computing initial density...
 time stepping begins...            1         1440
 done adjusting wind stress ...
 done flow b.c.
 done hvis...
 done backtracking
 done 1st preparation
 done 2nd preparation
 done solver; etatot=   3.1245774014922456E-002 ; average |eta|=  4.553235604713400E-005
 done solving momentum eq...
 done solving w
 done solving transport equation
 done recomputing levels...
 done density calculation...
TIME STEP=            1;  TIME=           300.000000
….
```

!!!notes
    The ‘average |eta|’ above can be used as a quick and easy way to check if the run is progressing smoothly; it is the average of the absolute value of surface elevation at all nodes. If it’s too large or NaN, you have a problem.

## Global output
SCHISM netcdf4 outputs are emitted in a directory called outputs/. This directory must exist or you 
will get an immediate crash from the model. 

Depending on whether or not you turned on OLDIO, the global netcdf outputs will look different.

1) Scribed I/O (OLDIO is OFF)
Under this mode, the netcdf outputs are global (combined) outputs, and you can visualize or process thm
 using latest FORTRAN (e.g., read_output10*), matlab or python scripts. 
For example the latest VisIT plugins can visualize these outputs directly.

All 2D variables (e.g. `elevation`, `sigWaveHeight` etc) as well as static information such as geometry
 and connectivity info are grouped into `out2d_*.nc`. On the other hand, each 3D variable has its own 
 output, and vector variables have X and  Y components in separate outputs (e.g. `horizontalVelX_*.nc` and
 `horizontalVelY_*.nc`).

2) Old I/O (OLDIO is ON)
Under this mode, each MPI process will dump its own output and a post-processing script (`combin_output*.f90`)
 will need to be used to combine these into global netcdf outputs.  

An example output file name is `outputs/schout_000000_2.nc`. More generally, the file name is: `schout_[processor_no]_[time_block /stack #].nc`

!!!notes "Processor number"
    The mpi_processor number starts at 0 and represents the MPI processor ID from the task that wrote the output. 

!!!notes "Time block" 
    The time blocks (‘stack’) start from 1 and are sequential. The model buffers and writes data occasionally. Every `ihfskip` time steps it opens a new stack. For instance, if the time step is 120 seconds and `ihfskip = 10080`, each stack will be 14-day long.

    - "Neat" time lengths that will make meaningful analysis (e.g. daily, 10 days etc) are usually easiest later when you post-process;
    - Some of the post-processing scripts will run a lot better if the length of your simulation is an even multiple of ihfskip. This can be done by altering ihfskip or the simulation length - at the risk of lengthening the simulation a bit, the latter often produces a neater result.
    - If your simulation length is not an even multiple of the time block length, the last time block will be truncated on the last block. This will cause some minor errors and warnings in the post-processing tools. In addition, if you then restart the run it is best to repeat and overwrite the truncated block - the post-processing tools do not work well with blocks that grow and shrink in the middle of the run.
    - If the output blocks match the end of the simulation very neatly, the model (at the time of writing) will open a new block that is very small in size. This is useful for the `autocombine_MPI_elfe.pl`, as the latter always waits until a new block to come out before starting to combine the previous block (and so it'd hang if the last empty block were not written out).

!!!notes "Variable names"
    The variables inside .nc correspond to 3D grid and state variables info at each output time step. The state variables (arrays) may have different centerings, e.g. at node/element/side. At the moment, most variables are centered around nodes (and most post-processing FORTRAN scripts also work on node-centered variables).

!!!notes "Combine outputs"
    The per-processor outputs need to be gathered into combined nc4 outputs first before you can visualize or post-process them. The script that does this is called `combine_output11.f90` (a simple perl script `autocombine_MPI_elfe.pl` exists to combine all available outputs transparently; you just need to update the path to the compiled `combine_output11` inside the script, and it can be launched before or after a run is done). See the header of `combine_output11.f90` on sample compilation commands and usage. Once you are done combining, you should have nc4 files called something like `schout_2.nc` etc. Note that the stack # remains but the MPI process number is gone. There is no utility for gathering the outputs in stack/time; instead most post-processing tools are able to work with multiple stacks.

    Note that SCHISM allows users to easily add more customized outputs, using the routine writeout_nc() inside schism_step. The combine scripts will automatically combine the additional outputs.

    You can visualize the combined nc4 outputs using VisIT (with SCHISM plug-ins). More info can be found in [Visualization](./../getting-started/visualization.md)

!!!notes "Other global outputs"
    The user may be interested in some maximum quantities. At the moment, SCHISM outputs two max files for elevation and depth-averaged velocity (`outputs/maxelev_*` and `outputs/maxdahv_*`). These files can be combined using `Utility/Combining_Scripts/combine_gr3.f90` to generate `maxelev.gr3` and `maxdahv.gr3`. Another type of global outputs are hotstart outputs, which should be combined (using `combine_hotstart7.f90`) to generate a restart input.

## Station outputs
These outputs are invoked with `iout_sta=1`, and are found in `outputs/staout_[1..,9]`, corresponding respectively to elev, air pressure, wind u, wind v, T, S, u, v, w. Each output has a simple ASCII format:

```
Time(sec), variable @ station 1,2,…. (specified in station.in)
```

## Warning and fatal messages
Warning message (`nonfatal_*`) contains non-fatal warnings, while fatal message file (`fatal.error`) contains fatal errors. In addition, you’d also check the system error outputs from your parallel job.
