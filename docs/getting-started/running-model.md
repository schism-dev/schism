To simulate the model, the model executable needs to be run inside the folder where all the model files reside. The following inputs are minimum - 

1. hgrid.gr3
2. vgrid.in
3. param.nml
4. [bottom_friction].gr3
5. bctides.in

The model is run through `mpirun -np NPROC pschism`, where NPROC is the number of process used for parallel computing.

!!!caution Scribe IO
    Shortly after v5.9.0 we have implemented a new I/O mode called scribed I/O. Under this mode, the user needs to specify at runtime how many 'scribe' cores they want to use. The # of scribes= # of 3D outputs (vectors counted as 2) plus 2. Each 3D output has its own netcdf output (e.g. salinity_*.nc) and all 2D outputs share same output out2d_*.nc. The code always outputs zCoordinates_*.nc for VisIT. Efficient message passing is done inside the code for I/O that minimizes latency.

    Users must specify # of scribes on cmd line as `mpirun -np NPROC pschism nscribe`; the specified number can be >= min required as in `param.nml` and explained above. If not you get an error. If you specify more than needed, you waste some cores but otherwise fine.