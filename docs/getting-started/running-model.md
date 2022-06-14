To simulate the model, the model executable needs to be run inside the folder where all the model files reside. The following inputs are minimum - 

1. hgrid.gr3
2. vgrid.in
3. param.nml
4. [bottom_friction].gr3
5. bctides.in

The model is usually run through a batch script, which in essence executes the code like:

- `mpirun -np NPROC pschism  <# scribes>` if OLDIO is OFF 
- `mpirun -np NPROC pschism ` if OLDIO is ON

where NPROC is the number of process used for parallel computing. Note that your system may require 
 other commands than `mpirun` or more arguments.

!!!caution "Scribed IO"
    Shortly after [v5.9.0](https://github.com/schism-dev/schism/commit/8efc374) we have implemented a new I/O mode called scribed I/O. 
    
Under scribed IO mode, the outputs are combined during the model simulation by dedicated cores for combining. Efficient asynchronous message passing is done inside the code for I/O that minimizes latency.

Some details for using scribed IO mode are following:
    
- The user needs to specify at runtime how many 'scribe' cores they want to use. The # of scribes= # of 3D outputs (vectors counted as 2) plus 1 (which is used for all 2D variables).
- Each 3D output has its own netcdf output (e.g. `salinity_*.nc`) and all 2D outputs share same output `out2d_*.nc`. Each vector output is splitted into X,Y components, e.g. `horizontalVel[X,Y]_*.nc`.
- The outputs `out2d_*.nc` and `zCoordinates_*.nc` are needed for VisIT. 
- Users must specify # of scribes on cmd line as `mpirun -np NPROC pschism nscribe`. The specified number can be >= min required as in `param.nml` and explained above. If not you'll get an error. If you specify more than needed, you waste some cores but otherwise fine.
