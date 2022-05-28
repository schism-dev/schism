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
 other command then `mpirun` or more arguments.

!!!caution "Scribe IO"
    Shortly after v5.9.0 we have implemented a new I/O mode called scribed I/O. Under this mode, 
the user needs to specify at runtime how many 'scribe' cores they want to use. The # 
of scribes= # of 3D outputs (vectors counted as 2) plus 1 (which is used for all 2D variables). 
Each 3D output has its own netcdf output (e.g. `salinity_*.nc`) and all 2D outputs share same output 
`out2d_*.nc`. Each vector output is splitted into X,Y components, e.g. `horizontalVel[X,Y]_*.nc`.
The output `zCoordinates_*.nc` is needed for VisIT. Efficient asynchronous message passing is 
done inside the code for I/O that minimizes latency.
    Users must specify # of scribes on cmd line as `mpirun -np NPROC pschism nscribe`; the 
specified number can be >= min required as in `param.nml` and explained above. If not 
you'll get an error. If you specify more than needed, you waste some cores but otherwise fine.
