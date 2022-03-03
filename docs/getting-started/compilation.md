Your system must have a FORTRAN and C compiler (and MPI wrapper like mpif90, mpicc), netcdf (version 4.4 or newer is required), python (version 2.7 and above) and perl installed.

Two build system is currenlty available - GNU Make and Cmake

## GNU Make
You need the following two directories for gnu make. `src/` where the sourcecode resides and the `Makefile` which is the main makefile. In general this `src/Makefile` should not be changed. In `mk`, there are other makefiles `Make.defs.*` which can be used to provide system specific information. 

You need a `mk/Make.defs.local` file for compilation. To help you design this file, we have put `Make.defs.*` which are from various clusters around the world as well as common linux operating systems. You should copy/link one of these to `Make.defs.local`. For example, if your MPI compiler is Intel based ifort, you may try - 

``` bash
# Copy
cp Make.defs.bora Make.defs.myown

# or link
ln -sf Make.defs.myown Make.defs.local
```

Then you need to edit `Make.defs.local` for e.g. the MPI compiler name, path names for `netcdf` library (v4.4 and above preferred) on your local cluster (sometimes you may also need to add `FORTRAN` support for netcdf, e.g. `–lnetcdff`), and name of the executable (`EXEC`). Also the graph partitioning lib ParMETIS is compiled each time together with the SCHISM code, and you may need to update the MPI C compiler names in `src/ParMetis-*` (consult also `INSTALL` inside the `ParMETIS` directory). Lastly, turn on/off modules in `include_modules` (note that `TVD_LIM` should always have a valid value). Make sure `mk/sfmakedepend.pl` and `mk/cull_depends.py` are executable (otherwise make them executable with `chmod +x`).

After all of these are done - 

``` bash
cd ../src
make clean
make pschism
```

You might get some errors if you did not use git to clone so the make cannot find the hash information; ignore it and proceed in making. The final name of executable has the cluster name and also additional suffixes if you turn on modules: pschism_* (note that you get a preview of this name when you do make clean so you know which modules you have turned on).

Note: the Makefile will automatically invoke `Core/gen_version.py`, which will generate `Core/schism_version.F90`, needed for querying the hash/version.

## Cmake
Alternatively, the cmake utility is a very powerful way of building the source code and utility script bundle. The main cmake files are found in cmake/. You’ll need two essential files in this directory.

1. `SCHISM.local.build`: used to toggle on/off optional modules (similar to include_modules);
2. `SCHISM.local.cluster_name`: similar to Make.defs.local, this file specifies the most important environment variables like name/path of compiler, netcdf library etc. In general, cmake is quite adept at inferring some of these variables but sometimes you might need to overwrite the `defaults` in this file. You can start from an existing file for a similar cluster e.g. `cp –L SCHISM.local.whirlwind SCHISM.local.myown`

Once these two files are set, run - 

```bash
mkdir ../build
cd ../build; rm -rf * # Clean old cache
cmake -C ../cmake/SCHISM.local.build -C ../cmake/SCHISM.local.myown ../src/
```

The cmake is essentially a pre-processor for make, and it creates cache files (e.g. build/ CMakeCache.txt, where you can inspect all env variables). After cmake is done, make can be executed in parallel or in serial mode.

```bash
make -j8 pschism # efficient parallel build, replace 8 with number of process you want
# or
make VERBOSE=1 pschism # serial build with a lot of messages
```

The executable (`pschism*`) is found in `build/bin/` and the compiled libraries are in `build/lib/`. If `pschism` is omitted above, the main executable and all utility scripts will be built in `bin/`.