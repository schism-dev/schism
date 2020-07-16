###W&M Whirlwind cluster
#  Set the base name of the executable.
#  The main reason for this is to include something like a cluster/architecture name.
#  Do not add the file extension (none for linux, .exe for Windows etc)
#  or the list of enabled modules, both of which will be automatically appended.
  set (SCHISM_EXE_BASENAME pschism_GNU_WW CACHE STRING "Base name (modules and file extension to be added of the executable. If you want a machine name, add it here")


###Relative paths won't work
set(CMAKE_Fortran_COMPILER gfortran CACHE PATH "Path to serial Fortran compiler")
set(CMAKE_C_COMPILER gcc CACHE PATH "Path to serial C compiler")
set(NetCDF_FORTRAN_DIR "$ENV{NETCDF_FORTRAN}" CACHE PATH "Path to NetCDF Fortran library")
set(NetCDF_C_DIR  "$ENV{NETCDF}"  CACHE PATH "Path to NetCDF C library")

###Compile flags
#CMAKE_EXE_LINKER_FLAGS did not work so I had to remove -static
set(CMAKE_Fortran_FLAGS_RELEASE "-O2 -ffree-line-length-none -static-libgfortran -finit-local-zero" CACHE STRING "Fortran flags" FORCE)
