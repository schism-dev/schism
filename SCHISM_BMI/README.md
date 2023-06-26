Building SCHISM and BMI Driver on the RDHPCS Hera supercomputer

# Import modules from the Hera supercomputer
module load cmake/3.20.1
module load intel/2022.3.0
module load impi/2022.3.0
module load netcdf-hdf5parallel/4.7.4

# Export Linux variables for intel compiler executables for CMake to find during its build cycle
export FC=mpiifort
export CXX=mpiicpc
export CC=mpiicc

# Go into SCHISM source code main directory and configure the SCHISM build for BMI shared libraries and for the Hera supercomputer
1. cd schism_NWM_BMI/cmake 
2. Insert/replace the following three lines in the file: "SCHISM.local.build"
  set (OLDIO ON CACHE BOOLEAN "Old nc output (each rank dumps its own data)")
  set (USE_ATMOS ON CACHE BOOLEAN "Coupling with atmospheric model via ESMF")
  set (USE_NWM_BMI ON CACHE BOOLEAN "Use NWM BMI for source and some b.c.")
3. Insert the following lines within CMakeLists.txt in ../src directory to enable BMI compliancy for SCHISM
  ########### SCHISM BMI ADDITIONAL FLAGS TO ENABLE PARMETIS/METIS SHARED LIBRARIES ##########
  set (CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -fPIC")
  set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fPIC")
  set (CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -fPIC")
4. Modiy the following lines within CMakeLists.txt in ../src directory to enable BMI compliancy for SCHISM
  define_opt(OLDIO "Old nc output option" ON)
  define_opt(USE_ATMOS "Atmospheric model output option" ON)
  define_opt(USE_NWM_BMI "Use NWM BMI for source and some b.c." ON)
6. To build SCHISM system using cmake version 3.12 or higher:
  mkdir ../build
  cd ../build; rm -rf *
  cmake -C ../cmake/SCHISM.local.build
  make VERBOSE=1 pschism
7. Copy over include and lib SCHISM libraries/modules directory to the SCHISM BMI directory where we will link these shared libraries respectively
  scp -r lib ../SCHISM_BMI/SCHISM_LIB_NWM_BMI
  scp -r include ../SCHISM_BMI/SCHISM_LIB_NWM_BMI

# Compile SCHISM BMI shared libraries and BMI driver that will allow the user to test the SCHISM BMI functionality for a given simulation
1. cd ../SCHISM_BMI/
2. mkdir build ; cd build
3. cmake ../
4. cmake --build . --target testbmifortranmodel
5. cmake --build . --target schism_driver
6. Once these steps are completed, you should see that a "schism_driver" executable has been created. This is linked to the SCHISM BMI shared libraries and should allow you to evaluate a SCHISM model using BMI functionality.
