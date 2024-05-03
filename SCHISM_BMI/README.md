# SCHISM Building and documentation
The manual may be found on the SCHISM wiki at http://ccrm.vims.edu/schismweb/. Build instructions are described in Chapter 1.

The online documentation can be accessed at https://schism-dev.github.io/schism.


# SCHISM modules required to compile model and BMI implementation 
* cmake >= 3.20
* intel >= 2018
* impi >= 2018
* netcdf/hdf5 >= 4.7

# Export Linux variables for intel compiler executables for CMake to find during its build cycle
* export FC=mpiifort
* export CXX=mpiicpc
* export CC=mpiicc

# Go into SCHISM source code main directory and configure the SCHISM build for BMI shared libraries
1. cd $root/cmake 
2. Insert/replace the following three lines in the file: "SCHISM.local.build"
  * set (OLDIO ON CACHE BOOLEAN "Old nc output (each rank dumps its own data)")
  * set (USE_ATMOS ON CACHE BOOLEAN "Coupling with atmospheric model via ESMF")
  * set (USE_NWM_BMI ON CACHE BOOLEAN "Use NWM BMI for source and some b.c.")
3. Insert the following lines within CMakeLists.txt in ../src directory to enable BMI compliancy for SCHISM
  ########### SCHISM BMI ADDITIONAL FLAGS TO ENABLE PARMETIS/METIS SHARED LIBRARIES ##########
  * set (CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -fPIC")
  * set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fPIC")
  * set (CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -fPIC")
4. Build the SCHISM model engine:
  * mkdir $root/build; cd $root/build
  * cmake -C ../cmake/SCHISM.local.build ../src
  * make VERBOSE=1 pschism
    Once the SCHISM model sucessfully compiles it's model engine, copy over include and lib SCHISM directories listed in your build directory to the SCHISM BMI directory where we will link these shared libraries respectively
  * scp -r lib ../SCHISM_BMI/SCHISM_LIB_NWM_BMI
  * scp -r include ../SCHISM_BMI/SCHISM_LIB_NWM_BMI

# Compile SCHISM BMI shared libraries and BMI driver that will allow the user to test the SCHISM BMI functionality for a given simulation
1. cd ../SCHISM_BMI/
2. tar -xvf ParMetis-4.0.3.tar.gz ; iso_c_fortran_bmi.tar.gz ; metis-5.1.0.tar.gz
3. mkdir build ; cd build
4. cmake ../
5. cmake --build . --target testbmifortranmodel
6. cmake --build . --target schism_driver
7. Once these steps are completed, you should see that a "schism_driver" executable has been created. This is linked to the SCHISM BMI shared libraries and should allow you to evaluate a SCHISM model using BMI functionality.
8. Copy over the "schism_driver" executable to your SCHISM model directory to evaluate the SCHISM BMI functionality of the code base in a standalone test. Assuming the program executes smoothly, the SCHISM BMI should be ready to be utilized within a BMI coupling framework.

# Current Limitations to SCHISM BMI and Future Work
1. MPI parallelization is currently not developed for the SCHISM BMI. We are underway currently with that enhancement.
2. Water level and source/sink connectors have not been fully validated yet for SCHISM.
3. SCHISM BMI source code here has not yet been linked to main SCHISM-dev repository. Will be linking this code branch in the coming months once development and testing of the BMI has been finsihed. 
