# -*- mode: cmake -*-

#
# MSTK NetCDF Find Module
# Shamelessly stolen from Amanzi open source code https://software.lanl.gov/ascem/trac
#
# Usage:
#    Control the search through NetCDF_DIR or setting environment variable
#    NetCDF_ROOT to the NetCDF installation prefix.
#
#    This module does not search default paths! 
#
#    Following variables are set:
#    NetCDF_FOUND            (BOOL)       Flag indicating if NetCDF was found
#    NetCDF_INCLUDE_DIR      (PATH)       Path to the NetCDF include file
#    NetCDF_INCLUDE_DIRS     (LIST)       List of all required include files
#    NetCDF_LIBRARY_DIR      (PATH)       Path to the NetCDF library
#    NetCDF_LIBRARY          (FILE)       NetCDF library
#    NetCDF_LIBRARIES        (LIST)       List of all required NetCDF libraries
#
#    Additional variables set
#    NetCDF_C_LIBRARY_DIR    (PATH)       Path to the NetCDF C library
#    NetCDF_C_LIBRARY        (FILE)       NetCDF C library
#    NetCDF_CXX_LIBRARY      (FILE)       NetCDF C++ library
#    NetCDF_LARGE_DIMS       (BOOL)       Checks the header files for size of 
#                                          NC_MAX_DIMS, NC_MAX_VARS and NC_MAX_VARS_DIMS
#                                          Returns TRUE if
#                                          NC_MAX_DIMS >= 655363
#                                          NC_MAX_VARS >= 524288
#                                          NC_MAX_VAR_DIMS >= 8
#
# #############################################################################

# Standard CMake modules see CMAKE_ROOT/Modules
include(FindPackageHandleStandardArgs)

# MSTK CMake functions see <root>/cmake/modules for source
#include(PrintVariable)
#include(AddPackageDependency)
if ( NetCDF_LIBRARIES AND NetCDF_INCLUDE_DIRS )

    # Do nothing. Variables are set. No need to search again
else(NetCDF_LIBRARIES AND NetCDF_INCLUDE_DIRS)
   
    # Cache variables
    if(NetCDF_DIR)
        set(NetCDF_DIR "${NetCDF_DIR}" CACHE PATH "Path to search for NetCDF include and library files")
    endif()

    if(NetCDF_INCLUDE_DIR)
        set(NetCDF_INCLUDE_DIR "${NetCDF_INCLUDE_DIR}" CACHE PATH "Path to search for NetCDF include files")
    endif()

    if(NetCDF_LIBRARY_DIR)
        set(NetCDF_LIBRARY_DIR "${NetCDF_LIBRARY_DIR}" CACHE PATH "Path to search for NetCDF library files")
    endif()

    
    # Search for include files
    # Search order preference:
    #  (1) NetCDF_INCLUDE_DIR - check existence of path AND if the include files exist
    #  (2) NetCDF_DIR/<include>
    #  (3) Default CMake paths See cmake --html-help=out.html file for more information.
    #
    set(netcdf_inc_names "netcdf.inc")
    if (NetCDF_INCLUDE_DIR)

        if (EXISTS "${NetCDF_INCLUDE_DIR}")

            find_path(cdf_test_include_path
                      NAMES ${netcdf_inc_names}
                      HINTS ${NetCDF_INCLUDE_DIR}
                      NO_DEFAULT_PATH)
            if(NOT cdf_test_include_path)
                message(SEND_ERROR "Can not locate ${netcdf_inc_names} in ${NetCDF_INCLUDE_DIR}")
            endif()
            set(NetCDF_INCLUDE_DIR "${cdf_test_include_path}")

        else()
            message(SEND_ERROR "NetCDF_INCLUDE_DIR=${NetCDF_INCLUDE_DIR} does not exist")
            set(NetCDF_INCLUDE_DIR "NetCDF_INCLUDE_DIR-NOTFOUND")
        endif()

    else() 

        set(netcdf_inc_suffixes "include")
        if(NetCDF_DIR)

            if (EXISTS "${NetCDF_DIR}" )

                find_path(NetCDF_INCLUDE_DIR
                          NAMES ${netcdf_inc_names}
                          HINTS ${NetCDF_DIR}
                          PATH_SUFFIXES ${netcdf_inc_suffixes}
                          NO_DEFAULT_PATH)

            else()
                 message(SEND_ERROR "NetCDF_DIR=${NetCDF_DIR} does not exist")
                 set(NetCDF_INCLUDE_DIR "NetCDF_INCLUDE_DIR-NOTFOUND")
            endif()    


        elseif(NetCDF_FORTRAN_DIR)

            if (EXISTS "${NetCDF_FORTRAN_DIR}" )

                find_path(NetCDF_INCLUDE_DIR
                          NAMES ${netcdf_inc_names}
                          HINTS ${NetCDF_FORTRAN_DIR}
                          PATH_SUFFIXES ${netcdf_inc_suffixes}
                          NO_DEFAULT_PATH)

            else()
                 message(SEND_ERROR "NetCDF_FORTRAN_DIR=${NetCDF_FORTRAN_DIR} does not exist")
                 set(NetCDF_INCLUDE_DIR "NetCDF_INCLUDE_DIR-NOTFOUND")
            endif()    

        else()

            find_path(NetCDF_INCLUDE_DIR
                      NAMES ${netcdf_inc_names}
                      PATH_SUFFIXES ${netcdf_inc_suffixes})

        endif()

    endif()


    if ( NOT NetCDF_INCLUDE_DIR )
        message(SEND_ERROR "Can not locate NetCDF include directory")
    endif()

    # Large dimension check here
    if ( NetCDF_INCLUDE_DIR ) 
       
        set(netcdf_h "${NetCDF_INCLUDE_DIR}/netcdf.inc" )
        message(STATUS "NetCDF include file ${netcdf_h} will be searched for define values")

        file(STRINGS "${netcdf_h}" netcdf_max_dims_string REGEX "^#define NC_MAX_DIMS")
        string(REGEX REPLACE "[^0-9]" "" netcdf_max_dims "${netcdf_max_dims_string}")

        file(STRINGS "${netcdf_h}" netcdf_max_vars_string REGEX "^#define NC_MAX_VARS")
        string(REGEX REPLACE "[^0-9]" "" netcdf_max_vars "${netcdf_max_vars_string}")

        file(STRINGS "${netcdf_h}" netcdf_max_var_dims_string REGEX "^#define NC_MAX_VAR_DIMS")
        string(REGEX REPLACE "[^0-9]" "" netcdf_max_var_dims "${netcdf_max_var_dims_string}")

        #if ( 
        #     ( (netcdf_max_dims EQUAL 65536)  OR (netcdf_max_dims GREATER 65536) ) AND
        #     ( (netcdf_max_vars EQUAL 524288) OR (netcdf_max_vars GREATER 524288) ) AND
        #     ( (netcdf_max_var_dims EQUAL 8)  OR  (netcdf_max_var_dims GREATER 8)  )

        #            )
        #    set(NetCDF_LARGE_DIMS TRUE)
        #else()
        #    message(WARNING "The NetCDF found in ${NetCDF_DIR} does not have the correct NC_MAX_DIMS, NC_MAX_VARS and NC_MAX_VAR_DIMS\n"
        #                     "It may not be compatible with other TPL libraries such MOAB and ExodusII\n" )
        #    set(NetCDF_LARGE_DIMS FALSE)
        #endif()

    endif()    

    # Search for libraries 
    # Search order preference:
    #  (1) NetCDF_LIBRARY_DIR - check existence of path AND if the include files exist
    #  (2) NetCDF_DIR/<lib,Lib>
    #  (3) Default CMake paths See cmake --html-help=out.html file for more information.
    #
    if (NetCDF_LIBRARY_DIR)

        if (EXISTS "${NetCDF_LIBRARY_DIR}")
            if (NOT EXISTS "${NetCDF_C_LIBRARY_DIR}")
              set(NetCDF_C_LIBRARY_DIR "${NetCDF_LIBRARY_DIR}")
            endif()
            find_library(NetCDF_C_LIBRARY
                         NAMES netcdf
                         HINTS ${NetCDF_C_LIBRARY_DIR}
                         NO_DEFAULT_PATH)

            find_library(NetCDF_Fortran_LIBRARY
                         NAMES netcdff
                         HINTS ${NetCDF_FORTRAN_LIBRARY_DIR}
                         NO_DEFAULT_PATH)

        else()
            message(SEND_ERROR "NetCDF_LIBRARY_DIR=${NetCDF_LIBRARY_DIR} does not exist")
            set(NetCDF_LIBRARY "NetCDF_C_LIBRARY-NOTFOUND")
#            set(NetCDF_LIBRARY "NetCDF_CXX_LIBRARY-NOTFOUND")
        endif()

    else() 

        if(NetCDF_DIR)

            if (EXISTS "${NetCDF_DIR}" )

                find_library(NetCDF_C_LIBRARY
                             NAMES netcdf
                             HINTS ${NetCDF_DIR}
                             PATH_SUFFIXES "lib" "Lib"
                             NO_DEFAULT_PATH)

                find_library(NetCDF_Fortran_LIBRARY
                             NAMES netcdff
                             HINTS ${NetCDF_DIR}
                             PATH_SUFFIXES "lib" "Lib"
                             NO_DEFAULT_PATH)



#                find_library(NetCDF_CXX_LIBRARY
#                             NAMES netcdf_c++
#                             HINTS ${NetCDF_DIR}
#                             PATH_SUFFIXES "lib" "Lib"
#                             NO_DEFAULT_PATH)

            else()
                 message(SEND_ERROR "NetCDF_DIR=${NetCDF_DIR} does not exist")
                 set(NetCDF_LIBRARY "NetCDF_C_LIBRARY-NOTFOUND")
                 set(NetCDF_LIBRARY "NetCDF_FORTRAN_LIBRARY-NOTFOUND")
            endif()    

        else()
          if(NetCDF_C_DIR)
            if (EXISTS "${NetCDF_C_DIR}")
                find_library(NetCDF_C_LIBRARY
                             NAMES netcdf
                             HINTS ${NetCDF_C_DIR}
                             PATH_SUFFIXES "lib" "Lib"
                             NO_DEFAULT_PATH)
            else()
                 message(SEND_ERROR "NetCDF_C_DIR=${NetCDF_C_DIR} does not exist")
                 set(NetCDF_LIBRARY "NetCDF_C_LIBRARY-NOTFOUND")
            endif()
          endif()
              
          if(NetCDF_FORTRAN_DIR)
            if (EXISTS "${NetCDF_FORTRAN_DIR}")
                find_library(NetCDF_Fortran_LIBRARY
                             NAMES netcdff
                             HINTS ${NetCDF_FORTRAN_DIR}
                             PATH_SUFFIXES "lib" "Lib"
                             NO_DEFAULT_PATH)
            else()
                 message(SEND_ERROR "NetCDF_FORTRAN_DIR=${NetCDF_FORTRAN_DIR} does not exist")
                 set(NetCDF_LIBRARY "NetCDF_FORTRAN_LIBRARY-NOTFOUND")
            endif()    
          endif()

           #  find_library(NetCDF_C_LIBRARY
           #               NAMES netcdf
           #               PATH_SUFFIXES ${netcdf_lib_suffixes})

           #  find_library(NetCDF_Fortran_LIBRARY
           #               NAMES netcdff
           #               PATH_SUFFIXES ${netcdf_lib_suffixes})
           # endif()

        endif()

    endif()

    if ( NOT NetCDF_C_LIBRARY )
        message(SEND_ERROR "Can not locate NetCDF C library")
    endif()    

    if ( NOT NetCDF_Fortran_LIBRARY )
        message(WARNING "\n********\nCan not locate separate NetCDF Fortran library (netcdff) in your NetCDF Installation. \nFor older versions of NetCDF the fortran library was not separate (everything was in the netcdf library) and this is not a problem. \nIf you experience a lot of linker errors for symbols starting with 'nf_' , lack of netcdff is a likely cause.")
    endif()    

    
#    if ( NOT NetCDF_CXX_LIBRARY )
#        message(SEND_ERROR "Can not locate NetCDF CXX library")
#    endif()    


   
    # Define the LIBRARIES and INCLUDE_DORS
    set(NetCDF_INCLUDE_DIRS ${NetCDF_INCLUDE_DIR})
    set(NetCDF_LIBRARIES    ${NetCDF_C_LIBRARY} ${NetCDF_Fortran_LIBRARY})

   

    # Need to find the NetCDF config script to check for HDF5
    if ( NetCDF_DIR OR NetCDF_BIN_DIR )
        
        find_program(netcdf_config nc-config
                       HINTS ${NetCDF_DIR} ${NetCDF_BIN_DIR}
                       PATH_SUFFIXES bin Bin
                       DOC "NetCDF configuration script")

        if (netcdf_config AND (NOT (${CMAKE_SYSTEM_NAME} MATCHES "Windows")))
            message(STATUS "Found NetCDF configuration script: ${netcdf_config}")
            execute_process(COMMAND "${netcdf_config}" "--has-hdf5"
                            RESULT_VARIABLE _ret_code
                            OUTPUT_VARIABLE _stdout
                            ERROR_VARIABLE  _stderr
                           )
            string(REGEX REPLACE "[\n\r ]" "" _hdf5_answer ${_stdout})
            message(STATUS "${netcdf_config} --has-hdf5 returned '${_hdf5_answer}'")
            string(COMPARE EQUAL "${_hdf5_answer}" "yes" _has_hdf5)
            if (${_has_hdf5} ) 
                set(NetCDF_NEEDS_HDF5 True)
            else()
                message(STATUS "NetCDF does not require HDF5")
            endif()   

            execute_process(COMMAND "${netcdf_config}" "--version"
                  RESULT_VARIABLE _ret_code
                  OUTPUT_VARIABLE _stdout
                  ERROR_VARIABLE  _stderr
                  )

            string(REGEX REPLACE "[\n\r ]" "" _netcdf_version ${_stdout})
            string(REGEX REPLACE "netCDF" "" _netcdf_version ${_netcdf_version})
            message(STATUS "netcdf version: ${_netcdf_version}")
            if(${_netcdf_version} VERSION_GREATER "4")
               set( NETCDF_4 TRUE)
                message(STATUS "support nc4:${NETCDF_4}")
            endif() 
        else()
            if (NOT DEFINED NetCDF_NEEDS_HDF5)
                message(SEND_ERROR "netcdf_config not available, NetCDF_NEEDS_HDF5 must be set manually")
            endif()
        endif()
    endif()    

    if(NetCDF_NEEDS_HDF5) 
        message(STATUS "NetCDF requires HDF5")
        #add_package_dependency(NetCDF DEPENDS_ON HDF5)
    endif()

endif(NetCDF_LIBRARIES AND NetCDF_INCLUDE_DIRS )    

# Send useful message if everything is found
find_package_handle_standard_args(NetCDF DEFAULT_MSG
                                           NetCDF_LIBRARIES
                                           NetCDF_INCLUDE_DIRS)

# find_package)handle)standard_args should set NetCDF_FOUND but it does not!
if ( NetCDF_LIBRARIES AND NetCDF_INCLUDE_DIRS)
    set(NetCDF_FOUND TRUE)
else()
    set(NetCDF_FOUND FALSE)
endif()

mark_as_advanced(
  NetCDF_INCLUDE_DIR
  NetCDF_INCLUDE_DIRS
  NetCDF_C_LIBRARY
  NetCDF_CXX_LIBRARY
  NetCDF_LIBRARIES
  NetCDF_LIBRARY_DIR
)

