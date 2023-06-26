# - Try to find GOTM
# Once done this will define
#  GOTM_FOUND - System has GOTM
#  GOTM_INCLUDE_DIRS - The GOTM include directories
#  GOTM_LIBRARIES - The libraries needed to use GOTM
#  GOTM_DEFINITIONS - Compiler switches required for using GOTM


# This is not a long term solution to matching the compiler to GOTM compiler, but it will work for now

message(STATUS "GOTM search in ${GOTM_DIR}")

find_path(GOTM_INCLUDE_DIR turbulence.mod
          HINTS "${GOTM_DIR}/include" )

find_library(GOTM_TURB_LIBRARY NAMES turbulence HINTS ${GOTM_DIR}/lib )
find_library(GOTM_UTIL_LIBRARY NAMES util HINTS ${GOTM_DIR}/lib )

set(GOTM_LIBRARY ${GOTM_TURB_LIBRARY} ${GOTM_UTIL_LIBRARY} )
set(GOTM_LIBRARIES ${GOTM_LIBRARY} )
set(GOTM_INCLUDE_DIRS ${GOTM_INCLUDE_DIR} )

include(FindPackageHandleStandardArgs)
# handle the QUIETLY and REQUIRED arguments and set GOTM_FOUND to TRUE
# if all listed variables are TRUE
find_package_handle_standard_args(GOTM  DEFAULT_MSG
                                  GOTM_LIBRARY GOTM_INCLUDE_DIR)

mark_as_advanced(GOTM_INCLUDE_DIR GOTM_LIBRARY )
