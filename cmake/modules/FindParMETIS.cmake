# This software is licensed under the terms of the Apache Licence Version 2.0
# which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.

# Author: Panagiotis Velissariou <panagiotis.velissariou@noaa.gov>, NOAA
#
# - Aug 25 2022 Panagiotis Velissariou <panagiotis.velissariou@noaa.gov>, NOAA
#     * Find ParMETIS libraries and include files in a basic
#       Linux/GNU/gfortran system
#

# Try to find ParMETIS includes and libraries.
# Supports static and shared libaries.
#
# This module defines
#
#   - ParMETIS_FOUND                - System has ParMETIS
#   - ParMETIS_EXTERNAL_FOUND       - System has ParMETIS
#   - ParMETIS_INCLUDE_DIRS         - the ParMETIS include directories
#   - ParMETIS_LIBRARIES            - the ParMETIS libraries
#   - ParMETIS_VERSION              - the version of ParMETIS
#   - METIS_VERSION                 - the version of METIS
#
#   - PARMETIS_FOUND                - System has ParMETIS
#   - PARMETIS_EXTERNAL_FOUND       - System has ParMETIS
#   - PARMETIS_INCLUDE_DIRS         - the ParMETIS include directories
#   - PARMETIS_LIBRARIES            - the ParMETIS libraries
#   - PARMETIS_VERSION              - the version of ParMETIS
#   - METIS_VERSION                 - the version of METIS
#
#   - Parmetis_FOUND                - System has ParMETIS
#   - Parmetis_EXTERNAL_FOUND       - System has ParMETIS
#   - Parmetis_INCLUDE_DIRS         - the ParMETIS include directories
#   - Parmetis_LIBRARIES            - the ParMETIS libraries
#   - Parmetis_VERSION              - the version of ParMETIS
#   - METIS_VERSION                 - the version of METIS
#
# The following paths will be searched in order if set in CMake (first priority) or environment (second priority)
#
#   - ParMETIS_ROOT                 - root of ParMETIS installation
#   - ParMETIS_PATH                 - root of ParMETIS installation
#   - ParMETIS_HOME                 - root of ParMETIS installation
#   - ParMETISROOT                  - root of ParMETIS installation
#   - ParMETISPATH                  - root of ParMETIS installation
#   - ParMETISHOME                  - root of ParMETIS installation
#
# The search process begins with locating ParMETIS Include headers.  If these are in a non-standard location,
# set one of the following CMake or environment variables to point to the location:
#
#  - ParMETIS_INCLUDE_DIR or PARMETIS_INCLUDE_DIR or Parmetis_INCLUDE_DIR
#  - ParMETIS_INCLUDE_DIRS or PARMETIS_INCLUDE_DIRS or Parmetis_INCLUDE_DIRS
#

list( APPEND pmetis_search_components ParMETIS METIS GKlib )

unset(_search_components)
unset(_new_search_components)
unset(_search_hints)
unset(pmetis_include_search_hints)

## Include names for each component
set( ParMETIS_INCLUDE_NAME  parmetis.h )
set( METIS_INCLUDE_NAME     metis.h )
set( GKlib_INCLUDE_NAME     GKlib.h )

## Library names for each component
set( ParMETIS_LIBRARY_NAME  parmetis )
set( METIS_LIBRARY_NAME     metis )
set( GKlib_LIBRARY_NAME     GKlib )


## Search hints for finding include directories and libraries
foreach( _comp IN ITEMS "_" "" )
  foreach( _name IN ITEMS PARMETIS ParMETIS Parmetis METIS Metis GKLIB GKlib)
    foreach( _var IN ITEMS HOME ROOT PATH)
        get_filename_component(_tmp_val "${${_name}${_comp}${_var}}" ABSOLUTE REALPATH)
      list(APPEND pmetis_search_hints ${_tmp_val})
        get_filename_component(_tmp_val "$ENV{${_name}${_comp}${_var}}" ABSOLUTE REALPATH)
      list(APPEND pmetis_search_hints ${_tmp_val} ${_tmp_val}/.. ${_tmp_val}/../.. )
      list(APPEND pmetis_include_search_hints
                ${${_name}${_comp}INCLUDE_DIR} $ENV{${_name}${_comp}INCLUDE_DIR}
                ${${_name}${_comp}INCLUDE_DIRS} $ENV{${_name}${_comp}INCLUDE_DIRS} )
    endforeach()
  endforeach()
endforeach()
unset(_comp)
unset(_name)
unset(_var)
unset(_tmp_val)
list( REMOVE_DUPLICATES pmetis_search_hints )
list( REMOVE_DUPLICATES pmetis_include_search_hints )


########################################
### BEG:: Find headers for each component
########################################
set(ParMETIS_INCLUDE_DIRS)
set(_new_search_components)
foreach( _comp IN LISTS pmetis_search_components )
  if(NOT ${PROJECT_NAME}_${_comp}_FOUND)
      list(APPEND _new_search_components ${_comp})
  endif()
  find_file(${_comp}_INCLUDE_FILE
    NAMES ${${_comp}_INCLUDE_NAME}
    DOC "ParMETIS ${_comp} include directory"
    HINTS ${pmetis_include_search_hints} ${pmetis_search_hints}
    PATH_SUFFIXES include include/parmetis include/metis
    NO_DEFAULT_PATH
  )
  mark_as_advanced(${_comp}_INCLUDE_FILE)
  message(DEBUG "${_comp}_INCLUDE_FILE: ${${_comp}_INCLUDE_FILE}")
  if( ${_comp}_INCLUDE_FILE )
    get_filename_component(${_comp}_INCLUDE_FILE ${${_comp}_INCLUDE_FILE} ABSOLUTE)
    get_filename_component(${_comp}_INCLUDE_DIR ${${_comp}_INCLUDE_FILE} DIRECTORY)
    list(APPEND ParMETIS_INCLUDE_DIRS ${${_comp}_INCLUDE_DIR})
  endif()
endforeach()

if(ParMETIS_INCLUDE_DIRS)
  list(REMOVE_DUPLICATES ParMETIS_INCLUDE_DIRS)
  set(ParMETIS_INCLUDE_DIRS "${ParMETIS_INCLUDE_DIRS}" CACHE STRING "ParMETIS include paths" FORCE)
endif()
########################################
### END:: Find headers for each component
########################################


########################################
### BEG:: Find libraries for each component
########################################
set( ParMETIS_LIBRARIES )
foreach( _comp IN LISTS pmetis_search_components )
  find_library( ${_comp}_LIBRARY
    NAMES ${${_comp}_LIBRARY_NAME}
    DOC "ParMETIS ${_comp} library"
    HINTS ${${_comp}_INCLUDE_DIRS} ${pmetis_search_hints}
    PATH_SUFFIXES lib64 lib ../lib64 ../lib ../../lib64 ../../lib NO_DEFAULT_PATH )
  mark_as_advanced( ${_comp}_LIBRARY )
  get_filename_component(${_comp}_LIBRARY ${${_comp}_LIBRARY} ABSOLUTE)
  set(${_comp}_LIBRARY ${${_comp}_LIBRARY} CACHE STRING "${_comp} library" FORCE)
  message(DEBUG "${_comp}_LIBRARY: ${${_comp}_LIBRARY}")

  if( ${_comp}_LIBRARY )
    list( APPEND ParMETIS_LIBRARIES ${${_comp}_LIBRARY} )
  endif()
endforeach()

if(ParMETIS_LIBRARIES)
  list(REMOVE_DUPLICATES ParMETIS_LIBRARIES)
  set(ParMETIS_LIBRARIES "${ParMETIS_LIBRARIES}" CACHE STRING "ParMETIS library paths" FORCE)
endif()
########################################
### END:: Find libraries for each component
########################################


########################################
### BEG:: Find version for each component
########################################
if(ParMETIS_INCLUDE_DIRS)
  foreach( _dir IN LISTS ParMETIS_INCLUDE_DIRS)
    if( EXISTS "${_dir}/parmetis.h" )
      file(STRINGS "${_dir}/parmetis.h" pmetis_version_lines
      REGEX "#define[ \t]+PARMETIS_(MAJOR|MINOR|SUBMINOR)_VERSION")
      string(REGEX REPLACE ".*PARMETIS_MAJOR_VERSION *\([0-9]*\).*" "\\1" pmetis_version_major "${pmetis_version_lines}")
      string(REGEX REPLACE ".*PARMETIS_MINOR_VERSION *\([0-9]*\).*" "\\1" pmetis_version_minor "${pmetis_version_lines}")
      string(REGEX REPLACE ".*PARMETIS_SUBMINOR_VERSION *\([0-9]*\).*" "\\1" pmetis_version_subminor "${pmetis_version_lines}")
      set(ParMETIS_VERSION "${pmetis_version_major}.${pmetis_version_minor}.${pmetis_version_subminor}")
      unset(pmetis_version_major)
      unset(pmetis_version_minor)
      unset(pmetis_version_subminor)
      unset(pmetis_version_lines)
    endif()
    if( EXISTS "${_dir}/metis.h" )
      file(STRINGS "${_dir}/metis.h" metis_version_lines
      REGEX "#define[ \t]+METIS_VER_(MAJOR|MINOR|SUBMINOR)")
      string(REGEX REPLACE ".*METIS_VER_MAJOR *\([0-9]*\).*" "\\1" metis_version_major "${metis_version_lines}")
      string(REGEX REPLACE ".*METIS_VER_MINOR *\([0-9]*\).*" "\\1" metis_version_minor "${metis_version_lines}")
      string(REGEX REPLACE ".*METIS_VER_SUBMINOR *\([0-9]*\).*" "\\1" metis_version_subminor "${metis_version_lines}")
      set(METIS_VERSION "${metis_version_major}.${metis_version_minor}.${metis_version_subminor}${metis_version_note}")
      unset(metis_version_major)
      unset(metis_version_minor)
      unset(metis_version_subminor)
      unset(metis_version_lines)
    endif()
  endforeach()
  
  unset(_dir)
endif()
########################################
### END:: Find version for each component
########################################

if(ParMETIS_INCLUDE_DIRS AND ParMETIS_LIBRARIES)
  set(ParMETIS_FOUND TRUE CACHE BOOL "ParMETIS Found" FORCE)
  set(${PROJECT_NAME}_ParMETIS_FOUND TRUE)
else()
  set(ParMETIS_FOUND FALSE CACHE BOOL "ParMETIS Found" FORCE)
  set(${PROJECT_NAME}_ParMETIS_FOUND FALSE)
  unset(ParMETIS_INCLUDE_DIRS CACHE)
  unset(ParMETIS_LIBRARIES CACHE)
endif()

foreach( _prefix PARMETIS ParMETIS Parmetis ${CMAKE_FIND_PACKAGE_NAME} )
  set( ${_prefix}_INCLUDE_DIRS   ${ParMETIS_INCLUDE_DIRS} )
  set( ${_prefix}_LIBRARIES      ${ParMETIS_LIBRARIES})
  set( ${_prefix}_VERSION        ${ParMETIS_VERSION} )
  set( ${_prefix}_FOUND          ${${CMAKE_FIND_PACKAGE_NAME}_FOUND} )
  set( ${_prefix}_EXTERNAL_FOUND ${${CMAKE_FIND_PACKAGE_NAME}_FOUND} )
endforeach()
