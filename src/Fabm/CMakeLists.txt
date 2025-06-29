set(schismmodlibs fabm_schism ${schismmodlibs} PARENT_SCOPE)
add_library(fabm_schism fabm_schism.F90)

# instructions taken from GOTM model, such that usage is similar
if(USE_FABM)
  find_path(FABM_BASE src/fabm.F90 DOC "Path to FABM source directory.")
  mark_as_advanced(CLEAR FABM_BASE)
  if(FABM_BASE)
    message(STATUS "Using FABM_BASE=${FABM_BASE}")
    if(NOT EXISTS ${FABM_BASE}/include/fabm_version.h)
      add_definitions(-D_FABM_API_VERSION_=0)
      message(STATUS "Added definition -D_FABM_API_VERSION_=0")
    endif()
    # Build FABM from source by including its CMakeLists.txt directory.
    set(FABM_FORCED_HOST schism CACHE STRING "Forcing FABM host to schism")
    message(STATUS "Forcing FABM_HOST=${FABM_FORCED_HOST}")
    # SCHISM realkind found in src/Core/schism_glbl.F90
    set(FABM_FORCED_REAL_KIND 8 CACHE STRING "Forcing FABM realkind to 8")
    message(STATUS "Forcing FABM_REAL_KIND=${FABM_FORCED_REAL_KIND}")
    #mark_as_advanced(FABM_FORCED_REAL_KIND)
    add_subdirectory(${FABM_BASE} fabm_src)
    set(FABM_LIBRARIES fabm)
    #message(STATUS "Adding FABM library -l${FABM_LIBRARIES}")
    set(FABM_INCLUDE_DIRS ${CMAKE_CURRENT_BINARY_DIR}/fabm_src/modules 
      ${FABM_BASE}/include ${FABM_BASE}/src/drivers/schism
    )
    #message(STATUS "Adding FABM includes ${FABM_INCLUDE_DIRS}")
    mark_as_advanced(FORCE FABM_PREFIX)
  endif()
  set_property(TARGET fabm_schism APPEND PROPERTY INCLUDE_DIRECTORIES "${FABM_INCLUDE_DIRS}")
  #set_property(TARGET fabm_schism PROPERTY Fortran_MODULE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/modules)
  #message(STATUS "fabm_schism dependencies: core ${FABM_LIBRARIES}")
  target_link_libraries(fabm_schism ${FABM_LIBRARIES} core)
endif(USE_FABM)

