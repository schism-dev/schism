
set(schismmodlibs fabm_schism ${schismmodlibs} PARENT_SCOPE)
add_library(fabm_schism fabm_schism.F90)

# instructions taken from GOTM model, such that usage is similar
if(USE_FABM)
  find_path(FABM_BASE src/fabm.F90 DOC "Path to FABM source directory.")
  mark_as_advanced(CLEAR FABM_BASE)
  if(FABM_BASE)
    # Build FABM from source by including its CMakeLists.txt directory.
    set(FABM_FORCED_HOST schism)
    add_subdirectory(${FABM_BASE}/src fabm_src)
    set(FABM_LIBRARIES fabm)
    set(FABM_INCLUDE_DIRS ${CMAKE_CURRENT_BINARY_DIR}/fabm_src/modules)
    mark_as_advanced(FORCE FABM_PREFIX)
  endif()
  set_property(TARGET fabm_schism APPEND PROPERTY INCLUDE_DIRECTORIES "${FABM_INCLUDE_DIRS}")
  #set_property(TARGET fabm_schism PROPERTY Fortran_MODULE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/modules)
  #message(STATUS "fabm_schism dependencies: core ${FABM_LIBRARIES}")
  target_link_libraries(fabm_schism ${FABM_LIBRARIES} core)
endif(USE_FABM)

