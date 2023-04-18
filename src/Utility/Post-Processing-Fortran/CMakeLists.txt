add_definitions(${C_PREPROCESS_FLAG})

# set list of files
set(post_processing_files
    compute_alongchannel_vel.f90
    compute_average3.f90
    compute_fluxes_ns.f90
    extract_subregion2.f90
    read_output8_allnodes.f90
    read_output9_transect.f90
    read_output9_xyt.f90
    read_output9_xyz.f90
    read_output9_xyzt.f90
    read_output10_allnodes.f90
    read_output10_xyz.f90
    read_output10_transect.f90
    read_output10_xyt.f90
    read_output10_xyzt.f90)

# define executables, their dependencies and libraries that will be used in linking
foreach(f ${post_processing_files})
  get_filename_component(fwe ${f} NAME_WE)
  add_executable(${fwe} ${f})
  target_link_libraries(${fwe} utillib ${NetCDFLIBS} ${HDF5_LIBRARIES})
  target_include_directories(${fwe} PRIVATE ${CMAKE_CURRENT_BINARY_DIR}/../UtilLib)
  add_dependencies(${fwe} utillib)
endforeach()
