
add_definitions(${C_PREPROCESS_FLAG})
add_executable( compute_alongchannel_vel compute_alongchannel_vel.f90)
add_executable( compute_average3 compute_average3.f90)
add_executable( compute_fluxes_ns compute_fluxes_ns.f90)
add_executable( extract_subregion2 extract_subregion2.f90)
add_executable( read_output8_allnodes read_output8_allnodes.f90) 
add_executable( read_output9_transect read_output9_transect.f90) 
add_executable( read_output9_xyt read_output9_xyt.f90) 
add_executable( read_output9_xyz read_output9_xyz.f90) 
add_executable( read_output9_xyzt read_output9_xyzt.f90) 


target_link_libraries (compute_alongchannel_vel utillib ${NetCDFLIBS} ${HDF5_LIBRARIES} )
target_link_libraries (compute_average3 utillib ${NetCDFLIBS} ${HDF5_LIBRARIES} )
target_link_libraries (compute_fluxes_ns utillib ${NetCDFLIBS} ${HDF5_LIBRARIES})
target_link_libraries (extract_subregion2 utillib ${NetCDFLIBS} ${HDF5_LIBRARIES})
target_link_libraries (read_output8_allnodes utillib ${NetCDFLIBS} ${HDF5_LIBRARIES})
target_link_libraries (read_output9_transect utillib ${NetCDFLIBS} ${HDF5_LIBRARIES})
target_link_libraries (read_output9_xyt utillib ${NetCDFLIBS} ${HDF5_LIBRARIES})
target_link_libraries (read_output9_xyz utillib ${NetCDFLIBS} ${HDF5_LIBRARIES})
target_link_libraries (read_output9_xyzt utillib ${NetCDFLIBS} ${HDF5_LIBRARIES})

add_dependencies(utility compute_alongchannel_vel compute_average3 compute_fluxes_ns extract_subregion2 read_output8_allnodes read_output9_transect read_output9_xyt read_output9_xyz read_output9_xyzt)

