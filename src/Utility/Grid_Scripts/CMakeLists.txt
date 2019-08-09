

add_executable( grid_subsampler grid_subsampler.f90 )
#target_link_libraries (combine_output6 ${NetCDFLIBS} ${HDF5_LIBRARIES})
#mpi_wrap( combine_output6 combine_output5 combine_output4 )
add_dependencies(utility grid_subsampler)

