message(STATUS "In /Utility NetCDF_LIBS ${NetCDFLIBS}")
add_definitions(${C_PREPROCESS_FLAG})
add_executable( combine_gr3 combine_gr3.f90 )
add_executable( combine_outHA combine_outHA.f )
add_executable( combine_hotstart7 combine_hotstart7.f90 )
add_executable( combine_output10 combine_output10.f90 netcdf_var_names.f90)
add_executable( combine_output11 combine_output11.f90 netcdf_var_names.f90)

target_link_libraries (combine_hotstart7 utillib ${NetCDFLIBS} ${HDF5_LIBRARIES})
target_link_libraries (combine_output10 utillib ${NetCDFLIBS} ${HDF5_LIBRARIES} )
target_link_libraries (combine_output11 utillib ${NetCDFLIBS} ${HDF5_LIBRARIES} )


add_dependencies(utility combine_gr3 combine_outHA combine_hotstart7 combine_output10 combine_output11)

