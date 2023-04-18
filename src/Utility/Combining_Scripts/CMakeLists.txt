message(STATUS "In /Utility NetCDF_LIBS ${NetCDFLIBS}")
add_definitions(${C_PREPROCESS_FLAG})

add_executable(combine_gr3 combine_gr3.f90)
add_executable(combine_hotstart7 combine_hotstart7.f90)
add_executable(combine_output11 combine_output11.f90 netcdf_var_names.f90)

target_link_libraries(combine_hotstart7 utillib ${NetCDFLIBS} ${HDF5_LIBRARIES})
target_link_libraries(combine_output11 utillib ${NetCDFLIBS} ${HDF5_LIBRARIES})

target_include_directories(combine_hotstart7 PRIVATE ${CMAKE_CURRENT_BINARY_DIR}/../UtilLib)
target_include_directories(combine_output11 PRIVATE ${CMAKE_CURRENT_BINARY_DIR}/../UtilLib)

add_dependencies(combine_gr3 utillib)
add_dependencies(combine_hotstart7 utillib)
add_dependencies(combine_output11 utillib)
