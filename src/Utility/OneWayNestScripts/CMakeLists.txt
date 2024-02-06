add_executable(interpolate_variables7 interpolate_variables7.f90)
add_executable(timeint_3Dth2 timeint_3Dth2.f90)
add_executable(interpolate_variables8 interpolate_variables8.f90)

add_dependencies(interpolate_variables7 utillib)
add_dependencies(timeint_3Dth2 utillib)
add_dependencies(interpolate_variables8 utillib)

target_link_libraries(interpolate_variables8 utillib ${NetCDFLIBS} ${HDF5_LIBRARIES})
target_include_directories(interpolate_variables8 PRIVATE ${CMAKE_CURRENT_BINARY_DIR}/../UtilLib)

target_link_libraries(interpolate_variables7 utillib ${NetCDFLIBS} ${HDF5_LIBRARIES})
target_include_directories(interpolate_variables7 PRIVATE ${CMAKE_CURRENT_BINARY_DIR}/../UtilLib)

