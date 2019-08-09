
add_executable(interpolate_variables7 interpolate_variables7.f90 )
add_executable(timeint_3Dth2 timeint_3Dth2.f90)

add_dependencies(utility interpolate_variables7 timeint_3Dth2)
target_link_libraries(interpolate_variables7 utillib ${NetCDFLIBS} ${HDF5_LIBRARIES} )
