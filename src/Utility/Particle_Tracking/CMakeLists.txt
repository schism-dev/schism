add_executable(ptrack3 ptrack3.f90)
add_executable(ptrack4 ptrack4.f90)

add_dependencies(utility ptrack3 ptrack4)
target_link_libraries(ptrack3 utillib ${NetCDFLIBS} ${HDF5_LIBRARIES})
target_link_libraries(ptrack4 utillib ${NetCDFLIBS} ${HDF5_LIBRARIES})


