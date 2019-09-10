add_executable(ptrack3 ptrack3.f90)

add_dependencies(utility ptrack3)
target_link_libraries(ptrack3 utillib ${NetCDFLIBS} ${HDF5_LIBRARIES})


