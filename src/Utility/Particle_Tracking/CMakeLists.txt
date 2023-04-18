add_executable(ptrack3 ptrack3.f90)
add_dependencies(ptrack3 utillib)
target_link_libraries(ptrack3 utillib ${NetCDFLIBS} ${HDF5_LIBRARIES})
target_include_directories(ptrack3 PRIVATE ${CMAKE_CURRENT_BINARY_DIR}/../UtilLib)
