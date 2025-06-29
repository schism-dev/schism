add_executable(ptrack3 ptrack3.f90)
#add_executable(ptrack4 ptrack4.f90)

add_dependencies(utility ptrack3) ## ptrack4)
target_link_libraries(ptrack3 utillib ${NetCDFLIBS} ${HDF5_LIBRARIES})
target_include_directories(ptrack3 PRIVATE ${CMAKE_CURRENT_BINARY_DIR}/../UtilLib)
#target_link_libraries(ptrack4 utillib ${NetCDFLIBS} ${HDF5_LIBRARIES})
#target_include_directories(ptrack4 PRIVATE ${CMAKE_CURRENT_BINARY_DIR}/../UtilLib)
