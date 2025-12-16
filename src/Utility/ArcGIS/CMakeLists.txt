add_executable( binary_2_ncstruc binary_2_ncstruc.f90)

#add_executable( read_output7b_group_zfs2 read_output7b_group_zfs2.f90)

# Link to utillib which contains stripesearch_unstr.f90 (called by binary_2_ncstruc)
target_link_libraries(binary_2_ncstruc ${NetCDFLIBS} ${HDF5_LIBRARIES} utillib)
target_include_directories(binary_2_ncstruc PRIVATE ${CMAKE_CURRENT_BINARY_DIR}/../UtilLib)
# Add dependency on utillib to ensure it's built before compilation
add_dependencies(binary_2_ncstruc utillib)
add_dependencies(utility binary_2_ncstruc)

