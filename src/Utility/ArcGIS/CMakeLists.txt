

add_executable( binary_2_ncstruc binary_2_ncstruc.f90)

#add_executable( read_output7b_group_zfs2 read_output7b_group_zfs2.f90)


target_link_libraries (binary_2_ncstruc ${NetCDFLIBS} ${HDF5_LIBRARIES} schism_geometry)
add_dependencies(utility binary_2_ncstruc read_output7b_group_zfs2)

