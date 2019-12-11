add_executable( gen_fg gen_fg.f90 )

#target_link_libraries (combine_hotstart7 utillib ${NetCDFLIBS} ${HDF5_LIBRARIES})

add_dependencies(utility gen_fg)

