message(STATUS "In /Utility NetCDF_LIBS ${NetCDFLIBS}")
add_executable( check_openbnd check_openbnd.f90 )
add_executable( convert_3Dth_nc convert_3Dth_nc.f90 )
add_executable( convert_nudge_nc convert_nudge_nc.f90 )
add_executable( edit_vgrid edit_vgrid.f90)
add_executable( change_vgrid change_vgrid.f90)
add_executable( fix_bad_quads fix_bad_quads.f90 )
add_executable( find_small_elem find_small_elem.f90 )
add_executable( gen_hdif gen_hdif.f90 )
add_executable( gen_nudge gen_nudge.f90 )
add_executable( gen_nudge2 gen_nudge2.f90 )
add_executable( gen_nudge_from_hycom gen_nudge_from_hycom.f90 )
add_executable( gen_slope_filter gen_slope_filter.f90 )
add_executable( gen_slope_filter2 gen_slope_filter2.f90 )
add_executable( gen_source gen_source.f90 )
add_executable( gen_vqs_1 gen_vqs_1.f90 )
add_executable( gen_vqs_2masters gen_vqs_2masters.f90 )
add_executable( list_bnd_pt list_bnd_pt.f90)
add_executable( manning manning.f90)
##add_executable( relocate_wet_source_elem relocate_wet_source_elem.f90)
add_executable( split_quads_wwm split_quads_wwm.f90)
add_executable( viz_source viz_source.f90 )

# set list of files
set(utility_list
    check_openbnd.f90
    convert_3Dth_nc.f90
    convert_nudge_nc.f90
    edit_vgrid.f90
    change_vgrid.f90
    fix_bad_quads.f90
    find_small_elem.f90
    gen_hdif.f90
    gen_nudge.f90
    gen_nudge2.f90
    gen_nudge_from_hycom.f90
    gen_slope_filter.f90
    gen_slope_filter2.f90
    gen_source.f90
    gen_vqs_1.f90
    gen_vqs_2masters.f90
    list_bnd_pt.f90
    manning.f90
    split_quads_wwm.f90
    viz_source.f90
    viz_source_sink.f90)

add_dependencies(utility check_openbnd convert_3Dth_nc convert_nudge_nc fix_bad_quads find_small_elem gen_hdif gen_vqs_1 gen_vqs_2masters gen_nudge gen_nudge2 gen_nudge_from_hycom gen_slope_filter gen_slope_filter2 gen_source list_bnd_pt manning split_quads_wwm viz_source )

# add extra libraries for linking
target_link_libraries(convert_3Dth_nc utillib ${NetCDFLIBS} ${HDF5_LIBRARIES})
target_link_libraries(convert_nudge_nc utillib ${NetCDFLIBS} ${HDF5_LIBRARIES})
target_link_libraries(gen_vqs_1 utillib)
target_link_libraries(gen_vqs_2masters utillib)
target_link_libraries(gen_nudge_from_hycom utillib ${NetCDFLIBS} ${HDF5_LIBRARIES})

# add include directories for linking
target_include_directories(gen_vqs_1 PRIVATE ${CMAKE_CURRENT_BINARY_DIR}/../UtilLib)
target_include_directories(gen_vqs_2masters PRIVATE ${CMAKE_CURRENT_BINARY_DIR}/../UtilLib)
target_include_directories(gen_nudge_from_hycom PRIVATE ${CMAKE_CURRENT_BINARY_DIR}/../UtilLib)
