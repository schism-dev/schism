message(STATUS "In Utility/Gen_Hotstart")
add_executable( convert_hotstart_nc convert_hotstart_nc.f90 )
add_executable( gen_hot_3Dth_from_hycom gen_hot_3Dth_from_hycom.f90)

target_link_libraries (convert_hotstart_nc utillib ${NetCDFLIBS} ${HDF5_LIBRARIES} )
target_link_libraries (gen_hot_3Dth_from_hycom utillib ${NetCDFLIBS} ${HDF5_LIBRARIES} )
