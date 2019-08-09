message(STATUS "In Utility/Gen_Hotstart")
add_executable( convert_hotstart_nc convert_hotstart_nc.f90 )

target_link_libraries (convert_hotstart_nc ${NetCDFLIBS} ${HDF5_LIBRARIES} )
