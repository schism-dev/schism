message(STATUS "In Utility/Gen_Hotstart")

add_executable( convert_hotstart_nc convert_hotstart_nc.f90 )

# Removed from cmake because it contains huge arrays and kind of forces mcmodel=medium and doesn't
# necessarily compile gracefully even with that. Below are some flags to add these things on a
# file specific basis, but that isn't working yet.
# add_executable( gen_hot_3Dth_from_hycom gen_hot_3Dth_from_hycom.f90)
# target_compile_options(gen_hot_3Dth_from_hycom PRIVATE "-O2 -CB -mcmodel=medium")
# target_link_options(gen_hot_3Dth_from_hycom PRIVATE "-mcmodel=medium -CB")

target_link_libraries (convert_hotstart_nc utillib ${NetCDFLIBS} ${HDF5_LIBRARIES} )
#target_link_libraries (gen_hot_3Dth_from_hycom utillib ${NetCDFLIBS} ${HDF5_LIBRARIES} )
