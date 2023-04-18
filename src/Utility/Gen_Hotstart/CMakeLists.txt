message(STATUS "In Utility/Gen_Hotstart")

# set list of files
set(gen_hotstart_file_list
    convert_hotstart_nc.f90
    gen_hot_3Dth_from_hycom.f90
    gen_hot_from_hycom.f90
    gen_3Dth_from_hycom.f90)

# define executables and their dependencies
foreach(f ${gen_hotstart_file_list})
  get_filename_component(fwe ${f} NAME_WE)
  add_executable(${fwe} ${f})
  target_include_directories(${fwe} PRIVATE ${CMAKE_CURRENT_BINARY_DIR}/../UtilLib)
  target_link_libraries(${fwe} utillib ${NetCDFLIBS} ${HDF5_LIBRARIES})
  add_dependencies(${fwe} utillib)
endforeach()
