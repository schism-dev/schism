message(STATUS "In Utility/Grid_Scripts")

add_executable( cpp_bp cpp_bp.f90)
add_executable( cpp cpp.f90)
add_executable( grid_spring grid_spring.f90)
add_executable( grid_subsampler grid_subsampler.f)
add_executable( interpolate_depth_structured2 interpolate_depth_structured2.f90)
add_executable( interpolate_unstructured interpolate_unstructured.f90)
add_executable( periodic_grid periodic_grid.f90)
add_executable( gen_boundary_global  gen_boundary_global.f90)
add_executable( check_lonlat_grid check_lonlat_grid.f90)
add_executable( mergegrid5 mergegrid5.f90)

target_link_libraries (mergegrid5 utillib)
#mpi_wrap( combine_output6 combine_output5 combine_output4 )
#add_dependencies(utility cpp_bp cpp grid_spring grid_subsampler interpolate_unstructured periodic_grid gen_boundary_global check_lonlat_grid mergegrid5)

