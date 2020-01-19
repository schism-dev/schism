
add_custom_target(utility)
add_subdirectory( UtilLib)
add_subdirectory( Combining_Scripts )
add_subdirectory( Post-Processing-Fortran )
add_subdirectory( Pre-Processing )
add_subdirectory( Tides )
add_subdirectory( OneWayNestScripts )
add_subdirectory( Particle_Tracking )
add_subdirectory( Gen_Hotstart )
add_subdirectory( Grid_Scripts)

add_dependencies(utility utillib)


#     Currently a link problem with NCOM
# add_subdirectory( NCOM )

#add_subdirectory( ArcGIS )
