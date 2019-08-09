
### Algorithm options. Everything done using make_define_options will be applied globally.
make_define_option("USE_HA"  "Use the Harmonic Analsis option" OFF "harm.F90")

if(DEFINED SIMPLE_WIND)
  add_definitions(-DUSE_SIMPLE_WIND=${SIMPLE_WIND})
endif()

set(versionfile "${CMAKE_CURRENT_BINARY_DIR}/schism_version.F90")

add_custom_target(
    sversion
    COMMAND python ${CMAKE_CURRENT_SOURCE_DIR}/gen_version.py ${versionfile}
)

set_source_files_properties(${versionfile} PROPERTIES GENERATED TRUE)

add_library(core ${local_extra_code}
                 schism_glbl.F90
                 schism_msgp.F90
                 schism_assert.F90
		 misc_modules.F90
                 ${versionfile}
           )

add_dependencies(core sversion)

list( APPEND locallibs core)
set(schismmodlibs ${locallibs} ${schismmodlibs} PARENT_SCOPE)




