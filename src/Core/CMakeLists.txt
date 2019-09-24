
## Algorithm options. Everything done using make_define_options will be applied globally.

set(versionfile "${CMAKE_CURRENT_BINARY_DIR}/schism_version.F90")

add_custom_target(
    sversion
    COMMAND python ${CMAKE_CURRENT_SOURCE_DIR}/gen_version.py ${versionfile}
)

set_source_files_properties(${versionfile} PROPERTIES GENERATED TRUE)

set(core_extra_src "")
if(${USE_HA})
  set (core_extra_src ${core_extra_src} harm.F90)
endif()

add_library(core ${local_extra_code}
                 schism_glbl.F90
                 schism_msgp.F90
                 schism_assert.F90
                 schism_io.F90
		 misc_modules.F90
                 hydraulic_structures.F90
                 ${core_extra_src}
                 ${versionfile}
           )

add_dependencies(core sversion)

list( APPEND locallibs core)
set(schismmodlibs ${locallibs} ${schismmodlibs} PARENT_SCOPE)




