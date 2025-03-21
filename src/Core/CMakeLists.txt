
## Algorithm options. Everything done using make_define_options will be applied globally.

set(versionfile "${CMAKE_CURRENT_BINARY_DIR}/schism_version.F90")

if (NOT Python_FOUND)
  find_package(Python COMPONENTS Interpreter)
  message(INFO " Python_EXECUTABLE = ${Python_EXECUTABLE}")
endif()

add_custom_target(
    sversion
    COMMAND ${Python_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/gen_version.py ${versionfile}
    WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
)

set_source_files_properties(${versionfile} PROPERTIES GENERATED TRUE)

# @todo this is a temporary fix for gcc10 and greater until this is
# resolved in later versions of gcc
if (CMAKE_Fortran_COMPILER_ID STREQUAL "GNU")
  if (CMAKE_Fortran_COMPILER_VERSION VERSION_GREATER_EQUAL "10.0.0" )

   add_compile_options("$<$<COMPILE_LANGUAGE:Fortran>:-fallow-argument-mismatch>")
	 #add_compile_options(-fallow-argument-mismatch)
  endif()
  #add_compile_options(-ffree-line-length-none)
  add_compile_options("$<$<COMPILE_LANGUAGE:Fortran>:-ffree-line-length-none>")
  endif()


set(core_extra_src "")
if(${USE_HA})
  set (core_extra_src ${core_extra_src} harm.F90)
endif()

if(${USE_PAHM})
  set (core_extra_src ${core_extra_src} PaHM/PaHM_Global.F90 PaHM/csv_module.F90 PaHM/csv_parameters.F90 PaHM/csv_utilities.F90 PaHM/PaHM_Sizes.F90 PaHM/PaHM_Messages.F90 PaHM/PaHM_Utilities.F90 PaHM/PaHM_Vortex.F90 PaHM/parwind.F90 PaHM/SortUtils.F90 PaHM/TimeDateUtils.F90)
endif()

add_library(core ${local_extra_code}
                 schism_glbl.F90
                 schism_msgp.F90
                 schism_assert.F90
                 schism_io.F90
                 scribe_io.F90
		 misc_modules.F90
                 hydraulic_structures.F90
                 ${core_extra_src}
                 ${versionfile}
           )

add_dependencies(core sversion)

list( APPEND locallibs core)
set(schismmodlibs ${locallibs} ${schismmodlibs} PARENT_SCOPE)




