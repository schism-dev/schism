## Algorithm options. Everything done using make_define_options will be applied globally.

set(versionfile "${CMAKE_CURRENT_BINARY_DIR}/schism_version.F90")

# CMake-native version generation using built-in Git support
find_package(Git QUIET)

if(Git_FOUND)
    # Get git describe output
    execute_process(
        COMMAND ${GIT_EXECUTABLE} describe --always --dirty
        WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
        OUTPUT_VARIABLE GIT_DESCRIBE_OUTPUT
        OUTPUT_STRIP_TRAILING_WHITESPACE
        ERROR_QUIET
    )
    
    # Get short hash as fallback
    execute_process(
        COMMAND ${GIT_EXECUTABLE} rev-parse --short HEAD
        WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
        OUTPUT_VARIABLE GIT_SHORT_HASH
        OUTPUT_STRIP_TRAILING_WHITESPACE
        ERROR_QUIET
    )
    
    # Parse git describe output
    if(GIT_DESCRIBE_OUTPUT)
        set(VERSION_RAW ${GIT_DESCRIBE_OUTPUT})
        
        # Check for dirty flag
        string(FIND ${VERSION_RAW} "-dirty" IS_DIRTY_POS)
        if(IS_DIRTY_POS GREATER -1)
            set(IS_DIRTY true)
            string(REPLACE "-dirty" "" VERSION_RAW ${VERSION_RAW})
        else()
            set(IS_DIRTY false)
        endif()
        
        # Parse semantic version pattern: v1.2.3-4-g5abc6def
        string(REGEX MATCH "^(v?[0-9]+\\.[0-9]+\\.[0-9]+)(-[0-9]+)?(-g[a-f0-9]+)?$" VERSION_MATCH ${VERSION_RAW})
        
        if(VERSION_MATCH)
            set(SCHISM_VERSION ${CMAKE_MATCH_1})
            set(GIT_VERSION ${GIT_SHORT_HASH})
            
            # Check if there are commits after the tag
            if(CMAKE_MATCH_2 AND CMAKE_MATCH_3)
                string(SUBSTRING ${CMAKE_MATCH_2} 1 -1 COMMITS_SINCE_TAG)
                string(SUBSTRING ${CMAKE_MATCH_3} 2 -1 GIT_VERSION)
                string(APPEND SCHISM_VERSION "mod")
                set(GIT_VERSION "${GIT_VERSION} (${COMMITS_SINCE_TAG} commits since semantic tag, edits=${IS_DIRTY})")
            elseif(IS_DIRTY)
                string(APPEND SCHISM_VERSION "mod") 
                set(GIT_VERSION "${GIT_VERSION} (0 commits since semantic tag, edits=${IS_DIRTY})")
            endif()
        else()
            set(SCHISM_VERSION "semantic version not determined")
            set(GIT_VERSION ${GIT_SHORT_HASH})
        endif()
    else()
        set(SCHISM_VERSION "semantic version not determined")
        set(GIT_VERSION ${GIT_SHORT_HASH})
    endif()
else()
    set(SCHISM_VERSION "semantic version not determined")
    set(GIT_VERSION "none")
endif()

# Check for user-defined version file as fallback
if(SCHISM_VERSION STREQUAL "semantic version not determined")
    set(USER_VERSION_FILE "${CMAKE_CURRENT_SOURCE_DIR}/schism_version_user.txt")
    if(EXISTS ${USER_VERSION_FILE})
        file(READ ${USER_VERSION_FILE} USER_VERSION)
        string(STRIP ${USER_VERSION} USER_VERSION)
        string(LENGTH ${USER_VERSION} USER_VERSION_LENGTH)
        if(USER_VERSION_LENGTH GREATER_EQUAL 3)
            set(SCHISM_VERSION ${USER_VERSION})
        else()
            set(SCHISM_VERSION "develop")
        endif()
    else()
        set(SCHISM_VERSION "develop")
    endif()
endif()

# Ensure we have a git version
if(NOT GIT_VERSION OR GIT_VERSION STREQUAL "")
    set(GIT_VERSION "none")
endif()

# Generate the version file
set(TEMPLATE_FILE "${CMAKE_CURRENT_SOURCE_DIR}/schism_version.F90.template")
if(EXISTS ${TEMPLATE_FILE})
    file(READ ${TEMPLATE_FILE} TEMPLATE_CONTENT)
    string(REPLACE "@{VERSION_GIT}" "${GIT_VERSION}" VERSION_CONTENT ${TEMPLATE_CONTENT})
    string(REPLACE "@{VERSION_SCHISM}" "${SCHISM_VERSION}" VERSION_CONTENT ${VERSION_CONTENT})
    file(WRITE ${versionfile} ${VERSION_CONTENT})
else()
    # Fallback: create minimal version file if template doesn't exist
    file(WRITE ${versionfile}
"      character(LEN=*),parameter :: schism_version = '${SCHISM_VERSION}', git_version = '${GIT_VERSION}'
")
endif()

# Print version info
message(STATUS " SCHISM version: ${SCHISM_VERSION}")
message(STATUS " GIT commit: ${GIT_VERSION}")

# Create custom target for dependencies
add_custom_target(sversion ALL DEPENDS ${versionfile})

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
