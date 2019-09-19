## Create an cmake option associated with a compiler option (-D on compile line)
## Preference is given to CMake variables, but an attempt is made to look at the environment
## in which case the variable is not cached. 

macro (make_define_option name doc default extra_code)
  if(DEFINED ${name})
    message(STATUS "Option ${name} set by user cache or command to ${${name}}, default ${default}")
  else (DEFINED ${name})
    set(val ${default} )
    if ($ENV{name} MATCHES ".+" )
      set(val $ENV{name})
      message(STATUS "Option ${name} not defined in cache but discovered in environment. Using environment value ${val}")      
    endif ()
    option( ${name} ${doc} ${val} )
  endif( DEFINED ${name} )

  if(${name})
    set_property(GLOBAL APPEND PROPERTY DEFINE-LIST ${name})
    list(APPEND local_extra_code ${extra_code} )
    list(APPEND local_define_options "-D${name}")
  endif()
endmacro(make_define_option)

