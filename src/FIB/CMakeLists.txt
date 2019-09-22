
set(schismmodlibs fib  ${schismmodlibs} PARENT_SCOPE)


##### Create library and add to the list of compiled libraries

add_library(fib fib_param.F90 fib.F90)


##### Dependencies
add_dependencies(fib core)
target_link_libraries(fib core)




