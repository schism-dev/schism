set(schismmodlibs gen ${schismmodlibs} PARENT_SCOPE)

##### Create library and add to the list of compiled libraries
add_library(gen GEN_model.F90 )

##### Dependencies
add_dependencies(gen core ${schismmodlibs})
target_link_libraries(gen core ${schismmodlibs})
