
set(schismmodlibs cosine ${schismmodlibs} PARENT_SCOPE)


##### Create library and add to the list of compiled libraries

add_library(cosine cosine_mod.F90 cosine_init.F90 cosine.F90)


##### Dependencies
add_dependencies(cosine core)
target_link_libraries(cosine core)




