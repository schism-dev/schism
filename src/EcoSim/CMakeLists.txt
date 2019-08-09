
set(schismmodlibs ecosim ${schismmodlibs} PARENT_SCOPE)
add_library(ecosim bio_param.F90 eclight.F90 biology.F90 bio_init.F90 ecosim.F90 read_ecoin.F90 specir.F90)
target_link_libraries (ecosim core)

