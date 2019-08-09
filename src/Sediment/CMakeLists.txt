

##### Algorithmic options.Everything done using make_define_options will be applied globally. #####
# Adds the options to a local list plus a cumulative $all_define_options variable 
make_define_option("BEDLOAD_VR" "Bedload sediment option" OFF "")
make_define_option("BEDLOAD_MPM" "MPM model bedload sediment option" OFF "")
make_define_option("SUSPLOAD" "Suspended load option sediment option" OFF "")
make_define_option("SED_MORPH" "Sediment-induced morphology sediment option" OFF "")
make_define_option("UV_LOGDRAG" "Log drag formulation sediment option" ON "")
make_define_option("UV_QDRAG" "Quadratic drag formulation sediment option" OFF "")
make_define_option("UV_LDRAG" "Linear drag formulation sediment option" OFF "")


##### Create library and add to the list of compiled libraries

add_library(sediment sed_mod.F90 sed_misc_subs.F90 
            sed_init.F90 read_sed_input.F90 sed_init.F90
            sed_friction.F90 sed_filter.F90 
            sediment.F90 sed_bedload.F90)

set(schismmodlibs sediment ${schismmodlibs} PARENT_SCOPE)

##### Dependencies
add_dependencies(sediment core)
target_link_libraries(sediment core)




