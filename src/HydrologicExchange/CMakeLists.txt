add_library(hydrologic_exchange hydrologic_exchange_soil.F90
                               hydrologic_exchange_nwm.F90
                               hydrologic_exchange_mapping.F90
                               hydrologic_exchange.F90
                               hydrologic_exchange_run.F90)

add_dependencies(hydrologic_exchange core)
target_link_libraries(hydrologic_exchange core ${NetCDFLIBS})
