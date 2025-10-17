add_library(fabm_models_uzweijens OBJECT
            model_library.F90
            qsim_age.F90
           )

target_link_libraries(fabm_models_uzweijens PRIVATE fabm_base)
