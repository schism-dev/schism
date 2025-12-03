######## QSim ########
# if (USE_QSIM)
#   message(STATUS "\n### going to USE_QSIM in src/QSim/CMakeList.txt")
# endif()

add_library(fabm_models_uzweijens OBJECT
            model_library.F90
            qsim_variables.F90
            qsim_age.F90
            qsim_tracerzone.F90
           )

# target_link_libraries(fabm_models_uzweijens PRIVATE fabm_base)
target_link_libraries(fabm_models_uzweijens PUBLIC fabm_base)

######## QSim ########
# if (USE_QSIM)
#   add_definitions("-DUSE_QSIM")
#   message(STATUS "\n### USEd_QSIM in /src/QSim/CMakeList.txt")
# endif()
