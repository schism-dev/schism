# Require 3.18 for CMAKE_Fortran_PREPROCESS support
cmake_minimum_required (VERSION 3.18)

set(CMAKE_Fortran_PREPROCESS ON)

set(BMI_SOURCES bmi.f90 bmischism.f90 schism_model_container.f90)

add_library(schism_bmi SHARED ${BMI_SOURCES})

target_compile_definitions(schism_bmi PRIVATE NGEN_ACTIVE)
target_link_libraries(schism_bmi PUBLIC hydro core)

add_executable(schism_bmi_driver schism_bmi_driver_test.f90)
target_link_libraries(schism_bmi_driver PUBLIC schism_bmi)
