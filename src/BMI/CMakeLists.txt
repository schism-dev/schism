# Require 3.18 for CMAKE_Fortran_PREPROCESS support
cmake_minimum_required (VERSION 3.18)

set(CMAKE_Fortran_PREPROCESS)

add_subdirectory(iso_c_fortran_bmi ${CMAKE_BINARY_DIR}/iso_c_bmi)

set(BMI_SOURCES bmi.f90 bmischism.f90 schism_model_container.f90)

add_library(schism_bmi SHARED ${BMI_SOURCES})

target_link_libraries(schism_bmi PUBLIC iso_c_bmi hydro core)

target_compile_definitions(schism_bmi PRIVATE NGEN_ACTIVE)


add_executable(schism_driver schism_bmi_driver_test.f90)
target_link_libraries(schism_driver PUBLIC schism_bmi)


