

include_directories(${PROJECT_SOURCE_DIR}/Core)
add_definitions(${C_PREPROCESS_FLAG} -DSELFE_ASSERTIVE)   # make sure the preprocessor is on, assertion is a macro
make_schism_test(test_mpi_assert 4)
add_subdirectory( Hydro)
