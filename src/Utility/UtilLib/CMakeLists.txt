##### Create library and add to the list of compiled libraries

add_library(utillib schism_geometry.f90 compute_zcor.f90 extract_mod.f90 stripsearch_unstr.f90 argparse.f90 pt_in_poly.f90 pt_in_poly_ray_method.f90)


add_dependencies(utility utillib)



