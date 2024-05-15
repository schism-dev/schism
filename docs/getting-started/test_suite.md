## Overview
SCHISM provides benchmark tests for verifying your installation or new code development.

Due to the large file size of these tests, they are distributed via svn.
You’ll need svn v1.8 or above ([svn manual](http://svnbook.red-bean.com/)).
Svn clients on linux/unix/windows/Mac should all work.

The command to checkout the tests is:

`svn co https://columbia.vims.edu/schism/schism_verification_tests`

(You can also simply use wget to download)

Note that the test suite is kept up to date with the latest master branch on the [SCHISM's git repo](https://github.com/schism-dev/schism/tree/master).
As a result, there may be some differences (some parameters may have been removed or added) between the param.nml you are using and the one from the verification tests (master branch).
 It's important to use a correct `param.nml` corresponding to the version of SCHISM you are using. A sample param.nml is always provided under your SCHISM source code directory:

`$your_schism_dir/sample_inputs/param.nml`

The latest master version for this input can also be viewed on the [Github page](https://github.com/schism-dev/schism/blob/master/sample_inputs/param.nml).
Useful info can be found in the source code bundle `src/Readme.beta_notes` (including change of format for input files and bug fixes) if you wish to hop among different versions.

To ease the burden of beginners who wish to verify their own build of different SCHISM
 versions, we have also included a subdirectory `Tags` in the test suite, in which users can
 test all major tag releases with a simple test case (Quarter Annulus).


## SCHISM Modules required in the test suite
Depending on which verification test you are conducting, you may need to enable certain modules when compiling SCHISM.

Here is a reference:

| Test      | Module needed |
| ----------- | ----------- |
| Test_Btrack_Cone | None |
| Test_Btrack_Gausshill | None |
| Test_Btrack_Gausshill_CPU | None |
| Test_CORIE | None |
| Test_CORIE_LSC2 | None |
| Test_COSINE_SFBay | None |
| Test_Chezy | None |
| Test_Convergence_Grid1 | None |
| Test_Convergence_Grid2 | None |
| Test_Convergence_Grid3 | None |
| Test_Convergence_Grid4 | None |
| Test_Convergence_Grid5 | None |
| Test_ECO_Toy | ECO |
| Test_FABM_COSINE_SFBay | FABM |
| Test_Flat_Adjust | None |
| Test_GEN_MassConsv | GEN |
| Test_GEN_MassConsv2 | GEN |
| Test_Geostrophic | None |
| Test_HeatConsv_TVD | None |
| Test_HeatConsv_Upwind | None |
| Test_HeatPool | None |
| Test_HydraulicStruct | None |
| Test_ICM_ChesBay | ICM |
| Test_ICM_UB | ICM |
| Test_Inun_CircularIsland_CaseB | None |
| Test_Inun_CircularIsland_CaseB_3D | None |
| Test_Inun_CircularIsland_CaseB_CPU | None |
| Test_Inun_NWaves_2D | None |
| Test_Inun_NWaves_3D | None |
| Test_MassSource | None |
| Test_Nonhydro_Flat_Adjust | Nonhydro (not active) |
| Test_Nonhydro_StandingWaves | Nonhydro (not active) |
| Test_ParaBowl | None |
| Test_QuarterAnnulus | None |
| Test_QuarterAnnulus_hvis | None |
| Test_SED_Trench_Migration | SED |
| Test_SED_meander_2 | SED |
| Test_Sed2d_Trench_Migration | SED |
| Test_TIMOR_rouse | None |
| Test_VolConsv_2D_1 | None |
| Test_VolConsv_2D_2 | None |
| Test_VolConsv_3D_1 | None |
| Test_VolConsv_3D_2 | None |
| Test_WWM_Analytical | WWM |
| Test_WWM_Duck | WWM |
| Test_WWM_L31_2A | WWM |
| Test_WWM_VF_adiabatic_case | WWM |
| Test_WWM_limon_NODIF | WWM |
| Test_Williamson5 | None |
