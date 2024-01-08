The Community Earth System Model (CESM) is a fully coupled global climate model developed in collaboration with colleagues in the research community. CESM provides state of the art computer simulations of Earth's past, present, and future climate states.

The CESM project is supported primarily by the National Science Foundation (NSF). Administration of the CESM is maintained by the Climate and Global Dynamics Laboratory (CGD) at the National Center for Atmospheric Research (NCAR).

The [SCHISM NUOPC](nuopc.html) can be integrated into the CESM as an ocean component in principle, once it is integrated into their Common Infrastructure for Modeling the Earth (CIME, see below).  

## Obtaining CESM

The version of CESM supporting SCHISM can be obtained as 

```
git clone  https://github.com/mvertens/cesm.git -b feature/add_schism
./manage_externals/checkout_externals -v -o
```

## Obtaining and building with CIME

CIME, pronounced "SEAM", primarily consists of a Case Control System that supports the configuration, compilation, execution, system testing and unit testing of an Earth System Model. The two main components of the Case Control System are:

1. Scripts to enable simple generation of model executables and associated input files for different scientific cases, component resolutions and combinations of full, data and stub components with a handful of commands.
2. Testing utilities to run defined system tests and report results for different configurations of the coupled system.

A test case with SCHISM can be built with CIME by executing 

```
cd cime/scripts
./create_newcase --case test_schism --res TL319_rsfb --compset CSCHISM --run-unsupported --machine conda
```

After the case is created, one needs to execute `./case.build` and `./case.submit`, but these currently do not yet work.  Tasks are
- fill in `cime_config/buildlib`  
- fill in `cime_config/buildnml`

The development process is documented in https://github.com/schism-dev/schism-esmf/issues/20.  Help is much appreciated.

