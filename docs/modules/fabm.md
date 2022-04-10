
The `fabm` module is the SCHISM host implementation for the separatly available Framework for Aquatic Biogeochemical models (FABM).  Within the FABM framework, a number of host-agnostic ecosystem models have been implemented, such that they work within many different hydrodynamic host models, including a `python` host and `Fortran` hosts from zero to three-dimensional, and now including SCHISM.  For these hosts, FABM's internal loops are arranged for optimal numerical performance. 

## Obtaining the FABM framework code

The FABM framework code is hosted separately, as are some of the ecosystem model codes. The SCHISM host provided here is compatible with both the deprecated `version 0` of the FABM framework and the current `version 1`. 

The required FABM code can be downloaded from 
1. the official repository https://github.com/fabm-model/fabm 
2. our development fork that closely follows the official repository, and includes modularization of some the existing SCHISM modules refactored within FABM.
3. any fork of FABM as of their  `SHA:048673` (2017-09-15).

## Compiling with FABM

The FABM framework and its contained models are compiled within the SCHISM build when using `cmake`

```
cmake  [...] -DUSE_FABM=ON -DFABM_BASEDIR=/your/path/to/fabm/ 
```

## Developing a FABM biogeochemical model 

Documentation on how to implement an ecosystem model in FABM is provided on the FABM wiki page https://github.com/fabm-model/fabm/wiki/Developing-a-new-biogeochemical-model

## Implementation of SCHISM as a FABM host 

### Calls to FABM infrastructure from Hydro/schism_init.F90

 1. **Initialization of FABM**
    `schism_init.F90` calls `fabm_schism_init_model` and gets back the number of tracers to be transported. This step reads the FABM configuration, specifically the files `fabm.yaml` and `fabm_schism.in`. The FABM package is initialized, which builds up its infrastructure tree. Now, FABM knows how many state variables, surface state variables, bottom state variables and diagnostic variables are required. |
 2. **Post-Initialization of FABM**
    At this stage, the necessary arrays are allocated within SCHISM by the routine `fabm_schism_init_stage2`. Pointers to the state variables and forcings are set for FABM, further arrays are allocated based on the grid layout.  An output netcdf file is created by `fabm_schism_create_output_netcdf` to store diagnostics and bottom state variables. This netcdf file can also used for hotstart, however this step can be merged into SCHISM's netcdf infrastructure now. |
 3. **Initialize concentrations** 
    `fabm_schism_init_concentrations` is called to set the initial concentrations from the FABM configuration.
 4. **Read 2D hotstart data from netcdf**
    With `fabm_schism_read_horizontal_state_from_netcdf`, the 2d hotstart information, e.g. for bottom state variable concentrations on elements, is read from a global netcdf. |


### Calls to FABM infrastructure from Hydro/schism_step.F90

During the integration, SCHISM calls FABM to get the rates of change of the state variables and the swimming and sinking speeds of moving state variables at each timestep based on the current environmental conditions. Most of the forcing (e.g. temperature) can be linked during the initialization. Currently, the bottom stress in the element center is calculated exclusively for FABM within `schism_step.F90`.

 1. **timestep call** The major call to FABMs ecosystem models within `schism_step.F90` using `fabm_schism_do`. Within this call, the light conditions are calculated, further forcing is updated in FABM, diagnostic variables are updated, the rates of change of state variables are collected from FABM and put into the `bdy_frc`, `flux_bt`, and `flux_sf` arrays. The current swimming and sinking speeds are taken from FABM and put into the `wsett` array. The changes of bottom and surface state variables are calculated with a simple Euler integration step, which is a potentially weak feature. The temporal integration of bottom and surface state variables can potentially improved, e.g. iterated if necessary. |

 2. **write output by SCHISM** All state variables are written by `writeout_nc` at each output time step by SCHISM's output infrastructure.

 3. **write specific output** The routine `fabm_schism_write_output_netcdf` is called for the output timesteps. The diagnostic variables to be written are collected before, which includes temporal averaging, if required. This step can be merged into the SCHISM output infrastructure.
