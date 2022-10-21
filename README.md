# schism_Nea

SCHISM, Original Code cloned from https://github.com/schism-dev/schism
with Modifications started on 2022 October 16 by J.Lefevre, IRD :

Updates : 17-Oct-2022
---------------------
- Add GAHM (Generalized Asymmetric Holland Model) in PaHM/parwind.F90 as a second alternative to Holland Model 
 (see details https://wiki.adcirc.org/Generalized_Asymmetric_Holland_Model and https://noaa-ocs-modeling.github.io/PaHM/pahm_manual.pdf),
- Add support for Hurricanes in both South and North Hemisphere (in HM and GAHM models),
- Test using Niran, SW Pacific: (See inside src/Core/PaHM/inputs/ : a track for Cyclone Niran "niran2021-bdeck.dat" and my comparison with HM and GAHM versus a weather model output)
- CAUTION // CAUTION : 
  To switch from HM (1) or GAHM (10), the user still need to change the value of "modelType" in Pahm_Utilities.F90 (line 3210) and recompile schism 
- CAUTION / CAUTION : Unlike in noaa-ocs-modeling.github.io/PaHM, there isnot Control File support in SCHISM/PaHM yet
---------------------
Updates End

# SCHISM

The **S**emi-implicit **C**ross-scale **H**ydroscience **I**ntegrated **S**ystem **M**odel (SCHISM) is an open-source community-supported modeling system based on unstructured grids and designed for the seamless simulation of 3D baroclinic circulation across creek-lake-river-estuary-shelf-ocean scales.

# Building and documentation

The manual may be found on the SCHISM wiki at http://ccrm.vims.edu/schismweb/. Build instructions are described in Chapter 1.

The online documentation can be accessed at https://schism-dev.github.io/schism.

# Developing and contributing

When using the development version, note changes in flags and features described in `src/Readme.beta_notes` and `sample_inputs/param.nml`, `sample_inputs/bctides.in`, etc.

Please refer to `CONTRIBUTING.md` for more information on contributing to SCHISM.
