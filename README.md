# Evan Vegetation Edits to SCHISM

In progress work of adding more functions to the marsh migration module within SCHISM. Currently taking user requests for additions. Not currently user friendly as inputs have not been carried into param.nml at his point. Changes:
-Mortality due to maximum and average bed shear stress 
-Species specific mortality due to inundation %
-Max erosion that causes veg mortality
-Vegetation parameters related to marsh age (exponential and linear growth, no physical parameters until mature) for multiple species
-Different tolerance to physical stressors between old and young plants
-Change marsh spread. Choose between marsh spreading to adjacent cells and spreading to any user deems suitable
-Spin up time for migration. Model won't create marsh until x days
-Lag for creating marsh, if dead, won't attempt repolutation until x days
-marsh_age.prop is now a required input. Must say how old a marsh is in days to inform physical parameters.

Future Changes:
-Addition of life stages. User can input different life stages and tolerance to stresses at this life stage
-Rooting binding of sediment
-Add a moving window option for all parameters. 
-Switch for turning off migration for certain species (i.e. there was a request for SAV to be used in conjunction with this module)
-Change mortalilty to indundation to where depth of inundation is included. Potentially relate to height of plant
-Randomness switch. Add % chance a marsh establishes in a cell. May enhance creek formation
-Incorporate flexible vegetation
-Add input changes to param.nml or make a veg.nml file. If anyone wants to use this now, changes must be made in src/Hydro/schism_init.f90
-Make critical and max shear species specific
-Look into simple rules for species succession. 
-Add rules for growth stressor. Not mortality, but slows growth. 


# SCHISM

The **S**emi-implicit **C**ross-scale **H**ydroscience **I**ntegrated **S**ystem **M**odel (SCHISM) is an open-source community-supported modeling system based on unstructured grids and designed for the seamless simulation of 3D baroclinic circulation across creek-lake-river-estuary-shelf-ocean scales.

# Building and documentation

The manual may be found on the SCHISM wiki at http://ccrm.vims.edu/schismweb/. Build instructions are described in Chapter 1.

The online documentation can be accessed at https://schism-dev.github.io/schism.

# Developing and contributing

When using the development version, note changes in flags and features described in `src/Readme.beta_notes` and `sample_inputs/param.nml`, `sample_inputs/bctides.in`, etc.

Please refer to `CONTRIBUTING.md` for more information on contributing to SCHISM.
