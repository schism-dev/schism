!!! attention
    Note that some modules are under active development and we will update the info as it becomes available.

SCHISM modules can be broadly divided into two categories: tracer and non-tracer modules. The main difference is that tracer modules share more infrastructure with the main hydro code base, e.g. using the transport solver, with I.C. and B.C.’s that resemble those for the temperature and salinity, and sharing the source inputs (`msource.th`). Most modules also require additional inputs of their own (e.g. `wwminput.nml` for WWM).

There are 12 tracer modules and they are (in order of appearance and precedence in `bctides.in`; the names in brackets are used in input names; e.g. `TEM_1.th` etc) - 

1. Temperature [`TEM`]
2. Salinity [`SAL`]
3. Generic tracer [`GEN`]: generic tracer module with a settle velocity `gen_wsett` in `param.nml`). The user can use this module as a template to add their own tracer behavior etc (by modifying the code sections bounded by `USE_GEN`);
4. AGE [`AGE`]: water age module of Delhez & Deleersnijder (2002) and Shen and Haas (2004);
5. SED3D [`SED`]: 3D non-cohesive sediment transport module;
6. EcoSim [`ECO`]: EcoSim of Paul Bissett;
7. ICM [`ICM`]: USACE’s water quality model of CE-QUAL-ICM
8. CoSINE [`COS`]: Carbon, Silicate, Nitrogen Ecosystem model of Prof. Fei Chai (U. of Maine)
9. Fecalbacteria [`FIB`]: fecal indicating bacteria model;
10. TIMOR: not active at the moment
11. FABM [`FBM`]: Framework for Aquatic Biogeochemical Models, a flexible biogeochemical model framework;
12. DVD [`DVD`]: numerical mixing analysis of Klingbeit et al. (2014)

The B.C. flags for each invoked tracer module are specified in `bctides.in`. For example, if you invoked `GEN` (1 class), `SED` (1 class), and `ICM` (21 state variables inside), the boundary condition at an open segment may look like - 

```
39 2 0 3 4 1 2 2 ![# of nodes], elev, vel, T,S, GEN, SED, ICM
0.5 !constant elev
0.1 !relax for T
0.1 !relax for S
0.5 !relax for GEN
0. !constant SED concentration
1.e-3 !relax for SED
0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. !ICM concentrations
1.e-3  !relax for ICM
[next segment...]
```

And in addition, you’ll need to prepare inputs: `SAL_3D.th.nc` and `GEN_*.th`. 

Similarly, msource.th should also include the tracer concentrations for all invoked modules; for the example shown above, a line in the msource.th should look like this (assuming 2 sources in source_sink.in):

```
86400. 10. 10. 0. 0. -9999. -9999. 0. 0. -9999. ...-9999. !time (sec), T, S, GEN, SED, ICM

                                         <------------->
                                       21x2 values for ICM
```

For some tracer modules the user also needs to specify number of tracer classes inside the module in `param.nml`:

```
ntracer_gen = 2 !user defined module (USE_GEN)
ntracer_age = 4 !age calculation (USE_AGE). Must be =2*N where N is # of age tracers
sed_class = 5 !SED3D (USE_SED)
eco_class = 27 !EcoSim
```

The output flags for all modules are `iof_[name]`, where name is the lower case of the module name; e.g. `iof_wwm()`. See the [sample `param.nml`](https://github.com/schism-dev/schism/blob/master/sample_inputs/param.nml) for a complete list of output flags as well as the variable names that appear in the netcdf outputs. Some modules have additional parameters specified in `param.nml`; e.g., `gen_wsett`, `flag_fib` etc. See the [sample `param.nml`](https://github.com/schism-dev/schism/blob/master/sample_inputs/param.nml) for explanation.
