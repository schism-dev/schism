Following is a curated changelog of the code. The IDs starts with `R` represents the `svn` era id, others are `git` hash.
To find all the changes between tag releases, search for 'Tag'.

## Bug fixes and major algorithmic changes
91. `R5178`: fix bug of saturate `DO` in [ICM](modules/icm.md)
92. `R5187`: incorporated dry bnd treatment from bndry (c/o Jens Wyrwa). However, to get good volume conservation it's best to keep all open flow bnd's wet all the time as before;
93. `R5191`: fixed a bug in `iwind_form` (accidentally re-init'ed after reading);
94. `R5202`: fixed a bug on negative ICM tracer concentration - check negative values before sending back to hydro.
95. `R5204`: a new option to control behavior of `btrack` when trajectory hots an open bnd (set `vel=0` and exit). Goal is to eventually use this as default;
96. [`64b7181`](https://github.com/schism-dev/schism/commit/64b7181): fixed an efficiency issue and added nc checks in `ptrack3` (c/o Marcel Rieker)
97. [`2138e77`](https://github.com/schism-dev/schism/commit/2138e77): Fixed a major bug introduced during reshuffling of `_init` for PDAF: need to call `nodavel` before sflux, as uu2,vv20 are needed there.
98. [`5ff2198`](https://github.com/schism-dev/schism/commit/5ff2198): revamped [AGE](modules/age.md) module. B.C. now all should use 0. Only 0 and 1 should be used in `AGE_hvar_[1:ntr/2].ic`, where ntr is @ of age tracers, and the code will 'hold' concentration at 1 at those prisms (level specified in param.nml) that have 1 as i.c. I.C. in `AGE_hvar_[ntr/2+1:ntr].ic` should `=0`. Injecting '1' at dry elem's is allowed but results may be harder to interpret.
99. [`6cec698`](https://github.com/schism-dev/schism/commit/6cec698): bug fix for spherical coord in [WWM](modules/wwm.md) (around dateline);
100. [`c46c7bd`](https://github.com/schism-dev/schism/commit/c46c7bd) (Feb 14, 2020): Fei tweaked WENO solver (wrt upwind);
101. [`d3adb2c`](https://github.com/schism-dev/schism/commit/d3adb2c) (Feb 14, 2020): changed to bilinear interp for HYCOM hot and nudging scripts (to speed up);
102. [`7506147`](https://github.com/schism-dev/schism/commit/7506147) (Mar 1, 2020): add limiter for settling vel to avoid char line out of bnd;
103. [`5665418`](https://github.com/schism-dev/schism/commit/5665418) (Mar 18, 2020): PR from Baptiste mengual for [SED3D](modules/sed3d.md) (`sed_frition.F90` merged after the rest).
104. [`eb00a00`](https://github.com/schism-dev/schism/commit/eb00a00) (Mar 31, 2020): added 12th tracer model (`USE_DVD` of Klingbeil)
105. [`aaf1306`](https://github.com/schism-dev/schism/commit/aaf1306) (April 4, 2020): added hybrid ELM transport for efficiency (via branch hybrid_ELM);
106. [`40955ab`](https://github.com/schism-dev/schism/commit/40955ab) (April 16, 2020): replaced the fatal error in `nadv=2` (aptivity trap) with exit;
107. [`e703189`](https://github.com/schism-dev/schism/commit/e703189) (May 6, 2020): reverted to old simple approach for depositional mass (c/o Baptiste Mengual);
108. [`14e201f`](https://github.com/schism-dev/schism/commit/14e201f) (May 7, 2020): switched to ParMETIS v4.0.3;
109. [`33fb37f`](https://github.com/schism-dev/schism/commit/33fb37f) (Jun 30, 2020): fixed a bug in `ielm_transport` (Fei Y.);
110. [`55910d1`](https://github.com/schism-dev/schism/commit/55910d1) (July 31, 2020; via branch mem_leak_nc4): reverted the change in schism_io (close/open nc output after/b4 each write, to accommodate PDAF's flexible mode), as this has caused a lot of mem leaking (discovered by Nicole C. when testing ICM).
111. [`eaf9093`](https://github.com/schism-dev/schism/commit/eaf9093) (Aug 2020): Baptiste M. modified the computation of the bottom shear stress under combined waves and currents in SED3D
112. [`4ac56df`](https://github.com/schism-dev/schism/commit/4ac56df) (Aug 2020): Paul Ryan and Claire Trenham (CSIRO) fixed some init array issues in hydro and WWM.
113. [`7e9ff0a`](https://github.com/schism-dev/schism/commit/7e9ff0a) (Aug 20, 2020): Fei Y fixed a msource bug (init of msource at 0 at elem's not in source_sink.in led to ice rain;
114. [`f311026`](https://github.com/schism-dev/schism/commit/f311026) (Aug. 27, 2020): reverted init of msource to 0 for tracers other than T,S, c/o Nicole Cai. Using ambient values for other tracers can lead to artificial and additional accumulation of nutrients;
115. [`00b6f16`](https://github.com/schism-dev/schism/commit/00b6f16) (Nov. 9, 2020): PR #21 from Kijin merged (SED): changed method to compute near bottom vel for $LSC^2$.
116. [`8083cfa`](https://github.com/schism-dev/schism/commit/8083cfa) (Nov. 10, 2020): added option for transport solver only (`itransport_only`)
117. [`6dc44d3`](https://github.com/schism-dev/schism/commit/6dc44d3) (Dec 14, 2020): Fixed a bug, thanks to Kevin Kartin, on vortex formulism of wave force (missing a factor of area() in $I_4$);
118. [`28ff9d1`](https://github.com/schism-dev/schism/commit/28ff9d1) (Dec 16, 2020): added optional self-attraction loading tides. The option shares some constants with tidal potential: freq names and cut-off depth;
119. [`b1bcaa0`](https://github.com/schism-dev/schism/commit/b1bcaa0) (April 20, 2021): removed most of `goto`. Remaining ones: `harm.F90`, `lap.F90`, `WWM`;
120. [`0cec024`](https://github.com/schism-dev/schism/commit/0cec024) (April 21, 2021): bug fix in `misc_subs` c/o of Fei (`inunfl=1`: nodeA in final extrap stage may be interface node);
121. [`fb30239`](https://github.com/schism-dev/schism/commit/fb30239) (April 22, 2021): bug fix for ellipsoidal earth (tensor frames).
122. [`fb79cdc`](https://github.com/schism-dev/schism/commit/fb79cdc) (May 26, 2021): in `interpolate_depth_structured2*`, added an option to shift 1/2 cell for ll corner (c/o Charles Seaton);
123. [`5e87c24`](https://github.com/schism-dev/schism/commit/5e87c24) (June 9, 2021): more changes in `interpolate_depth_structured2*` to extrap also into right/upper sides.
124. [`843c40f`](https://github.com/schism-dev/schism/commit/843c40f) (19 June, 2021): fixed a bug in `ptrack3` (pt_in_poly3; c/o Jilian Xiong) that affects quads;
125. [`b2cf92b`](https://github.com/schism-dev/schism/commit/b2cf92b) (29 June, 2021): fixed a bug in station outputs in basin scale cases `ics=2` (local proj is not accurate if the station is far away from the local frame);

- Tag v5.9.0: 657157a (1 July, 2021)

126. `a9f8f6c` (8 July, 2021): fixed a bug in iloadtide=1 (index error). Also revamped reading of input files
                              for large core counts.
127. `ab46def` (20 July, 2021): merged multi ice model (CICE) from Qian Wang.

128. `70db927` (11 Aug 2021): bug fixes from Qian Wang on ice matrix (ghost nodes)
129. `8efc374` (30 Aug 2021): first version of scribe dictated I/O (as an option)
130. `53689f4` (12 Sept 2021): fixed a bug in weno transport introduced in e794d8c (23 Aug 2021).
131. `c7dd516` (27 Sept 2021): AGE module tweaked (clamped @ i.c. at injection points for both concentrations).
132. (19 Oct 2021): filtered isolated wet nodes during outputs (should not change other results);
133. `7455d5d` (30 Dec 2021): added a new option nws=-1 (PaHM, Parametric Hurricane Model);
134. `7455d5d` (30 Dec 2021): added Fariall air-sea exchange (USE_BULK_FAIRALL) c/o Jerome Lefevre;
135. `e288f9a` (31 Jan 2022): removed the restriction of counter-clockwise orientation in sflux grid;
136. `dde345f` (7 April 2022): fixed a bug in sflux, invoked only if nws=2, ihconsv=1, and iwind_form=0:
                              uninit'ed stress values.
137. `6524e19` (27 April 2022): merged with La Rochelle group on WWM  (VOR);


## Changes in input and output format
The info below can also be found in src/Readme.beta_notes. Most changes are made in param.in (now renamed as [param.nml](input-output/param.md)).

- `f8e15ee` (11 Oct 2022): added lev_tr_source(:)=0 option to allow source injecton to the entire water column;
- `bc113eb2` (30 Sep 2022): added a new flag 'iflux_out_format' to control flux.out format;
- `24a2723` (27 Sep 2022): [merged PR from LRU] new options related to vegetation effects on WWM. Added
                       iof_wwm(31) as dissipation due to vegetation. Also
                       nws=-3: surface stress formulation from Donelan et al. (1993) based on the wave age;
- `c94b2ca` (1 Sep 2022): added 'iof_ugrid' to output UGRID metadata for 3D outputs;
- `a0960ed` (25 Aug 2022): added 'iprecip_off_bnd' to optionally turn off precip near land bnd (to avoid some spurious plumes
               near tropical islands etc);
- `ed26f60` (25 Aug 2022): removed 'impose_net_flux' and replaced it with IMPOSE_NET_FLUX (so this won't be actively supported);

Tag stofs3d-atl.v1.1.0 (identical to v5.10.0): for NOAA STOFS3D Atlantic operational forecast release

Tag v5.10.0: (13 July, 2022)

- `d03830e` (16 May 2022): removed the 2 extra lines required in `partition.prop` (so it's identical to `global_to_local.prop` now);
- `a71ea69` (12 May 2022): removed `ioffline_partition` and replaced it with CPP (`NO_PARMETIS`) to allow build without ParMETIS lib; renamed `global_to_local.in` as `partition.prop`;
- `da90a33` (2 May 2022): removed `itr_met=1,2`;
- `2165224` (29 April 2022): added optional offline partitioning `ioffline_partition` to bypass ParMETIS;
                         if on, needs `global_to_local.in`;
- `567b612` (28 April 2022): added new WWM outputs and re-ordered after `iof_wwm(26)`; `zCoordinate` outputs made optional (for performance test);
- `6524e19` (27 April 2022): a few more WWM related parameters under VOR (`turbinjds,alpha,wfwvor_streaming`)
                         and new WWM outputs, after merging with `lrdev_updates`;
- `d44dbdf` (1 April, 2022): added an optional threshold inflation ratio for mass conversation for ICM: `rinflation_icm`
- `6f367d0` (30 Mar 2022): added optional vertical levels to inject source tracer conventration: `lev_tr_source`;
- `7cc7f27` (30 Mar 2022): added an optional threshold depth for mass conversation for ICM: `h_massconsv`;
- `1788867` (16 Mar 2022): added optional parameters for S exchange (split from heat exchange): `i_hmin_salt_ex,hmin_salt_ex`;
- `dcd5356` (31 Jan 2022): removed ramp flags except `nramp_elev` (use `dramp*`<=0 to turn off ramp):
      `nramp, nrampbc, nrampwind, nrampwafo,nramp_ss`; removed the optional input `shapiro_min.gr3`;
- `80f1caa` (26 Jan 2022): added 4 parameters in `ice.nml`: `albsn,albsnm, albi,albm`
- `d792ed0` (18 Jan 2022): added 2 new ice parameters: `ice_atmos_stress_form` and `cdwin0`;
- `b4d2214` (7 Jan 2022): added a new parameter `nu_sum_mult` to allow sum/product of final relax constants for tracer nudging;
- `f9f0844` (3 Jan 2022): ICE module, added 2 new parameters: `lead_closing, Saterm`. Also changed ice-wind stress formula;
- `7455d5d` (30 Dec 2021): added a new option `nws=-1` (PaHM, Parametric Hurricane Model);
- `91c91a7` (17 Dec 2021): ICE module, added new parameters: `ncyc_fct, depth_ice_fct`. Also requires `ice_fct.gr3` now;
- `0e2f90e` (17 Dec 2021): added a new parameter in `ice.nml` `mevp_coef` and 2 related ones: `mevp_alpha[34]`
- `2a7647d` (24&29 Nov 2021): added a new parmater `loadtide_coef` for `iloadtide=2,3`;
- `a7a0b2b` (17 Oct 2021): restored `iwindoff`;
- `afc7b1a` (11 Oct 2021): changed `vgrid.in` format for `ivcor=1` to facilitate parallel read (use
                       `change_vgrid.f90` to convert)
- `e9aec89` (1 Oct 2021): added most module outputs; order of flags changed in: ICE, ANALYSIS, SED;
- `85b585a` (1 Sept 2021): add zcor output flag as `iof_hydro(25)`; rest of hydro output flags (`iof_hydro`) shifted by 1;
- `8efc374` (30 Aug 2021): first version of scribe dictated I/O (as an option);
- `1902c1e` (16 July 2021): changed `hotstart.nc`: added `nsteps_from_cold` and `cumsum_eta` for DA;

- Tag v5.9.0: 657157a (1 July, 2021)
- [`e281d94`](https://github.com/schism-dev/schism/commit/e281d94) (25 June 2021): added a new option for SAL (Stepanov & Hughes 2004): `iloadtide=3`;
- [`bfb4afc`](https://github.com/schism-dev/schism/commit/bfb4afc) (18 June, 2021): added an optional input `shapiro_min.gr3` to be used with `ishapiro=2`;
- [`d07d75d`](https://github.com/schism-dev/schism/commit/d07d75d) (April 21, 2021): added T,S in required inputs for offline transport (to use hydro only results);
- [`b66e554`](https://github.com/schism-dev/schism/commit/b66e554) (April 14, 2021): added `ishapiro=2` - Smagorinsky like filter option. In this option, `shapiro0` is the coefficient;
- [`b43afea`](https://github.com/schism-dev/schism/commit/b43afea) (April 2, 2021): changed x,y to double in nc outputs (for newer visIT);
- [`9150d86`](https://github.com/schism-dev/schism/commit/9150d86) (Mar 31, 2021): added a new SAL option (`iloadtide=2`) using a simple scaling;
- [`084e149`](https://github.com/schism-dev/schism/commit/084e149) (Feb 25, 2021): Removed parameters: `ibtrack_openbnd`, `iwindoff`, `dzb_decay` (with hardwire); also removed option for negative roughness. Fixed race condition for marsh module.
- [`3576c9c`](https://github.com/schism-dev/schism/commit/3576c9c) (Jan 24, 2021): added new option for source/sink input: `if_source=-1` requires
`source.nc` (which includes elem list inside; allows different time steps and # of records for volume/mass source/sinks. Also, now the source/sink values in `.th` must be single precision (not double);
- [`bcdfce6`](https://github.com/schism-dev/schism/commit/bcdfce6) (Jan 6, 2021): added main switch for nc output `nc_out` (useful for other programs to control outputs);
- [`28ff9d1`](https://github.com/schism-dev/schism/commit/28ff9d1) (Dec 16, 2020): added optional self-attraction loading tides; if `iloadtide=1`, need amp/phases in `loadtide_[FREQ].gr3` (freq's shared with tidal potential);
- [`f8ba470`](https://github.com/schism-dev/schism/commit/f8ba470) (Dec 3, 2020): Added a new option (`meth_sink`) to treat net sink: if an elem is dry with a net sink, vsource is reset to 0;
- [`f9043e2`](https://github.com/schism-dev/schism/commit/f9043e2) (Nov. 12, 2020): merged with LRU branch `wwm_lr`; removed `sav_cd` (`isav=1` now requires `sav_cd.gr3`). 
    - New parameters related to WWM: `fwvor_advxy_stokes`,`fwvor_advz_stokes`, `fwvor_advz_stokes`, `fwvor_breaking`, `cur_wwm`, `wafo_obcramp`; 
    - Changes in WWM .nml: `ALPBJ ->B_ALP`, `BRHD->BRCR`.
- [`8083cfa`](https://github.com/schism-dev/schism/commit/8083cfa) (Nov. 10, 2020): added option for offline transport solver only ('itransport_only');
- [`7e9ff0a`](https://github.com/schism-dev/schism/commit/7e9ff0a) (Aug 20, 2020): added two options (selected by `i_hmin_airsea_ex`) for locally turning off heat/salt exchange.
    - `i_hmin_airsea_ex=1`: exchange turned off if local grid depth<hmin_airsea_ex;
    - `i_hmin_airsea_ex=2`: exchange turned off if local water depth<hmin_airsea_ex
    - This replaces the change made by [`adde1aa`](https://github.com/schism-dev/schism/commit/adde1aa) (July 8, 2020): added a new parameter `hmin_airsea_ex` (min total water depth for heat/salt exchange);
- [`2f7bab1`](https://github.com/schism-dev/schism/commit/2f7bab1) (Jun 28,2020): reorganized ICM process (rates/fluxes) and added output flag in `param.nml`.
- [`1477780`](https://github.com/schism-dev/schism/commit/1477780) (May 8, 2020): changed `sediment.in` to sediment.nml. 
    - `BEDLOAD_COEFF` $\rightarrow$ `bedload_coeff`, `NEWLAYER_THICK` $\rightarrow$ `newlayer_thick`, `IMETH_BED_EVOL` $\rightarrow$ `imeth_bed_evol`
    - `SAND_SD50` $\rightarrow$ `Sd50`, `SAND_ERATE` $\rightarrow$ `Erate`, `SED_TYPE` $\rightarrow$ `iSedtype`, `SAND_SRHO` $\rightarrow$ `Srho`, `SAND_WSED` $\rightarrow$ `Wsed`, `SAND_TAU_CE` $\rightarrow$ `tau_ce`, `sed_morph_fac` $\rightarrow$ `morph_fac`, `poro_cst` $\rightarrow$ `porosity`);
