Please refer to sample [param.nml](https://github.com/schism-dev/schism/blob/master/sample_inputs/param.nml) in the source code directory (sample_inputs/) while you read the following.

The file uses the FORTRAN namelist format. The order of input parameters is not important. Governing rules for this file are:

- lines beginning with `!` are comments; blank lines are ignored;
- the format for each parameter is: `keywords=value`; keywords are case sensitive; spaces allowed between `keywords` and `=` and `value`; comments starting with `!` after value are ignored;
- `value` is an integer, double, or 2-char string (use ' ' (single quotes) for this); for double, any of the format is acceptable: 40 or 40. or 4.e1 but the last 2 are preferred. Use of decimal point for integers is discouraged;
- if multiple entries for a parameter are found, the last one wins - please avoid this
- array inputs follow column major (like FORTRAN) and can spill to multiple line.

The namelist file is divided into 3 major sections: [CORE](#core), [OPT](#opt) and [SCHOUT](#schout). [CORE](#core) lists out all core parameters that _must_ be specified by the user, i.e., no defaults are provided by the code. [OPT](#opt) and [SCHOUT](#schout) sections contain optional parameters and I/O flags, all of which have default values so the user does not have to specify any of these (the values shown in the sample are defaults unless otherwise stated). SCHISM will also echo the input values in the output file `param_out.nml`.

Most parameters (and their keywords) are explained as follows; some are ‘developers handles’ that should 
not be tweaked usually. Also the sample has suggested values for many parameters. Note that you do not 
have to follow the order below (cf. `param.nml`). In many cases we have grouped related parameters 
for easier explanation, but you should specify them on separate lines. Also you'll find additional useful
 info in the comments of the sample `param.nml`. The parameters are listed out below in alphabetic order.

## CORE block
The following parameters have to be specified by the user; otherwise you’ll get a fatal error.

### ipre (int)
Pre-processing flag (1: on; 0: off). `ipre=0`: normal run.

Pre-processing flag is very useful for checking integrity of the horizontal grid and some inputs. `ipre=1`: code will output centers.bp, sidecenters.bp, (centers build point, sidcenters build point), and `mirror.out` and stop. Check errors in `fatal.error`. 

!!!important 
    `ipre/=0` only works for single CPU! If you use scribed I/O, make sure the number of 'computes' is 1, plus extra for scribes. Also under this option, the pre-processing run may finish without a clean exist from all cores. Check `outputs/` (`mirror.out` and `fatal.error`); if the run is finished you can manually kill it.

### ibc (int), ibtp (int)
Barotropic/baroclinic flags.

If `ibc=0`, a baroclinic model is used and regardless of the value for `ibtp`, the transport equation is solved. If `ibc=1`, a barotropic model is used, and the transport equation may (when `ibtp=1`) or may not (when `ibtp=0`) be solved; in the former case, S and T are treated as passive tracers.

### rnday (double)
Total simulation time in days.

### dt (double)
Time step in seconds.

### msc2 (int), mdc2 (int)
These two parameters are only used if the wave module WWM is invoked (`USE_WWM` is on and `icou_elfe_wwm=1`).
 The values represent the spectral resolution used in WWM and must match those in 
[wwminput.nml](https://github.com/schism-dev/schism/blob/master/sample_inputs/wwminput.nml.WW3);

### eco_class, ntracer_gen, ntracer_age, sed_class (int)
These parameters set the # of tracer ‘classes’ for each tracer module (EcoSim, GEN, AGE and SED3D), and are required if these modules are invoked in makefile. Note that other tracers modules (ICM, CoSiNE) set their own # of classes.

### nspool, ihfskip (int)
These two flags control the global netcdf outputs. Output is done every `nspool` steps, and a new output stack is opened every `ihfskip` steps.

## OPT block
The optional parameters below are explained in alphabetical order. The default values can be seen below and also in the sample file (sample_inputs/).

### dtb_min=10, dtb_max=30 (double)
Min/max sub-steps allowed in btrack; actual sub-steps are calculated based on local gradients.

### flag_ic(:)=1 (int array)
Options for specifying initial tracer field for cold start, where each array entry corresponds to individual tracer model (e.g. TEM, SAL, SED etc). If `flag_ic=1`, a vertically homogeneous but horizontally varying initial tracer field is specified in inputs like `temp.ic`, `salt.ic`, `[MOD]_hvar_[class #].ic` etc. If `flag_ic=2`, a horizontally homogeneous but vertically varying initial tracer field, prescribed in a series of z-levels, is specified in inputs like `ts.ic`, `[MOD]_vvar_[class #].ic`. For more general 3D initial tracer fields, use the hot start option.

### h0=0.01 (double)
Minimum depth (in m) for wetting and drying (recommended value: `1cm`). When the total depth is less than `h0`, the corresponding nodes/sides/elements are considered dry. It should always be positive.

### h[1,2]_bcc=50,100 (double)
Option on how the baroclinic gradient is calculated below bottom. The 'below-bottom' gradient is zeroed out if `h>=h2_bcc` (i.e. like Z) or uses constant extrapolation (i.e. like terrain-following) if `h<=h1_bcc(<h2_bcc)`. A linear transition is used if the local depth `h1_bcc<h<h2_bcc`.

### ibcc_mean=0 (int)
Mean T,S profile option. If `ibcc_mean=1` (or `ihot=0` and `icst=2`), mean T/S profile is read in from `ts.ic`, and will be removed when calculating baroclinic force. No `ts.ic` is needed if `ibcc_mean=0`.

### ic_elev=0 (int), nramp_elev=0 (int) (int)
Elevation initial condition flag for cold start only (`ihot=0`). If `ic_elev=1`, `elev.ic` (in `*.gr3` format) 
is needed to specify the I.C. Otherwise elevation is initialized to 0 everywhere (cold start only). 
If `ic_elev=1`, the user can ramp up the elevation
 smoothly starting from the specified `elev.ic` or `hotstart.nc` (if $ihot\neq 0$) by setting `nramp_elev=1` (and the ramp-up
 period in this case is `dramp`).

### icou_elfe_wwm=0, iwbl=0 (int)
Coupler flag with WWM; needed if USE_WWM pre-processor is enabled. `icou_elfe_wwm = 0`: no feedback from WWM to SCHISM (decoupled); `1`: coupled SCHISM-WWM.

`iwbl=1`: modified Grant-Madsen formulation for wave enhanced boundary layer; =2: Soulsby (1997) formulation; `=0`: off.

If `icou_elfe_wwm=1`, additional parameters are:

- `nstep_wwm=1` (int): call WWM every this many time steps;
- `hmin_radstress=1.0` (double): min. total water depth used only in radiation stress calculation; the radiation stress is zero if `local depth <hmin_radstress`.

### ics=1 (int)
Coordinate frame flag. If `ics=1`, Cartesian coordinates are used; if `ics=2`, both `hgrid.ll` and `hgrid.gr3` use degrees latitude/longitude (and they should be identical to each other in this case).

### i_hmin_airsea_ex (int), hmin_airsea_ex (double)
Option to locally turn off heat exchange.

- `i_hmin_airsea_ex=1`: exchange turned off if `local grid depth<hmin_airsea_ex`
- `i_hmin_airsea_ex=2`: exchange turned off if `local water depth<hmin_airsea_ex`

### i_hmin_salt_ex (int), hmin_salt_ex (double)
Simialr tpo `i_hmin_airsea_ex` and `hmin_airsea_ex`.

### ielm_transport = 0, max_subcyc = 10 (int)
Hybrid ELM-FV transport for performance; used only with `itr_met>=3`. If `ielm_transport=1`, 
the hybrid scheme is invoked and `max_subcyc` represents the max # of subcycling per time step in 
transport allowed; if the actual # of subcycling in a prism at a time step exceeds this threshold, 
more efficient ELM transport is used locally (at the expense of mass conservation, 
so make sure this option is used sparingly).

### ieos_type=0, ieos_pres=0 (int)
By default, use the nonlinear equation of state: `ieos_type==0`. If the potential temperature is used, the pressure effect has been accounted for: `ieos_pres=0.`

### if_source=0 (int), dramp_ss=2. (double), lev_tr_source(:)=-9 (int array)
Point sources/sinks option (0: no; 1: on). If `if_source=1`, needs `source_sink.in`, `vsource.th`, `vsink.th`, 
and `msource.th` (see sample files in the source code directory src/ for their formats). 
If `if_source=-1`, the input is `source.nc`, which includes element list inside and allows for different time 
steps and # of records for volume/mass source/sinks. If `if_source/=0`, specify ramp-up period (in days) with `dramp_ss` (no
 ramp-up if <=0). 

The tracers are injected into an element at a particular level, as specified by `lev_tr_source(1:ntr)` (where `ntr` is
 # of tracer modules, i.e. 1 input for each module). The code will extrapolate below bottom/above surface if necessary, 
 so e.g., '-9' means bottom.

### iflux=0 (int)
Parameter for checking volume and salt conservation. If turned on (`=1`), the conservation will be checked in regions specified by `fluxflag.prop`.

### iharind=0 (int)
Harmonic analysis flag. If `iharind≠0`, an input `harm.in` is needed.

### ihconsv=0, isconsv=0 (int)
Heat budget and salt conservation models flags. If `ihconsv=0`, the heat budget model is not used. If `ihconsv=1`, a heat budget model is invoked, and a number of netcdf files for radiation flux input are read in from `sflux/sflux_rad*.nc`. If `isconsv=1`, the evaporation and precipitation model is evoked but the user needs to turn on the pre-processing flag `PREC_EVAP` in makefile and recompile. In this case, `ihconsv` must be `1`, and additional netcdf inputs for precipitation (`sflux/sflux_prc*.nc`) are required.

### ihdif=0 (int)
Flag to use non-zero horizontal diffusivity. If `ihdif=0`, it is not used. If `ihdif≠0`, input `hdif.gr3` is needed.

### ihot=0 (int)
Hot start flag. If `ihot=0`, cold start; if `ihot≠0`, hot start from `hotstart.nc`. If `ihot=1`, the time and time step are reset to zero, and outputs start from `t=0` accordingly (and you need to adjust other inputs like `.th` etc). If `ihot=2`, the run (and outputs) will continue from the time specified in `hotstart.nc`, and in this case, you do not need to adjust other inputs.

### ihydraulics=0 (int)
Hydraulic model option. If `ihydraulics≠0`, `hydraulics.in` is required.


### iloadtide=0 (int)  loadtide_coef (double)
Option to specify Self Attraction and Loading (SAL) tide, usually used for basin-scale applications. 
If `iloadtide=0`, SAL is off. If `iloadtide=1`, the SAL input is interpolated values from a tide database,
 e.g., FES2014, given in `loadtide_[FREQ].gr3`, where `[FREQ]` are frequency names (shared with 
tidal potential, in upper cases like M2) and the two 'depths' inside are amplitude (m) 
and phases (degrees behind GMT). In this option, SAL is
 lumped into tidal potential so it shares some parameters with tidal potential
 in bctides.in (cut-off depth, frequencies).

If iloadtide=2 or 3, use a simple scaling for gravity approach (in this option,
 SAL is applied everywhere and does not share parameters with tidal potential).
If `iloadtide=2`, a simple scaling is used to reduce 
the gravity. If `iloadtide=3`, the scaling is dependent on the local depth _a la_ Stepanov & Hughes (2004),
 with a maximum value of `loadtide_coef`.

### imm=0, ibdef=10 (int)
Bed deformation option. Default: `0` (no bed deformation); `1`: with bed deformation (needs `ibdef` (# of steps during which deformation occurs), and `bdef.gr3`); 2: 3D bottom deformation (need to interact with code).

### indvel=0 (int), ihorcon=0 (int), hvis_coef0=0.025 (double), ishapiro=1, niter_shap=1 (int), shapiro0=0.5 (double)
These parameters (and `inter_mom` below) control the numerical dissipation in momentum solver; see SCHISM paper (Zhang et al. 2016) for details.

`indvel` determines the method of converting side velocity to node velocity. If `indvel=0`, the node velocity is allowed to be discontinuous across elements and additional viscosity/filter is needed to filter out grid-scale noises (spurious 'modes'). If `indvel=1`, an inverse-distance interpolation procedure is used instead and the node velocity is continuous across elements; this method requires no additional viscosity/filter unless kriging ELM is used (`inter_mom >0`). In general, `indvel=0` leads to smaller numerical dissipation and better accuracy, but does generally require a velocity B.C.

Due to spurious modes or dispersion (oscillation), viscosity/filter should be applied. `ihorcon=0`:no horizontal viscosity; `=1`: Laplacian (implemented as a filter); `=2`: bi-harmonic. For `ihorcon/=0`, `hvis_coef0` specifies the non-dimensional viscosity. In addition to the viscosity, one can add the Shapiro filter, which is specified by `ishapiro` =0,±1, 2 (turn off/on Shapiro filter). If `ishapiro=1`, `shapiro0` specifies the Shapiro filter strength. If `ishapiro=-1`, an input called `shapiro.gr3` is required which specifies the filter strength at each node. If `ishapiro=2`, a Smagorinsky-like filter is applied and `shapiro0` is a coefficient $(\gamma_0)$, which is on the order of $10^3$:

\begin{equation}
\label{eq01}
\begin{aligned}
\gamma &= 0.5\tanh (\gamma_0 \Delta t \hat \mu)\\
\hat\mu &= \sqrt{u_x^2 + v_y^2 + \frac{(u_y + v_x)^2}{2}}
\end{aligned}
\end{equation}

If `ishapiro/=0`, `niter_shap` specifies the number of times the filter is applied. Note that for non-eddying regime applications (nearshore, estuary, river), an easiest option is: `indvel=0`, `ishapiro=1` (`shapiro0=0.5`), `ihorcon= inter_mom=0`.

For applications that include the eddying regime, grid resolution in the eddying regime needs to vary smoothly (Zhang et al. 2016), and the user needs to tweak dissipation carefully. A starting point can be: `indvel=ishapiro=inter_mom=0`, `ihorcon=2`, `hvis_coef0=0.025`. If the amount of dissipation is insufficient in the non-eddying regime, consider using `ishapiro=-1`, with an appropriate `shapiro.gr3` to turn on Shapiro filter locally to add dissipation, or use `ishapiro=2` and `shapiro0=1000`.

### inter_mom=0, kr_co=1 (int)
Interpolation method at foot of characteristic line during ELM. `inter_mom=0`: linear interpolation; `=1`: dual kriging method. If `inter_mom=-1`, the depth in `krvel.gr3` (0 or 1) will determine the order of interpolation (linear or kriging). If the kriging ELM is used, the general covariance function is specified in `kr_co`: 1: linear $f(h)=-h$; 2: $(h^2*log(h)$; 3: cubic $(h^3)$; 4: $(-h^5)$.

In general, `indvel=0` should be used with `inter_mom=0` or `1` to avoid large dispersion (with additional viscosity/filter also). `indvel=1` can be used with any covariance function without viscosity/filter.

### inu_elev=0, inu_uv=0 (int)
Sponge layer for elevation and velocity. In SCHISM, relaxation/nudging of a generic variable is implemented as:

\begin{equation}
\label{eq02}
\widetilde \varphi = (1-\gamma)\varphi + \gamma\varphi_{target}
\end{equation}

which is a discrete analogue of the restoration equation:

\begin{equation}
\label{eq03}
\frac{\partial \varphi}{\partial t} = \frac{\gamma}{\Delta t} \left( \varphi_{target} - \varphi \right)
\end{equation}

If `inu_elev=0`, no relaxation is applied to elevation. If `inu_elev=1`, relaxation constants are specified in `elev_nudge.gr3` (depth=0 means no relaxation, depth=1 means strongest nudging) and the elevations are relaxed toward 0. Similarly for `inu_uv` (with input `uv_nudge.gr3`).

### inu_tr(:)=0 (int array), nu_sum_mult(int), step_nu_tr=86400. (double)
Nudging flag for tracer models (e.g. temperature), and nudging step (in sec). When `inu_tr=0`, no nudging is done.

When `inu_tr=1`, relax back to initial conditions.

When `inu_tr=2`, nudge to values specified in `[MOD]_nu.nc`, which has a time step of `step_nu_tr`.

If `inu_tr≠0`, the horizontal relaxation factors are specified in `[MOD]_nudge.gr3` (as depths info), 
 and the vertical relaxation factors are specified as a linear function of depths with: `vnh[1,2]` 
(transitional depths) and `vnf[1,2]` (relaxation constants at the 2 depths). The final relaxation 
constant is either the sum (if `nu_sum_mult=1`) or product (if `nu_sum_mult=2`) of the two, 
i.e. (horizontal `+ or *` vertical relaxation factors) times `dt`.

### inunfl=0 (int)
Choice of inundation algorithm. `inunfl=1` can be used if the horizontal resolution is fine enough, and this is critical for tsunami simulations. Otherwise use `inunfl=0`.


### isav=0 (int)
Parameters for submerged or emergent vegetation. If `isav=1` (module on), you need to supply 4 extra inputs: `sav_cd.gr3` (form drag coefficient), `sav_D.gr3` (depth is stem diameter in meters); `sav_N.gr3` (depth is # of stems per m2); and `sav_h.gr3` (height of canopy in meters).

### itr_met=3 (int), h_tvd=5. (double)
Transport option for all tracers. `itr_met=3` for $TVD^2$, and `itr_met=4` for 3rd order WENO. **Note that `itr_met=1,2` (explicit) is only kept for comparison and benchmark purpose and so should NOT be used normally- use ‘3’ or ‘4’ instead**. If `itr_met>=3`, 1 extra parameter is needed: `h_tvd` which specifies the transition depth (in meters) between upwind and TVD schemes; i.e. more efficient upwind is used when the `local depth < h_tvd`. Also in this case, additional toggle between upwind and TVD is specified in `tvd.prop`. The TVD limiter function is specified in `TVD_ LIM` in `mk/include_modules` (for efficiency purpose). If itr_met=3, 2 tolerances are also required (use recommended values).

### itur=0 (int)
Turbulence closure model selection.

If `itur=0`, constant diffusivities are used for momentum and transport, and the diffusivities are specified in `dfv0`, `dfh0`.

If `itur=-2`, vertically homogeneous but horizontally varying diffusivities are used, which are read in from `hvd.mom` and `hvd.tran`.

If `itur=-1`, horizontally homogeneous but vertically varying diffusivities are used, which are read in from `vvd.dat`.

If `itur=2`, the zero-equation Pacanowski and Philander closure is used. In this case, a few extra parameters are required: `h1_pp`, `vdmax_pp1`, `vdmin_pp1`, `tdmin_pp1`, `h2_pp`, `vdmax_pp2`, `vdmin_pp2`, `tdmin_pp2`. Eddy viscosity is computed as: $\text{vdiff}=\text{vdiff_max}/(1+\text{rich})^2+\text{vdiff_min}$, and diffusivity $\text{tdiff}=\text{vdiff_max}/(1+\text{rich})^2+\text{tdiff_min}$, where $\text{rich}$ is a Richardson number. The limits (`vdiff_max`, `vdiff_min` and `tdiff_min`) vary linearly with depth between depths `h1_pp` and `h2_pp`.

If `itur=3`, then the two-equation closure schemes from the GLS model of Umlauf and Burchard (2003) are used. In this case, 2 additional parameters are required: `mid`, `stab`, which specify the closure scheme and stability function used: `mid=` `MY` is Mellor & Yamada; `=KL` is GLS as k-kl; KE is GLS as $k-\varepsilon$ =KW is GLS as $k-\omega$ =UB is Umlauf & Burchard's optimal; `stab=` GA is Galperin's clipping (only for MY); =KC is Kantha & Clayson's stability function) Also the user needs to specify max/min diffusivity/viscosity in `diffmax.gr3` and `diffmin.gr3`, as well as a surface mixing length scale constant `xlsc0`.

If `itur=4`, GOTM turbulence model is invoked; the user needs to compile the GOTM libraries first 
(see README inside GOTM/ for instructions), and turn on pre-processing flag `USE_GOTM` in makefile 
and recompile. In this case, the minimum and maximum viscosity/diffusivity are still specified in `diffmin.gr3` 
and `diffmax.gr3`. In addition, GOTM also requires an input called `gotmturb.inp`. There are some 
ready-made samples for this input in the source code bundle. If you wish to tune some parameters 
inside, you may consult [gotm.net](https://gotm.net) for more details.

### level_age(:)=-999 (int array)
If `USE_AGE` is on, this array specifies the vertical level indices used to inject age tracers. Use -999 to inject the tracer at all levels.

### meth_sink=0 (int)
Option for sinks. If `meth_sink =1`, the sink value is reset to `0` if an element is dry with 
a net sink value locally to prevent further drawdown of groundwater.

### nadv=1 (int)
Advection on/off option. If `nadv=0`, advection is selectively turned off based on the input file `adv.gr3`. 
If `nadv=1` or `2`, advection is on for the whole domain, and backtracking is done using either 
Euler or 2nd-order Runge-Kutta scheme.

### nchi=0 (int)
Bottom friction option. If `nchi=-1`, and Manning's n is specified in `manning.gr3`. If `nchi=0`, spatially varying drag coefficients are read in from `drag.gr3` (as depth info). For `nchi=1`, bottom roughnesses (in meters) are read in from `rough.gr3`.

If `nchi=-1`, an additional parameter is required: `hmin_man` (in meters) which sets the minimum depth used in the Manning formulation.

If `nchi=1`, one additional parameter is required: `dzb_min` (in meters). In this case the drag coefficients are calculated using the log drag law when the bottom cell thickness $\delta_b>=\text{dzb_min}$. when $\delta_b<\text{dzb_min}$, $\text{Cd}=\text{Cdmax}$, where $\text{Cdmax}=\text{Cd}(\delta_b=\text{dzb_min})$. This is to avoid exaggeration of $\text{Cd}$ in very shallow water.

### ncor=0 (int)
Coriolis option. If `ncor=0` or `-1`, a constant Coriolis parameter is specified. If `ncor=0`, `coricoef` specifies the Coriolis factor. If `ncor=-1`, `rlatitude` specifies the mean latitude used to calculate the Coriolis factor.

If `ncor=1`, a variable Coriolis parameter, based either on a beta-plane approximation (`ics=1`) or on the latitude-dependent Coriolis (`ics=2`), is used, with the lat/lon coordinates read in from `hgrid.ll`. For `ics=1`, the center of beta-plane approximation must be correctly specified in `sfea0`.

###  dramp=1. (double), drampbc=1. (double)
Ramp periods in days for the tides, B.C. or baroclincity. 
If `ibc=0`, the ramp-up for baroclinicity is specified with `drampbc` (in days). 
Turn off ramp-up by setting the ramp-up period <=0. 
The ramp function is a hyperbolic tangent function: $f(t) = \tanh(2t/86400/\text{drampbc})$.

### nws=0 (int), drampwind=1. (double), iwind_form=-1 (int), iwindoff(int), wtiminc=dt (double)
Wind forcing options and the interval (in seconds) with which the wind input is read in. If `nws=0`, no 
wind is applied (and `wtiminc` becomes unused). If `nws=1`, constant wind is applied to the whole domain 
at any given time, and the time history of wind is read in from `wind.th`. If `nws=2`, spatially 
and temporally variable wind is applied and the input consists of a number of netcdf files in the directory
 `sflux/`. The option `nws=3` is reserved for coupling with atmospheric model via ESMF caps. If `nws=4`, 
the required input `wind.th` specifies wind and pressure at each node and at time of multiple of `wtiminc`. 
If `nws=-1`, use Holland parametric wind model (barotropic only with wind and atmos. pressure).
 In this case, the Holland model is called every step so wtiminc is not used. An extra
 input is needed: `hurricane-track.dat`.


If `nws>0`, the ramp-up period (in days) is specified with `drampwind`. Also
 the user has the option to scale the wind speed using `iwindoff`=1 (which requires an additional input `windfactor.gr3`).

The wind stress formulation is selected with `iwind_form`. If `nws=1` or `4`, or `nws=2 && ihconsv=0`, 
or `nws=2 && iwind_form=-1`, the stress is calculated from Pond & Pichard formulation (originally from Garret). 
If `nws=1` or `4`, or `nws=2 && ihconsv=0`, or `nws=2 && iwind_form=1`, the stress is calculated from 
Hwang (2018) formulation. If `nws=2, ihconsv=1 && iwind_form=0`, the stress is calculated 
from heat exchange routine. If `WWM` is enabled and `icou_elfe_wwm > 0` and `iwind_form=-2`, 
stress is calculated by WWM; otherwise the formulations above are used.

### rho0=1000, shw=4184. (double)
Reference water density for Boussinesq approximation and specific heat of water in J/kg/K.

### rmaxvel=10. (double)
Maximum velocity. This is needed mainly for the air-water exchange as the model may blow up if the water velocity is above 20m/s.

### s1_mxbnt=0.5, s2_mxnbnt=3.5 (double)
Dimensioning parameters used in inter-subdomain backtracking. Start from `s[12]_mxnbt=0.53`, and increase them (gradually) if you get a fatal error like “btrack: overflow”. Accuracy is not affected by the choice of these two parameters.

### slam0=-124, sfea0=45 (double)
Centers of projection used to convert lat/lon to Cartesian coordinates. These are used if a variable Coriolis parameter is employed (ncor=1).

### start_year=2000, start_month=1, start_day=1 (int), start_hour=0, utc_start=8 (double)
Starting time for simulation. `utc_start` is hours behind the GMT, and is used to adjust time zone. For example, `utc_start=5` is US Eastern Time, and `utc_start= -8` is Beijing Time.

### thetai=0.6 (double)
Implicitness parameter (between 0.5 and 1). Recommended value: 0.6.

## SCHOUT block
### iout_sta=0, nspool_sta=10 (int)
Station output flag. If `iout_sta≠1`, an input `station.in` is needed. In addition, `nspool_sta` specifies the spool for station output.

### nc_out =1(int)
Main switch to turn on/off netcdf outputs, useful for other programs to control outputs.

### nhot=0, nhot_write=8640 (int)
Hot start output control parameters. If `nhot=0`, no hot start output is generated. If `nhot=1`, 
hotstart output is named `outputs/hotstart_[process_id]_[time_step].nc` every `nhot_write` 
steps, where `it` is the corresponding time iteration number. `nhot_write` must be a multiple of 
`ihfskip`. If you want to hotstart a run from step `it`, you need to combine all process-specific 
hotstart outputs into a `hotstart.nc` using `combine_hotstart7.f90` (`./combine_hotstart7 –h` for help).

### iof_* (int)
Global output (in netcdf4 format) options, where `*` stands for module name (e.g. "hydro", "wwm" etc). 
The frequency of global outputs is controlled by 2 parameters in [CORE](#core): 
[`nspool`](#nspool-ihfskip-int) and [`ihfskip`](#nspool-ihfskip-int). Output is done every
 `nspool` steps, and a new output stack is opened every `ihfskip` steps. 
Under OLDIO, the outputs are named as `outputs/schout_[MPI process id]_[1,2,3,...].nc` etc. The combine scripts 
are then used to gather each output variable across all MPI processes into a single output, e.g., `schout_[1,2,3…].nc`.
With new scribed I/O, outputs look like `out2d_1,2,3…].nc` etc (and no combining is necessary).

Each output variable is controlled by `1` flag in `param.nml`. We only show a few examples below; 
the rest are similar. Note that variables may be centered at nodes/sides/elements horizontally and 
whole/half levels vertically. However, at the moment most variables are centered at nodes and whole levels,
 and most post-processing FORTRAN scripts can only handle this type of outputs.

```
iof_hydro(1) = 1 !global elevation output control. If iof_hydro(1)=0, no global elevation is recorded. 
                 !If iof_hydro(1)= 1, global elevation for each node is recorded.
                 !The output is either starting from scratch or appended to existing ones depending
                 !on ihot.
```

Some outputs are conditional upon you turn on certain module; e.g. `iof_sed(7) = 1` won’t output 
the bottom depth change unless you turn on `USE_SED` in makefile.

Some ‘native’ variables (e.g., element- or side-centered) are:

```
iof_hydro(27) = 1 !horizontal velocity defined at side [m/s]. These are the original velocity inside SCHISM
```
