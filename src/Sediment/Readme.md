In order of invocation in Hydro code:

sed_mod - parameters for sediment model

read_sed_input - subroutine to read the sediment model inputs written in sediment.in file; optionally
                 calculate settling vel and critical shear stress;

sed_init - subroutine to set initial conditions for sediment tracer fields and initial bed conditions

sed_roughness (inside sed_friction.F90) - calculate bottom roughness

sediment - routine to compute the sediment source and sinks. It includes the vertical settling of 
           the sediment, erosional and depositional flux, transport of multiple grain sizes, bedload 
           based on Meyer Peter and Mueller (not active yet), van Rijn, Soulsby and Damgaard (2005)
           or Wu and Lin (2014). Specific bedload transport caused by wave asymmetry (acceleration
           skewness) can be included.
           Morphological model is found inside.
           Vertical settling is handled with semi-Lagrangian PPM and eventually cast into a body force for
           the transport solver.

[other routines:
sed_bedload - bedload routines (including call to solve_jcg() for morphol. model)
sed_friction - bottom formulations
sed_filter - morphological filters
sed_misc_subs - misc routines
sed_poro - compute sediment porosity
]

----------------------------------------------------------------------------------------------------
Joseph's notes on doing SED simulations:
(1) using a very small dzb_min (i.e. unlimit the Cd) will lead to reduced near -bottom vel and stress -> little erosion!
(2) From erosion & deposition tests, the suspended load is sensitive to vertical diffusion. Try itur=3: k-kl and diffmin=1.e-4.
