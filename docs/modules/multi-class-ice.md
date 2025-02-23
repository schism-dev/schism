#Theory
The column physics package of the Los Alamos sea ice model (CICE), the Icepack (v1.3.4) has been couple with SCHISM. Besides the zero-layer thermodynamics, two more sophisticated thermodynamic formulations, the Bitz and Lipscomb (1999; BL99) thermodynamics formulation for constant salinity profiles, and the mushy layer thermodynamics formulation for evolving salinity (Turner et al., 2013), are also implemented. At the sub-grid scale, thin ice and thick ice coexist, and therefore an ice thickness distribution (ITD; Lipscomb, 2001; Bitz et al., 2001) has been implemented in order to describe the unresolved spatial heterogeneity of the thickness field. The ITD offers a prognostic statistical description of the sea ice thickness, which it divides into multiple categories, along with the ice area fraction corresponding to each category – a more detailed approach than the single fraction used in the previous implementation. More tracers and more ice processes are added in this coupled model by Icepack, including multiple melt pond parameterizations (Hunke et al., 2013) and a mechanical redistribution parameterization (Lipscomb et al., 2007) that responds to sea ice convergence by piling up thin sea ice and therefore mimicking ridging and rafting events. The interaction between the shortwave radiation and the sea ice in Icepack is addressed using two formulations: the Community Climate System Model (CCSM3) formulation, which relates the surface albedo to the surface sea ice temperature, and the delta–Eddington formulation (Briegleb et al., 2007), which relates the albedo to inherent optical properties of sea ice and snow. The dynamic solver is not included in Icepack and we used two approaches: (1) the classic elastic–viscous–plastic method (EVP; Hunke and Dukowicz, 1997) and (2) the modified elastic–viscous–plastic method (mEVP; Kimmritz et al., 2015). Both methods are inherited from the previous single-class ice and snow formulation (Zhang et al., 2023).

#Usage
1. Compile with USE_MICE and USE_EVAP on;
2. Provide snow flux data in sflux/sflux_prc* files (named as srate), the format is the same as rainfall data (prate);
3. The main parameter input for this module is mice.nml and namelist.icepack (you can find a sample in sample_input/);
4. For the parameters in the mice.nml, ice_advection=6 (hybrid TVD-upwind) is recommended. For the parameters in namelist.icepack, you can find them all in the [Icepack manual](https://github.com/CICE-Consortium/Icepack/wiki/Icepack-Release-Table).
5. For the output, all variables for iof_mice in param.nml are enabled under new scribe IO, but only the variables in namelist.icepack are enabled under OLDIO; this needs to be updated in newer versions.
6. Optional ihot_mice in mice.nml to allow different hotstart modes between ice and hydro. In general, ihot_mice works like ihot. ihot_mice=0 is the cold start. Under ihot_mice=1 the initial file (hotstart.nc) should contain ice concentration, thickness and velocity and from external model like HYCOM (scripts is gen_hot_from_hycom_ice.f90). ihot_mice=2 is a restart function like ihot=2, where the ice variables are stored in outputs/hotstart*. The routine for this is inside icedrv_io.f90 (restart_icepack), but this mode is not stable if there are too many variables, like floe size distribution (tr_fsd).


#References
Bitz, C. M. and Lipscomb, W. H.: An energy-conserving thermodynamic model of sea ice, J. Geophys. Res.-Oceans, 104, 15669–15677, https://doi.org/10.1029/1999JC900100, 1999. 

Bitz, C. M., Holland, M. M., Weaver, A. J., and Eby, M.: Simulating the ice-thickness distribution in a coupled climate model, J. Geophys. Res.-Oceans, 106, 2441–2463, https://doi.org/10.1029/1999JC000113, 2001. 

Briegleb, B. P. and Light, B.: A Delta-Eddington multiple scattering parameterization for solar radiation in the sea ice component of the Community Climate System Model, Tech. Rep. NCAR/TN 472+STR, National Center for Atmospheric Research, Boulder, Colorado USA, https://doi.org/10.5065/D6B27S71, 2007. 

Hunke, E., Allard, R., Bailey, D. A., Blain, P., Craig, A., Dupont, F., DuVivier, A., Grumbine, R., Hebert, D., Holland, M., Jeffery, N., Lemieux, J.-F., Osinski, R., Rasmussen, T., Ribergaard, M., Roach, L., Roberts, A., Turner, M., and Winton, M.: CICE-Consortium/Icepack: Icepack 1.3.4 (1.3.4), Zenodo [code], https://doi.org/10.5281/zenodo.8336034, 2023. 

Hunke, E. C. and Dukowicz, J. K.: An elasticviscous-plastic model for sea ice dynamics, J. Phys. Oceanogr., 27, 1849–1867, https://doi.org/10.1175/1520-0485(1997)027<1849:AEVPMF>2.0.CO;2, 1997. 

Kimmritz, M., Danilov, S., and Losch, M.: On the convergence of the modified elastic-viscous-plastic method for solving the sea ice momentum equation, J. Comput. Phys., 296, 90–100, https://doi.org/10.1016/j.jcp.2015.04.051, 2015. 

Lipscomb, W. H.: Remapping the thickness distribution in sea ice models, J. Geophys. Res., 106, 13989–14000, https://doi.org/10.1029/2000JC000518, 2001. 

Lipscomb, W. H. and Hunke, E. C.: Modeling sea ice transport using incremental remapping, Mon. Weather Rev., 132, 1341–1354, https://doi.org/10.1175/1520-0493(2004)132<1341:MSITUI>2.0.CO;2, 2004. 

Turner, A. K., Hunke, E. C., and Bitz, C. M.: Two modes of sea-ice gravity drainage: A parameterization for largescale modeling, J. Geophys. Res.-Oceans, 118, 2279–2294, https://doi.org/10.1002/jgrc.20171, 2013. 

Zhang, Y. J., Wu, C., Anderson, J., Danilov, S., Wang, Q., Liu, Y., and Wang, Q.: Lake ice simulation using a 3D unstructured grid model, Ocean Dynam., 73, 219–230, https://doi.org/10.1007/s10236-023-01549-9, 2023. 

