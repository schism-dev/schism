SCHISM (Semi-implicit Cross-scale Hydroscience Integrated System Model) is a 3D seamless cross-scale model grounded on unstructured grids (UG) for hydrodynamics and ecosystem dynamics. SCHISM is a derivative work from the original SELFE model (v3.1dc as of Dec. 13, 2014). SCHISM has been implemented by Dr. Joseph Zhang (College of William & Mary) and other developers around the world, and licensed under Apache. SELFE was developed at the Oregon Health Sciences University. However, there are now significant differences between the two models.

# Major characteristics of SCHISM
- Finite-element/finite-volume formulation
- Unstructured mixed triangular/quadrangular grid in the horizontal dimension
- Hybrid SZ coordinates or new $LSC^2$ in the vertical dimension
- Polymorphism: a single grid can mimic 1D/2DV/2DH/3D configurations
- Semi-implicit time stepping (no mode splitting): no CFL stability constraints → numerical efficiency
- Robust matrix solver
- Higher-order Eulerian-Lagrangian treatment of momentum advection
- Natural treatment of wetting and drying suitable for inundation studies
- Mass conservative and monotonicity-enforced transport: TVD2
- No bathymetry smoothing necessary
- Very tolerant of bad-quality meshes in the non-eddying regime

# Modeling system and application area
- 3D baroclinic cross-scale lake-river-estuary-plume-shelf-ocean circulations
- Tsunami hazards
- Storm surge
- Sediment transport
- Biogeochemistry/ecology/water quality
- Oil spill
- Short wave-current interaction

# Main differences between SCHISM and SELFE
1) Apache license
2) Vertical grid system ($LSC^2$, Zhang et al. 2015);
3) Mixed triangular-quadrangular horizontal grid;
4) Implicit vertical advection scheme for transport ($TVD^2$);
5) Higher-order transport in the horizontal dimension (3rd-order WENO);
6) Advection scheme for momentum: optional higher-order kriging with ELAD filter;
7) A new horizontal viscosity scheme (including bi-harmonic viscosity) to effectively filter out inertial spurious modes without introducing excessive dissipation.

# Citation
We suggest the following language for citing the model:
SCHISM is a derivative product built from the original SELFE (v3.1dc; Zhang and Baptista 2008), with many enhancements and upgrades including new extension to large-scale eddying regime and a seamless cross-scale capability from creek to ocean (Zhang et al. 2016).
- Zhang, Y. and Baptista, A.M. (2008) SELFE: A semi-implicit Eulerian-Lagrangian finiteelement model for cross-scale ocean circulation", Ocean Modelling, 21(3-4), 71-96.
- Zhang, Y., Ye, F., Stanev, E.V., Grashorn, S. (2016). Seamless cross-scale modeling with SCHISM, Ocean Modelling, 102, 64-81. doi:10.1016/j.ocemod.2016.05.002

# How to read this manual
This manual contains detailed information on physical and numerical formulations, as well as usage for SCHISM. For beginners, we suggest you familiarize yourself with the basic notations in Chapters 2 and 3 but skip some details in those two chapters; there is also a ‘cheat sheet’ near the end of this chapter. Chapters 4&5 describe how to set up the model, including grid generation, and so should be read carefully, in consultation with appropriate sections in the previous 2 chapters.
Since SCHISM is quite a sophisticated package, we strongly recommend you start from the simple test cases in Chapter 5.5 and gradually progress toward more complex 3D baroclinic or coupled applications.

# Notation used in this manual
We will use bold characters to denote vectors and matrices, and unbold characters to denote scalars in mathematical equations. In addition, superscripts usually denote time step and subscripts denote spatial locations. E.g., $T_{i,k}^{n+1} may mean the temperature at step n+1 (i.e., new time step) and prism (i,k), where i is the element number and k is the (whole) vertical index. We will use red words to denote input file names (e.g. param.nml), purple words to denote output file names (e.g. mirror.out), and green words to denote parameters specified in param.nml; e.g. nchi (the bottom friction option flag). Code/pseudo-code fragments are written in italic. 

Below are some notations used in this manual:

$N_p$: # (number) of nodes
$N_s$: # $ of sides
$N_e$: # of elements
$N_z$: maximum # of vertical levels
$i34(j)$: type of an element j (3: triangle; 4: quad)
$Nb(i)$: # of surrounding elements of a node i;
$kbp(i)$: bottom index as seen by a node i
$kbs(i)$: bottom index as seen by a side i
$kbe(i)$: bottom index as seen by an element i
$A$: area of an element
$\Delta z$: layer thickness (at a node, side or elem.)
$\delta_{ij}$: Dirac’s Delta function (=1 when i=j; 0 otherwise)