# How to read this manual
This manual contains detailed information on physical and numerical formulations, as well as usage for 
SCHISM. For beginners, we suggest you familiarize yourself with the basic notations in 
[Physical formulation](physical-formulation.md) and [Numerical formulation](geometry-discretization.md) 
but skip some details in those two chapters; there is also a ‘cheat sheet’ in 'Typical workflow'. [Getting started](#) describe how to set up the model, including grid generation, and 
so should be read carefully, in consultation with appropriate sections in the previous 2 chapters.
Since SCHISM is quite a sophisticated package, we strongly recommend you start from the simple 
[test cases](../verification-tests.md) and gradually progress toward more complex 3D baroclinic or coupled applications.

# Notation used in this manual
We will use bold characters to denote vectors and matrices, and unbold characters to denote scalars in mathematical equations. In addition, superscripts usually denote time step and subscripts denote spatial locations. E.g., $T_{i,k}^{n+1}$ may mean the temperature at step $n+1$ (i.e., new time step) and prism $(i,k)$, where $i$ is the element number and $k$ is the (whole) vertical index. We will use inline code blocks to denote input file names (e.g. `param.nml`), purple words to denote output file names (e.g. mirror.out), and green words to denote parameters specified in param.nml; e.g. nchi (the bottom friction option flag). Code/pseudo-code fragments are written in `code-blocks`. 

Below are some notations used in this manual:

$N_p$: number of nodes

$N_s$: number of sides

$N_e$: number of elements

$N_z$: maximum number of vertical levels

$i34(j)$: type of an element $j$ (3: triangle; 4: quad)

$Nb(i)$: number of surrounding elements of a node $i$;

$kbp(i)$: bottom index as seen by a node $i$

$kbs(i)$: bottom index as seen by a side $i$

$kbe(i)$: bottom index as seen by an element $i$

$A$: area of an element

$\Delta z$: layer thickness (at a node, side or elem.)

$\delta_{ij}$: Dirac’s Delta function ($=1$ when $i=j$; $0$ otherwise)
