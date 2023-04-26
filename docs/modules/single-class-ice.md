#Theory

The single-class ice/snow module is taken from FESOM2, courtesy of Dr. Sergey Danilov's group. The info below is from the FESOM2 manual with minor modifications 
 to show some new additions.

Dynamical part solves the ice motion equation:

\begin{equation}
\label{ice01}
 m(\partial_t + \pmb{f} \times)=a\tau - aC_d\rho_o (\pmb{u}-\pmb{u}_0)|\pmb{u}-\pmb{u}_0| +\pmb{F} -mg\nabla H
\end{equation}

Here $m$ is the ice plus snow mass per unit area, $C_d$ the ice-ocean drag coefficient, $\rho_o$ the
water density, $a$ the ice compactness, $u$ and $u_o$ the ice and ocean velocities, $\tau$ the wind
stress, $H$ the sea surface elevation, $g$ the acceleration due to gravity and $F_j = \partial_i σ_{ij}$ is
the force from stresses within the ice. Here we use Cartesian coordinates for brevity, with
$i, j = 1, 2$ implying $x$ and $y$ directions, and the implementation of spherical coordinates will
be discussed further. Summation over repeating coordinate indices is implied. The stress
tensor is symmetric. The mass $m$ is the combination of ice and snow contributions

with ρice and ρs, respectively, the densities of ice and snow and hice and hs their mean
thicknesses (volumes per unit area).
For the visco-plastic (VP) rheology (Hibler, 1979) one writes

where

is the deformation rate tensor, η and ζ are the moduli (‘viscosities’) and P the pressure.
Both the stress and deformation rate tensors are symmetric, so they are characterized by
three independent components. The standard VP rheology adopts the following scheme of
computing pressure P and moduli η and ζ:


 are default values, they are adjusted in practice. In this scheme, 
∆min serves for viscous regularization of plastic behavior in areas where ∆ is very small. We note that recent multi-category
ice implementations (such as CICE, see Hunke and Lipscomb 2008) use different parameterization for P0, which takes into account the distribution of ice over thickness categories. This
does not change the basic equations (1, 2).
In our case we deal with three tracers, the area coverage a, ice volume hice and snow
volume hs. They are advected with ice velocities and modified through thermodynamical
forcing

with Sa, Sice and Ss the sources due to the exchange with atmosphere and ocean.

The system (1), (2) and (3), augmented with appropriate model of sources and boundary conditions, defines the sea ice model. We use the boundary conditions of no slip for
momentum and no flux for tracers at lateral walls. The well known difficulty in solving
these equation is related to the rheology term in the momentum equation, which makes this
equation very stiff and would require time steps of fractions of second if stepped explicitly.
There are two ways of handling this difficulty. The first one treats the rheology part in an
implicit way, with linearization for the moduli, as suggested by Zhang and Hibler (1997). As
mentioned elsewhere (see, e. g., Lemieux and Tremblay, 2009), it does not warrant full convergence, which requires a full nonlinear solver (for example, a Jacobian-free Newton-Krylov
one, see Lemieux et al. 2012). The latter is still too expensive computationally, so the VP
option used by us is similar in spirit to that of Zhang and Hibler (see section 2.4). The second
way is to reformulate the VP approach by adding pseudo-elasticity. It raises the time order
of the system (1, 2) to the second, which makes the CFL limitation on the explicit time step
much less severe than in the original VP framework.

Below is the formulation for a modified elastic-visco-plastic (mEVP) due to Bouillon et al. (2013).
It can be considered as a
pseudo-time solver for VP rheology. In this case one writes

for stresses and

for velocity. Here α and β are some large constants. The subscript p is related to pseudotime
iterations, replacing subcycling of the standard EVP, and n is that of external time stepping.
Fields are initialized with values at time step n for p = 1, and their values for the last
iteration p = NEV P are taken as solutions for time step n+ 1. Clearly, if NEV P is sufficiently
large, the VP solutions have to be recovered. In order that CFL limitations be satisfied,
the product αβ should be sufficiently large (see Bouillon et al. (2013)). The regime of the
standard EVP scheme (NEV P = 120 and TEV P = ∆t/3) will be approximately recovered
for α = β = 40 and NEV P = 120, but much larger values have to be used on fine meshes to
warrant absence of noise.
One expects that if this scheme is stable and converged, it should produce solutions
identical to those of converging VP solver. In reality, it will not be run for full convergence,
which is still too expensive, and some difference will be preserved.
The sea-ice model described here implements the three approaches mentioned above,
which will be referred further as VP, EVP and mEVP. The reason of keeping all them is
two-fold. First, it facilitates the comparison of results with other models which may be
using one of these approaches. Second, their numerical efficiency and performance depend
on applications, and one may wish to select the most appropriate approach.

#Usage

#References

