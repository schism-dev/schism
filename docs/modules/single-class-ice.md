#Theory

The single-class ice/snow module is taken from FESOM2, courtesy of Dr. Sergey Danilov's group. The info below is from the FESOM2 manual with minor modifications 
 to show some new additions.

Dynamical part solves the ice motion equation:

\begin{equation}
\label{ice01}
 m(\partial_t + \pmb{f} \times) \pmb{u}=a\tau - aC_d\rho_o (\pmb{u}-\pmb{u}_o)|\pmb{u}-\pmb{u}_o| +\pmb{F} -mg\nabla H
\end{equation}

Here $m$ is the ice plus snow mass per unit area, $C_d$ the ice-ocean drag coefficient, $\rho_o$ the
water density, $a$ the ice compactness, $\pmb{u}$ and $\pmb{u}_o$ the ice and ocean velocities, $\tau$ the wind
stress, $H$ the sea surface elevation, $g$ the acceleration due to gravity and $F_j = \partial_i \sigma_{ij}$ is
the force from stresses within the ice. Here we use Cartesian coordinates for brevity, with
$i, j = 1, 2$ implying $x$ and $y$ directions, and the implementation of spherical coordinates will
be discussed further. Summation over repeating coordinate indices is implied. The stress
tensor is symmetric. The mass $m$ is the combination of ice and snow contributions

\begin{equation}
 m=\rho_{ice}h_{ice}+\rho_sh_s
\end{equation}

with $\rho_{ice}$ and $\rho_s$, respectively, the densities of ice and snow and $h_{ice}$ and $h_s$ their mean
thicknesses (volumes per unit area).

For the visco-plastic (VP) rheology (Hibler, 1979) one writes

\begin{equation}
\label{ice02}
\sigma_{ij}=2\eta (\dot{\epsilon}_{ij} -0.5 \delta_{ij} \dot{\epsilon}_{kk}) + \zeta \delta_{ij} \dot{\epsilon}_{kk} -0.5\delta_{ij} P
\end{equation}

where

\begin{equation}
 \dot{\epsilon}_{ij}=0.5(\partial u_i / \partial x_j + \partial u_j / \partial x_i)
\end{equation}


is the deformation rate tensor, $\eta$ and $\zeta$ are the moduli (‘viscosities’) and $P$ the pressure.
Both the stress and deformation rate tensors are symmetric, so they are characterized by
three independent components. The standard VP rheology adopts the following scheme of
computing pressure $P$ and moduli $\eta$ and $\zeta$: 

\begin{equation}
 P_0=h_{ice}p^* exp[-C(1-a)], P=P_0\frac{\Delta}(\Delta+\Delta_{min}}, \zeta=\frac{0.5P_0}{\Delta+\Delta_{min}}, \eta=\frac{\zeta}{e^2}
\end{equation}


$e=2$  (the ellipticity parameter), $C=20$, $\Delta_{min}=2.e-9 s^{-1}$, and $p^*=15000 Pa$ 
are default values, they are adjusted in practice. In this scheme, 
$\Delta_{min}$ serves for viscous regularization of plastic behavior in areas where $\Delat$ is very small. We note that recent multi-category
ice implementations (such as CICE, see Hunke and Lipscomb 2008) use different parameterization for $P_0$, which takes into account the distribution of ice over thickness categories. This
does not change the basic equations $\ref{ice01}$, $\ref{ice02}$.

In our case we deal with three tracers, the area coverage $a$, ice volume $h_{ice}$ and snow
volume $h_s$. They are advected with ice velocities and modified through thermodynamical
forcing

\begin{equation}
\label{ice03}
\partial_t a+\nabla(a\pmb{u})=S_a, \partial_t h_{ice}+\nabla(h_{ice}\pmb{u})=S_{ice}, \partial_t H-s+\nabla(h_s\pmb{u})=S_s
\end{equation}


with $S_a$, $S_{ice}$ and $S_s$ the sources due to the exchange with atmosphere and ocean.

The system $\ref{ice01}, \ref{ice02}$, \ref{ice03}$, augmented with appropriate model of sources and boundary conditions, defines the sea ice model. We use the boundary conditions of no slip for
momentum and no flux for tracers at lateral walls. The well known difficulty in solving
these equation is related to the rheology term in the momentum equation, which makes this
equation very stiff and would require time steps of fractions of second if stepped explicitly.
There are two ways of handling this difficulty. The first one treats the rheology part in an
implicit way, with linearization for the moduli, as suggested by Zhang and Hibler (1997). As
mentioned elsewhere (see, e. g., Lemieux and Tremblay, 2009), it does not warrant full convergence, which requires a full nonlinear solver (for example, a Jacobian-free Newton-Krylov
one, see Lemieux et al. 2012). The latter is still too expensive computationally, so the VP
option used by us is similar in spirit to that of Zhang and Hibler (see section 2.4). The second
way is to reformulate the VP approach by adding pseudo-elasticity. It raises the time order
of the system $\ref{ice01}, \ref{ice02}$ to the second, which makes the CFL limitation on the explicit time step
much less severe than in the original VP framework.

Below is the formulation for a modified elastic-visco-plastic (mEVP) due to Bouillon et al. (2013).
It can be considered as a pseudo-time solver for VP rheology. In this case one writes

\begin{equation}
\begin{aligned}
 \alpha (\sigma_1^{p+1} - \sigma_1^{p})=\frac{P_0}{\Delta^p+\Delta_{min}} (\dot{\epsilon}_1^p -\Delta^p)-\sigma_1^{p} \\
 \alpha (\sigma_2^{p+1} - \sigma_2^{p})=\frac{P_0}{(\Delta^p+\Delta_{min})e^2} \dot{\epsilon}_2p -\sigma_2^{p} \\
 \alpha (\sigma_{12}^{p+1} - \sigma_{12}^{p})=\frac{P_0}{(\Delta^p+\Delta_{min})e^2} \dot{\epsilon}_{12}^p -\sigma_{12}^{p} 
\end{aligned}
\end{equation}


for stresses and

\begin{equation}
\beta(\pmb{u}^{p+1}-\pmb{u}^p)=-\pmb{u}^{p+1}+\pmb{u}^n-\Delta t \pmb{f}\times \pmb{u}^{p+1}+\frac{\Delta t}{m}[\pmb{F}^{p+1}+a\tau+C_da\rho_o(\pmb{u}_o^n-\pmb{u}^{p+1})|\pmb{u}_o^n-\pmb{u}^p| -mg\nabla H^n]
\end{equation}

for velocity. Here $\alpha$ and $\beta$ are some large constants. The subscript $p$ is related to pseudotime
iterations, replacing subcycling of the standard EVP, and $n$ is the previous time step.
Fields are initialized with values at time step $n$ for $p = 1$, and their values for the last
iteration $p = N_{EVP}$ are taken as solutions for time step $n+ 1$. Clearly, if $N_{EVP}$ is sufficiently
large, the VP solutions have to be recovered. In order that CFL limitations be satisfied,
the product $\alpha\beta$ should be sufficiently large (see Bouillon et al. (2013)). The regime of the
standard EVP scheme ($N_{EVP} = 120$ and $T_{EVP} = \Delta t/3$) will be approximately recovered
for $\alpha =\beta  = 40$ and $N_{EVP} = 120$, but much larger values have to be used on fine meshes to
warrant absence of noise.

One expects that if this scheme is stable and converged, it should produce solutions
identical to those of converging VP solver. In reality, it will not be run for full convergence,
which is still too expensive, and some difference will be preserved.
The sea-ice model described here implements the three approaches mentioned above,
which will be referred further as VP, EVP and mEVP. The reason of keeping all them is
two-fold. First, it facilitates the comparison of results with other models which may be
using one of these approaches. Second, their numerical efficiency and performance depend
on applications, and one may wish to select the most appropriate approach.

The tracer transport equation $\ref{ice03}$ is solved using a FCT scheme. 

#Usage

#References

