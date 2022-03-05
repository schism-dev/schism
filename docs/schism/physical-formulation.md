## Governing equations
We will focus only on the hydrostatic solver in side SCHISM. Under this mode, we solve the standard Navier-Stokes equations with hydrostatic and Boussinesq approximations, including the effects of vegetation.

Momentum equation: 

\begin{equation}
\begin{aligned}
\frac{Du}{dt} = \mathbf{f} - g \nabla \eta + \mathbf{m}_z - \alpha \left| \mathbf{u} \right| \mathbf{u} L(x, y, z)\\
\mathbf{f} = f(v, -u) - \frac{g}{\rho_0} \int_z^{\eta} \nabla \rho d\zeta - \frac{\nabla p_A}{\rho_0} + \alpha a \nabla \Psi + \mathbf{F}_m + other
\end{aligned}
\label{eq:momentum}
\end{equation}

Continuity equation in 3D and 2D depth-integrated forms:

\begin{align}
\nabla \cdot \mathbf{u} + \frac{\partial w}{\partial z} = 0\\
\frac{\partial \eta}{\partial t} + \nabla \cdot \int_{-h}^{\eta} \mathbf{u} dz = 0
\end{align}

Transport equations:

\begin{equation}
\frac{\partial C}{\partial t} + \nabla \cdot (\mathbf{u}C) = \frac{\partial}{\partial z} \left( \kappa \frac{\partial C}{\partial z} \right) + F_h 
\end{equation}

Equation of state:

\begin{equation*}
\rho = \rho(S, T, p)
\end{equation*}

Where, 

- $\nabla$ : $\left( \frac{\partial}{\partial x}, \frac{\partial}{\partial y} \right)$
- $\frac{D}{Dt}$ : material derivative
- $(x, y)$: horizontal Cartesian coordinates
- $z$: vertical coordinate, positive upward
- $t$: time
- $\eta(x, y, t)$: free-surface elevation in meters [$m$]
- $h(x, y)$: bathymetric depth (measured from a fixed datum) [$m$]
- $\mathbf{u}(x, y, z, t)$: horizontal velocity, with Cartesian components $(u, v)$ [$m/s$]
- $w(x, y, z, t)$: vertical velocity [$m/s$]
- $p$: hydrostatic pressure [$Pa$]
- $p_A$: atmospheric pressure reduced to mean sea level (MSL) [$Pa$]
- $\rho, \rho_0$: water density and reference water density [$kg/m^3$]
- $\mathbf{f}$: other forcing terms in momentum (baroclinic gradient, horizontal viscosity, Coriolis, earth tidal potential, atmospheric pressure, radiation stress). These are usually treated explicitly in the numerical formulation
- $g$: acceleration of gravity, in [$m/s^2$]
- $C$: tracer concentration (e.g., salinity, temperature, sediment etc)
- $\nu$: vertical eddy viscosity, in [$m^2/s$]
- $\kappa$: vertical eddy diffusivity, for tracers, in [$m^2/s$]
- $\mathbf{F}_m$: horizontal viscosity [$m^2/s$]
- $F_h$: horizontal diffusion and mass sources/sinks [$m^2/s$]

Vegetation effects have been accounted for in Eq. 1. The main vegetation parameter is $\alpha(x, y) = D_v N_v C_{Dv}/2$ is a vegetation related density variable in [$m^{-1}]$, where $D_v$ is the stem diameter, $N_v$ is the vegetation density (number of stems per $m^2$), and $C_{Dv}$ is the bulk form drag coefficient. Selection of $C_{Dv}$ is the topic of other studies with values between 0 and 3 (Nepf and Vivoni 2000; Tanino and Nepf 2008), and is validated against reported lab study values. The underlying assumption used here is to treat the vegetation as arrays of solid cylinders, which is only a first-order approximation of the problem. Flexibility of the vegetation, sheltering effects within a cluster of vegetation can lead to one to two orders of reduction in the drag forces, and Gaylord et al. (2008) showed that the drag formulation is also species dependent. These additional complexities are outside the scope of the current study. In this paper, we assume $C_{Dv}$ is a constant, but a vertically varying $C_{Dv}$ (as suggested by Nepf and Vivoni 2000 and others) can be easily added as well; the latter can be used to approximate flexible stems (Nepf and Vivoni 2000; Luhar and Nepf 2011).

Since SCHISM allows ‘polymorphism’ with mixed 2D and 3D cells in a single grid (Zhang et al. 2016), we have different forms for the vertical eddy viscosity term $\mathbf{m}_z$ and vegetation term $L(x, y, z)$. 

\begin{equation}
\begin{aligned}
    \mathbf{m}_z= 
\begin{cases}
    \frac{\partial}{\partial z}\left( \nu \frac{\partial \mathbf{u}}{\partial z} \right),& \text{3D cells}\\
    \frac{\mathbf{\tau}_w - \chi \mathbf{u}}{H}, & \text{2D cells}
\end{cases}\\
    L(x, y, z)= 
\begin{cases}
    \mathcal{H}(z_v - z),& \text{3D}\\
    1,              & \text{2D}
\end{cases}
\end{aligned}
\end{equation}

where, $\nu$ is the eddy viscosity, $\mathbf{\tau}_w$ is the surface wind stress, $H=h+\eta$ is the total water depth (with $h$ being the depth measured from a fixed datum), $\chi = C_d \left| \mathbf{u} \right|$, $C_D$ is the bottom drag coefficient, $z_v$  is the z-coordinate of the canopy, and $\mathcal{H}()$ is the Heaviside step function - 

\begin{equation*}
\mathcal{H} = 
\begin{cases}
    1, & x \geq 0\\
    0, & x \lt 0
\end{cases}
\end{equation*}

Note that $\mathbf{u}$ denotes the depth-averaged velocity in a 2D region.

## Boundary conditions (B.C.)
The differential equations above need initial condition (I.C.) and B.C. In general, all state variables ($\eta$, $\mathbf{u}$, $C$) are specified at $t=0$ as I.C. and these are also specified at all open _lateral_ boundary segments (open ocean, rivers etc). However, not all variables need to be specified at all boundary segments and we’ll revisit this in the input-output section, i.e., [bctides](../input-output/bctides.md).

The vertical B.C. for (Eq 1-4) are described in detail below as these impact the numerical scheme. Note that these only apply to 3D cells; for 2D cells, Eq. 1a has taken into account the B.C.

At the sea surface, SCHISM enforces the balance between the internal Reynolds stress and the applied shear stress.

\begin{equation}
\nu \frac{\partial \mathbf{u}}{\partial z} = \mathbf{\tau}_w, \text{ at } z = \eta
\end{equation}

where the stress $\mathbf{\tau}_z$ can be parameterized using the approach of Zeng et al. (1998) or the simpler approach of Pond and Pickard (1998). If the Wind Wave Model is invoked, it can also be calculated from the wave model.

Because the bottom boundary layer is usually not well resolved in ocean models, the no-slip condition at the sea or river bottom ($\mathbf{u} = w = 0$) is replaced by a balance between the internal Reynolds stress and the bottom frictional stress. 

\begin{equation}
\nu \frac{\partial \mathbf{u}}{\partial z} = \mathbf{\tau}_b, \text{ at } z=-h
\end{equation}

The specific form of the bottom stress $\mathbf{\tau}_b$ depends on the type of boundary layer used and here we will only discuss the turbulent boundary layer below (Blumberg and Mellor 1987), given its prevalent usage in ocean modeling. The bottom stress is then - 

\begin{equation}
\mathbf{\tau}_b = C_D \left| \mathbf{u}_b \right| \mathbf{u}_b \equiv \chi \mathbf{u}_b
\end{equation}

The velocity profile in the interior of the bottom boundary layer obeys the logarithmic law, which is smoothly matched to the exterior flow at the top of the boundary layer.

\begin{equation}
\mathbf{u} = \frac{ln[(z+h)/z_0]}{ln(\delta_b/z_0)}\mathbf{u}_b, z_0-h \leq z \leq \delta_b -h
\end{equation}

Here, $\delta_b$ is the thickness of the bottom computational cell (assuming that the bottom is sufficiently resolved in SCHISM that the bottom cell is inside the boundary layer), $z_0$ is the bottom roughness, and $\mathbf{u}_b$ is the velocity measured at the top of the bottom computational cell. Therefore the Reynolds stress inside the boundary layer is derived as - 

\begin{equation}
\nu \frac{\partial \mathbf{u}}{\partial z} = \frac{\nu}{(z+h)ln(\delta_b/z_0)}  \mathbf{u}_b
\end{equation}

Utilizing the turbulence closure theory discussed below, we can show that the Reynolds stress is constant inside the boundary layer - 

\begin{equation}
\nu \frac{\partial \mathbf{u}}{\partial z} = \frac{\kappa_0}{ln(\delta_b/z_0)} C_D^{1/2} \left| \mathbf{u}_b \right| \mathbf{u}_b, z_0-h \leq z \leq \delta_b -h
\end{equation}

and the drag coefficient is calculated from Eq. 7, 8, and 11 as - 

\begin{equation}
C_D = \left( \frac{1}{\kappa_0} ln(\delta_b/z_0) \right) ^{-2}
\end{equation}

which is the drag formula as discussed in Blumberg and Mellor (1987). Eq. 11 also shows that the vertical viscosity term in the momentum equation vanishes inside the boundary layer. This fact will be utilized in the numerical formulation.

## Turbulence closure

## Air-sea exchange

