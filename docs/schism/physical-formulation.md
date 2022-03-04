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

## Turbulence closure

## Air-sea exchange

