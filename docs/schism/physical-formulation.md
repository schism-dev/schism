## Governing equations
We will focus only on the hydrostatic solver in side SCHISM. Under this mode, we solve the standard Navier-Stokes equations with hydrostatic and Boussinesq approximations, including the effects of vegetation.

Momentum equation: 

\begin{equation}
\begin{aligned}
\frac{Du}{dt} = \pmb{f} - g \nabla \eta + \pmb{m}_z - \alpha \left| \pmb{u} \right| \pmb{u} L(x, y, z)\\
\pmb{f} = f(v, -u) - \frac{g}{\rho_0} \int_z^{\eta} \nabla \rho d\zeta - \frac{\nabla p_A}{\rho_0} + \alpha a \nabla \Psi + \pmb{F}_m + other
\end{aligned}
\label{eq:momentum}
\end{equation}

Continuity equation in 3D and 2D depth-integrated forms:

\begin{align}
\nabla \cdot \pmb{u} + \frac{\partial w}{\partial z} = 0\\
\frac{\partial \eta}{\partial t} + \nabla \cdot \int_{-h}^{\eta} \pmb{u} dz = 0
\end{align}

Transport equations:

\begin{equation}
\frac{\partial C}{\partial t} + \nabla \cdot (\pmb{u}C) = \frac{\partial}{\partial z} \left( \kappa \frac{\partial C}{\partial z} \right) + F_h 
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
- $\pmb{u}(x, y, z, t)$: horizontal velocity, with Cartesian components $(u, v)$ [$m/s$]
- $w(x, y, z, t)$: vertical velocity [$m/s$]
- $p$: hydrostatic pressure [$Pa$]
- $p_A$: atmospheric pressure reduced to mean sea level (MSL) [$Pa$]
- $\rho, \rho_0$: water density and reference water density [$kg/m^3$]
- $\pmb{f}$: other forcing terms in momentum (baroclinic gradient, horizontal viscosity, Coriolis, earth tidal potential, atmospheric pressure, radiation stress). These are usually treated explicitly in the numerical formulation
- $g$: acceleration of gravity, in [$m/s^2$]
- $C$: tracer concentration (e.g., salinity, temperature, sediment etc)
- $\nu$: vertical eddy viscosity, in [$m^2/s$]
- $\kappa$: vertical eddy diffusivity, for tracers, in [$m^2/s$]
- $\pmb{F}_m$: horizontal viscosity [$m^2/s$]
- $F_h$: horizontal diffusion and mass sources/sinks [$m^2/s$]

Vegetation effects have been accounted for in Eq. 1. The main vegetation parameter is $\alpha(x, y) = D_v N_v C_{Dv}/2$ is a vegetation related density variable in [$m^{-1}]$, where $D_v$ is the stem diameter, $N_v$ is the vegetation density (number of stems per $m^2$), and $C_{Dv}$ is the bulk form drag coefficient. Selection of $C_{Dv}$ is the topic of other studies with values between 0 and 3 ([Nepf and Vivoni 2000](#); [Tanino and Nepf 2008](#)), and is validated against reported lab study values. The underlying assumption used here is to treat the vegetation as arrays of solid cylinders, which is only a first-order approximation of the problem. Flexibility of the vegetation, sheltering effects within a cluster of vegetation can lead to one to two orders of reduction in the drag forces, and [Gaylord et al. (2008)](#) showed that the drag formulation is also species dependent. These additional complexities are outside the scope of the current study. In this paper, we assume $C_{Dv}$ is a constant, but a vertically varying $C_{Dv}$ (as suggested by [Nepf and Vivoni 2000](#) and others) can be easily added as well; the latter can be used to approximate flexible stems ([Nepf and Vivoni 2000](#); Luhar and Nepf 2011).

Since SCHISM allows ‘polymorphism’ with mixed 2D and 3D cells in a single grid ([Zhang et al. 2016](#zhang2016)), we have different forms for the vertical eddy viscosity term $\pmb{m}_z$ and vegetation term $L(x, y, z)$. 

\begin{equation}
\begin{aligned}
    \pmb{m}_z= 
\begin{cases}
    \frac{\partial}{\partial z}\left( \nu \frac{\partial \pmb{u}}{\partial z} \right),& \text{3D cells}\\
    \frac{\pmb{\tau}_w - \chi \pmb{u}}{H}, & \text{2D cells}
\end{cases}\\
    L(x, y, z)= 
\begin{cases}
    \mathcal{H}(z_v - z),& \text{3D}\\
    1,              & \text{2D}
\end{cases}
\end{aligned}
\end{equation}

where, $\nu$ is the eddy viscosity, $\pmb{\tau}_w$ is the surface wind stress, $H=h+\eta$ is the total water depth (with $h$ being the depth measured from a fixed datum), $\chi = C_d \left| \pmb{u} \right|$, $C_D$ is the bottom drag coefficient, $z_v$  is the z-coordinate of the canopy, and $\mathcal{H}()$ is the Heaviside step function - 

\begin{equation*}
\mathcal{H} = 
\begin{cases}
    1, & x \geq 0\\
    0, & x \lt 0
\end{cases}
\end{equation*}

Note that $\pmb{u}$ denotes the depth-averaged velocity in a 2D region.

## Boundary conditions (B.C.)
The differential equations above need initial condition (I.C.) and B.C. In general, all state variables ($\eta$, $\pmb{u}$, $C$) are specified at $t=0$ as I.C. and these are also specified at all open _lateral_ boundary segments (open ocean, rivers etc). However, not all variables need to be specified at all boundary segments and we’ll revisit this in the input-output section, i.e., [bctides](../input-output/bctides.md).

The vertical B.C. for (Eq 1-4) are described in detail below as these impact the numerical scheme. Note that these only apply to 3D cells; for 2D cells, Eq. 1a has taken into account the B.C.

At the sea surface, SCHISM enforces the balance between the internal Reynolds stress and the applied shear stress.

\begin{equation}
\nu \frac{\partial \pmb{u}}{\partial z} = \pmb{\tau}_w, \text{ at } z = \eta
\end{equation}

where the stress $\pmb{\tau}_z$ can be parameterized using the approach of [Zeng et al. (1998)](#zeng1998) or the simpler approach of [Pond and Pickard (1998)](#pond1998). If the Wind Wave Model is invoked, it can also be calculated from the wave model.

Because the bottom boundary layer is usually not well resolved in ocean models, the no-slip condition at the sea or river bottom ($\pmb{u} = w = 0$) is replaced by a balance between the internal Reynolds stress and the bottom frictional stress. 

\begin{equation}
\nu \frac{\partial \pmb{u}}{\partial z} = \pmb{\tau}_b, \text{ at } z=-h
\end{equation}

The specific form of the bottom stress $\pmb{\tau}_b$ depends on the type of boundary layer used and here we will only discuss the turbulent boundary layer below (Blumberg and Mellor 1987), given its prevalent usage in ocean modeling. The bottom stress is then - 

\begin{equation}
\pmb{\tau}_b = C_D \left| \pmb{u}_b \right| \pmb{u}_b \equiv \chi \pmb{u}_b
\end{equation}

The velocity profile in the interior of the bottom boundary layer obeys the logarithmic law, which is smoothly matched to the exterior flow at the top of the boundary layer.

\begin{equation}
\pmb{u} = \frac{ln[(z+h)/z_0]}{ln(\delta_b/z_0)}\pmb{u}_b, z_0-h \leq z \leq \delta_b -h
\end{equation}

Here, $\delta_b$ is the thickness of the bottom computational cell (assuming that the bottom is sufficiently resolved in SCHISM that the bottom cell is inside the boundary layer), $z_0$ is the bottom roughness, and $\pmb{u}_b$ is the velocity measured at the top of the bottom computational cell. Therefore the Reynolds stress inside the boundary layer is derived as - 

\begin{equation}
\nu \frac{\partial \pmb{u}}{\partial z} = \frac{\nu}{(z+h)ln(\delta_b/z_0)}  \pmb{u}_b
\end{equation}

Utilizing the turbulence closure theory discussed below, we can show that the Reynolds stress is constant inside the boundary layer - 

\begin{equation}
\nu \frac{\partial \pmb{u}}{\partial z} = \frac{\kappa_0}{ln(\delta_b/z_0)} C_D^{1/2} \left| \pmb{u}_b \right| \pmb{u}_b, z_0-h \leq z \leq \delta_b -h
\end{equation}

and the drag coefficient is calculated from Eq. 7, 8, and 11 as - 

\begin{equation}
C_D = \left( \frac{1}{\kappa_0} ln(\delta_b/z_0) \right) ^{-2}
\end{equation}

which is the drag formula as discussed in [Blumberg and Mellor (1987)](#blumberg1987). Eq. 11 also shows that the vertical viscosity term in the momentum equation vanishes inside the boundary layer. This fact will be utilized in the numerical formulation.

## Turbulence closure
Eq 1-4 are not closed and must be supplemented by turbulence closure equations for the viscosity/diffusivity. We use the Generic Length-scale (GLS) model of [Umlauf and Burchard (2003)](#umlauf2003), which has the advantage of encompassing most of the Eq 6 closure models $k-\varepsilon$ ([Rodi 1984](#rodi1984)); $k-\omega$ ([Wilcox 1998](#wilcox1998); [Mellor and Yamada, 1982](#mellor1982)). In this framework, the transport, production, and dissipation of the turbulent kinetic energy ($K$) and of a generic length-scale variable ($\psi$) are governed by - 

\begin{equation}
\frac{Dk}{Dt} = \frac{\partial}{\partial z} \left( \nu_k^{\psi} \frac{\partial k}{\partial z} \right) + \nu M^2 + \kappa N^2 - \varepsilon + c_{fk} \alpha \left| \pmb{u} \right|^3 \mathcal{H} (z_v - z)
\end{equation}

\begin{equation}
\frac{D \psi}{Dt} = \frac{\partial}{\partial z} \left( \nu_{\psi} \frac{\partial \psi}{\partial z} \right) + \frac{\psi}{k}\left[ c_{\psi 1} \nu M^2 + c_{\psi 3} \kappa N^2 - c_{\psi 2}\varepsilon F_{wall} + c_{f\psi} \alpha \left| \pmb{u} \right|^3 \mathcal{H}(z_v-z) \right]
\end{equation}

where $\nu_k^{\psi}$ and $\nu_{\psi}$ are vertical turbulent diffusivities, $c_{\psi 1}$, $c_{\psi 2}$, and $c_{\psi 3}$ are model-specific constants ([Umlauf and Burchard 2003](#umlauf2003)), $F_{wall}$ is a wall proximity function, $M$ and $N$ are shear and buoyancy frequencies, and $\varepsilon$ is a dissipation rate. The generic length-scale is defined as - 

\begin{equation}
\psi = (c_{\mu}^0)^p K^m \ell^n
\end{equation}

where $c_{\mu}^0 = 0.3^{1/2}$ and $\ell$ is the turbulence mixing length. The specific choices of the constants $p$, $m$, and $n$ lead to the different closure models mentioned above. Finally, vertical viscosities and diffusivities as appeared in Eq 1-4 are related to $K$, $ell$, and stability functions - 

\begin{equation}
\begin{aligned}
\nu = \sqrt{2} s_m K^{1/2} \ell \\
\kappa = \sqrt{2} s_h K^{1/2} \ell \\
\nu_k^{\psi} = \frac{\nu}{\sigma_k^{\psi}} \\
\nu_{\psi} = \frac{\nu}{\sigma_{\psi}}
\end{aligned}
\end{equation}

where the Schmidt numbers $\sigma_k^{\psi}$ and $\sigma_{\psi}$ are model-specific constants. The stability functions ($s_m$ and $s_h$) are given by an Algebraic Stress Model (e.g.: [Kantha and Clayson 1994](#kantha1994), [Canuto et al. 2001](#canuto2001), or [Galperin et al. 1988](#galperin1988)). Following [Shimizu and Tsujimoto (1994; ST94 hereafter)](#), we set $c_{fk} = 0.07$ and $c_{f\psi} = 0.16$.

At the free surface and at the bottom of rivers and oceans, the turbulent kinetic energy and the mixing length are specified as Direchlet boundary conditions - 

\begin{equation}
K = \frac{1}{2} B_1^{2/3} \left| \pmb{\tau_b} \right|, \text{ or } \frac{1}{2} B_1^{2/3} \left| \pmb{\tau_w} \right| 
\end{equation}

\begin{equation}
\ell = \kappa_o d_b \text{ or } \kappa_0 d_s
\end{equation}

where $\pmb{\tau_b}$ is a bottom frictional stress, $\kappa_0 = 0.4$ is the von Karman’s constant, $B_1$ is a constant, and $d_b$ and $d_s$ are the distances to the bottom and the free surface, respectively.

## Air-sea exchange
We use the bulk aerodynamic module of [Zeng et al. (1998)](#zeng1998), which can be viewed [here](http://ccrm.vims.edu/yinglong/SVN_large_files/Zeng_etal_JClimate_1998-BulkAerodynamic-Model.pdf).


**References**

<span id="blumberg1987">Blumberg, A.F. and G.L. Mellor (1987) A description of a three-dimensional coastal ocean circulation model. In: Three-Dimensional Coastal Ocean Models, vol. 4, Coastal and Estuarine Studies, N. Heaps, editor, Washington, D.C.: AGU, pp. 1-16.</span>
<span id="canuto2001">Canuto, V.M., A. Howard, Y. Cheng and M.S. Dubovikov (2001) Ocean turbulence I: one-point closure model. Momentum and heat vertical diffusivities. J. Phys. Oceano., 31, pp. 1413-1426.</span>
<span id="galperin1988">Galperin, B., L. H. Kantha, S. Hassid and A. Rosati (1988) A quasi-equilibrium turbulent energy model for geophysical flows. J. Atmos. Sci., 45, pp. 55-62.</span>
<span id="kantha1994">Kantha, L.H. and C.A. Clayson (1994) An improved mixed layer model for geophysical applications. J. Geophy. Res, 99(25), pp. 235-266.</span>
<span id="mellor1982">Mellor, G.L. and T. Yamada (1982) Development of a turbulence closure model for geophysical fluid problems. Rev. Geophys., 20, pp. 851-875.</span>
<span id="pond1998">Pond, S. and G.L. Pickard (1998) Introductory Dynamical Oceanography, Butterworth-Heinmann.</span>
<span id="rodi1984">Rodi, W. (1984) Turbulence models and their applications in hydraulics: a state of the art review. Delft, The Netherlands, International Association for Hydraulics Research.</span>
<span id="umlauf2003">Umlauf, L. and H. Burchard (2003) A generic length-scale equation for geophysical turbulence models. J. Mar. Res., 6, pp. 235-265.</span>
<span id="wilcox1998">Wilcox, D.C. (1998) Reassessment of scale determining equation for advance turbulence models. AIAA J., 26, pp. 1299-1310.</span>
<span id="zeng1998">Zeng, X., M. Zhao and R.E. Dickinson (1998) Intercomparison of bulk aerodynamic algorithms for the computation of sea surface fluxes using TOGA COARE and TAO data. J. Clim., 11, pp. 2628-2644.</span>
<span id="zhang2016">Zhang, Y., Ye, F., Stanev, E.V., Grashorn, S. (2016). Seamless cross-scale modeling with SCHISM, Ocean Modelling, 102, 64-81. doi:10.1016/j.ocemod.2016.05.002</span>
