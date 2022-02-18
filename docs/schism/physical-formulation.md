# Governing equations
We will focus only on the hydrostatic solver in side SCHISM. Under this mode, we solve the standard Navier-Stokes equations with hydrostatic and Boussinesq approximations, including the effects of vegetation.

Momentum equation:

\begin{equation}
\frac{Du}{dt} = f - g \nabla \eta + m_z - \alpha \left| u \right| u L(x, y, z)
\end{equation}

where, 
\begin{equation}
f = f(v, -u) - \frac{g}{\rho_0} \int_z^{\eta} \nabla \rho d\zeta - \frac{\nabla p_A}{\rho_0} + \alpha a \nabla \psi + F_m + other
\end{equation}

# Boundary conditions (B.C.)

# Turbulence closure

# Air-sea exchange
