After the elevations are found, SCHISM solves the momentum Eqs. [1](physical-formulation.md)  along each vertical column at side centers. A semi-implicit Galerkin finite-element method is used, with the barotropic pressure gradient and the vertical viscosity terms being treated implicitly, and other terms treated explicitly. For 3D cells, we have - 

\begin{equation}
\label{eq01}
\int_{\delta_b-h}^{\eta} \psi_l \left[ \pmb{u}^{n+1} - \Delta t \frac{\partial}{\partial z} \left( \nu \frac{\partial \pmb{u}^{n+1}}{\partial z} \right) \right] dz = \int_{\delta_b - h}^{\eta} \pmb{g} \psi_l dz , (l = kbs + 1, \cdots, N_z)
\end{equation}

where $\psi$ is the hat function in the vertical dimension, $\delta_b$ is the bottom cell thickness, and 

\begin{equation}
\label{eq02}
\pmb{g} = \pmb{u}^* + \Delta t \left[ \pmb{f} - g\theta\nabla\eta^{n+1} - g(1-\theta)\nabla\eta^n \right]
\end{equation}

The two terms that are treated implicitly would have imposed the most severe stability constraints if treated explicitly. The explicit treatment of the baroclinic pressure gradient and the horizontal viscosity terms, however, does impose mild stability constraints.

The final FEM equations are - 

\begin{equation}
\label{eq03}
\begin{aligned}
\frac{\Delta z_{l+1}}{6} \left( 2 \pmb{u}_{l}^{n+1} + \pmb{u}_{l+1}^{n+1} \right) + \frac{\Delta z_l}{6} \left( 2 \pmb{u}_l^{n+1} + \pmb{u}_{l-1}^{n+1} \right) - \nu_{l+1/2} \Delta t \frac{\pmb{u}_{l+1}^{n+1} - \pmb{u}_l^{n+1}}{\Delta z_{l+1}} + \nu_{l-1/2} \Delta t \frac{\pmb{u}_l^{n+1} - \pmb{u}_{l-1}^{n+1}}{\Delta z_l} &= \frac{\Delta z_{l+1}}{6} \left( 2 \pmb{g}_l + \pmb{g}_{l+1} \right) + \frac{\Delta z_l}{6} \left( 2\pmb{g}_l + \pmb{g}_{l-1} \right), (l = kbs + 2, \cdots, N_z - 1)\\
\frac{\Delta z_{l+1}}{6} \left( 2\pmb{u}_{l}^{n+1} + \pmb{u}_{l+1}^{n+1} \right) - \nu_{l+1/2} \Delta t \frac{\pmb{u}_{l+1}^{n+1} - \pmb{u}_{l}^{n+1}}{\Delta z_{l+1}} \chi\Delta t \pmb{u}_{kbs+1}^{n+1} &= \frac{\Delta z_{l+1}}{6} \left( 2\pmb{g}_l + \pmb{g}_{l+1} \right), (l = kbs + 1)\\
\frac{\Delta z_l}{6} \left( 2\pmb{u}_l^{n+1} + \pmb{u}_{l-1}^{n+1} \right) + \nu_{l-1/2} \Delta t \frac{\pmb{u}_l^{n+1} - \pmb{u}_{l-1}^{n+1}}{\Delta z_l} &= \pmb{\tau}_{w}^{n+1} \Delta t + \frac{\Delta z_l}{6} \left( 2\pmb{g}_l + \pmb{g}_{l-1} \right), (l=N_z)
\end{aligned}
\end{equation}

The bottom velocity is - 

\begin{equation}
\label{eq04}
\begin{aligned}
\pmb{u}_{kbs}^{n+1} &= 0, \text{ if } \chi \neq 0\\
\pmb{u}_{kbs}^{n+1} &= \pmb{u}_{kbs + 1}^{n+1}, \text{ if } \chi = 0
\end{aligned}
\end{equation}

which is consistent with the bottom BL formulation we used.

After the velocities at all sides are found, the velocity at a node, which is needed in ELM, is evaluated using scheme MA or MB as discussed above.

If a cell is 2D locally, the velocity is simply solved as - 

\begin{equation}
\label{eq05}
\pmb{u}^{n+1} = \frac{\breve{H}}{H} \left[ \pmb{u}^* + \left( \pmb{f} + \pmb{\tau}_w/H \right) \Delta t - g\theta\nabla\eta^{n+1} - g(1-\theta)\nabla\eta^n \right]
\end{equation}