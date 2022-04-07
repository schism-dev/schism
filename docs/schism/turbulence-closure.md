[Umlauf and Burchard](#umlauf2003)’s Generic Length Scale (GLS) model is 

\begin{equation}
\label{eq01}
\begin{aligned}
\frac{DK}{Dt} &= \frac{\partial}{\partial z} ( \nu_k^\psi \frac{\partial K}{\partial z} ) + \nu M^2 + \mu N^2 - \varepsilon\\
\frac{D\psi}{Dt} &= \frac{\partial}{\partial z} ( \nu_\psi \frac{\partial \psi}{\partial z} ) + \frac{\psi}{K} ( c_{\psi 1} \nu M^2 + c_{\psi 3} \mu N^2 - c_{\psi 2} F_w \varepsilon)
\end{aligned}
\end{equation}

with natural B.C.

\begin{equation}
\label{eq02}
\begin{cases}
\nu_k^\psi \frac{\partial K}{\partial z} = 0, z = -h \text{ or } \eta\\
\nu_\psi \frac{\partial\psi}{\partial z} = \kappa_0 n \nu_\psi \frac{\psi}{l}, z = -h\\
\nu_\psi \frac{\partial\psi}{\partial z} = -\kappa_0 n \nu_\psi \frac{\psi}{l}, z = \eta
\end{cases}
\end{equation}

and essential B.C.:

\begin{equation}
\label{eq03}
\begin{cases}
K = ( c_\mu^0 )^{-2} \nu \left| \frac{\partial u}{\partial z} \right|\\
l = \kappa_0\Delta, \text{ at } z = -h \text{ or } \eta\\
\psi = ( c_\mu^0 )^p K^m ( \kappa_0 \Delta )^n
\end{cases}
\end{equation}

where $K$ is the TKE, $l$ is the mixing length, $c_{\psi *}$ are constants, $\psi = ( c_\mu^0 )^p K^m l^n$ is a generic length-scale variable, and $\Delta$ is the distance to the boundary (surface or bottom). The turbulence production and dissipation terms are:

\begin{equation}
\label{eq04}
M^2 = (\frac{\partial u}{\partial z})^2 + (\frac{\partial v}{\partial z})^2
\end{equation}

\begin{equation}
\label{eq05}
N^2 = \frac{g}{\rho_0}\frac{\partial \rho}{\partial z}
\end{equation}

\begin{equation}
\label{eq06}
\varepsilon = (c_\mu^0)^3 K^{1.5} l^{-1}
\end{equation}

In the code, the natural B.C. is applied first (see the FEM formulation below), and the essential B.C. is then used to overwrite the boundary values of the unknown, as suggested by GOTM. In the FEM formulation, $K$, $\psi$ are defined at nodes and whole levels. Furthermore, the sums of $M^2$ and $N^2$ are treated explicitly/implicitly depending on the sign. The final equations look similar to those for the momentum equation: 

\begin{equation}
\label{eq07}
\begin{aligned}
&\mathcal{H} (N_z -m) \left[ \frac{\Delta z_{m+1}}{6} (2K_m^{n+1} + K_{m+1}^{n+1} - 2K_m^n - K_{m+1}^n) - \Delta t (\nu_k^\psi)_{m+\frac{1}{2}}^n \frac{K_{m+1}^{n+1} - K_m^{n+1}}{\Delta z_{m+1}} \right]\\
&+ \mathcal{H}(m-kbp) \left[ \frac{\Delta z_m}{6} (2K_m^{n+1} + K_{m-1}^{n+1} - 2K_m^n - K_{m-1}^n) - \Delta t (\nu_k^\psi)_{m-\frac{1}{2}}^n \frac{K_m^{n+1} - K_{m-1}^{n+1}}{\Delta z_m} \right]\\
&= \mathcal{H} (N_z - m) \Delta t \left[ \begin{Bmatrix} \frac{\Delta z_m + 1}{2} (\nu_t M^2 + \nu_t^\theta N^2)_{m+\frac{1}{2}}^{n}\\ \frac{\Delta z_{m+1}}{2} (\nu_t M^2 + \nu_t^\theta N^2)_{m+\frac{1}{2}}^{n} \frac{2K_m^{n+1} + K_{m+1}^{n+1}}{K_{m+\frac{1}{2}}^{n}} \end{Bmatrix} - (c_\mu^0)^3 (K^{0.5} l^{-1})_{m+\frac{1}{2}}^{n} \frac{\Delta z_{m+1}}{6} (2K_m^{n+1} + K_{m+1}^{n+1}) \right]\\
&+ \mathcal{H}(m-kbp)\Delta t \left[ \begin{Bmatrix} \frac{\Delta z_m}{2} (\nu_t M^2 + \nu_t^\theta N^2)_{m-\frac{1}{2}}^{n}\\ \frac{\Delta z_m}{2} (\nu_t M^2 + \nu_t^\theta N^2)_{m-\frac{1}{2}}^{n} \frac{2K_m^{n+1} + K_{m-1}^{n+1}}{K_{m-\frac{1}{2}}^{n}}\end{Bmatrix} - (c_\mu^0)^3 (K^{0.5}l^{-1})_{m-\frac{1}{2}}^{n} \frac{\Delta z_m}{6} (2K_m^{n+1} + K_{m-1}^{n+1}) \right] , (l = kbp, \cdots, N_z)
\end{aligned}
\end{equation}

\begin{equation}
\label{eq08}
\begin{aligned}
&\mathcal{H}(N_z - m) \left[ \frac{\Delta z_{m+1}}{6} (2\psi^{n+1}_{m} + \psi_{m+1}^{n+1} - 2\psi_m^n - \psi_{m+1}^{n}) - \Delta t (\nu_\psi)_{m+\frac{1}{2}}^{n} \frac{\psi_{m+1}^{n+1} - \psi_{m}^{n+1}}{\Delta z_{m+1}} \right]\\
&+ \mathcal{H}(m-kbp) \left[ \frac{\Delta z_{m+1}}{6} (2\psi_m^{n+1} + \psi_{m-1}^{n+1} -2\psi_m^n - \psi_{m-1}^n) - \Delta t (\nu_\psi)_{m-\frac{1}{2}}^n \frac{\psi_m^{n+1} - \psi_{m-1}^{n+1}}{\Delta z_m} \right]\\
&+ \kappa_0 n \Delta t \left[ \delta_{m, N_z} \left(\frac{\nu_\psi}{l}\right)_{N_z}^n \psi_{N_z}^{n+1} + \delta_{m, kbp} \left( \frac{\nu_\psi}{l} \right)_{kbp}^n \psi_{kbp}^{n+1} \right]\\
&= \mathcal{H}(N_z - m) \Delta t \left[ \begin{Bmatrix} \frac{\Delta z_{m+1}}{2} (c_{\psi 1} \nu_t M^2 + c_{\psi 3} \nu_t^\theta N^2)_{m+\frac{1}{2}}^{n} \left(\frac{\psi}{K}\right)_{m+\frac{1}{2}}^{n} \\ \frac{\Delta z_{m+1}}{2} \left( \nu_t M^2 + \nu_t^\theta N^2 \right)_{m+\frac{1}{2}}^n \frac{2\psi_m^{n+1} + \psi_{m+1}^{n+1}}{2K_{m+\frac{1}{2}}^n} \end{Bmatrix} - c_{\psi 2} \left( c_\mu^0 \right)^3 \left(K^{0.5}l^{-1}F_w\right)_{m+\frac{1}{2}}^{n} \frac{\Delta z_{m+1}}{6} \left( 2\psi_m^{n+1} + \psi_{m+1}^{n+1} \right) \right] \\
&+\mathcal{H}(m-kbp) \Delta t \left[ \begin{Bmatrix} \frac{\Delta z_m}{2} \left( c_{\psi 1} \nu_t M^2 + c_{\psi 3} \nu_t^\theta N^2 \right)_{m-\frac{1}{2}}^{n} \left( \frac{\psi}{K} \right)_{m-\frac{1}{2}}^{n}\\ \frac{\Delta z_m}{2} \left( \nu_t M^2 + \nu_t^\theta N^2 \right)_{m-\frac{1}{2}}^{n} \frac{2\psi_{m}^{n+1} + \psi_{m-1}^{n+1}}{2K_{m-\frac{1}{2}}^{n}} \end{Bmatrix} - (c_\mu^0)^3 \left(K^{0.5}l^{-1}F_w\right)_{m-\frac{1}{2}}^{n} \frac{\Delta z_m}{6} \left( 2\psi_{m}^{n+1} + \psi_{m-1}^{n+1} \right) \right], (l = kbp, \cdots, N_z) 
\end{aligned}
\end{equation}

where $\begin{Bmatrix}\end{Bmatrix}$ indicates the alternative explicit/implicit schemes mentioned above, and $\mathcal{H}$ is a step function. We have applied the natural B.C. (Eqs. $\ref{eq02}$) in these equations, and after $K$ and $\psi$ are solved, the essential B.C. (Eqs $\ref{eq03}$) is then used to overwrite the boundary values.

**References**

<span id="umlauf2003">Umlauf, L. and H. Burchard (2003) A generic length-scale equation for geophysical turbulence models. J. Mar. Res., 6, pp. 235-265.</span>