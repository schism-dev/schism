The vertical velocity serves as a diagnostic variable for local volume conservation[^1], but is a physically important quantity, especially when a steep slope is present ([Zhang et al. 2004](#zhang2004)). To solve the vertical velocity, we apply a finite-volume method to a typical prism, as depicted in Figure [6](geometry-discretization.md#figure06), assuming that $w$ is constant within an element $i$, and obtain - 

[^1]: Although other definitions of volume/mass exist, we define volume/mass in the finite-volume sense throughout this paper and measure conservation based on this definition.

\begin{equation}
\label{eq01}
\hat{S}_{k+1} \left( \overline{u}_{k+1}^{n+1} n_{k+1}^{x} + \overline{v}_{k+1}^{n+1} n_{k+1}^{y} + w_{i, k+1}^{n+1} n_{k+1}^{z} \right) - \hat{S}_{k} \left( \overline{u}_{k}^{n+1} n_{k}^{x} + \overline{v}_{k}^{n+1} n_{k}^{y} + w_{i, k}^{n+1} n_{k}^{z} \right) + \sum_{m=1}^{3} \hat{P}_{js(i, m)} \left( \hat{q}_{js(i, m), k}^{n+1} + \hat{q}_{js(i, m), k+1}^{n+1}  \right)/2 = 0, (k=k^b, \cdots, N_z - 1)
\end{equation}

where $\hat{S}$ and $\hat{P}$ are the areas of the prism surfaces (Figure [6](geometry-discretization.md#figure06)), ($n^x, n^y, n^z$), are the normal vector (pointing upward), $\overline{u}$ and $\overline{v}$ the averaged horizontal velocities at the top and bottom surfaces, and $\hat{q}$ is the outward normal velocity at each side center. The vertical velocity is then solved from the bottom to the surface, in conjunction with the bottom boundary condition $(u, v, q)\cdot\pmb{n}=0$. In the case of earthquake module (`imm≠0`), the bed velocity is prescribed. A compact form for Eqs $\ref{eq01}$ is $\sum_{j\epsilon S^+} \left| Q_j\right| = \sum_{j\epsilon S^-} \left| Q_j\right|$, where $Q_j$ is the facial fluxes outward of a prism $i$. This conservation will be utilized in the transport equation as the foundation for mass conservation and constancy.

**References**
<span id="zhang2004">Zhang, Y., Baptista, A.M. and Myers, E.P. (2004) "A cross-scale model for 3D baroclinic circulation in estuary-plume-shelf systems: I. Formulation and skill assessment". Cont. Shelf Res., 24: 2187-2214.</span>