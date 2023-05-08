SCHISM supports a few FV solvers for the transport equation. All of the tracers, including T,S, sediment (if invoked) etc are solved simultaneously for efficiency. 

The transport equation for a generic tracer C is given by:

\begin{equation}
\label{eq01}
\frac{\partial C}{\partial t} + \nabla \cdot (\pmb{u}C) = \frac{\partial}{\partial z} \left( \kappa \frac{\partial C}{\partial z} \right) + F_h
\end{equation}

where $F_h$ includes vertical settling term (see [Vertical movement](#vertical-movement)), source/sink and also horizontal viscosity terms. The vertical B.C. is:

\begin{equation}
\label{eq02}
\begin{aligned}
\kappa\frac{\partial C}{\partial z} &= \hat{C}, z = \eta\\
\kappa\frac{\partial C}{\partial z} &= \hat{C}_b, z = -h
\end{aligned}
\end{equation}

Note that the 3D continuity equation ensures the constancy condition for the transport equation, i.e. $C=\text{const}$ initially will remain so in the absence of sinks/source.

In the following, we describe the numerical algorithm starting from 
the simplest 1st-order upwind scheme to the more complex 3rd-order WENO. In the newer versions, the pure upwind and explicit 
 TVD schemes have been deprecated.

## Upwind
Since most of the variables below are defined at prism center, we will use shorthand like $i$ etc to denote a prism at level $k + 1/2$ when there is no confusion. Also we often omit the superscript $n$ in the explicit terms for brevity.

A FV discretization of Eq. $\ref{eq01}$ for prism $i$ is :

\begin{equation}
\label{eq03}
C_i^{m+1} = C_i^m - \frac{\Delta t'}{V_i}\sum_{j\in S} Q_j C_{j*} + (F_h)_i^n \Delta t' + \frac{A_i \Delta t'}{V_i} \left[ \kappa_{i, k} \frac{C_{i, k+1}^{m+1} - C_{i,k}^{m+1}}{\Delta z_{i, k+\frac{1}{2}}} - \kappa_{i, k-1} \frac{C_{i,k}^{m+1} - C_{i,k-1}^{m+1}}{\Delta z_{i, k-\frac{1}{2}}} \right], (k = kbe + 1, \cdots, N_z)
\end{equation}

where $\Delta t' \neq \Delta t$  is the transport time step (subject to Courant condition below), $V_i$ is the volume of the prism, $C_i$ s a shorthand for $C_{i, k}$ (i.e., concentration at prism $(i,k)$), and $Q_j$ is the flux at face $j$ outward of the prism. Note that we have treated the diffusion term implicitly. For the sake of brevity we’ll drop the source and diffusion terms from now on and focus on the advection term. With the upwind scheme, the face concentration is defined as:

\begin{equation}
\label{eq04}
C_{j*} = \begin{cases}
C_i, j \in S^+\\
C_j, j \in S^-
\end{cases} \equiv C_{up}
\end{equation}

where we have used shorthand for concentration at prism $j$ (i.e. the prism adjacent to $(i,k)$ from face j), and $S^+$ and $S^-$ are outflow and inflow faces respectively. The face concentration take different forms with higher-order schemes.

!!!note "Mass conservation"
    Eq. $\ref{eq03}$ is the starting point of all FV solvers in SCHISM, from which a conservation statement can be derived. Assuming no zero fluxes at surface, bottom and lateral boundary, summing up over all prisms leads to:

    \begin{equation*}
    \sum_i V_i C_i^{m+1} = \sum_i V_i C_i^m - \Delta t' \sum_{j\in FS} Q_j C_j
    \end{equation*}

    where $FS$ stands for free surface. Note that the volume $V_i$ is evaluated at previous step $m$. The 2nd term in Eq. $\ref{eq03}$ cancels out at all faces (or vanish at lateral boundary) except at the free surface. The 2nd term in above equation represents the contribution from the surface movement and is supposed to account for the movement from $m$ to $m+1$. However, this balance is not precise (time truncation error). In the case of $TVD^2$ or $WENO$ (`itr_met>2`), there is also additional splitting error. Therefore mass conservation is only good up to time truncation error.

Retaining only the advection term, Eq. $\ref{eq03}$ then becomes:

\begin{equation}
\label{eq05}
C_i^{m+1} = C_i \left( 1 - \frac{\Delta t'}{V_i} \sum_{j\in S^+} \left| Q_j \right| \right) + \frac{\Delta t'}{V_i}\sum_{j\in S^+} \left| Q_j \right| C_j
\end{equation}

We have utilized the volume conservation:

\begin{equation}
\label{eq06}
\sum_{j\in S^+} \left| Q_j \right| = \sum_{j\in S^-} \left| Q_j \right|
\end{equation}

Therefore the Courant condition is:

\begin{equation}
\label{eq07}
1 - \frac{\Delta t'}{V_i} \sum_{j\in S^+} \left| Q_j \right| \geq 0
\end{equation}

SCHISM uses this eq. to dynamically adjust the time step for transport for each step. Moreover, to improve efficiency, the vertical flux terms in Eq. $\ref{eq05}$ are treated implicitly, and the corresponding terms are then removed in the Courant condition - Eq. $\ref{eq07}$. This is allowable because upwind is a linear method.

## TVD (explicit)
The only difference between TVD and upwind schemes lies in the evaluation of the interfacial concentration:

\begin{equation}
\label{eq08}
C_{j*} = C_{up} + \frac{\varphi_j}{2}\left( C_{jD} - C_{up} \right)
\end{equation}

where $C_{up}$ is given in Eq. $\ref{eq04}$, $C_{jD}$ is the downstream concentration, and $0 \leq \varphi_j \leq 2$ is a limiter function. TVD scheme nominally approaches 2nd order accuracy due to the anti-diffusion term. 

After some algebraic manipulation, the final eq. for TVD is:

\begin{equation}
\label{eq09}
C_i^{m+1} = C_i + \frac{\Delta t'}{V_i} \sum_{j\in S^-} \left| Q_j \right| (C_j - C_i) + \frac{\Delta t'}{V_i} \sum_{j\in S} \left| Q_j \right| \frac{\varphi_j}{2} (C_i - C_j) + \text{source} + \text{diffusion}
\end{equation}

and the Courant condition is:

\begin{equation}
\label{eq10}
\Delta t' \leq \frac{V_i}{\sum_{j\in S^-}\left| Q_j \right| \left( 1-\frac{\varphi_j}{2} + \delta_i\right)}
\end{equation}

where:

\begin{equation}
\label{eq11}
\delta_i = \sum_{p\in S^+} \frac{\varphi(r_p)}{2r_p}
\end{equation}

and the upwind ratio, which involves upwind of upwind neighboring prism, is given by:

\begin{equation}
\label{eq12}
r_p = \frac{\sum_{m\in S^-} \left| Q_m \right| (C_m - C_i)}{\left| Q_p \right| (C_i - C_p)}, p \in S^+
\end{equation}

In Eqs. $\ref{eq09}$ and $\ref{eq10}$, the faces $S$, $S^+$, and $S^-$ need to exclude the locations where upwind is applied: all horizontal and vertical boundaries, and interfaces between wetting and drying. In those places, $\varphi_j = \delta_i = 0$. Again SCHISM automatically calculates the time step according to the Courant condition (Eq. $\ref{eq10}$); the sub-time step used is the minimum of all prisms. The choices for the limiter function include: *MINMOD*, *OSHER*, *van Leer*, *Super Bee* etc.

Since TVD scheme here is a nonlinear method, we cannot treat the vertical fluxes implicitly, and so all fluxes have 
to be treated explicitly. TVD method is therefore more expensive than upwind. A more efficient [$TVD^2$](#tvd2) using a fractional
 time step method should be used. 

## $TVD^2$
The TVD scheme shown above is explicit in 3D space and thus subject to the Courant condition, which comprises of horizontal and vertical fluxes across each of the prism faces ([Casulli and Zanolli 2005](#casulli2005)). The restriction related to the vertical fluxes is especially severe due to smaller grid size used in the vertical dimension, and therefore a large number of sub-cycles within each time step are usually required. To partially mitigate the issue, a hybrid upwind-TVD approach can be used in which the more efficient upwind scheme, with an implicit treatment of the vertical fluxes, is used when the flow depth falls below a given threshold (with the assumption that stratification is usually much smaller in the shallows). However, this approach does not work in deeper depths of eddying regime, as large vertical velocities are not uncommon along steep bathymetric slopes. Together with the fact that a large number of vertical levels are usually required in the eddying regime, the explicit scheme leads to subpar computational performance and usually takes over 90% of the total CPU time.

We therefore develop an implicit TVD scheme in the vertical dimension in SCHISM. We start from the FVM formulation of the 3D transport equation at a prism $i$:

\begin{equation}
\label{eq13}
C_i^{n+1} = C_i^n - \frac{\Delta t}{V_i} \sum_{j\in S^-} \left| Q_j \right| (C_i - C_j) - \frac{\Delta t}{V_i} \sum_{j\in S} Q_j C_{jr} + \frac{A_i \Delta t}{V_i} \left[ \left(\kappa \frac{\partial C}{\partial z} \right)_{i,k} - \left(\kappa \frac{\partial C}{\partial z} \right)_{i, k-1} \right] + \frac{\Delta t}{V_i} \int_{V_i} F_h dV
\end{equation}

where $C_j$ is the concentration at the neighboring prism of $i$ across a prism face $j\in S = S^+ \cup S^-$, with $S^+/S^-$ denoting outflow/inflow faces (which can be horizontal or vertical) respectively, $V_i$ is the prism volume, $A_i$ is the area of the associated surficial element, and $Q_j$ is the flux at a face. In Eq. $\ref{eq13}$ we have utilized the volume conservation in a prism (which is enforced by the solution of the vertical velocity): $\sum_{j\in S^-}\left| Q_j \right| = \sum_{j\in S^+} \left| Q_j \right|$. We have also approximated the concentration at a face as the sum of an upwind and a correction part as:

\begin{equation}
\label{eq14}
C\Biggr|_j = C_{jup} + C_{jr}
\end{equation}

Note that in the 2nd term of RHS of Eq. $\ref{eq13}$, we have $C_j = C_{jup}$ as $j$ is an inflow face. In addition, we have intentionally left out the time level in some terms in Eq. $\ref{eq13}$ as they will be treated explicitly or implicitly in the following.

We split the solution of Eq. $\ref{eq13}$ into 3 sub-steps:

\begin{equation}
\label{eq15}
C_i^{m+1} = C_i^n + \frac{\Delta t_m}{V_i} \sum_{j\in S_H^-} \left| Q_j \right| (C_j^m - C_i^m) - \frac{\Delta t_m}{V_i} \sum_{j\in S_H} Q_j \hat\psi_j^m, (m = 1, \cdots, M)
\end{equation}

\begin{equation}
\label{eq16}
\widetilde C_i = C_i^{M+1} + \frac{\Delta t}{V_i} \sum_{j\in S_V^-} \left| Q_j \right| (\widetilde C_j - \widetilde C_i) - \frac{\Delta t}{V_i} \sum_{j\in S_V} Q_j (\Phi_j + \Psi_j), (j = kbe+1, \cdots, N_z)
\end{equation}

\begin{equation}
\label{eq17}
C_i^{n+1} = \widetilde C_i + \frac{A_i \Delta t}{V_i} \left[ \left( \kappa \frac{\partial C}{\partial z} \right)_{i,k}^{n+1} - \left( \kappa \frac{\partial C}{\partial z} \right)_{i, k-1}^{n+1} \right] + \frac{\Delta t}{V_i} \int_{V_i} F_h^n dV, (k = kbe+1, \cdots, N_z)
\end{equation}

The 1st step solves the horizontal advection part (for all 3D prisms $i$), the 2nd step deals with the vertical advection part (where $k_b$ is the bottom level index and $N_z$ is the surface level index), and the last step tackles the remaining terms. We could have combined the 1st and 3rd steps into a single step at the expense of efficiency, because sub-cycling is used in the 1st step. In Eq. $\ref{eq15}$, sub-cylcing in $M$ sub-steps is required because of the horizontal Courant number condition, $\Delta t_m$ is the sub-time step used, and $\hat\psi_j^m$ is a standard TVD limiter function. Eq. $\ref{eq15}$ is then solved with a standard TVD method. The last step (Eq. $\ref{eq17}$) requires the solution of a simple tri-diagonal matrix. So we will only focus on the 2nd step.

Following [Duraisamy and Baeder (2007, hereafter DB07)](#duraisamy2007), we use two limiter functions in Eq. $\ref{eq16}$: $\Phi_j$ is the space limiter and $\Psi_j$ is the time limiter - thus the name $TVD^2$. The origin of these two limiters is the approximation Eq. $\ref{eq14}$ via a Taylor expansion in both space and time ([DB07](#duraisamy2007)):

\begin{equation}
\label{eq18}
\begin{aligned}
C_j^{n+\frac{1}{2}} &= C_{jup}^{n+1} + \Phi_j + \Psi_j\\
&= C_{jup}^{n+1} + \pmb{r}\cdot\Bigr[ \nabla C \Bigr]_{jup}^{n+1} - \frac{\Delta t}{2} \Bigr[ \frac{\partial C}{\partial t} \Bigr]_{jup}^{n+1}
\end{aligned}
\end{equation}

Note that the interface value is taken at time level $n+\frac{1}{2}$ to gain 2nd-order accuracy in time. The vector $\pmb{r}$ points from prism center $jup$ to face center $j$. Due to the operator splitting method, $C^{n+1}$ now actually corresponds to $\widetilde{C}$. Customary in a TVD method, we then replace the last 2 terms with limiter functions:

\begin{equation}
\label{eq19}
C_j^{n+\frac{1}{2}} = \widetilde C_{jup} + \frac{\phi_j}{2} (\widetilde C_{jD} - \widetilde C_{jup}) - \frac{\psi_j}{2} (\widetilde C_{jup} - C_{jup}^{M+1})
\end{equation}

and so:

\begin{equation}
\label{eq20}
\begin{aligned}
\Phi_j = \frac{\phi_j}{2} (\widetilde C_{jD} - \widetilde C_{jup})\\
\Psi_j = \frac{\psi_j}{2} (\widetilde C_{jup} - C_{jup}^{M+1})
\end{aligned}
\end{equation}

where ‘jD’ stands for the downwind prism of $i$ along the face $j$, and $\phi_j$ and $\psi_j$ are 2 limiter functions in space and time respectively. Note that $\phi_j = \psi_j = 1$ leads to 2nd-order accuracy in both space and time.

Substituting Eq. $\ref{eq20}$ in to Eq. $\ref{eq16}$ and after some algebra we obtain a nonlinear equation for the unknown concentration:

\begin{equation}
\label{eq21}
\widetilde C_i + \frac{\frac{\Delta t}{V_i} \sum_{j\in S_V^-} \left| Q_j \right| \left[ 1 + \frac{1}{2} \left( \sum_{p\in S_V^+} \frac{\phi_p}{r_p} - \phi_j \right) \right] (\widetilde C_i - \widetilde C_j)}{1 + \frac{\Delta t}{2V_i} \sum_{j\in S_V^+} \left| Q_j \right| \left( \sum_{q\in S_V^-} \frac{\psi_q}{s_q} - \psi_j \right)} = C_i^{M+1}
\end{equation}

where $r_p$ and $s_q$ are upwind and downwind ratios respectively:

\begin{equation}
\label{eq22}
\begin{aligned}
r_p &= \frac{\sum_{q\in S_V^-}\left| Q_q \right| (\widetilde C_q - \widetilde C_i)}{\left| Q_p \right| (\widetilde C_i - \widetilde C_p)}, p \in S_V^+\\
s_q &= \frac{(\widetilde C_i - C_i^{M+1}) \sum_{p\in S_V^+} \left| Q_p \right|}{\left| Q_p \right| (\widetilde C_q - C_q^{M+1})}, q\in S_V^-
\end{aligned}
\end{equation}

[DB07](#duraisamy2007) showed that a sufficient TVD condition for Eq. $\ref{eq21}$ is that the coefficient of the 2nd LHS term be non-negative, i.e.:

\begin{equation}
\label{eq23}
1 + \frac{1}{2}\left(\sum_{p\in S_V^+} \frac{\phi_p}{r_p} - \phi_j \right) \geq 0
\end{equation}

\begin{equation}
\label{eq24}
1 + \frac{\Delta t}{2V_i} \sum_{j\in S_V^+} \left| Q_j \right| \left( \sum_{q\in S_V^-} \frac{\psi_q}{s_q} - \psi_j \right) \geq \delta \gt 0
\end{equation}

where $\delta$ is a small positive number. Eq. $\ref{eq23}$ can be satisfied with any choice of standard limiter functions in space, and Eq. $\ref{eq24}$ must be solved together with Eq. $\ref{eq21}$ iteratively, because $\psi$ and $s_q$ are functions of $\widetilde C$. We need to discuss 3 scenarios for prism $i$:

**Scenario 1.**  vertically convergent flow: in this case, the outer sum in Eq. $\ref{eq24}$ is 0, so the inequality is always true;

**Scenario 2.** divergent flow: the numerator of the 2nd LHS term in Eq. $\ref{eq21}$ is 0, and so $\widetilde C_i = C_i^{M+1}$;

**Scenario 3.** uni-directional flow (either upward or downward): in this case, prism $i$ has exactly 1 inflow and 1 outflow face vertically, so a sufficient condition for Eq. $\ref{eq24}$ is:

\begin{equation}
\label{eq25}
1 - \frac{\Delta t}{2 V_i} \left| Q_j \right| \psi_j \geq \delta \gt 0, j\in S_V^+
\end{equation}

Therefore we choose the following form for the limiter:

\begin{equation}
\label{eq26}
\psi_j = \text{max}\left[0, \text{min}\left[1, \frac{2(1-\delta)V_i}{\left| Q_j \right|\Delta t} \right] \right], j\in S_V^+
\end{equation}

where we have imposed a maximum of 1 in an attempt to obtain 2nd-order accuracy in time. Note that the limiter is a function of the vertical Courant number: it decreases as the Courant number increases. Eqs. $\ref{eq21}$ and $\ref{eq26}$ are then solved using a simple Picard iteration method starting from $\psi = 0$ everywhere, and fast convergence within a few iterations is usually observed.

Simple benchmark tests indicate that $TVD^2$ is accurate for a wide range of Courant numbers as found in typical geophysical flows. It works equally well in eddying and non-eddying regimes, from very shallow to very deep depths, and is thus ideal for cross-scale applications. You are encouraged to use this option as much as possible.

## Third-order WENO scheme

This option starts from the same Eqs ($\ref{eq15}$ - $\ref{eq17}$), but solves Eq. ($\ref{eq15}$) using a third-order WENO scheme.
Essentially we use a higher-order reconstruction method to approximate the numerical flux and the details can be found in 
 [Ye et al. (2019)](http://ccrm.vims.edu/yinglong/wiki_files/Ye_etal_OM_2019-SCHISM-WENO.pdf).

## Hybridization with Eulerian-Lagrangian Method

The FV schemes described above all have explicit component (in the horizontal dimension) that is subject to stability 
 constraint (Courant condition). This constraint can become severe, e.g. with high mesh resolution in the shallows in watershed
 applications. To alleviate this constaint, SCHISM allows local hybridization between any of the FV schemes and the Eulerian-Lagrangian Method (ELM). Like the ELM scheme used for momentum advection, the ELM scheme for the transport advection is unconditionally
 stable and monotone (if a linear interpolation is used at the foot of the characteristic line). The only downside of the ELM
 is that it does not conserve mass in Eulerian sense (but it does in Lagrangian sense), and therefore should be used sparingly,
 i.e., only to locally speed up the transport solver.

The user can invoke this scheme by setting `ielm_transport=1` and prescribing a maximum allowed number of sub-sucyclings per
 time step `max_subcyc`. To limit the use of ELM only for extreme cases, it's important to set a proper `max_subcyc`. A
 rule of thumb is that a sub time step (for the explicit FV sovler) should be around 10 sec, so max_subcyc should be `dt/10`.

Another issue with this hybridized scheme is that the combination of WENO and ELM in shallows can sometimes lead to 
large numerical dispersion; see [Known issues](../known_issues.md#numerical-dispersion-with-weno) for details. A simple
 solution for this issue is to make the affected regions upwind via `tvd.prop`.

## Vertical movement
Many tracers have ‘behaviors’ in the form of vertical migration (upward or downward) in the water column. This is modeled with a ‘settling’ term: 

\begin{equation}
\label{eq27}
\frac{\partial C}{\partial t} + \nabla \cdot (\pmb{u}C) = \frac{\partial}{\partial z} \left( \kappa \frac{\partial C}{\partial z} \right) + \frac{\partial (w_s C)}{\partial z}
\end{equation}

where $w_s$ is the settling velocity (__positive downward__). This term is treated implicitly to avoid stability issues; in particular, it’s solved in the 3rd step together with the diffusion term in Eq. $\ref{eq17}$. The benefit of this approach is that often the settling term balances the diffusion at boundary (e.g., sediment).

## Horizontal B.C. for transport
In either upwind or TVD schemes, the concentration at the neighboring prism $T_j$ at the open boundary is known. For outflow, $T_j=T_i$ and the signal is advected out of the domain without hindrance. For incoming flow, $T_j$ is specified by the B.C. (either in `bctides.in` or `*.th`), and SCHISM nudges to this value with a relaxation constant (specified in `bctides.in`), in order to prevent sharp gradient there. For a complete list of horizontal B.C. supported by SCHISM, see [bctides.in](../input-output/bctides.md).

**References**

<span id="casulli2005">Casulli, V. and P. Zanolli (2005) High resolution methods for multidimensional advection–diffusion problems in free-surface hydrodynamics. Ocean Modelling, 10, pp.137-151.</span>

<span id="duraisamy2007">Duraisamy, K. and J.D. Baeder (2007), Implicit Scheme For Hyperbolic Conservation Laws Using Nonoscillatory Reconstruction In Space And Time, Siam J. Sci. Comput. 29(6), 2607–2620.</span>
