!$Id: lengthscaleeq.F90,v 1.8 2007-01-06 11:49:15 kbk Exp $
#include"cppdefs.h"
!-----------------------------------------------------------------------
!BOP
!
! !ROUTINE: The dynamic q2l-equation\label{sec:lengthscaleeq}
!
! !INTERFACE:
   subroutine lengthscaleeq(nlev,dt,depth,u_taus,u_taub,z0s,z0b,h,NN,SS)

!
! !DESCRIPTION:
! Following suggestions of \cite{Rotta51a}, \cite{MellorYamada82}
! proposed an equation for the product $q^2 l$ expressed by
! \begin{equation}
!   \label{MY}
!   \dot{\overline{q^2 l}}
!   = {\cal D}_l + l ( E_1  P + E_3 G - E_2  F \epsilon )
!   \comma
! \end{equation}
! where $\dot{\overline{q^2 l}}$ denotes the material derivative of $q^2 l$.
! The production terms $P$ and $G$ follow from \eq{PandG}, and $\epsilon$
! can be computed either directly from \eq{epsilonMY}, or from \eq{epsilon}
! with the help \eq{B1}.
!
! The so-called wall function, $F$, appearing in \eq{MY} is defined by
! \begin{equation}
!   \label{F}
!   F = 1 + E_2 \left( \dfrac{l}{\kappa {\cal L}_z} \right)^2
!   \comma
! \end{equation}
! $\kappa$ being the von K{\'a}rm{\'a}n constant and ${\cal L}_z$ some
! measure for the distance from the wall. Different possiblities
! for  ${\cal L}_z$ are implemented in GOTM, which can be activated
! be setting the parameter {\tt MY\_length} in {\tt gotmturb.nml} to
! appropriate values. Close to the wall, however, one always has
! ${\cal L}_z= \overline{z}$, where $\overline{z}$ is the distance from
! the wall.
!
! For horizontally homogeneous flows, the transport term ${\cal D}_l$
! appearing in \eq{MY} is expressed by a simple gradient formulation,
! \begin{equation}
!   \label{diffusionMYlength}
!   {\cal D}_l = \frstder{z} \left( q l S_l \partder{q^2 l}{z} \right)
!  \comma
! \end{equation}
! where $S_l$ is a constant of the model. The values for the model
! constants recommended by \cite{MellorYamada82} are displayed in
! \tab{tab:MY_constants}. They can be set in {\tt gotmturb.nml}. Note,
! that the parameter $E_3$ in stably stratifed flows is in principle
! a function of the so-called steady state Richardson-number,
! as discussed by \cite{Burchard2001c}, see discussion in the context
! of \eq{Ri_st}.
! \begin{table}[ht]
!   \begin{center}
! \begin{tabular}{ccccccc}
!                           & $B_1$  & $S_q$ & $S_l$ & $E_1$ & $E_2$ & $E_3$    \\[1mm]
!      \hline
!     \cite{MellorYamada82} & $16.6$ & $0.2$ & $0.2$ & $1.8$ & $1.33$ & $1.8$\\
!   \end{tabular}
!   \caption{\label{tab:MY_constants} Constants appearing in \eq{MY}
!     and \eq{epsilonMY}}
!   \end{center}
! \end{table}
!
! At the end of this routine the length-scale can be constrained according to a
! suggestion of \cite{Galperinetal88}. This feature is optional and can be activated
! by setting {\tt length\_lim = .true.} in {\tt gotmturb.nml}.
!
! !USES:
   use turbulence, only: P,B
   use turbulence, only: tke,tkeo,k_min,eps,eps_min,L
   use turbulence, only: kappa,e1,e2,e3,b1
   use turbulence, only: MY_length,cm0,cde,galp,length_lim
   use turbulence, only: q2l_bc, psi_ubc, psi_lbc, ubc_type, lbc_type
   use turbulence, only: sl
   use util,       only: Dirichlet,Neumann

   IMPLICIT NONE
!
! !INPUT PARAMETERS:

!  number of vertical layers
   integer,  intent(in)                :: nlev

!  time step (s)
   REALTYPE, intent(in)                :: dt

!  local water depth (m)
   REALTYPE, intent(in)                :: depth

!  surface and bottom
!  friction velocity (m/s)
   REALTYPE, intent(in)                :: u_taus,u_taub

!  surface and bottom
!  roughness length (m)
   REALTYPE, intent(in)                :: z0s,z0b

!  layer thickness (m)
   REALTYPE, intent(in)                :: h(0:nlev)

!  square of shear and buoyancy
!  frequency (1/s^2)
   REALTYPE, intent(in)                :: NN(0:nlev),SS(0:nlev)
!
! !REVISION HISTORY:
!  Original author(s): Lars Umlauf
!                     (re-write after first version of
!                      H. Burchard and K. Bolding
!
!  $Log: lengthscaleeq.F90,v $
!  Revision 1.8  2007-01-06 11:49:15  kbk
!  namelist file extension changed .inp --> .nml
!
!  Revision 1.7  2005/11/15 11:35:02  lars
!  documentation finish for print
!
!  Revision 1.6  2005/11/03 20:53:37  hb
!  Patankar trick reverted to older versions for 
!  stabilising 3D computations
!
!  Revision 1.5  2005/06/27 13:44:07  kbk
!  modified + removed traling blanks
!
!  Revision 1.4  2003/03/28 09:20:35  kbk
!  added new copyright to files
!
!  Revision 1.3  2003/03/10 09:02:05  gotm
!  Added new Generic Turbulence Model + 
!  improved documentation and cleaned up code
!
!
!EOP
!------------------------------------------------------------------------
!
! !LOCAL VARIABLES:
   REALTYPE                  :: DiffQ2lup,DiffQ2ldw,pos_bc
   REALTYPE                  :: prod,buoyan,diss
   REALTYPE                  :: prod_pos,prod_neg,buoyan_pos,buoyan_neg
   REALTYPE                  :: ki,epslim,NN_pos
   REALTYPE                  :: ds,db,Lcrit
   REALTYPE                  :: cnpar=_ONE_
   REALTYPE                  :: q2l(0:nlev),q3(0:nlev)
   REALTYPE                  :: avh(0:nlev)
   REALTYPE                  :: Lz(0:nlev)
   REALTYPE                  :: Lsour(0:nlev),Qsour(0:nlev)

   REALTYPE                  :: l_min

   integer                   :: i
!
!------------------------------------------------------------------------
!BOC

! compute lower bound for length scale
  l_min = cde*k_min**1.5/eps_min

!  some quantities in Mellor-Yamada notation
   do i=1,nlev-1
      q2l(i)=2.*tkeo(i)*L(i)
      q3 (i)=sqrt(8.*tke(i)*tke(i)*tke(i))
   end do

!  diagnostic length scale for wall function
   db=_ZERO_
   ds=_ZERO_
   do i=1,nlev-1
      db=db+h(i)
      ds=depth-db
      ! Parabola shape
      if (MY_length.eq.1) Lz(i)=kappa*(ds+z0s)*(db+z0b)/(ds+z0s+db+z0b)
      ! Triangle shape
      if (MY_length.eq.2) Lz(i)=kappa*min(ds+z0s,db+z0b)
      ! For infinite depth
      if (MY_length.eq.3) Lz(i)=kappa*(ds+z0s)
   end do

! prepare the production terms
   do i=1,nlev-1

!     compute diffusivity
      avh(i)      =  sl*sqrt(2.*tke(i))*L(i)

!     compute production terms in q^2 l - equation
      prod        =  e1*L(i)*P(i)
      buoyan      =  e3*L(i)*B(i)
      diss        =  q3(i)/b1*(1.+e2*(L(i)/Lz(i))*(L(i)/Lz(i)))

!     compute positive and negative parts of RHS
      if (prod+buoyan .gt. 0) then
         Qsour(i) =  prod + buoyan
         Lsour(i) = -diss/q2l(i)
      else
         Qsour(i) =  prod 
         Lsour(i) = -(diss-buoyan)/q2l(i)
      end if

   end do

!  TKE and position for upper BC
   if (psi_ubc.eq.Neumann) then
!     tke at center "nlev"
      ki = tke(nlev-1)

!     flux at center "nlev"
      pos_bc = 0.5*h(nlev)
   else
!     tke at face "nlev-1"
      ki = tke(nlev-1)

!     value at face "nlev-1"
      pos_bc = h(nlev)
   end if

!  obtain BC for upper boundary of type "ubc_type"
   DiffQ2lup  = q2l_bc(psi_ubc,ubc_type,pos_bc,ki,z0s,u_taus)


!  TKE and position for lower BC
   if (psi_lbc.eq.Neumann) then
!     tke at center "1"
      ki = tke(1)

!     flux at center "1"
      pos_bc = 0.5*h(1)
   else
!     tke at face "1"
      ki = tke(1)

!     value at face "1"
      pos_bc = h(1)
   end if

!  obtain BC for lower boundary of type "lbc_type"
   DiffQ2ldw  = q2l_bc(psi_lbc,lbc_type,pos_bc,ki,z0b,u_taub)


!  do diffusion step
   call diff_face(nlev,dt,cnpar,h,psi_ubc,psi_lbc,                          &
                  DiffQ2lup,DiffQ2ldw,avh,Lsour,Qsour,q2l)


!  fill top and bottom value with something nice
!  (only for output)
   q2l(nlev)  = q2l_bc(Dirichlet,ubc_type,z0s,tke(nlev),z0s,u_taus)
   q2l(0   )  = q2l_bc(Dirichlet,lbc_type,z0b,tke(0   ),z0b,u_taub)



! compute L and epsilon
  do i=0,nlev
     L(i)=q2l(i)/(2.*tke(i))

!    apply the length-scale clipping of Galperin et al. (1988)
     if ((NN(i).gt.0).and.(length_lim)) then
        Lcrit=sqrt(2*galp*galp*tke(i)/NN(i))
        if (L(i).gt.Lcrit) L(i)=Lcrit
     end if

!    compute dissipation rate
     eps(i) = cde*sqrt(tke(i)*tke(i)*tke(i))/L(i)

!    check for very small lengh scale
     if (L(i).lt.l_min) L(i)=l_min

!    substitute minimum value
     if (eps(i).lt.eps_min) then
        eps(i) = eps_min
          L(i) = cde*sqrt(tke(i)*tke(i)*tke(i))/eps_min
     endif
  end do


  return
  end subroutine lengthscaleeq
!EOC

!-----------------------------------------------------------------------
! Copyright by the GOTM-team under the GNU Public License - www.gnu.org
!-----------------------------------------------------------------------
