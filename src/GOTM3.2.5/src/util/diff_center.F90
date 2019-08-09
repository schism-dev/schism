!$Id: diff_center.F90,v 1.1 2005-06-27 10:54:33 kbk Exp $
#include"cppdefs.h"
!-----------------------------------------------------------------------
!BOP
!
! !ROUTINE: Diffusion schemes --- grid centers\label{sec:diffusionMean}
!
! !INTERFACE:
   subroutine diff_center(N,dt,cnpar,h,Bcup,Bcdw,                       &
                          Yup,Ydw,nuY,Lsour,Qsour,Taur,Yobs,Y)
!
! !DESCRIPTION:
! This subroutine solves the one-dimensional diffusion equation
! including source terms,
!  \begin{equation}
!   \label{YdiffCenter}
!    \partder{Y}{t}
!    = \partder{}{z} \left( \nu_Y \partder{Y}{z} \right)
!    - \frac{1}{\tau_R}(Y-Y_{obs})
!    + Y L_{\text{sour}} + Q_{\text{sour}}
!    \comma
!  \end{equation}
! for al variables defined at the centers of the grid cells, and
! a diffusion coefficient $\nu_Y$ defined at the faces.
! Relaxation with time scale $\tau_R$ towards observed values
! $Y_{\text{obs}}$ is possible. $L_{\text{sour}}$ specifies a
! linear source term, and $Q_{\text{sour}}$ a constant source term.
! Central differences are used to discretize the problem
! as discussed in \sect{SectionNumericsMean}. The diffusion term,
! the linear source term, and the linear part arising from the
! relaxation term are treated
! with an implicit method, whereas the constant source term is treated
! fully explicit.
!
! The input parameters {\tt Bcup} and {\tt Bcdw} specify the type
! of the upper and lower boundary conditions, which can be either
! Dirichlet or Neumann-type. {\tt Bcup} and {\tt Bcdw} must have integer
! values corresponding to the parameters {\tt Dirichlet} and {\tt Neumann}
! defined in the module {\tt util}, see \sect{sec:utils}.
! {\tt Yup} and {\tt Ydw} are the values of the boundary conditions at
! the surface and the bottom. Depending on the values of {\tt Bcup} and
! {\tt Bcdw}, they represent either fluxes or prescribed values.
!
! Note that fluxes \emph{entering} a boundary cell are counted positive
! by convention. The lower and upper position for prescribing these fluxes
! are located at the lowest und uppermost grid faces with index "0" and
! index "N", respectively. If values are prescribed, they are located at
! the centers with index "1" and index "N", respectivly.
!
! !USES:
   use util,          only  : Dirichlet, Neumann
   use mtridiagonal
!
! !INPUT PARAMETERS:

!  number of vertical layers
   integer,  intent(in)                :: N

!  time step (s)
   REALTYPE, intent(in)                :: dt

!  "implicitness" parameter
   REALTYPE, intent(in)                :: cnpar

!  layer thickness (m)
   REALTYPE, intent(in)                :: h(0:N)

!  type of upper BC
   integer,  intent(in)                :: Bcup

!  type of lower BC
   integer,  intent(in)                :: Bcdw

!  value of upper BC
   REALTYPE, intent(in)                :: Yup

!  value of lower BC
   REALTYPE, intent(in)                :: Ydw

!  diffusivity of Y
   REALTYPE, intent(in)                :: nuY(0:N)

!  linear source term
!  (treated implicitly)
   REALTYPE, intent(in)                :: Lsour(0:N)

!  constant source term
!  (treated explicitly)
   REALTYPE, intent(in)                :: Qsour(0:N)

!  relaxation time (s)
   REALTYPE, intent(in)                :: Taur(0:N)

!  observed value of Y
   REALTYPE, intent(in)                :: Yobs(0:N)
!
! !INPUT/OUTPUT PARAMETERS:
   REALTYPE                            :: Y(0:N)
!
! !REVISION HISTORY:
!  Original author(s): Lars Umlauf
!
!  $Log: diff_center.F90,v $
!  Revision 1.1  2005-06-27 10:54:33  kbk
!  new files needed
!
!
!EOP
!
! !LOCAL VARIABLES:
   integer                   :: i
   REALTYPE                  :: a,c,l
!
!-----------------------------------------------------------------------
!BOC
!
!  set up matrix
   do i=2,N-1
      c     = 2.0*dt*nuY(i)  /(h(i)+h(i+1))/h(i)
      a     = 2.0*dt*nuY(i-1)/(h(i)+h(i-1))/h(i)
      l     =     dt*Lsour(i)

      cu(i) =-cnpar*c
      au(i) =-cnpar*a
      bu(i) = _ONE_ + cnpar*( a + c - l)
      du(i) = (_ONE_ - (_ONE_-cnpar)*(a + c - l))*Y(i)                  &
            + (_ONE_ - cnpar)*( a*Y(i-1) + c*Y(i+1) ) + dt*Qsour(i)
   end do

!   set up upper boundary condition
   select case(Bcup)
   case(Neumann)
      a     = 2.0*dt*nuY(N-1)/(h(N)+h(N-1))/h(N)
      l     = dt*Lsour(N)

      au(N) =-cnpar*a
      bu(N) =  _ONE_ - au(N) - cnpar*l
      du(N) = (_ONE_ + (_ONE_-cnpar)*l)*Y(N) + dt*(Qsour(N)+Yup/h(N))   &
            + (_ONE_ - cnpar)*a*(Y(N-1)-Y(N))
   case(Dirichlet)
      au(N) = _ZERO_
      bu(N) = _ONE_
      du(N) = Yup
   case default
      FATAL 'invalid boundary condition type for upper boundary'
      stop  'diff_center.F90'
   end select

!   set up lower boundary condition
   select case(Bcdw)
   case(Neumann)
      c     = 2.0*dt*nuY(1)/(h(1)+h(2))/h(1)
      l     = dt*Lsour(1)

      cu(1) =-cnpar*c
      bu(1) = _ONE_ - cu(1) - cnpar*l
      du(1) = (_ONE_ + (_ONE_-cnpar)*l)*Y(1) + dt*(Qsour(1)+Ydw/h(1))   &
            + (_ONE_ - cnpar)*c*(Y(2)-Y(1))
   case(Dirichlet)
      cu(1) = _ZERO_
      bu(1) = _ONE_
      du(1) = Ydw
   case default
      FATAL 'invalid boundary condition type for lower boundary'
      stop  'diff_center.F90'
   end select

   ! relaxation to observed value
   if (minval(Taur).lt.1.E10) then
      do i=1,N
         bu(i)=bu(i)+dt/Taur(i)
         du(i)=du(i)+dt/Taur(i)*Yobs(i)
      end do
   end if

!  solve linear system
   call tridiagonal(N,1,N,Y)

   return
   end subroutine diff_center
!EOC

!-----------------------------------------------------------------------
! Copyright by the GOTM-team under the GNU Public License - www.gnu.org
!-----------------------------------------------------------------------
