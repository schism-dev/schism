!$Id: algebraiclength.F90,v 1.7 2007-01-06 11:49:15 kbk Exp $
#include"cppdefs.h"
!-----------------------------------------------------------------------
!BOP
!
! !ROUTINE: Some algebraic length-scale relations \label{sec:algebraiclength}
!
! !INTERFACE:
   subroutine algebraiclength(method,nlev,z0b,z0s,depth,h,NN)
!
! !DESCRIPTION:
! This subroutine computes the vertical profile of the turbulent
! scale $l$ from different types of analytical expressions. These
! range from simple geometrical forms to more complicated expressions
! taking into account the effects of stratification and shear. The
! users can select their method in the input file {\tt gotmturb.nml}.
! For convenience, we define here $d_b$ and $d_s$ as the distance
! from the bottom and the surface, respectively. The water
! depth is then given by $H=d_b+d_s$, and $z_0^b$ and
! $z_0^s$ are the repective roughness lengths. With these
! abbreviations, the expressions implemented in GOTM are as follows.
! \begin{enumerate}
!  \item The parabolic profile is defined according to
!    \begin{equation}
!      l=\kappa \frac{(d_s+z_0^s) (d_b+z_0^b)}
!                    {d_s+d_b+z_0^b+z_0^s}
!    \comma
!    \end{equation}
!    where it should be noted that only for large water depth
!    this equation converges to $\kappa(z+z_0)$ near the bottom
!   or near the surface.
!  \item The triangular profile is defined according to
!    \begin{equation}
!       l = \kappa \, \min(d_s+z_0^s,d_b+z_0^b)
!    \comma
!    \end{equation}
!    which converges always to $\kappa(z+z_0)$ near the bottom
!   or near the surface.
!  \item A distorted parabola can be constructed by
!     using a slightly modified form of the equation
!     used by \cite{XingDavies95},
!     \begin{equation}
!        l = \kappa \frac{(d_s+z_0^s)(d_b^\text{Xing}+z_0^b)}
!                      {d_s+d_b^\text{Xing}+z_0^s+z_0^b}
!        \comma
!        d_b^\text{Xing} =
!        d_b \exp{\left(-\beta \frac{d_b}{H} \right)}
!    \comma
!    \end{equation}
!    where it should be noted that only for large water depth
!    this equation converges to $\kappa(z+z_0)$ near the bottom
!   or near the surface. The constant $\beta$ is a form parameter
!   determining the distortion of the profile. Currently we use
!   $\beta = 2$ in GOTM.
!  \item A distorted parabola can be constructed by
!     using a slightly modified form of the equation
!     used by \cite{RobertOuellet87},
!    \begin{equation}
!       l = \kappa (d_b+z_0^b)
!           \sqrt{1-\frac{d_b-z_0^s}{H}}
!    \comma
!   \end{equation}
!    where it should be noted that only for large water depth
!    this equation converges to $\kappa(z+z_0)$ near the bottom.
!    Near the surface, the slope of $l$ is always different from
!    the law of the wall, a fact that becomes important when model
!    solutions for the case of breaking waves are computed, see
!    \sect{sec:analyse}.
!  \item Also the famous formula of \cite{Blackadar62} is based on
!     a parabolic shape, extended by an extra length--scale $l_a$.
!    Using the form of \cite{Luytenetal96}, the algebraic relation
!    is expressed by
!     \begin{equation}
!        l = \left( \frac{1}{\kappa (d_s+z_0^s)}
!                  +\frac{1}{\kappa (d_b+z_0^b)}
!                  +\dfrac{1}{l_a} \right)
!    \comma
!     \end{equation}
!    where
!   \begin{equation}
!        l_a = \gamma_0 \frac{\int_{-H}^\eta k^\frac{1}{2} z dz}
!                            {\int_{-H}^\eta k^\frac{1}{2} dz}
! \end{equation}
!    is the natural kinetic energy scale resulting from the
!    first moment of the rms turbulent velocity. The constant
!    $\gamma_0$ usually takes the value $\gamma_0 = 0.2$.
!    It should be noted that this expression for $l$
!    converges to $\kappa(z+z_0)$ at the surface and the bottom
!    only for large water depth, and when $l_a$ plays only a
!    minor role.
!  \item The so--called ISPRAMIX method to compute the length--scale
!   is described in detail in \sect{sec:ispramix}.
! \end{enumerate}
! After the length--scale has been computed, it is optionally
! limited by the method suggested by \cite{Galperinetal88}. This
! option can be activated in {\tt gotmturb.nml} by setting
! {\tt length\_lim = .true.} The rate of dissipation is computed
! according to \eq{epsilon}.
!
! !USES:
   use turbulence, only: L,eps,tke,k_min,eps_min
   use turbulence, only: cde,galp,kappa,length_lim
   IMPLICIT NONE
!
! !INPUT PARAMETERS:

!  type of length scale
   integer,  intent(in)                :: method

!  number of vertical layers
   integer,  intent(in)                :: nlev

!  surface and bottom roughness (m)
   REALTYPE, intent(in)                :: z0b,z0s

!  local depth (m)
   REALTYPE, intent(in)                :: depth

!  layer thicknesses (m)
   REALTYPE, intent(in)                :: h(0:nlev)

!  buoyancy frequency (1/s^2)
   REALTYPE, intent(in)                :: NN(0:nlev)
!
! !DEFINED PARAMETERS:
   integer, parameter                  :: Parabola=1
   integer, parameter                  :: Triangle=2
   integer, parameter                  :: Xing=3
   integer, parameter                  :: RobertOuellet=4
   integer, parameter                  :: Blackadar=5
   integer, parameter                  :: ispra_length=7
!
! !REVISION HISTORY:
!  Original author(s):  Manuel Ruiz Villarreal, Hans Burchard
!
!  $Log: algebraiclength.F90,v $
!  Revision 1.7  2007-01-06 11:49:15  kbk
!  namelist file extension changed .inp --> .nml
!
!  Revision 1.6  2005/11/15 11:35:02  lars
!  documentation finish for print
!
!  Revision 1.5  2005/06/27 13:44:07  kbk
!  modified + removed traling blanks
!
!  Revision 1.4  2003/03/28 09:20:35  kbk
!  added new copyright to files
!
!  Revision 1.3  2003/03/10 09:02:03  gotm
!  Added new Generic Turbulence Model + 
!  improved documentation and cleaned up code
!
!  Revision 1.2  2002/02/08 08:59:58  gotm
!
!  Revision 1.1.1.1  2001/02/12 15:55:58  gotm
!  initial import into CVS
!
!EOP
!-----------------------------------------------------------------------
! !LOCAL VARIABLES:
   integer                 :: i
   REALTYPE                :: ds,db,dbxing
   REALTYPE                :: beta,gamma,La,int_qz,int_q
   REALTYPE                :: Lcrit,L_min
!
!-----------------------------------------------------------------------
!BOC
!
! distance from bottom and surface initialised
   db=_ZERO_
   ds=_ZERO_
!
! parabolic shape
   select case (method)

      case(parabola)

         do i=1,nlev-1
            db=db+h(i)
            ds=depth-db
            L(i)=kappa*(ds+z0s)*(db+z0b)/(ds+db+z0b+z0s)
         end do

         L(0)    = kappa*z0b
         L(nlev) = kappa*z0s
!
! triangular shape

      case(triangle)
         do i=1,nlev-1
            db=db+h(i)
            ds=depth-db
            L(i)=kappa*min(ds+z0s,db+z0b)
         end do

         L(0)    = kappa*z0b
         L(nlev) = kappa*z0s
!
! modified Xing and Davies (1995)
! modification of parabolic mixing length
      case(Xing)
         beta = 2. ! a tuning parameter
         do i=1,nlev-1
           db=db+h(i)
           ds=depth-db
           ! lu changed a bug in here
           dbxing=db*exp(-beta*db/depth)
           L(i)=kappa*(ds+z0s)*(dbxing+z0b)/(ds+dbxing+z0s+z0b)
         end do

         L(0)    = kappa*z0b
         L(nlev) = kappa*z0s
!
! modified Robert and Ouellet(1987)
! modification of parabolic mixing length
      case(RobertOuellet)
         do i=1,nlev-1
           db=db+h(i)
           ds=depth-db
           ! lu changed a bug in here
           L(i)= kappa*(db+z0b)*sqrt(1-(db-z0s)/depth)
         end do

         L(0)    = kappa*z0b
         L(nlev) = kappa*(depth+z0b)*sqrt(z0s/depth) ! no log-law (!)
!
! Blackadar (1962).
! In the form suggested by Luyten et al. (1996) for two boundary layers.
      case(Blackadar)
         int_qz   = 0.
         int_q    = 0.

         do i=1,nlev-1
            db=db+h(i)
            ! compute the first moment turbulent velocity
            int_qz = int_qz + sqrt(tke(i))*(db+z0b)*h(i)
            ! compute vertically averaged turbulent velocity
            int_q  = int_q  + sqrt(tke(i))*h(i)
         end do

         gamma=0.2               ! an empirical factor
         La=gamma*int_qz/int_q   ! free turbulence length-scale

         db=0.0
         do i=1,nlev-1
            db=db+h(i)
            ds=depth-db
             L(i)=1/(1/(kappa*(ds+z0s))+1/(kappa*(db+z0b))+1/La)
         end do

         L(0)    = kappa*z0b
         L(nlev) = kappa*z0s
!
!  Ispramix
      case(ispra_length)
         call ispralength(nlev,NN,h,depth)

         L(0)    = kappa*z0b
         L(nlev) = kappa*z0s

      case default
   end select

   do i=0,nlev

!     clip the length-scale at the Galperin et al. (1988) value
!     under stable stratifcitation
      if ((NN(i).gt.0).and.(length_lim)) then
         Lcrit=sqrt(2*galp*galp*tke(i)/NN(i))
         if (L(i).gt.Lcrit) L(i)=Lcrit
      end if

!     compute the dissipation rate
      eps(i)=cde*sqrt(tke(i)*tke(i)*tke(i))/L(i)

     ! substitute minimum value
      if (eps(i).lt.eps_min) then
         eps(i) = eps_min
         L(i) = cde*sqrt(tke(i)*tke(i)*tke(i))/eps_min
      endif

   end do

   return
   end subroutine algebraiclength
!EOC

!-----------------------------------------------------------------------
! Copyright by the GOTM-team under the GNU Public License - www.gnu.org
!-----------------------------------------------------------------------
