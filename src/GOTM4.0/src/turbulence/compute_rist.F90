!$Id: compute_rist.F90,v 1.3 2007-01-06 11:49:15 kbk Exp $
#include"cppdefs.h"
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!BOP
!
! !ROUTINE: Calculate steady-state Richardson number from c3\label{sec:Rist}
!
! !INTERFACE:
   REALTYPE function  compute_rist(c1,c2,c3)
!
! !DESCRIPTION:
! Numerically computes the steady-state Richardson-number $Ri_{st}$
! for two-equations models from the given
! $c_{\psi 3}$ and the parameters
! $c_{\psi 1}$ and $c_{\psi 2}$ according to \eq{Ri_st}.
! A (very tricky) double Newton-iteration is used to solve the resulting
! implicit non-linear equation.
!
! !USES:
   use turbulence, only:           as,an,cmue1,cmue2
   use turbulence, only:           cm0
   use turbulence, only:           turb_method,stab_method
   use turbulence, only:           cm0_fix,Prandtl0_fix
   use turbulence, only:           Constant
   use turbulence, only:           MunkAnderson
   use turbulence, only:           SchumGerz
   use turbulence, only:           EiflerSchrimpf
   IMPLICIT NONE
!
! !INPUT PARAMETERS:
   REALTYPE, intent(in)           :: c1,c2,c3
!
! !REVISION HISTORY:
!  Original author(s): Hans Burchard, Lars Umlauf
!
!EOP
!-----------------------------------------------------------------------
! !LOCAL VARIABLES:
     integer                      :: i,j,imax=100
     REALTYPE                     :: cc3,fc,fp,e=1.e-9,step,ann
     REALTYPE                     :: ffc,ffp,Ri,Rii,ee=1.e-4
     REALTYPE                     :: NN(0:2),SS(0:2)
     logical                      :: converged=.true.
!
!-----------------------------------------------------------------------
!BOC
   NN=0.
   SS=0.
   Ri=0.18  ! Initial guess for Rist

   do j=0,imax
      Rii=Ri
      ann=0.1
      do i=0,imax
         an(1)=ann
         as(1)=an(1)/Rii
         if (turb_method.eq.2) then
            select case(stab_method)
               case(Constant)
                  cmue1=cm0_fix
                  cmue2=cm0_fix/Prandtl0_fix
               case(MunkAnderson)
                  call cmue_ma(2)
               case(SchumGerz)
                  call cmue_sg(2)
               case(EiflerSchrimpf)
                  call cmue_rf(2)
            end select
         else
            call cmue_d(2)
         end if
         fc=cmue1(1)*an(1)/Rii-cmue2(1)*an(1)-cm0**(-3)
         an(1)=ann+e
         as(1)=an(1)/Rii
         if (turb_method.eq.2) then
            select case(stab_method)
               case(Constant)
                  cmue1=cm0_fix
                  cmue2=cm0_fix/Prandtl0_fix
               case(MunkAnderson)
                  call cmue_ma(2)
               case(SchumGerz)
                  call cmue_sg(2)
               case(EiflerSchrimpf)
                  call cmue_rf(2)
            end select
         else
            call cmue_d(2)
         end if
         fp=cmue1(1)*an(1)/Rii-cmue2(1)*an(1)-cm0**(-3)
         step=-fc/((fp-fc)/e)
         ann=ann+0.5*step

         if (abs(step).gt.100.) then
            STDERR '                                                 '
            STDERR 'Method for calculating the steady-state Richardson'
            STDERR 'number Rist does not converge.'
            STDERR 'Probably, the prescribed value for c3'
            STDERR 'is outside the range of the chosen stability function.'
            STDERR 'Please change gotmturb.nml accordingly.'
            STDERR 'You have chosen the stability function no.',stab_method
            STDERR 'If the problem persists, then use another'
            STDERR 'stability function.'
            STDERR 'I will continue anyway. Good luck!'
            STDERR '                                                 '
            converged=.false.
            goto 333
         endif

         if (abs(step).lt.1.e-10) goto 111
      end do
111   an(1)=ann
      as(1)=an(1)/Rii
      if (turb_method.eq.2) then
         select case(stab_method)
            case(Constant)
               cmue1=cm0_fix
               cmue2=cm0_fix/Prandtl0_fix
            case(MunkAnderson)
               call cmue_ma(2)
            case(SchumGerz)
               call cmue_sg(2)
            case(EiflerSchrimpf)
               call cmue_rf(2)
         end select
      else
         call cmue_d(2)
      end if
      cc3=c2+(c1-c2)/Rii*cmue1(1)/cmue2(1)
      ffc=cc3-c3

      Rii=Ri+ee
      ann=0.1
      do i=0,imax
         an(1)=ann
         as(1)=an(1)/Rii
         if (turb_method.eq.2) then
            select case(stab_method)
               case(Constant)
                  cmue1=cm0_fix
                  cmue2=cm0_fix/Prandtl0_fix
               case(MunkAnderson)
                  call cmue_ma(2)
               case(SchumGerz)
                  call cmue_sg(2)
               case(EiflerSchrimpf)
                  call cmue_rf(2)
            end select
         else
            call cmue_d(2)
         end if
         fc=cmue1(1)*an(1)/Rii-cmue2(1)*an(1)-cm0**(-3)
         an(1)=ann+e
         as(1)=an(1)/Rii
         if (turb_method.eq.2) then
            select case(stab_method)
               case(Constant)
                  cmue1=cm0_fix
                  cmue2=cm0_fix/Prandtl0_fix
               case(MunkAnderson)
                  call cmue_ma(2)
               case(SchumGerz)
                  call cmue_sg(2)
               case(EiflerSchrimpf)
                  call cmue_rf(2)
            end select
         else
            call cmue_d(2)
         end if
         fp=cmue1(1)*an(1)/Rii-cmue2(1)*an(1)-cm0**(-3)
         step=-fc/((fp-fc)/e)
         ann=ann+0.5*step
         if (abs(step).gt.100.) then
            STDERR 'Method for calculating the steady-state Richardson'
            STDERR 'number Rist does not converge.'
            STDERR 'Probably, the prescribed value for c3'
            STDERR 'is outside the range of the chosen stability function.'
            STDERR 'Please change gotmturb.nml accordingly.'
            STDERR 'You have chosen the stability function no.',stab_method
            STDERR 'If the problem persists, then use another'
            STDERR 'stability function.'
            STDERR 'I will continue anyway. Good luck!'
            STDERR '                                                 '
            converged=.false.
            goto 333
         endif
         if (abs(step).lt.1.e-10) goto 222
      end do
222   an(1)=ann
      as(1)=an(1)/Rii
      if (turb_method.eq.2) then
         select case(stab_method)
            case(Constant)
               cmue1=cm0_fix
               cmue2=cm0_fix/Prandtl0_fix
            case(MunkAnderson)
               call cmue_ma(2)
            case(SchumGerz)
               call cmue_sg(2)
            case(EiflerSchrimpf)
               call cmue_rf(2)
         end select
      else
         call cmue_d(2)
      end if
      cc3=c2+(c1-c2)/Rii*cmue1(1)/cmue2(1)
      ffp=cc3-c3

      step=-ffc/((ffp-ffc)/ee)
      Ri=Ri+0.25*step
      if (abs(step).gt.100.) then
         STDERR '                                                 '
         STDERR 'Method for calculating the steady-state Richardson'
         STDERR 'number Rist does not converge.'
         STDERR 'Probably, the prescribed value for c3'
         STDERR 'is outside the range of the chosen stability function.'
         STDERR 'Please change gotmturb.nml accordingly.'
         STDERR 'You have chosen the stability function no.',stab_method
         STDERR 'If the problem persists, then use another'
         STDERR 'stability function.'
         STDERR 'I will continue anyway. Good luck!'
         STDERR '                                                 '
         converged=.false.
         goto 333
    endif
      if (abs(step).lt.1.e-10) goto 333
   end do

333   compute_rist=Ri

   if (.not.converged) compute_rist=-999.

   return
   end function compute_rist
!EOC
!-----------------------------------------------------------------------
! Copyright by the GOTM-team under the GNU Public License - www.gnu.org
!-----------------------------------------------------------------------
