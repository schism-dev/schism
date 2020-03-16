!   Copyright 2014 College of William and Mary
!
!   Licensed under the Apache License, Version 2.0 (the "License");
!   you may not use this file except in compliance with the License.
!   You may obtain a copy of the License at
!
!     http://www.apache.org/licenses/LICENSE-2.0
!
!   Unless required by applicable law or agreed to in writing, software
!   distributed under the License is distributed on an "AS IS" BASIS,
!   WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
!   See the License for the specific language governing permissions and
!   limitations under the License.

!===============================================================================
!===============================================================================
! MORSELFE MISCELLANEOUS SUBROUTINES
!
! subroutine sed_settleveloc
! subroutine sed_taucrit
! subroutine sed_write_debug
! subroutine sed_comp_poro_noncoh
!
!===============================================================================
!===============================================================================

      SUBROUTINE sed_settleveloc(ised)
!--------------------------------------------------------------------!
! This subroutine computes settling velocity from Soulsby (1997)     !
!                                                                    !
! Author: florian ganthy (fganthy@lnec.pt ; florian.ganthy@gmail.com)!
! Date:   21/12/2012                                                 !
!                                                                    !
! History:                                                           !
!                                                                    !
!--------------------------------------------------------------------!

      USE schism_glbl, ONLY : rkind,errmsg
      USE schism_msgp, ONLY: parallel_abort
      USE sed_mod,   ONLY : rhom,g,Srho,Sd50,Wsed

      IMPLICIT NONE
      SAVE

!- Local variables --------------------------------------------------!
      
      INTEGER,INTENT(IN)    :: ised

      REAL(rkind),PARAMETER :: nu=1.36d-6 ! Cinematic viscosity
      REAL(rkind)           :: ratio,dstar,tmp

!- Start Statement --------------------------------------------------!

! - Preliminary parameters
      ratio    = Srho(ised)/rhom
      dstar    = Sd50(ised) * (g*(ratio-1.d0)/nu**2.d0)**(1.d0/3d0)

! - Settling velocity (in m.s-1)
      tmp=10.36d0**2.d0+1.049d0*dstar**3.d0
      if(tmp<0.d0) then
        write(errmsg,*)'SED3D; sed_misc_subs: tmp<0;',tmp
        call parallel_abort(errmsg)
      endif
      !Wsed(ised) = (nu/Sd50(ised)) * ((10.36d0**2.d0 + 1.049d0*      &
      !&                                dstar**3.d0)**0.5d0 - 10.36d0)
      Wsed(ised) = (nu/Sd50(ised)) * ( sqrt(tmp) - 10.36d0)
      if(Wsed(ised)<=0.d0) then
        write(errmsg,*)'SED3D; sed_misc_subs: Wsed<=0; ',Wsed(ised)
        call parallel_abort(errmsg)
      endif

!--------------------------------------------------------------------!
      END SUBROUTINE sed_settleveloc    

!===============================================================================
!===============================================================================

      SUBROUTINE sed_taucrit(ised)
!--------------------------------------------------------------------!
! This subroutine computes critical shear stres for erosion from     !
! Soulsby (1997)                                                     !
!                                                                    !
! Author: florian ganthy (fganthy@lnec.pt ; florian.ganthy@gmail.com)!
! Date:   21/12/2012                                                 !
!                                                                    !
! History:                                                           !
!                                                                    !
!--------------------------------------------------------------------!

      USE schism_glbl, ONLY : rkind,errmsg
      use schism_msgp, only: parallel_abort
      USE sed_mod,   ONLY : rhom,g,Srho,Sd50,tau_ce

      IMPLICIT NONE
      SAVE

!- Local variables --------------------------------------------------!
      
      INTEGER,INTENT(IN)    :: ised

      REAL(rkind),PARAMETER :: nu=1.36d-6 ! Cinematic viscosity
      REAL(rkind)           :: ratio,dstar,theta_cr

!- Start Statement --------------------------------------------------!

! - Preliminary parameters
      ratio    = Srho(ised)/rhom
      dstar    = Sd50(ised) * (g*(ratio-1.d0)/nu**2.d0)**(1.d0/3d0)
      if(1.d0+1.2d0*dstar==0.d0) call parallel_abort('SED3D; sed_misc_subs: dev. by 0')
      theta_cr = 0.3d0/(1.d0+1.2d0*dstar) + 0.055d0 *                &
      &         (1.d0-EXP(-0.02d0*dstar))

!'- Critical shear stress (N.m-2) - will be scaled to m^2/s/s in read_sed_input
      tau_ce(ised) = theta_cr*g*(Srho(ised)-rhom)*Sd50(ised)

      if(tau_ce(ised)<=0.d0) then
        write(errmsg,*)'SED3D; sed_misc_subs: tau_ce(ised)<=0; ',tau_ce(ised)
        call parallel_abort(errmsg)      
      endif

!--------------------------------------------------------------------!
      END SUBROUTINE sed_taucrit

!===============================================================================
!===============================================================================

      SUBROUTINE sed_write_debug(it)
!--------------------------------------------------------------------!
! This subroutine writes debug or additional information to a file   !
! e.g. mirror.out                                                    !
!                                                                    !
! Author: Ligia Pinto                                                !
! Date: xx/xx/xxxx                                                   !
!                                                                    !
! History: 2012/12 - F.Ganthy : form homogenisation of sediments     !
!          routines                                                  !
!                                                                    !
!--------------------------------------------------------------------!

      USE sed_mod
      USE schism_glbl, ONLY : nea,dt
      USE schism_msgp, ONLY : myrank

      IMPLICIT NONE

      INTEGER,INTENT(IN)     :: it
!- Local variables --------------------------------------------------!

      INTEGER :: i,j,k

!- Start Statement --------------------------------------------------!


      IF(myrank==0) THEN
        WRITE(12,*)'SED: sed_write_debug'
        WRITE(12,*)'it: ',it
        WRITE(12,*)'dt: ',dt
        WRITE(12,*)'nstp: ',nstp
        WRITE(12,*)'nnew: ',nnew

        DO i=1,nea
            WRITE(12,*)'nea:',i,' bustr:',real(bustr(i)),      &
     &            ' bvstr:',real(bvstr(i)),' tau_w:',real(tau_w(i))
        ENDDO !nea

        DO i=1,Nbed
          DO j=1,nea
            DO k=1,ntr_l
              WRITE(12,*)'Nbed:',i,      &
     &                ' nea:',j,' class:',k,' bed_frac:',real(bed_frac(i,j,k)),   &
     &                ' bed_mass(1):',real(bed_mass(i,j,1,k)),                       &
     &                ' bed_mass(2):',real(bed_mass(i,j,2,k))
            ENDDO !k
          ENDDO !j
        ENDDO !i

        DO i=1,Nbed
          DO j=1,nea
            WRITE(12,*) 'Nbed:',i,' nea:',j,                                       &
     &              ' bed_thick:',real(bed_thick(j)),' bed(ithck):',real(bed(i,j,ithck)),  &
     &              ' bed(iaged):',real(bed(i,j,iaged)),' bed(iporo):',real(bed(i,j,iporo))
          ENDDO !j
        ENDDO !i

        DO i=1,nea
          WRITE(12,*)'nea:',i,' bottom(itauc):',real(bottom(i,itauc)),                  &
     &            ' bottom(isd50):',real(bottom(i,isd50)),' bottom(iwsed):',real(bottom(i,iwsed)),&
     &            ' bottom(idens):',real(bottom(i,idens)),' bottom(iactv):',real(bottom(i,iactv))
        ENDDO !i

!!           DO i=1,nea
!!             WRITE(12,*)'nea:',i,                    &
!!     &              ' bottom(nea,itauc):',bottom(i,itauc),                     &
!!     &              ' bottom(nea,isd50):',bottom(i,isd50),                     &
!!     &              ' bottom(nea,iwsed):',bottom(i,iwsed),                     &
!!     &              ' bottom(nea,idens):',bottom(i,idens)
!!           ENDDO !nea

      ENDIF !myrank

!--------------------------------------------------------------------!
      END SUBROUTINE sed_write_debug

!!==============================================================================
  SUBROUTINE sed_comp_poro_noncoh(frac_sed, poro_gravsan)

   !--------------------------------------------------------------------------
   !                 ***  SUBROUTINE sed_comp_poro_noncoh  ***
   !
   ! ** Purpose : Compute the porosity representative of a non-cohesive matrix
   !
   ! ** Called by : sed_init and sediment
   !
   ! ** References :
   !    - Wooster et al. (2008) : Sediment supply and relative size
   !                              distribution effects on fine sediment
   !                              infiltration into immobile gravels
   ! ** History :
   !       !  2020-02  B. Mengual
   !
   !--------------------------------------------------------------------------

    USE sed_mod
    USE schism_glbl, ONLY : rkind

    IMPLICIT NONE
    include 'mpif.h'

    !! * Arguments
    REAL(rkind), INTENT(IN)    :: frac_sed(ntr_l)
    REAL(rkind), INTENT(OUT)   :: poro_gravsan

    !! * Local declarations
    INTEGER         ::  ised
    REAL(rkind), PARAMETER :: eps = 1.0d-4
    REAL(rkind)   ::    psim,siggeo2,siggeo,stdgeo, &
                        frac_gravsan,diam_gravsan
    REAL(rkind)   ::    frac_noncoh(ntr_l), psi_sed(ntr_l)

   !!---------------------------------------------------------------------------
   !! * Executable part

    frac_gravsan=0.0d0
    DO ised=1,ntr_l
      IF (iSedtype(ised) > 0) THEN
        frac_gravsan=frac_gravsan+frac_sed(ised)
      END IF
    END DO

    IF (frac_gravsan .GT. eps) THEN

      frac_noncoh(:)=0.0d0
      DO ised=1,ntr_l
        IF (iSedtype(ised) > 0) THEN
          frac_noncoh(ised)=frac_sed(ised)/frac_gravsan
        END IF
      END DO

      psim=0.0d0
      psi_sed(:)=0.0d0
      DO ised=1,ntr_l
        IF (iSedtype(ised) > 0) THEN
          psi_sed(ised)=LOG(Sd50(ised)*1000.0d0)/LOG(2.0d0)
          psim=psim+(psi_sed(ised)*frac_noncoh(ised))
        END IF
      END DO

      siggeo2=0.0d0
      DO ised=1,ntr_l
        IF (iSedtype(ised) > 0) THEN
          siggeo2=siggeo2+( frac_noncoh(ised)*(psi_sed(ised)-psim)**2.0d0 )
        END IF
      END DO

      siggeo  = siggeo2**(0.5d0)
      stdgeo  = MAX(1.0d0,MIN(2.0d0**(siggeo),4.0d0)) ! to remain in the range of values described by Wooster et al. (2008)

      diam_gravsan = (2.0d0**(psim))/1000.0d0 ! psi_sed is computed with diam_sed in mm --> /1e3
      poro_gravsan = Awooster*(stdgeo**(Bwooster))

    ELSE

      poro_gravsan = Awooster

    END IF

  END SUBROUTINE sed_comp_poro_noncoh
!!==============================================================================


