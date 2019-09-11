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
!   WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or
!   implied.
!   See the License for the specific language governing permissions and
!   limitations under the License.
 
          SUBROUTINE fib
!
!!= Marta Rodrigues=====================================================
!!= National Laboratory for Civil Enginnering                          !
!!======================================================================
!! March/2010 - Original code                                          ! 
!! June/2010 - Extension to include the settling due to aggregation    !
!!             with the sediments in the water column                  !
!!             (adapted from Lígia Pinto - MORSELFE, which was         !
!!             originally from ROMS-Licensed under a MIT/X style       !
!!             license. See License_ROMS.txt)                          !
!! February/2015 - Extension to include another option to compute      !
!!                 the decay of the FIB, based on Servais et al.,2007  ! 
!! November/2015 - Coupling with SCHISM                                !
!!======================================================Marta Rodrigues=
!!                                                                     !
!! This subroutine computes sources and sinks terms for the Fecal      !                 	
!! Indicator Bacteria (FIB) model                                      !
!!                                                                     !
!! Two tracers are considered, usually                                 !
!! 1. E. coli (or Fecal Coliforms)                                     !
!! 2. Enterococci                                                      !
!! but these tracers can be user-defined, based on the decay rates     !
!! due to mortality and sinking defined by the user                    ! 
!! (e.g. 1. Fecal Coliforms, 2. E. coli)                               !
!!======================================================================
!
      USE schism_glbl, only : bdy_frc,flx_sf,flx_bt,nea,nvrt,tr_el,& !tsel,
                            &idry_e,kbe,ze,npa,dt,eta2,dp,&
                            &ielg,errmsg,irange_tr,&
                            &iwater_type,srad,elnode,i34,flag_fib
      USE schism_msgp, only : myrank,parallel_abort
      USE fib_param

      IMPLICIT NONE
!      SAVE

      !real(r8),intent(IN) :: srad(npa)
      !integer,intent(IN) :: iwater_type(npa)
      real(r8), save :: solar_rad,rr,d_1,d_2
      real(r8), save :: sdp1,sdp2,solar_rad1,solar_rad2
      integer, save:: i,k,water_type,ifec,itmp1,itmp2
      
! Sinking variables
      integer, save :: ks 
      real(r8), save :: cffR,cffL,dltR,dltL,cu,cff
      real(r8), save :: Hz_inv3

      integer, save, allocatable :: ksource(:)
      real(r8), save, allocatable :: FC(:)
      real(r8), save, allocatable :: Hz_inv(:)
      real(r8), save, allocatable :: Hz_inv2(:)
      real(r8), save, allocatable :: qc(:)
      real(r8), save, allocatable :: qR(:)
      real(r8), save, allocatable :: qL(:)
      real(r8), save, allocatable :: WR(:)
      real(r8), save, allocatable :: WL(:)
      real(r8), save, allocatable :: Hz(:)
      

      if(.not.allocated(ksource)) allocate(ksource(nvrt))
      if(.not.allocated(FC)) allocate(FC(nvrt))
      if(.not.allocated(Hz_inv)) allocate(Hz_inv(nvrt))
      if(.not.allocated(Hz_inv2)) allocate(Hz_inv2(nvrt))
      if(.not.allocated(qc)) allocate(qc(nvrt))
      if(.not.allocated(qR)) allocate(qR(nvrt))
      if(.not.allocated(qL)) allocate(qL(nvrt))
      if(.not.allocated(WR)) allocate(WR(nvrt))
      if(.not.allocated(WL)) allocate(WL(nvrt))
      if(.not.allocated(Hz)) allocate(Hz(nvrt))

! Initialize some variables      
      Hz(:)=0.d0
      qL(:)=0.d0
      qR(:)=0.d0
      qc(:)=0.d0
      ks=0
      ksource(:)=0
      FC(:)=0.d0      
         
! Define the start-index and end-index for the tracers 
      itmp1=irange_tr(1,9)
      itmp2=irange_tr(2,9)      

! Calculates sources and sinks for SELFE
!----------------------------------------------------------------------------------------
! Decay rate (due to mortality)
!----------------------------------------------------------------------------------------
      IF(flag_fib==1)THEN ! User-defined decay rate 

        DO i=1,nea
          IF(idry_e(i)==1) CYCLE
            DO k=kbe(i)+1,nvrt              
               
              bdy_frc(itmp1,k,i)=-(kk_fib(i,1)/86400.d0)*tr_el(itmp1,k,i)
              bdy_frc(itmp2,k,i)=-(kk_fib(i,2)/86400.d0)*tr_el(itmp2,k,i)
            END DO                           
        END DO	 
			 			 
      ELSEIF(flag_fib==2)THEN  ! Canteras et al., 1995

        DO i=1,nea
          IF(idry_e(i)==1) CYCLE
            DO k=kbe(i)+1,nvrt    
         

! Compute solar radiation at elements (consider the light extiction in the
! water column)

               solar_rad=sum(srad(elnode(1:i34(i),i)))/i34(i)
               water_type=maxval(iwater_type(elnode(1:i34(i),i)))
!  Water type
               SELECT CASE(water_type)
                  CASE(1)
                    rr=0.58; d_1=0.35; d_2=23
                  CASE(2)
                    rr=0.62; d_1=0.60; d_2=20
                  CASE(3)
                    rr=0.67; d_1=1.00; d_2=17
                  CASE(4)
                    rr=0.77; d_1=1.50; d_2=14
                  CASE(5)
                    rr=0.78; d_1=1.40; d_2=7.9
                  CASE(6)
                    rr=0.62; d_1=1.50; d_2=20
                  CASE(7)
                    rr=0.80; d_1=0.90; d_2=2.1
                  CASE DEFAULT
                  call parallel_abort('Unknown water type (2)')
                END SELECT
                !MFR - Solar radiation at level
                !solar_rad=solar_rad*(rr*exp(ze(k,i)/d_1)+(1-rr)*exp(ze(k,i)/d_2))
                !solar_rad=max(solar_rad, 0.d0)
        
                !MFR - Solar radiation at prism center
                sdp1=min(ze(nvrt,i)-ze(k-1,i),500.d0) !to prevent underflow
                sdp2=min(ze(nvrt,i)-ze(k,i),500.d0) !to prevent underflow
                solar_rad1=solar_rad*(rr*exp(-sdp1/d_1)+(1-rr)*exp(-sdp1/d_2))
                solar_rad2=solar_rad*(rr*exp(-sdp2/d_1)+(1-rr)*exp(-sdp2/d_2))
                solar_rad=max(solar_rad2-solar_rad1,0.d0)/(ze(k,i)-ze(k-1,i))

! Compute decay rate                           
             
                kk_fib(i,:)=(2.533*1.040**(tr_el(1,k,i)-20.d0)*&
                           &1.012**tr_el(2,k,i)+0.113*solar_rad)
                bdy_frc(itmp1,k,i)=-(kk_fib(i,1)/86400.d0)*tr_el(itmp1,k,i)
                bdy_frc(itmp2,k,i)=-(kk_fib(i,2)/86400.d0)*tr_el(itmp2,k,i)
               
           ENDDO
        ENDDO

      ELSEIF(flag_fib==3)THEN ! Servais et al.,2007

        DO i=1,nea
          IF(idry_e(i)==1) CYCLE
            DO k=kbe(i)+1,nvrt

              kk_fib(i,:)=1.08*(exp(-(tr_el(1,k,i)-25.d0)**2/400.d0)/&
                            &exp(-25.d0/400.d0))

              bdy_frc(itmp1,k,i)=-(kk_fib(i,1)/86400.d0)*tr_el(itmp1,k,i)
              bdy_frc(itmp2,k,i)=-(kk_fib(i,2)/86400.d0)*tr_el(itmp2,k,i)
            END DO
        END DO

      ENDIF

!----------------------------------------------------------------------------------------
! Settling (due to aggregation with sediments in the water column)
! At the moment sinking is only available when using the 3D version of
! the code

      IF(nvrt>2)THEN !Using 3D model 
!----------------------------------------------------------------------------------------
! MFR - Adapted from MORSELFE (Author: Lígia Pinto)
! ---------------------------------------------------------------------------------------
! Compute thickness and actual depths in the middle of the volume of 
! control, Hz-layer thickness, z_w==ze - layer depth at RHO
! points and w points.
! ---------------------------------------------------------------------------------------
      DO i=1,nea    ! MFR - one cycle in the elements
        IF (idry_e(i)==1) CYCLE
	
        DO k=kbe(i)+1,nvrt
           Hz(k)=ze(k,i)-ze(k-1,i)
           IF(Hz(k)<=0) CALL parallel_abort('Hz<=0!')
        END DO

!
!...  Copy concentration of suspended sediment into scratch array "qc"
!...  (q-central, restrict it to be positive) which is hereafter
!...  interpreted as a set of grid-box averaged values for sediment
!...  concentration.
!
        SINK_LOOP: DO ifec=1,2        !MFR - ifec
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

          DO k=kbe(i)+1,nvrt     
            IF(tr_el(irange_tr(ifec,9),k,i)<1E-15_r8) tr_el(irange_tr(ifec,9),k,i)=0_r8  
            qc(k)=tr_el(irange_tr(ifec,9),k,i)
          END DO

!
!----------------------------------------------------------------------------------------
!  Vertical sinking of bacteria attached to suspended sediment.
!----------------------------------------------------------------------------------------
!
!...  Reconstruct vertical profile of suspended sediment "qc" in terms
!...  of a set of parabolic segments within each grid box. Then, compute
!...  semi-Lagrangian flux due to sinking.
!
!          IF(myrank==0) write(16,*)'reconstruct vertical'

          Hz_inv2=0.d0
          DO k=kbe(i)+1,nvrt-1
             Hz_inv2(k)=1.0_r8/(Hz(k)+Hz(k+1))
          END DO

          DO k=nvrt-1,kbe(i)+1,-1   !LLP
            FC(k)=(qc(k+1)-qc(k))*(Hz_inv2(k))
          END DO

          cffR=0.d0
          cffL=0.d0
          dltR=0.d0
          dltL=0.d0
          DO k=kbe(i)+2,nvrt-1       
            dltR=Hz(k)*FC(k)
            dltL=Hz(k)*FC(k-1)
            cff=Hz(k-1)+2.0_r8*Hz(k)+Hz(k+1)
            cffR=cff*FC(k)
            cffL=cff*FC(k-1)
!
!...  Apply PPM monotonicity constraint to prevent oscillations within the
!...  grid box.
!
            IF (dltR*dltL.le.0.0_r8) THEN
              dltR=0.0_r8
              dltL=0.0_r8
            ELSE IF (ABS(dltR).gt.ABS(cffL)) THEN
              dltR=cffL
            ELSE IF (ABS(dltL).gt.ABS(cffR)) THEN
              dltL=cffR
            END IF

!...  Compute right and left side values (qR,qL) of parabolic segments
!...  within grid box Hz(k); (WR,WL) are measures of quadratic variations. 
!
!...  NOTE: Although each parabolic segment is monotonic within its grid
!...        box, monotonicity of the whole profile is not guaranteed,
!...        because qL(k+1)-qR(k) may still have different sign than
!...        qc(k+1)-qc(k).  This possibility is excluded, after qL and qR
!...        are reconciled using WENO procedure.
!
            Hz_inv3=0.d0
            Hz_inv3=1.0_r8/(Hz(k-1)+Hz(k)+Hz(k+1))
            cff=(dltR-dltL)*Hz_inv3
            dltR=dltR-cff*Hz(k+1)
            dltL=dltL+cff*Hz(k-1)
            qR(k)=qc(k)+dltR
            qL(k)=qc(k)-dltL
            WR(k)=(2.0_r8*dltR-dltL)**2
            WL(k)=(dltR-2.0_r8*dltL)**2
          END DO !k=kbe(i)+2,nvrt-1

          cff=1.0E-14_r8

          DO k=kbe(i)+2,nvrt-2             
            dltL=MAX(cff,WL(k))
            dltR=MAX(cff,WR(k+1))
            qR(k)=(dltR*qR(k)+dltL*qL(k+1))/(dltR+dltL)
            qL(k+1)=qR(k)
          END DO

          FC(nvrt)=0.0_r8              ! no-flux boundary condition
!#  if defined LINEAR_CONTINUATION !MFR - use default
!          qL(nvrt,i)=qR(nvrt-1,i) !MFR - use default
!          qR(nvrt,i)=2.0_r8*qc(nvrt,i)-qL(nvrt,i) !MFR - use default
!#  elif defined NEUMANN !MFR - use default
!          qL(nvrt,i)=qR(nvrt-1,i) !MFR - use default
!          qR(nvrt,i)=1.5*qc(nvrt,i)-0.5_r8*qL(nvrt,i) !MFR - use default
!#  else !MFR - use default
          qR(nvrt)=qc(nvrt)         ! default strictly monotonic
          qL(nvrt)=qc(nvrt)         ! conditions
          qR(nvrt-1)=qc(nvrt)
!#  endif !MFR - use default
!#  if defined LINEAR_CONTINUATION  !MFR - use default
!            qR(kbe(i)+1,i)=qL(kbe(i)+2,i)   !MFR - use default  
!            qL(kbe(i)+1,i)=2.0_r8*qc(kbe(i)+1,i)-qR(kbe(i)+1,i)  !MFR - use default 
!#  elif defined NEUMANN  !MFR - use default
!            qR(kbe(i)+1,i)=qL(kbe(i)+2,i)  !MFR - use default   
!            qL(kbe(i)+1,i)=1.5_r8*qc(kbe(i)+1,i)-0.5_r8*qR(kbe(i)+1,i)  !MFR - use default   
!#  else  !MFR - use default 
          qL(kbe(i)+2)=qc(kbe(i)+1)                 ! bottom grid boxes are 
          qR(kbe(i)+1)=qc(kbe(i)+1)                 ! re-assumed to be    
          qL(kbe(i)+1)=qc(kbe(i)+1)                 ! piecewise constant. 
!#  endif !MFR - use default

!...  Apply monotonicity constraint again, since the reconciled interfacial
!...  values may cause a non-monotonic behavior of the parabolic segments
!...  inside the grid box.
!
          DO k=kbe(i)+1,nvrt    
            dltR=qR(k)-qc(k)
            dltL=qc(k)-qL(k)
            cffR=2.0_r8*dltR
            cffL=2.0_r8*dltL
            IF (dltR*dltL.lt.0.0_r8) THEN
              dltR=0.0_r8
              dltL=0.0_r8
            ELSE IF (ABS(dltR).gt.ABS(cffL)) THEN
              dltR=cffL
            ELSE IF (ABS(dltL).gt.ABS(cffR)) THEN
              dltL=cffR
            END IF
            qR(k)=qc(k)+dltR
            qL(k)=qc(k)-dltL
          END DO !k=kbe(i)+1,nvrt

!...  After this moment reconstruction is considered complete. The next
!...  stage is to compute vertical advective fluxes, FC. It is expected
!...  that sinking may occur relatively fast, the algorithm is designed
!...  to be free of CFL criterion, which is achieved by allowing
!...  integration bounds for semi-Lagrangian advective flux to use as
!...  many grid boxes in upstream direction as necessary.
!
!...  In the two code segments below, WL is the z-coordinate of the
!...  departure point for grid box interface z_w (ze)  with the same indices;
!...  FC is the finite volume flux; ksource(k) is index of vertical
!...  grid box which contains the departure point (restricted by N(ng)). 
!...  During the search: also add in content of whole grid boxes
!...  participating in FC.
!
!        if(myrank==0) write(16,*)dt,Wsed(ised)
   
        !cff=dt*ABS(Wsed(ised)) !MFR

	  cff=dt*ABS(sink_fib(i))    !MFR
          DO k=kbe(i)+1,nvrt
            FC(k-1)=0.0_r8
            WL(k)=ze(k-1,i)+cff      !layer depth+sinking(dt*wsed)
            WR(k)=Hz(k)*qc(k)    !thickness*sediment concentration
            ksource(k)=k             !index of vertical grid box of departurepoint
          END DO !k=kbe(i)+1,nvrt

          DO k=kbe(i)+1,nvrt    
            DO ks=k,nvrt-1
              IF (WL(k).gt.ze(ks,i)) THEN
                ksource(k)=ks+1
                FC(k-1)=FC(k-1)+WR(ks)
              END IF
            END DO !ks=k,nvrt-1
          END DO !k=kbe(i)+1,nvrt

!  Finalize computation of flux: add fractional part.
!
!...  Compute inverse thicknessed to avoid repeated divisions.
          Hz_inv=0.d0 
          cu=0.d0
          DO k=kbe(i)+1,nvrt
            Hz_inv(k)=1.0_r8/Hz(k)
          END DO
            
          DO k=kbe(i)+1,nvrt   
            ks=ksource(k)
            IF(ks<kbe(i)+1.or.ks>nvrt) THEN
              write(errmsg,*)'SINKING: out of bound, ',ks,ielg(i),k
              CALL parallel_abort(errmsg)
            ENDIF
            cu=MIN(1.0_r8,(WL(k)-ze(ks-1,i))*Hz_inv(ks))
            FC(k-1)=FC(k-1)+Hz(ks)*cu* &
     &(qL(ks)+cu*(0.5_r8*(qR(ks)-qL(ks))-(1.5_r8-cu)*(qR(ks)+qL(ks)-2.0_r8*qc(ks))))
          END DO !k

          DO k=kbe(i)+1,nvrt    
             qc(k)=(FC(k)-FC(k-1))*Hz_inv(k)
          END DO

          
!----------------------------------------------------------------------------------------
!...  Update global tracer variables (m Tunits).
!----------------------------------------------------------------------------------------
!... divide by dt 
          DO k=kbe(i)+1,nvrt     !new
            bdy_frc(irange_tr(ifec,9),k,i)=bdy_frc(irange_tr(ifec,9),k,i)+&
	                      qc(k)*fraction_fib(i)/dt  !LLP !+ MFR
          END DO
        
        END DO SINK_LOOP
      
      END DO !i=1,nea    !MFR
      
      END IF !using 3D model 	  
!----------------------------------------------------------------------------------------
! MFR - End bacteria attached to suspend sediment sinking	  
!----------------------------------------------------------------------------------------	  


! To avoid that the bdy_frc is larger that the tracer concentration
! (another method can be explored later ....)
      DO i=1,nea
          IF(idry_e(i)==1) CYCLE
          DO k=kbe(i)+1,nvrt
              IF(tr_el(itmp1,k,i)+bdy_frc(itmp1,k,i)*dt<1E-7_r8) bdy_frc(itmp1,k,i)=-tr_el(itmp1,k,i)/dt
              IF(tr_el(itmp2,k,i)+bdy_frc(itmp2,k,i)*dt<1E-7_r8) bdy_frc(itmp2,k,i)=-tr_el(itmp2,k,i)/dt 
          ENDDO
      ENDDO  

! ---------------------------------------------------------------------------------------
! Calculates bottom flux for SELFE
! ---------------------------------------------------------------------------------------
!       flx_bt(mntr,nea)
        flx_bt(itmp1:itmp2,:) = 0.d0
!
! ---------------------------------------------------------------------------------------
! Calculates surface flux for SELFE
! ---------------------------------------------------------------------------------------
!       flx_sf(mntr,nea)
        flx_sf(itmp1:itmp2,:) = 0.d0
! ---------------------------------------------------------------------------------------

      RETURN
      END SUBROUTINE fib
