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


    SUBROUTINE bio_init

        
!!======================================================================
!! April/May, 2007 - Original code                                     ! 
!!======================================================Marta Rodrigues=
!!                                                                     !
!! This subroutine is from ROMS (where is called ana_biology.h):       !
!!                                                                     !
!! February, 2009 - Extended for oxygen cycle                          !
!!======================================================================	 
!!
!!======================================================================
!! Copyright (c) 2002-2007 The ROMS/TOMS Group                         !
!!   Licensed under a MIT/X style license                              !
!!   See License_ROMS.txt                                              !
!!                                                                     !
!=======================================================================
!                                                                      !
!  This routine sets initial conditions for biological tracer fields   !
!  using analytical expressions.                                       !
!                                                                      !
!=======================================================================
!
      USE bio_param
!      USE schism_glbl, only : trel0,nvrt,nea,tr_el,idry_e
      USE schism_glbl, only : tr_nd,nvrt,npa,irange_tr
      USE biology 

      IMPLICIT NONE

!  Imported variable declarations.
!      real(r8), dimension(ntracers,nvrt,nea), intent(out) :: trel0
!
!  Local variable declarations.
!
      integer :: i, iseco, itrc, j, k, tmp00, tmp1, tmp2, tmp3

      real(r8) :: cff1, cff2, cff3, cff4, cff5, cff6, cff7, cff8, cff9
      real(r8) :: cff10, cff11, cff12, cff13, cff14, cff15
      real(r8) :: salt, sftm, temp

!
!---------------------------------------------------------------------
!  EcoSim initial fields.
!---------------------------------------------------------------------
!
! Assumed maximum temperature gradient.
!
      cff3=1.0_r8/14.0_r8
      cff4=1.0_r8/16.0_r8
      cff5=32.0_r8
      cff7=1.0_r8/0.0157_r8
      cff8=1.0_r8/6.625_r8
      cff9=1.0_r8/16.0_r8
      cff10=1.0_r8/15.0_r8
      cff11=1.0_r8/8.0_r8
      cff12=1.0_r8/128.0_r8
      cff13=1.0_r8/1000.0_r8
      cff14=1.0_r8/12.0_r8
      cff15=cff5*cff8*cff14                  ! mole N : gram Chl

!Marta Rodrigues - For vertical and node      
      DO k= nvrt, 1, -1
        DO i=1, npa
          !if(idry_e(i)==1) cycle

! Initialization of surface chlorophyll.
!
!Marta Rodrigues - temperature and salinity values
            sftm=tr_nd(1,nvrt,i) !tsel(1,nvrt,i)	!sea surface temperature
            temp=tr_nd(1,k,i) !tsel(1,k,i)
            salt=tr_nd(2,k,i) !tsel(2,k,i)

	    	    
            cff1=-0.0827_r8*sftm+2.6386_r8
            cff2=MAX(0.00001_r8,cff1*(1.0_r8-(sftm-temp)*cff3))
!
! Initialization of nutrients.
!
            tmp00 = irange_tr(1,6)-1 
            
            tmp1 = tmp00+iNH4_
            tr_nd(tmp1,k,i)=0.053_r8*temp+0.7990_r8 !NH4
            tmp2 = tmp00+iNO3_
            tr_nd(tmp2,k,i)=8.5_r8-cff2*cff15-tr_nd(tmp1,k,i) !NO3
            tmp3 = tmp00+iPO4_
            tr_nd(tmp3,k,i)=(tr_nd(tmp1,k,i)+tr_nd(tmp2,k,i))*cff4 !PO4
            
            IF(IRON==1)THEN
              tmp1 = tmp00+iFeO_
              tr_nd(tmp1,k,i)=1.0_r8  !MFR
            ENDIF
!
! Assuming diatoms are 75% of initialized chlorophyll.
!
            tmp1 = tmp00+iSiO_
            tr_nd(tmp1,k,i)=5.5_r8-(cff2*0.75_r8)*cff15*1.20_r8
            tmp1 = tmp00+iDIC_
            tr_nd(tmp1,k,i)=2000.0_r8

! Marta Rodrigues - Initialization of DO and COD (oxygen)

            tmp1 = tmp00+iDO_
            tr_nd(tmp1,k,i)=150.0_r8
            tmp1 = tmp00+iCOD_
            tr_nd(tmp1,k,i)=0.0_r8

!
! Bacteria Initialization.
!    
              DO iseco=1,Nbac
                tmp1 = tmp00+iBacC(iseco)
                tr_nd(tmp1,k,i)=0.85_r8
                tmp2 = tmp00+iBacN(iseco)
                tr_nd(tmp2,k,i)=tr_nd(tmp1,k,i)*N2cBAC
                tmp2 = tmp00+iBacP(iseco)
                tr_nd(tmp2,k,i)=tr_nd(tmp1,k,i)*P2cBAC
                IF(IRON==1)THEN
                  tmp2 = tmp00+iBacF(iseco)
                  tr_nd(tmp2,k,i)=tr_nd(tmp1,k,i)*Fe2cBAC  !MFR
                ENDIF
              END DO
!

! Marta Rodrigues
! Zooplankton initialization.

            DO iseco=1,Nzoo
              tmp1 = tmp00+iZooC(iseco)
	      tr_nd(tmp1,k,i)=2_r8
              tmp1 = tmp00+iZooN(iseco)
	      tr_nd(tmp1,k,i)=0.2_r8
              tmp1 = tmp00+iZooP(iseco)
	      tr_nd(tmp1,k,i)=0.02_r8
	    END DO  

! Initialize phytoplankton populations.
!
!            tmp1 = tmp00+iPhyC(1)
!            tr_nd(tmp1,k,i)=MAX(0.02_r8,                            &
!     &                              0.75_r8*0.75_r8*cff5*cff2*cff14)
            DO iseco=1,Nphy
              tmp1 = tmp00+iPhyC(iseco)
              tr_nd(tmp1,k,i)=MAX(0.02_r8,                            &
     &                              0.75_r8*0.75_r8*cff5*cff2*cff14)
              tmp2 = tmp00+iPhyN(iseco) 
              tr_nd(tmp2,k,i)=tr_nd(tmp1,k,i)*cff8
              tmp3 = tmp00+iPhyP(iseco)
              tr_nd(tmp3,k,i)=tr_nd(tmp2,k,i)*cff4

              IF(IRON==1)THEN
                 tmp3 = tmp00+iPhyF(iseco)
                 tr_nd(tmp3,k,i)=tr_nd(tmp1,k,i)*cff13  !MFR
              ENDIF               
  
              IF (iPhyS(iseco).gt.0) THEN
                tmp3 = tmp00+iPhyS(iseco)
                tr_nd(tmp3,k,i)=tr_nd(tmp2,k,i)*1.20_r8
              END IF

!  Initialize Pigments in ugrams/liter (not umole/liter).
!
              cff6=12.0_r8/cff5
              tmp2 = tmp00+iPigs(iseco,1)
              tr_nd(tmp2,k,i)=cff6*tr_nd(tmp1,k,i)
!
!  Chlorophyll-b.
!
              cff6=cff5-b_C2Cl(iseco)
              IF (iPigs(iseco,2).gt.0) THEN
                 tmp3 = tmp00+iPigs(iseco,2) 
                 tr_nd(tmp3,k,i)=tr_nd(tmp2,k,i)*         &
     &                                  (b_ChlB(iseco)+                 &
     &                                   mxChlB(iseco)*cff6)
              END IF
!
!  Chlorophyll-c.
!
              IF (iPigs(iseco,3).gt.0) THEN
                 tmp3 = tmp00+iPigs(iseco,3)
                 tr_nd(tmp3,k,i)=tr_nd(tmp2,k,i)*         &
     &                                  (b_ChlC(iseco)+                 &
     &                                   mxChlC(iseco)*cff6)
              END IF
!
!  Photosynthetic Carotenoids.
!
              IF (iPigs(iseco,4).gt.0) THEN
                 tmp3 = tmp00+iPigs(iseco,4)
                 tr_nd(tmp3,k,i)=tr_nd(tmp2,k,i)*         &
     &                                  (b_PSC(iseco)+                  &
     &                                   mxPSC(iseco)*cff6)
              END IF
!
!  Photoprotective Carotenoids.
!
              IF (iPigs(iseco,5).gt.0) THEN
                 tmp3 = tmp00+iPigs(iseco,5)
                 tr_nd(tmp3,k,i)=tr_nd(tmp2,k,i)*         &
     &                                  (b_PPC(iseco)+                  &
     &                                   mxPPC(iseco)*cff6)
              END IF
!
!  Low Urobilin Phycoeurythin Carotenoids.
!
              IF (iPigs(iseco,6).gt.0) THEN
                 tmp3 = tmp00+iPigs(iseco,6)
                 tr_nd(tmp3,k,i)=tr_nd(tmp2,k,i)*         &
     &                                  (b_LPUb(iseco)+                 &
     &                                   mxLPUb(iseco)*cff6)
              END IF
!
!  High Urobilin Phycoeurythin Carotenoids.
!
              IF (iPigs(iseco,7).gt.0) THEN
                 tmp3 = tmp00+iPigs(iseco,7)
                 tr_nd(tmp3,k,i)=tr_nd(tmp2,k,i)*         &
     &                                  (b_HPUb(iseco)+                 &
     &                                   mxHPUb(iseco)*cff6)
              END IF
            END DO
!
! DOC initialization.
!
            cff6=MAX(0.001_r8,-0.9833_r8*salt+33.411_r8)
            tmp1 = tmp00+iDOMC(1)
            tr_nd(tmp1,k,i)=0.1_r8
            tmp2 = tmp00+iDOMN(1)
            tr_nd(tmp2,k,i)=tr_nd(tmp1,k,i)*cff8
            tmp3 = tmp00+iDOMP(1)
            tr_nd(tmp3,k,i)=tr_nd(tmp2,k,i)*cff9
            IF(CDOC==1) THEN
              tmp2 = tmp00+iCDMC(1)
              tr_nd(tmp2,k,i)=tr_nd(tmp1,k,i)*cDOCfrac_c(1) !MFR
            ENDIF
            IF(Ndom==2)THEN  !MFR
              tmp1 = tmp00+iDOMC(2)
              tr_nd(tmp1,k,i)=15.254_r8*cff6+70.0_r8
              tmp2 = tmp00+iDOMN(2)
              tr_nd(tmp2,k,i)=tr_nd(tmp1,k,i)*cff10
              tmp1 = tmp00+iDOMP(2)
              tr_nd(tmp1,k,i)=0.0_r8
              IF(CDOC==1) THEN
                tmp1 = tmp00+iCDMC(2)
                tr_nd(tmp1,k,i)=(0.243_r8*cff6+0.055_r8)*cff7   !MFR
              ENDIF
            ENDIF 
!
! Fecal Initialization.
!
            DO iseco=1,Nfec
              tmp1 = tmp00 + iFecC(iseco)
              tr_nd(tmp1,k,i)=0.002_r8
              tmp2 = tmp00 + iFecN(iseco)
              tr_nd(tmp2,k,i)=tr_nd(tmp1,k,i)*cff11
              tmp2 = tmp00 + iFecP(iseco)  
              tr_nd(tmp2,k,i)=tr_nd(tmp1,k,i)*cff12
              IF(IRON==1) THEN
                tmp2 = tmp00 + iFecF(iseco)
                tr_nd(tmp2,k,i)=tr_nd(tmp1,k,i)*cff13 !MFR
              ENDIF
              tmp2 = tmp00 + iFecS(iseco)
              tr_nd(tmp2,k,i)=tr_nd(tmp1,k,i)*cff11
            END DO
        END DO
      END DO          
 
      RETURN
      END SUBROUTINE bio_init


      
      
      
