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

    
          SUBROUTINE ecosim(wtratio)
!
!!======================================================================
!! April/May, 2007 - Original code                                     ! 
!!======================================================Marta Rodrigues=
!!                                                                     !
!! This subroutine is from ROMS                                        !                 	
!! The model was extend for zooplankton                                !
!! February, 2009 - Extended for oxygen cycle                          ! 
!!======================================================================
!************************************************** Hernan G. Arango ***
!  Copyright (c) 2002-2007 The ROMS/TOMS Group                         !
!    Licensed under a MIT/X style license                              !
!    See License_ROMS.txt                                              !
!*************************************************** W. Paul Bissett ***
!  Copyright (c) 1997 W. Paul Bissett, FERI                            !
!***********************************************************************
!                                                                      !
!  The EcoSim code has been developed for research purposes only. It   !
!  consists of unpublished, proprietary formulations protected under   !
!  U.S. copyright law. It is freely available on request from the      !
!  Florida Environmental Research Institute (FERI). Commercial usage   !
!  of these formulations is forbidden without express written          !
!  permission from FERI. All rights reserved.                          !
!                                                                      !
!************************************************** Hernan G. Arango ***
!                                                                      !
!  This routine computes the EcoSim sources and sinks and adds them    !
!  to the global biological fields.                                    !
!                                                                      !
!  Reference:                                                          !
!                                                                      !
!    Bissett, W.P., J.J. Walsh, D.A. Dieterle, K.L. Carder, 1999:      !
!      Carbon cycling in the upper waters of the Sargasso Sea: I.      !
!      Numerical  simulation of  differential carbon and nitrogen      !
!      fluxes,  Deep-Sea Res., 46, 205-269.                            !
!                                                                      !
!    Bissett, W.P., K.L. Carder, J.J. Walsh, D.A. Dieterle, 1999:      !
!      Carbon cycling in the upper waters of the Sargasso Sea: II.     !
!      Numerical  simulation  of  apparent  and  inherent optical      !
!      properties, Deep-Sea Res., 46, 271-317                          !
!                                                                      !
!  NOTES to EcoSim                                                     !
!  This version uses a descending index for depth that is different    !
!  than the original coding.                                           !
!                                                                      !
!***********************************************************************
!
!
      USE bio_param
      USE biology
      USE schism_glbl, only : nea,nvrt,tr_el,&
                              bdy_frc,flx_sf,flx_bt,&
                              idry_e,kbe,ze,dt,&
                              irange_tr,& !,ntracers2,ntracers
                              windx1,windx2,windy1,windy2,&
                              elnode,i34 
      USE schism_msgp, only : myrank,parallel_abort

      IMPLICIT NONE
      SAVE
!
!  Imported variable declarations.
!
!      integer, intent(in) :: nea, nvrt
      
!      real(r8), intent(in) :: Hz(nvrt,nea)
!      real(r8), intent(in) :: zr(nvrt,nea)
!      real(r8), intent(in) :: z_w(nvrt,nea)      
!      real(r8), intent(in) :: dt

!      real(r8), intent(in) :: specir(nea,NBands)
!      real(r8), intent(in) :: avcos(nea,NBands)
     
!
!  Local variable declarations.
!
      integer :: Iter, Tindex, ibio, id, itrc, i, j, k
      integer :: ibac, iband, idom, ifec, iphy, ipig, ntr_local
!Marta Rodrigues
      integer :: izoo  
      
      real(r8), parameter :: MinVal = 0.0_r8

      real(r8) :: FV1, FV2, FV3, FV4, FV5, FV6, FV7, dtbio
      real(r8) :: DOC_lab, Ed_tot, Nup_max, aph442, aPHYN_wa
      real(r8) :: avgcos_min, par_b, par_bb, photo_DIC, photo_DOC
      real(r8) :: photo_decay, slope_AC, tChl, theta_m, total_photo
      real(r8) :: tot_ab, tot_b, tot_bb

      real(r8) :: Het_BAC
      real(r8) :: N_quota, RelDOC1, RelDON1, RelDOP1, RelFe

! Marta Rodrigues
      real(r8), dimension(nvrt) :: Hz, zr

      real(r8), dimension(4) :: Bac_G

      real(r8), dimension(NBands) :: dATT_sum

      real(r8), dimension(nvrt,NBands) :: specir_d
      real(r8), dimension(nvrt,NBands) :: avgcos, dATT

      real(r8), dimension(nvrt,Nphy) :: C2CHL, C2CHL_w
      real(r8), dimension(nvrt,Nphy) :: Gt_fl, Gt_ll, Gt_nl
      real(r8), dimension(nvrt,Nphy) :: Gt_sl, Gt_pl 
      real(r8), dimension(nvrt,Nphy) :: alfa
      real(r8), dimension(nvrt,Nphy) :: pac_eff

      real(r8), dimension(nvrt,Nphy,Npig) :: Pigs_w

      integer :: Keuphotic

      real(r8), dimension(nvrt) :: E0_nz
      real(r8), dimension(nvrt) :: Ed_nz
      real(r8), dimension(nvrt) :: DOC_frac
      real(r8), dimension(nvrt) :: Hz_inv
      real(r8), dimension(nvrt) :: NitrBAC
      real(r8), dimension(nvrt) :: NH4toNO3
      real(r8), dimension(nvrt) :: NtoNBAC
      real(r8), dimension(nvrt) :: NtoPBAC
      real(r8), dimension(nvrt) :: NtoFeBAC
      real(r8), dimension(nvrt) :: totDOC_d
      real(r8), dimension(nvrt) :: totDON_d
      real(r8), dimension(nvrt) :: totDOP_d
      real(r8), dimension(nvrt) :: totFe_d
      real(r8), dimension(nvrt) :: totNH4_d
      real(r8), dimension(nvrt) :: totNO3_d
      real(r8), dimension(nvrt) :: totPO4_d
      real(r8), dimension(nvrt) :: totSiO_d

      real(r8), dimension(nvrt,Nbac) :: GtBAC
      real(r8), dimension(nvrt,Nbac) :: NupDOC_ba
      real(r8), dimension(nvrt,Nbac) :: NupDON_ba
      real(r8), dimension(nvrt,Nbac) :: NupDOP_ba
      real(r8), dimension(nvrt,Nbac) :: NupFe_ba
      real(r8), dimension(nvrt,Nbac) :: NupNH4_ba
      real(r8), dimension(nvrt,Nbac) :: NupPO4_ba

      real(r8), dimension(nvrt,Nphy) :: C2fALG
      real(r8), dimension(nvrt,Nphy) :: C2nALG
      real(r8), dimension(nvrt,Nphy) :: C2pALG
      real(r8), dimension(nvrt,Nphy) :: C2sALG
      real(r8), dimension(nvrt,Nphy) :: GtALG
      real(r8), dimension(nvrt,Nphy) :: GtALG_r
      real(r8), dimension(nvrt,Nphy) :: NupDOP
      real(r8), dimension(nvrt,Nphy) :: NupDON
      real(r8), dimension(nvrt,Nphy) :: NupFe
      real(r8), dimension(nvrt,Nphy) :: NupNH4
      real(r8), dimension(nvrt,Nphy) :: NupNO3
      real(r8), dimension(nvrt,Nphy) :: NupPO4
      real(r8), dimension(nvrt,Nphy) :: NupSiO
      real(r8), dimension(nvrt,Nphy) :: graz_act
      real(r8), dimension(nvrt,Nphy) :: mu_bar_f
      real(r8), dimension(nvrt,Nphy) :: mu_bar_n
      real(r8), dimension(nvrt,Nphy) :: mu_bar_p
      real(r8), dimension(nvrt,Nphy) :: mu_bar_s
      real(r8), dimension(nvrt,Nphy) :: refuge

      real(r8), dimension(nvrt,Nfec) :: Regen_C
      real(r8), dimension(nvrt,Nfec) :: Regen_F
      real(r8), dimension(nvrt,Nfec) :: Regen_N
      real(r8), dimension(nvrt,Nfec) :: Regen_P
      real(r8), dimension(nvrt,Nfec) :: Regen_S

      real(r8), dimension(nvrt,NBands) :: specir_scal
      real(r8), dimension(nvrt,Nphy,NBands) :: aPHYN_al
      real(r8), dimension(nvrt,Nphy,NBands) :: aPHYN_at

      real(r8), dimension(nvrt,NBIT2) :: Bio
      real(r8), dimension(nvrt,NBIT) :: Bio_old
      real(r8), dimension(nvrt,NBIT) :: Bio_new
!      real(r8),save,allocatable :: Bio(:,:),Bio_old(:,:),Bio_new(:,:)

! Marta Rodrigues
      real(r8), dimension(nvrt,Nzoo) :: GtZOO
      real(r8), dimension(nvrt,Nzoo) :: Gt_z
      real(r8), dimension(nvrt,Nzoo) :: foodC
      real(r8), dimension(nvrt,Nzoo) :: foodN
      real(r8), dimension(nvrt,Nzoo) :: foodP
      real(r8), dimension(nvrt,Nzoo,Nphy) :: Gt_zoo_phy
      
      real(r8), dimension(Nphy) :: respPhy
      real(r8), dimension(Nbac) :: respBac
      real(r8), dimension(Nbac) :: fBac
      real(r8), dimension(Nzoo) :: respZoo
      real(r8), dimension(nvrt) :: Denit
      real(r8), dimension(nea)  :: Kreaer,K20,wind_speed,OW,temp_scale
      real(r8), dimension(nea)  :: DOW,DOsat,solub,solubCO2
      real(r8), dimension(nea)  :: reaerDO,reaerCO2
      real(r8), dimension(nea)  :: ScO2,ScCO2,pCO2w
      real(r8), dimension(nea)  :: kdiss1,kdiss2,CO2star
      real(r8), dimension(nea)  :: temp,salt
      real(r8) :: reox,Hconc,eqstate

      real(r8) :: u_wind1,u_wind2,v_wind1,v_wind2

      !real(r8), intent(in) :: uwind(nea)
      !real(r8), intent(in) :: vwind(nea)
      real(r8), intent(in) :: wtratio 

!
!=======================================================================
!  Add EcoSim Source/Sink terms.
!=======================================================================
!
      !! MFR - not used   J_LOOP : DO j=Jstr,Jend
! For each element
      DO i=1,nea 
        IF(idry_e(i)==1) CYCLE      

! Marta Rodrigues
! ----------------------------------------------------------------------
! Compute thickness and actual depths at RHO-points (in the middle of
! the volume of control)
! ----------------------------------------------------------------------
                      
!        DO i=1, nea
!          if(idry_e(i)==1) cycle

          DO k=kbe(i)+1, nvrt
              Hz(k)=ze(k,i)-ze(k-1,i)
          END DO
!        END DO ! element

!        DO i=1, nea
!          if(idry_e(i)==1) cycle
          DO k=kbe(i)+1, nvrt
            zr(k)=((ze(k-1,i)-ze(k,i))/2)+ze(k,i)
          END DO
!        END DO ! element
!
!-----------------------------------------------------------------------
!  Compute inverse thickness to avoid repeated divisions.
!-----------------------------------------------------------------------
!
!        DO i=1,nea
!          if(idry_e(i)==1) cycle
          DO k=kbe(i)+1,nvrt 
            Hz_inv(k)=1.0_r8/Hz(k)
          END DO
!        END DO

!
!-----------------------------------------------------------------------
!  Extract biological variables from tracer arrays, place them into
!  scratch arrays, and restrict their values to be positive definite.
!-----------------------------------------------------------------------
!

! MFR - Feb/2008

        Bio=0.0_r8
        Bio_old=0.0_r8
        Bio_new=0.0_r8
        
        DO ibio=1,NBIT
          itrc=idbio(ibio)
!          DO i=1,nea
!            if(idry_e(i)==1) cycle
            DO k=kbe(i)+1,nvrt               
!MFR              Bio(k,itrc)=MAX(MinVal,tr_el(itrc,k,i)*Hz_inv(k))
              Bio(k,itrc)=MAX(MinVal,tr_el(itrc+irange_tr(1,6)-1,k,i))
              Bio_old(k,itrc)=Bio(k,itrc)
              
!! HGA - The new tendency terms were not initialized.  This gives
!!       unexpected behavior on different computers since a variable
!!       was used before it was assigned.  This may explain earlier
!!       problems with the algorithm.  Perhaps, this time-stepping can
!!       be modified latter to avoid unnecessary storage between
!!       Bio, Bio_old, and Bio_new.
!!
              Bio_new(k,itrc)=0.0_r8
            END DO
!          END DO
        END DO
!      
!  Extract temperature and salinity - MFR.
!
!        DO i=1,nea
!          if(idry_e(i)==1) cycle
          DO k=kbe(i)+1,nvrt
            Bio(k,itemp)=tr_el(1,k,i) !tsel(1,k,i) !MFR *Hz_inv(k,i)
            Bio(k,isalt)=tr_el(2,k,i) !tsel(2,k,i) !MFR *Hz_inv(k,i)
          END DO
!        END DO
!
!-----------------------------------------------------------------------
!  Compute temperature and salinity dependent variables.
!-----------------------------------------------------------------------
!
!  Refuge depth calculation.
!
        DO iphy=1,Nphy
!          DO i=1,nea
!            if(idry_e(i)==1) cycle
            DO k=kbe(i)+1,nvrt
              refuge(k,iphy)=MinRefuge(iphy)
            END DO
!          END DO
        END DO
!
!  Initialize fecal regeneration arrays (N, P, and Fe from Moore et al.,
!  DSRII 2001; silica is given by values from Bidle and Azam, Nature,
!  1999).
!
        IF (Regen_flag) THEN
          DO ifec=1,Nfec
!            DO i=1,nea
!              if(idry_e(i)==1) cycle
              DO k=kbe(i)+1,nvrt
                FV1=EXP(RegTfac(ifec)*(Bio(k,itemp)-               &
     &                  RegTbase(ifec)))
                Regen_C(k,ifec)=RegCR(ifec)*FV1
                Regen_N(k,ifec)=RegNR(ifec)*FV1
                Regen_P(k,ifec)=RegPR(ifec)*FV1
                Regen_F(k,ifec)=RegFR(ifec)*FV1
                Regen_S(k,ifec)=RegSR(ifec)*FV1
              END DO
!            END DO
          END DO
        END IF
!
!  Calculate temperature dependent growth rate.
!
        DO iphy=1,Nphy
!          DO i=1,nea
!            if(idry_e(i)==1) cycle
            DO k=kbe(i)+1,nvrt
              GtALG(k,iphy)=GtALG_max(iphy)*                       &
     &                        EXP(PhyTfac(iphy)*                     &
     &                            (Bio(k,itemp)-PhyTbase(iphy)))

!
!  Calculate mu_bar for droop equation.
!
              FV1=maxC2nALG(iphy)*(1.0_r8+GtALG(k,iphy))
              mu_bar_n(k,iphy)=GtALG(k,iphy)*                       &
     &                           FV1/(FV1-minC2nALG(iphy))
              IF (HsSiO(iphy).lt.LARGE) THEN
                FV1=maxC2SiALG(iphy)*(1.0_r8+GtALG(k,iphy))
                mu_bar_s(k,iphy)=GtALG(k,iphy)*                     &
     &                             FV1/(FV1-minC2SiALG(iphy))
              ELSE
                mu_bar_s(k,iphy)=LARGE
              END IF
              IF (HsPO4(iphy).lt.LARGE) THEN
                FV1=maxC2pALG(iphy)*(1.0_r8+GtALG(k,iphy))
                mu_bar_p(k,iphy)=GtALG(k,iphy)*                     &
     &                             FV1/(FV1-minC2pALG(iphy))
              ELSE
                mu_bar_p(k,iphy)=LARGE
              END IF
              IF (HsFe(iphy).lt.LARGE) THEN
                FV1=maxC2FeALG(iphy)*(1.0_r8+GtALG(k,iphy))
                mu_bar_f(k,iphy)=GtALG(k,iphy)*                     &
     &                             FV1/(FV1-minC2FeALG(iphy))
              ELSE
                mu_bar_f(k,iphy)=LARGE
              END IF
            END DO
!          END DO
        END DO

!
!  Bacterial growth rate from Fasham et al., 1990.
!
        DO ibac=1,Nbac
!          DO i=1,nea
!            if(idry_e(i)==1) cycle
            DO k=kbe(i)+1,nvrt            
              GtBAC(k,ibac)=GtBAC_max(ibac)*                       &
     &                        EXP(BacTfac(ibac)*                     &
     &                            (Bio(k,itemp)-BacTbase(ibac)))
            END DO
!          END DO
        END DO
!
!  Grazing rate calculation.
!  NOTE: ES1 included separation calculations for grazing beneath the
!        zone of refuge (250 m). This has been removed and may
!        result in differences in deeper waters.
!! Revisions, WPB 10/20/02.  New grazing formulation that is better
!! representation of basal loss rates and biomass accumulations.
!

! Marta Rodrigues
! This term has been changed to account for the new variable zooplankton
! Now HsGRZ represents only the natural mortality (lises) 
! of the the phytoplankton
        
!	DO iphy=1,Nphy
!          DO k=kbe(i)+1,N(ng)
!            DO i=Istr,Iend
!              FV1=MAX(1.0_r8,(Bio(k,i,iPhyC(iphy))/refuge(k,i,iphy)))
!              graz_act(k,i,iphy)=HsGRZ(iphy,ng)*LOG(FV1)
!            END DO
!          END DO
!        END DO
! 

! Marta Rodrigues
! Zooplankton growth rate
        
        DO izoo=1, Nzoo
!          DO i=1,nea
!            if(idry_e(i)==1) cycle
            DO k=kbe(i)+1, nvrt
              IF(zoo_sp(izoo)==1)THEN !Acartia tonsa, naupilus
	        GtZOO(k,izoo)=0.0914_r8*EXP(0.0701_r8*Bio(k,itemp))&
                                *sec2day        
              ELSE IF(zoo_sp(izoo)==2)THEN !Acartia tonsa, copepodites
	        GtZOO(k,izoo)=0.0517_r8*EXP(0.130_r8*Bio(k,itemp))&
                                *sec2day
	      ELSE IF(zoo_sp(izoo)==3)THEN !Acartia clausi, naupilus
	        GtZOO(k,izoo)=0.0591_r8*EXP(0.0775_r8*Bio(k,itemp))&
                                *sec2day 
	      ELSE IF(zoo_sp(izoo)==4)THEN !Acartia clausi, copepodites
	        GtZOO(k,izoo)=0.0435_r8*EXP(0.114_r8*Bio(k,itemp))&
                                *sec2day
	      ELSE
	         call parallel_abort('Bad flag for zooplankton maximum growth rate calculation!')
!'
	      END IF !zoo_sp
	    END DO
!	  END DO
	END DO      
	    
    
!-----------------------------------------------------------------------
!  Iterate biology source and sink terms.
!-----------------------------------------------------------------------
!
!MFR,Apr/2017
!        dtbio=dt/REAL(BioIter,r8)
         dtbio=dt 

!        ITER_LOOP : DO Iter=1,BioIter
 
!          DO i=1,nea
!            if(idry_e(i)==1) cycle
            DO k=kbe(i)+1,nvrt                       
              totNH4_d(k)=0.0_r8
              totNO3_d(k)=0.0_r8
              totPO4_d(k)=0.0_r8
              totSiO_d(k)=0.0_r8
              totFe_d (k)=0.0_r8
              totDOC_d(k)=0.0_r8
              totDON_d(k)=0.0_r8
              totDOP_d(k)=0.0_r8
            END DO
!          END DO
          DO iphy=1,Nphy
!            DO i=1,nea
!              if(idry_e(i)==1) cycle
              DO k=kbe(i)+1,nvrt              
                NupNH4(k,iphy)=0.0_r8
                NupNO3(k,iphy)=0.0_r8
                NupPO4(k,iphy)=0.0_r8
                NupSiO(k,iphy)=0.0_r8
                NupFe (k,iphy)=0.0_r8
                NupDON(k,iphy)=0.0_r8
                NupDOP(k,iphy)=0.0_r8
              END DO
!            END DO
          END DO
!
!  Compute Ratio Arrays.
!  (Calculating only those that are accessed more than once.)
!
          DO iphy=1,Nphy
!            DO i=1,nea
!              if(idry_e(i)==1) cycle
              DO k=kbe(i)+1,nvrt
                C2nALG(k,iphy)=0.0_r8
                IF (Bio(k,iPhyN(iphy)).gt.0.0_r8) THEN
                  C2nALG(k,iphy)=Bio(k,iPhyC(iphy))/                &
     &                             Bio(k,iPhyN(iphy))
                END IF
                C2pALG(k,iphy)=0.0_r8
                IF (Bio(k,iPhyP(iphy)).gt.0.0_r8) THEN
                  C2pALG(k,iphy)=Bio(k,iPhyC(iphy))/                &
     &                             Bio(k,iPhyP(iphy))
                END IF 
                C2sALG(k,iphy)=0.0_r8
                IF (iPhyS(iphy).gt.0) THEN
                  IF (Bio(k,iPhyS(iphy)).gt.0.0_r8) THEN
                    C2sALG(k,iphy)=Bio(k,iPhyC(iphy))/              &
     &                               Bio(k,iPhyS(iphy))
                  END IF
                END IF 

! MFR - Only if Fe is being calculated...

                IF(IRON==1)THEN
                  C2fALG(k,iphy)=0.0_r8
                  IF (Bio(k,iPhyF(iphy)).gt.0.0_r8) THEN
                    C2fALG(k,iphy)=Bio(k,iPhyC(iphy))/                &
     &                               Bio(k,iPhyF(iphy))
                  END IF
                END IF  !IRON 
             
              END DO
!            END DO
          END DO
!


!  Marta Rodrigues
!-----------------------------------------------------------------------
!  Compute zooplankton growth rate
!-----------------------------------------------------------------------	
		    
!  Compute food available for zooplankton
         
          DO izoo=1,Nzoo
!            DO i=1,nea
!              if(idry_e(i)==1) cycle
              DO k=kbe(i)+1, nvrt            
                foodC(k,izoo)=0.0_r8
	        foodN(k,izoo)=0.0_r8
	        foodP(k,izoo)=0.0_r8
	        Gt_z(k,izoo)=0.0_r8
	      END DO
!	    END DO
	  END DO
   
          DO izoo=1,Nzoo
!            DO i=1,nea
!              if(idry_e(i)==1) cycle
              DO k=kbe(i)+1, nvrt                
                DO iphy=1,Nphy
                  foodC(k,izoo)=foodC(k,izoo)+DeltaZoo(izoo,iphy)*EfcCap(izoo,iphy) &
		                                  *Bio(k,iphyC(iphy))
                  foodN(k,izoo)=foodN(k,izoo)+DeltaZoo(izoo,iphy)*EfcCap(izoo,iphy) &
		                                  *Bio(k,iphyN(iphy))
                  foodP(k,izoo)=foodP(k,izoo)+DeltaZoo(izoo,iphy)*EfcCap(izoo,iphy) &
		                                  *Bio(k,iphyP(iphy))
                END DO
              END DO  
!            END DO
          END DO  
       
!  Growth rate
          DO izoo=1,Nzoo
!            DO i=1,nea
!              if(idry_e(i)==1) cycle
              DO k=kbe(i)+1,nvrt
                DO iphy=1,Nphy
                  IF(foodC(k,izoo)==0.d0)THEN
                     FV1=0.d0
                     FV2=0.d0
                  ELSE 
		     FV1=GtZOO(k,izoo)*DeltaZoo(izoo,iphy)*EfcCap(izoo,iphy) &
		         *Bio(k,iphyC(iphy))/foodC(k,izoo)
		     FV2=foodC(k,izoo)/(foodC(k,izoo)+HsZoo(izoo))
                  ENDIF
                  Gt_zoo_phy(k,izoo,iphy)=FV1*FV2
                  Gt_z(k,izoo)=Gt_z(k,izoo)+FV1*FV2
                END DO
              END DO
!            END DO   
          END DO

!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!  Daylight Computations.
!-----------------------------------------------------------------------
!
! Initialize.
!
!          DO i=1,nea
!            if(idry_e(i)==1) cycle
            Ed_nz(nvrt)=0.0_r8
            E0_nz(nvrt)=0.0_r8
            Keuphotic=nvrt+1
            IF (specir(i,21).gt.VSMALL) THEN
!MFR - SELFE 2D/3D
              IF((kbe(i)+1).ne.nvrt)THEN
                DO k=kbe(i)+1,nvrt-1
                  Ed_nz(k)=0.0_r8
                  E0_nz(k)=0.0_r8
                END DO
              ENDIF
!MFR
              DO iband=1,NBands
                dATT_sum(iband)=0.0_r8
                DO k=kbe(i)+1,nvrt
                  dATT(k,iband)=0.0_r8
                END DO
                DO iphy=1,Nphy
                  DO k=kbe(i)+1,nvrt
                    aPHYN_at(k,iphy,iband)=0.0_r8
                    aPHYN_al(k,iphy,iband)=0.0_r8
                  END DO
                END DO
              END DO
!
!  Calculate average cosine zenith angle at surface.
!  (See equation 14 Morel, 1991 Prog. Ocean.)
!
              Ed_tot=0.0_r8
              DO iband=1,NBands
                Ed_tot=Ed_tot+specir(i,iband)*DLAM
                avgcos(nvrt,iband)=avcos(i,iband)
              END DO
!
!  Total aph(442). adp(442) is set to 50% of aph(442).
!  NOTE: choosing sbands=9 which is band 442 using v8r16
!        sbands formulation. If spectral resolution changes, this
!        value must change!
!
              DO k=nvrt,kbe(i)+1,-1
                IF (Ed_tot.ge.1.0_r8) THEN
                  aph442=0.0_r8
                  tChl=0.0_r8
                  DO iphy=1,Nphy
                    IF (Bio(k,iPhyC(iphy)).gt.0.0_r8) THEN
                      tChl=tChl+Bio(k,iPigs(iphy,ichl))
                      pac_eff(k,iphy)=1.0_r8
                      IF (b_PacEff(iphy).gt.SMALL) THEN
                        FV2=Bio(k,iPigs(iphy,ichl))/                  &
     &                      (Bio(k,iPhyC(iphy))*12.0_r8)
                        pac_eff(k,iphy)=MAX(0.5_r8,                     &
     &                                      (MIN(1.0_r8,                &
     &                                           b_PacEff(iphy)+     &
     &                                           mxPacEff(iphy)*     &
     &                                           (FV2-                  &
     &                                            b_C2Cl(iphy)))))
                      END IF
                      iband=9
                      DO ipig=1,Npig
                        IF (iPigs(iphy,ipig).gt.0) THEN
                          aph442=aph442+                                &
     &                           Bio(k,iPigs(iphy,ipig))*             &
     &                           apigs(ipig,iband)*pac_eff(k,iphy)
                        END IF
                      END DO
                    END IF
                  END DO
!
!  Calculate absorption.
!  Calculating phytoplankton absorption for attentuation calculation.
!  NOTE: 12 factor to convert to ugrams (mg m-3)
!
                  aph442=0.5_r8*aph442
                  DO iband=1,NBands
                    tot_ab=0.0_r8
                    DO iphy=1,Nphy
                      DO ipig=1,Npig 
                        IF (iPigs(iphy,ipig).gt.0) THEN
                          aPHYN_at(k,iphy,iband)=                     &
     &                                      aPHYN_at(k,iphy,iband)+   &
     &                                      Bio(k,iPigs(iphy,ipig))*  &
     &                                      apigs(ipig,iband)*          &
     &                                      pac_eff(k,iphy)
                        END IF
                      END DO
                      tot_ab=tot_ab+aPHYN_at(k,iphy,iband)
!
!  Removing absorption due to PPC for "alfa" calculation.
!
                      ipig=5
                      IF (iPigs(iphy,ipig).gt.0) THEN
                        aPHYN_al(k,iphy,iband)=                       &
     &                                    aPHYN_at(k,iphy,iband)-     &
     &                                    Bio(k,iPigs(iphy,ipig))*    &
     &                                    apigs(ipig,iband)*            &
     &                                    pac_eff(k,iphy)

! MFR - If only Clf a is calculated
                      ELSE IF (iPigs(iphy,ipig).eq.0) THEN
                         aPHYN_al(k,iphy,iband)=                       &
     &                                    aPHYN_at(k,iphy,iband) 
                      END IF
                    END DO
!
!  Adding detrital absorption.
!
                    tot_ab=tot_ab+                                      &
     &                     aph442*EXP(0.011_r8*(442.0_r8-               &
     &                                (397.0_r8+REAL(iband,r8)*DLAM)))
!
! MFR - Only if CDOC is being calculated

                    IF(CDOC==1)THEN 
!
!  Calculate CDOC absorption.
!  NOTE: 12 factor is to convert ugrams per liter, and 0.001 converts
!        to mg/liter.  Specific absorption
!        coefficients were calculated as m-1 / (mg DOC/liters sw).
!        net factor = (12*0.001) = 0.012
!
                     tot_ab=tot_ab+ &
     &                      0.012_r8*(Bio(k,iCDMC(ilab))* &
     &                                aDOC(ilab,iband)+ &
     &                                Bio(k,iCDMC(irct))* &
                                      aDOC(irct,iband))+awater(iband)
!     &                      awater(iband)

                     END IF  !CDOC=1
!
!  Calculate scattering and backscattering (see equation 19 Morel, 1991,
!  Prog. Ocean). Morel, 1988 puts spectral dependency in backscattering.
!  Since Morel (1991) does not have a backscattering equation, use 1988
!  paper. Morel 2001 has slight adjustment 0.01, rather than 0.02. 
!  This was altered, but never tested in ROMS 1.8 on 03/08/03.
!
                    par_b =0.3_r8*(tChl**0.62_r8)
                    par_bb=0.0_r8
                    IF (tChl.gt.0.0_r8) THEN
                      par_bb=par_b*(0.002_r8+0.01_r8*                   &
     &                              (0.5_r8-0.25_r8*LOG10(tChl))*       &
     &                              wavedp(iband))
                    END IF
                    par_bb=MAX(par_bb,0.0_r8)
!
!  However, for omega0 calculation, par_b must be spectral, so use
!  dependency from Sathy and Platt 1988
!
                    tot_b=bwater(iband)+par_b*wavedp(iband)
!
!  Morel, 1988 instead of 1991. See methods
!
                    tot_bb=0.5_r8*bwater(iband)+par_bb
!
!  Sathy and Platt JGR 1988.  This is set with the average cosine of
!  the box above, and used to calculate a new avgcos for this level.
!  This new average cosine is then used to recalculate the attenuation
!  coefficient
!
                    dATT(k,iband)=(tot_ab+tot_bb)/avgcos(k,iband)
!
!  See Mobley, 1995 for graphical depiction of this equation.
!
                    avgcos_min=avgcos(k,iband)+                         &
     &                         (0.5_r8-avgcos(k,iband))*                &
     &                         (tot_b/(tot_ab+tot_b))
!
!  Calculate average cosine. Linear fit to average cosine versus optical
!  depth relationship. The FV1 calculation keeps the denominator of the
!  slope calculation from going negative and above 1.
!
                    FV1=MAX(1.0_r8,                                     &
     &                      7.0_r8-dATT(k,iband)*ABS(zr(k)))
                    slope_AC =MIN(0.0_r8,                               &
     &                            (avgcos_min-avgcos(k,iband))/FV1)
                    avgcos(k,iband)=avgcos(k,iband)+                    &
     &                             slope_AC*dATT(k,iband)*Hz(k)
                    dATT(k,iband)=(tot_ab+tot_bb)/avgcos(k,iband)
!
!  Set avgcos for next level.
!
                    IF (k.ne.(kbe(i)+1)) THEN
                      avgcos(k-1,iband)=avgcos(k,iband)
                    END IF
!
!  Calculate spectral irradiance with depth.
!
                    FV1=dATT(k,iband)*Hz(k)
                    FV2=dATT_sum(iband)+0.5_r8*FV1
                    dATT_sum(iband)=dATT_sum(iband)+FV1
                    specir_d(k,iband)=specir(i,iband)*                &
     &                                EXP(-FV2)*DLAM
!
!  Calculate spectral scalar irradiance.  Morel, 1991 Prog. Ocean.
!
                    specir_scal(k,iband)=specir_d(k,iband)*           &
     &                                     (dATT(k,iband)/tot_ab)
                    E0_nz(k)=E0_nz(k)+specir_scal(k,iband)
!
!  Calculate Ed_nz.
!
                    Ed_nz(k)=Ed_nz(k)+specir_d(k,iband)
                  END DO
                  Ed_tot=E0_nz(k)
!
!  Set bottom of the euphotic zone.
!
                  Keuphotic=k
                END IF
              END DO
            END IF
!          END DO
!
!-----------------------------------------------------------------------
!  Bacterial nutrient uptake.
!-----------------------------------------------------------------------
!
          DO ibac=1,Nbac
!            DO i=1,nea
!              if(idry_e(i)==1) cycle
              DO k=kbe(i)+1,nvrt              
!
!  DOM uptake.
!
                IF ((Bio(k,iDOMC(ilab)).gt.0.0_r8).and.               &
     &              (Bio(k,iDOMN(ilab)).gt.0.0_r8).and.               &
     &              (Bio(k,iDOMP(ilab)).gt.0.0_r8)) THEN
                  NupDOC_ba(k,ibac)=GtBAC(k,ibac)*                  &
     &                                Bio(k,iBacC(ibac))*             &
     &                                I_Bac_Ceff*                       &
     &                                (Bio(k,iDOMC(ilab))/            &
     &                                (HsDOC_ba(ibac)+                  &
     &                                 Bio(k,iDOMC(ilab))))
                  NupDON_ba(k,ibac)=NupDOC_ba(k,ibac)*              &
     &                                Bio(k,iDOMN(ilab))/             &
     &                                Bio(k,iDOMC(ilab))
                  NupDOP_ba(k,ibac)=NupDOC_ba(k,ibac)*              &
     &                                Bio(k,iDOMP(ilab))/             &
     &                                Bio(k,iDOMC(ilab))
                ELSE
                  NupDOC_ba(k,ibac)=0.0_r8
                  NupDON_ba(k,ibac)=0.0_r8
                  NupDOP_ba(k,ibac)=0.0_r8
                END IF
                totDOC_d(k)=totDOC_d(k)+NupDOC_ba(k,ibac)
                totDON_d(k)=totDON_d(k)+NupDON_ba(k,ibac)
                totDOP_d(k)=totDOP_d(k)+NupDOP_ba(k,ibac)
!
!  NH4 uptake.
!
                NupNH4_ba(k,ibac)=GtBAC(k,ibac)*                    &
     &                              Bio(k,iBacN(ibac))*               &
     &                              Bio(k,iNH4_)/                     &
     &                              (HsNH4_ba(ibac)+Bio(k,iNH4_))
                totNH4_d(k)=totNH4_d(k)+NupNH4_ba(k,ibac)
!
!  PO4 uptake.
!
                NupPO4_ba(k,ibac)=GtBAC(k,ibac)*                    &
     &                              Bio(k,iBacP(ibac))*               &
     &                              Bio(k,iPO4_)/                     &
     &                              (HsPO4_ba(ibac)+Bio(k,iPO4_))
                totPO4_d(k)=totPO4_d(k)+NupPO4_ba(k,ibac)
!

!  MFR - Only if Fe is being calculated
               
                IF(IRON==1)THEN

!  Fe uptake.
!
                  NupFe_ba(k,ibac)=GtBAC(k,ibac)*                     &
     &                               Bio(k,iBacF(ibac))*                &
     &                               Bio(k,iFeO_)/                      &
     &                               (HsFe_ba(ibac)+Bio(k,iFeO_))
                  totFe_d(k)=totFe_d(k)+NupFe_ba(k,ibac)

                END IF  !IRON
 
              END DO
!            END DO
          END DO
!
!-----------------------------------------------------------------------
!  Phytoplankton dark nutrient uptake.
!-----------------------------------------------------------------------
!
          DO iphy=1,Nphy
!            DO i=1,nea
!              if(idry_e(i)==1) cycle
              DO k=kbe(i)+1,nvrt              
                IF (C2nALG(k,iphy).gt.C2nALGminABS(iphy)) THEN
!
!  NOTE: these are being saved to test for total nutrient uptake.
!        If nutrient uptake is greater than maximum nutrient, then
!        each of the uptakes are reduced by their fractional contri-
!        bution to the total.
!
                  Nup_max=GtALG(k,iphy)
                  NupNO3(k,iphy)=(Bio(k,iNO3_)/                     &
     &                             (HsNO3(iphy)+Bio(k,iNO3_))*        &
     &                             EXP(-BET_(iphy)*Bio(k,iNH4_)))
                  NupNH4(k,iphy)=Bio(k,iNH4_)/                      &
     &                             (HsNH4(iphy)+Bio(k,iNH4_))
!
!  Test that Wroblewski equation does not exceed 1.0.
!
                  FV1=NupNO3(k,iphy)+NupNH4(k,iphy)
                  IF (FV1.gt.1.0_r8) THEN
                    FV1=1.0_r8/FV1
                    NupNO3(k,iphy)=NupNO3(k,iphy)*FV1
                    NupNH4(k,iphy)=NupNH4(k,iphy)*FV1
                  END IF
!
!  Change from percentage of maximum to mass per second.
!
                  FV1=Nup_max*Bio(k,iPhyN(iphy))
                  NupNO3(k,iphy)=NupNO3(k,iphy)*FV1
                  NupNH4(k,iphy)=NupNH4(k,iphy)*FV1
!
!  Test for DON uptake.
!
                  IF (C2nALG(k,iphy).gt.C2nNupDON(iphy)) THEN
                    NupDON(k,iphy)=FV1*                               &
     &                               Bio(k,iDOMN(ilab))/              &
     &                               (HsDON(iphy)+                      &
     &                                Bio(k,iDOMN(ilab)))
                  END IF
!
!  Accumulate total demand for nutrients.
!
                  totNO3_d(k)=totNO3_d(k)+NupNO3(k,iphy)
                  totNH4_d(k)=totNH4_d(k)+NupNH4(k,iphy)
                  totDON_d(k)=totDON_d(k)+NupDON(k,iphy)
                END IF
!
!  Dark silica uptake, min C2Si test.
!  The LARGE test can be removed after testing phase.
!
                IF (HsSiO(iphy).lt.LARGE) THEN
                  IF (C2sALG(k,iphy).gt.C2SiALGminABS(iphy)) THEN
                    Nup_max=GtALG(k,iphy)
                    NupSiO(k,iphy)=Bio(k,iSiO_)/                    &
     &                               (HsSiO(iphy)+Bio(k,iSiO_))
!
!  Change from percentage of maximum to mass per second.
!
                    IF (iPhyS(iphy).gt.0) THEN
                      FV1=Nup_max*Bio(k,iPhyS(iphy))
                      NupSiO(k,iphy)=NupSiO(k,iphy)*FV1
                    ELSE
                      NupSiO(k,iphy)=0.0_r8
                    END IF
!
!  Accumulate total demand for nutrients.
!
                    totSiO_d(k)=totSiO_d(k)+NupSiO(k,iphy)
                  END IF
                END IF
!
!  Dark phophorus uptake, min C2P test.
!  The LARGE test can be removed after testing phase.
!
                IF (HsPO4(iphy).lt.LARGE) THEN
                  IF (C2pALG(k,iphy).gt.C2pALGminABS(iphy)) THEN
                    Nup_max=GtALG(k,iphy)
                    NupPO4(k,iphy)=Bio(k,iPO4_)/                    &
     &                               (HsPO4(iphy)+Bio(k,iPO4_))
!
!  Change from percentage of maximum to mass per second.
!
                    FV1=Nup_max*Bio(k,iPhyP(iphy))
                    NupPO4(k,iphy)=NupPO4(k,iphy)*FV1
!
!  Test for alk. phosphatase
!
                    IF (C2pALG(k,iphy).gt.C2pALKPHOS(iphy)) THEN
                      NupDOP(k,iphy)=FV1*                             &
                                       Bio(k,iDOMP(ilab))/            &
     &                                 (HsDOP(iphy)+                 &
     &                                  Bio(k,iDOMP(ilab)))
                    END IF
!
!  Accumulate total demand for nutrients.
!
                    totPO4_d(k)=totPO4_d(k)+NupPO4(k,iphy)
                    totDOP_d(k)=totDOP_d(k)+NupDOP(k,iphy)
                  END IF
                END IF
!
!  MFR - Only if iron is being calculated...
                 
                IF(IRON==1)THEN

!  Dark iron uptake, min C2Fe test.
!  The LARGE test can be removed after testing phase.
!
                 IF (HsFe(iphy).lt.LARGE) THEN
                   IF (C2fALG(k,iphy).gt.C2FeALGminABS(iphy)) THEN
                     Nup_max=GtALG(k,iphy)
                     NupFe(k,iphy)=Bio(k,iFeO_)/                     &
     &                               (HsFe(iphy)+Bio(k,iFeO_))
!
!  Change from percentage of maximum to mass per second.
!
                     FV1=Nup_max*Bio(k,iPhyF(iphy))
                     NupFe(k,iphy)=NupFe(k,iphy)*FV1
!
!  Accumulate total demand for nutrients.
!
                     totFe_d(k)=totFe_d(k)+NupFe(k,iphy)
                   END IF
                 END IF
               
               END IF  !IRON
              END DO
!            END DO
          END DO
!
!  Calculate bacterial nitrification as a Michaelis-Menton function
!  of ambient NH4 concentration, beneath the euphotic zone (light
!  inhibits nitrification).
!
!          DO i=1,nea
!            if(idry_e(i)==1) cycle
            DO k=kbe(i)+1,nvrt
              NitrBAC(k)=0.0_r8
              NH4toNO3(k)=0.0_r8
              NtoNBAC(k)=0.0_r8
              NtoPBAC(k)=0.0_r8
              NtoFeBAC(k)=0.0_r8
              IF (k.lt.Keuphotic) THEN
	      
!  Marta Rodrigues	      
!  Inclusion of the oxygen - an alternative formulation was added to calculate
!  nitrification

              IF(NIT_flag==1)THEN	      
                NH4toNO3(k)=RtNIT*                                &
     &                        Bio(k,iNH4_)/(HsNIT+Bio(k,iNH4_))
              ELSE     !Vichi et al., 2007
	        NH4toNO3(k)=RtNIT*QN**((Bio(k,itemp)-10.d0)/10.d0)*Bio(k,iNH4_)*&
     &                        Bio(k,iDO_)/(HsNIT+Bio(k,iDO_))
              END IF
!
!  Nitrification fixes DIC into POC.
!  Conversion factor of 7.0 from Kaplan 1983 "Nitrogen in the Sea"
!  factor equals (1.0 / (7.0 * C2nBAC)). Adds NH4 uptake as biomass.
!
                NitrBAC(k)=NH4toNO3(k)/7.0_r8
                NtoNBAC(k)=NitrBAC(k)*N2cBAC
                NtoPBAC(k)=NitrBAC(k)*P2cBAC
                NtoFeBAC(k)=NitrBAC(k)*Fe2cBAC
                totNH4_d(k)=totNH4_d(k)+NH4toNO3(k)+NtoNBAC(k)
                totPO4_d(k)=totPO4_d(k)+NtoPBAC(k)
                totFe_d (k)=totFe_d (k)+NtoFeBAC(k)
              END IF
            END DO
!          END DO
!
!-----------------------------------------------------------------------
!  Test that total nutrient demand does not exceed supply.  If it does
!  total demand is normalized to the total supply. Each species demand
!  is reduced to its weighted average percentage of the supply.
!-----------------------------------------------------------------------
!
!          DO i=1,nea
!            if(idry_e(i)==1) cycle
            DO k=kbe(i)+1,nvrt            
              FV2=totNO3_d(k)*dtbio
              IF (FV2.gt.Bio(k,iNO3_)) THEN
                FV1=(Bio(k,iNO3_)-VSMALL)/FV2
                DO iphy=1,Nphy
                  NupNO3(k,iphy)=NupNO3(k,iphy)*FV1
                END DO
              END IF
!
              FV2=totNH4_d(k)*dtbio
              IF (FV2.gt.Bio(k,iNH4_)) THEN
                FV1=(Bio(k,iNH4_)-VSMALL)/FV2
                DO iphy=1,Nphy
                  NupNH4(k,iphy)=NupNH4(k,iphy)*FV1
                END DO
                DO ibac=1,Nbac
                  NupNH4_ba(k,ibac)=NupNH4_ba(k,ibac)*FV1
                END DO
                NH4toNO3(k)=NH4toNO3(k)*FV1
                NtoNBAC(k)=NtoNBAC(k)*FV1
              END IF
!
              FV2=totSiO_d(k)*dtbio
              IF (FV2.gt.Bio(k,iSiO_)) THEN
                FV1=(Bio(k,iSiO_)-VSMALL)/FV2
                DO iphy=1,Nphy
                  NupSiO(k,iphy)=NupSiO(k,iphy)*FV1
                END DO
              END IF
!
              FV2=totPO4_d(k)*dtbio
              IF (FV2.gt.Bio(k,iPO4_)) THEN
                FV1=(Bio(k,iPO4_)-VSMALL)/FV2
                DO iphy=1,Nphy
                  NupPO4(k,iphy)=NupPO4(k,iphy)*FV1
                END DO
                DO ibac=1,Nbac
                  NupPO4_ba(k,ibac)=NupPO4_ba(k,ibac)*FV1
                END DO
                NtoPBAC(k)=NtoPBAC(k)*FV1
              END IF
!

! MFR - Only if Fe is being calculated

              IF(IRON==1)THEN
               FV2=totFe_d(k)*dtbio
               IF (FV2.gt.Bio(k,iFeO_)) THEN
                 FV1=(Bio(k,iFeO_)-VSMALL)/FV2
                 DO iphy=1,Nphy
                   NupFe(k,iphy)=NupFe(k,iphy)*FV1
                 END DO
                 DO ibac=1,Nbac
                   NupFe_ba(k,ibac)=NupFe_ba(k,ibac)*FV1
                 END DO
                 NtoFeBAC(k)=NtoFeBAC(k)*FV1
               END IF
              END IF  !IRON
!
!  Bacteria are the only group to take up DOC.  Remove BAC DON and
!  BAC DOP uptake from total uptake; adjust uptake and add back.
!
              FV2=totDOC_d(k)*dtbio
              IF (FV2.gt.Bio(k,iDOMC(ilab))) THEN
                FV1=(Bio(k,iDOMC(ilab))-VSMALL)/FV2
                totDOC_d(k)=totDOC_d(k)*FV1
                DO ibac=1,Nbac
                  NupDOC_ba(k,ibac)=NupDOC_ba(k,ibac)*FV1
                  totDON_d(k)=totDON_d(k)-NupDON_ba(k,ibac)
                  NupDON_ba(k,ibac)=NupDON_ba(k,ibac)*FV1
                  totDON_d(k)=totDON_d(k)+NupDON_ba(k,ibac)
                  totDOP_d(k)=totDOP_d(k)-NupDOP_ba(k,ibac)
                  NupDOP_ba(k,ibac)=NupDOP_ba(k,ibac)*FV1
                  totDOP_d(k)=totDOP_d(k)+NupDOP_ba(k,ibac)
                END DO
              END IF
!
!  Remove BAC DON uptake from total uptake; adjust uptake and add back.
!
              FV2=totDON_d(k)*dtbio
              IF (FV2.gt.Bio(k,iDOMN(ilab))) THEN
                FV1=(Bio(k,iDOMN(ilab))-VSMALL)/FV2
                totDON_d(k)=totDON_d(k)*FV1
                totDOC_d(k)=totDOC_d(k)*FV1
                DO iphy=1,Nphy
                  NupDON(k,iphy)=NupDON(k,iphy)*FV1
                END DO
                DO ibac=1,Nbac
                  NupDON_ba(k,ibac)=NupDON_ba(k,ibac)*FV1
                  NupDOC_ba(k,ibac)=NupDOC_ba(k,ibac)*FV1
                  totDOP_d(k)=totDOP_d(k)-NupDOP_ba(k,ibac)
                  NupDOP_ba(k,ibac)=NupDOP_ba(k,ibac)*FV1
                  totDOP_d(k)=totDOP_d(k)+NupDOP_ba(k,ibac)
                END DO
              END IF
!
!  Remove BAC DOP uptake from total uptake; adjust uptake and add back.
!
              FV2=totDOP_d(k)*dtbio
              IF (FV2.gt.Bio(k,iDOMP(ilab))) THEN
                FV1=(Bio(k,iDOMP(ilab))-VSMALL)/FV2
                totDOP_d(k)=totDOP_d(k)*FV1
                totDOC_d(k)=totDOC_d(k)*FV1
                DO iphy=1,Nphy
                  NupDOP(k,iphy)=NupDOP(k,iphy)*FV1
                END DO
                DO ibac=1,Nbac
                  NupDOP_ba(k,ibac)=NupDOP_ba(k,ibac)*FV1
                  totDON_d(k)=totDON_d(k)-NupDON_ba(k,ibac)
                  NupDON_ba(k,ibac)=NupDON_ba(k,ibac)*FV1
                  totDON_d(k)=totDON_d(k)+NupDON_ba(k,ibac)
                  NupDOC_ba(k,ibac)=NupDOC_ba(k,ibac)*FV1
                END DO
              END IF
            END DO
!          END DO
!
!  Increase particulate nutrients by the amount of the uptake.
!
          DO iphy=1,Nphy
!            DO i=1,nea
!              if(idry_e(i)==1) cycle
              DO k=kbe(i)+1,nvrt
                Bio_new(k,iPhyN(iphy))=Bio_new(k,iPhyN(iphy))+      &
     &                                   NupNO3(k,iphy)+              &
     &                                   NupNH4(k,iphy)+              &
     &                                   NupDON(k,iphy)
                Bio_new(k,iPhyP(iphy))=Bio_new(k,iPhyP(iphy))+      &
     &                                   NupPO4(k,iphy)+              &
     &                                   NupDOP(k,iphy)
                IF(IRON==1) Bio_new(k,iPhyF(iphy))=                   &
     &                                   Bio_new(k,iPhyF(iphy))+      &
     &                                   NupFe(k,iphy)
                IF (iPhyS(iphy).gt.0) THEN
                  Bio_new(k,iPhyS(iphy))=Bio_new(k,iPhyS(iphy))+    &
     &                                     NupSiO(k,iphy)
                END IF
!
!  Update nutrient arrays for growth and budgets. Bacterial uptake
!  included below.
!
                Bio_new(k,iNO3_)=Bio_new(k,iNO3_)-                  &
     &                             NupNO3(k,iphy)
                Bio_new(k,iNH4_)=Bio_new(k,iNH4_)-                  &
     &                             NupNH4(k,iphy)
                Bio_new(k,iSiO_)=Bio_new(k,iSiO_)-                  &
     &                             NupSiO(k,iphy)
                Bio_new(k,iPO4_)=Bio_new(k,iPO4_)-                  &
     &                             NupPO4(k,iphy)
                IF(IRON==1) Bio_new(k,iFeO_)=                         &
     &                             Bio_new(k,iFeO_)-                  &
     &                             NupFe (k,iphy)
                Bio_new(k,iDOMN(ilab))=Bio_new(k,iDOMN(ilab))-      &
     &                                   NupDON(k,iphy)
                Bio_new(k,iDOMP(ilab))=Bio_new(k,iDOMP(ilab))-      &
     &                                   NupDOP(k,iphy)
              END DO
!            END DO
          END DO
!
!  Nitrification fixes DIC into DOC.
!
!          DO i=1,nea
!            if(idry_e(i)==1) cycle
            DO k=kbe(i)+1,nvrt           
              Bio_new(k,iDIC_)=Bio_new(k,iDIC_)-                    &
     &                           NitrBAC(k)
            END DO
!          END DO
!
!  Add nitrifying bacteria biomass to heterotrophic bacteria biomass.
!  Adding PON, POP, POFe to BacC arrays at current C2_BAC ratios.
!
          DO ibac=1,Nbac
!            DO i=1,nea
!              if(idry_e(i)==1) cycle
              DO k=kbe(i)+1,nvrt
                Bio_new(k,iBacC(ibac))=Bio_new(k,iBacC(ibac))+      &
     &                                   NitrBAC(k)
                Bio_new(k,iBacN(ibac))=Bio_new(k,iBacN(ibac))+      &
     &                                   NtoNBAC(k)
                Bio_new(k,iBacP(ibac))=Bio_new(k,iBacP(ibac))+      &
     &                                   NtoPBAC(k)
                IF(IRON==1) Bio_new(k,iBacF(ibac))=                   &
     &                                   Bio_new(k,iBacF(ibac))+      &
     &                                   NtoFeBAC(k)
              END DO
!            END DO
          END DO
!
!  Update nutrient arrays for nitrification.
!
!          DO i=1,nea
!            if(idry_e(i)==1) cycle
            DO k=kbe(i)+1,nvrt            
              Bio_new(k,iNO3_)=Bio_new(k,iNO3_)+                    &
     &                           NH4toNO3(k)
              Bio_new(k,iNH4_)=Bio_new(k,iNH4_)-                    &
     &                           (NH4toNO3(k)+NtoNBAC(k))
              Bio_new(k,iPO4_)=Bio_new(k,iPO4_)-                    &
     &                           NtoPBAC(k)
              IF(IRON==1) Bio_new(k,iFeO_)=Bio_new(k,iFeO_)-        &
     &                           NtoFeBAC(k)
            END DO
!          END DO
!
!-----------------------------------------------------------------------
!  Light mediated carbon growth.
!-----------------------------------------------------------------------
!
!         May use rank specific outputs
!          if(myrank==0) WRITE(600,*)'Keuphotic',Keuphotic(2161)
!          DO i=1,nea
!            if(idry_e(i)==1) cycle
            GtALG_r=0.d0
            Gt_fl=0.d0
            Gt_ll=0.d0
            Gt_nl=0.d0
            Gt_sl=0.d0
            Gt_pl=0.d0
            
            DO k=nvrt,Keuphotic,-1
              DO iphy=1,Nphy
                IF (Bio(k,iPhyC(iphy)).gt.0.0_r8) THEN
!
!  Calculate weighted average spectral absorption.
!
                  aPHYN_wa=0.0_r8
                  DO iband=1,NBands
                    aPHYN_wa=aPHYN_wa+(aPHYN_al(k,iphy,iband)*        &
     &                                 specir_scal(k,iband))
                  END DO
!
!  If Keuphotic(i) < N+1, and E0_nz(k,i)=0, this will cause pigments to
!  blow up. This should never happen, unless Keuphotic is not calcuated
!  properly. WPB
!
                  aPHYN_wa=aPHYN_wa/E0_nz(k)
!
!  Calculate "alfa" for HTAN function of P vs. I.
!  (conversion:  Ein/microEin * 10e3)
!
                  alfa(k,iphy)=(aPHYN_wa/Bio(k,iPhyC(iphy)))*         &
     &                          qu_yld(iphy)*0.001_r8
!
!  Light limited growth rate.
!
                  FV1=MAX(0.0_r8,E0_nz(k)-E0_comp(iphy))
                  FV2=E0_nz(k)-E0_inhib(iphy)
                  IF (FV2.gt.0.0_r8) THEN
                    Gt_ll(k,iphy)=GtALG(k,iphy)*                      &
     &                            TANH(alfa(k,iphy)*FV1/                &
     &                            GtALG(k,iphy))*                     &
     &                            EXP(-inhib_fac(iphy)*FV2)
                  ELSE
                    Gt_ll(k,iphy)=GtALG(k,iphy)*                      &
     &                            TANH(alfa(k,iphy)*FV1/                &
     &                            GtALG(k,iphy))
                  END IF
!
!  Nutrient limited growth rates.
!
!  REMEMBER that sinking speed to be set by gradient of limiting
!       nutrient, allowing for negative sinking. Try storing growth
!       rate terms in an array and using MAXLOC for if test.
!
!  Nitrogen limited growth rate.
!
                  IF (Bio(k,iPhyN(iphy)).gt.0.0_r8) THEN
                    FV1=Bio(k,iPhyC(iphy))/                           &
     &                  (Bio(k,iPhyN(iphy))+Bio_new(k,iPhyN(iphy)))
                    Gt_nl(k,iphy)=mu_bar_n(k,iphy)*                   &
     &                            (1.0_r8-ImaxC2nALG(iphy)*FV1)
                    Gt_nl(k,iphy)=MAX(0.0_r8,                           &
     &                                MIN(Gt_nl(k,iphy),                &
     &                                    GtALG(k,iphy)))
                  END IF
!
!  Silica limited growth rate.
!  Testing for silica incorporation.
!
                  IF (iPhyS(iphy).gt.0) THEN
                    IF ((HsSiO(iphy).lt.LARGE).and.                  &
     &                  (Bio(k,iPhyS(iphy)).gt.0.0_r8)) THEN
                      FV1=Bio(k,iPhyC(iphy))/                         &
     &                    (Bio(k,iPhyS(iphy))+                        &
     &                     Bio_new(k,iPhyS(iphy)))
                      Gt_sl(k,iphy)=mu_bar_s(k,iphy)*                 &
     &                              (1.0_r8-ImaxC2SiALG(iphy)*FV1)
                      Gt_sl(k,iphy)=MAX(0.0_r8,                         &
     &                                  MIN(Gt_sl(k,iphy),              &
     &                                      GtALG(k,iphy)))
                    ELSE
                      Gt_sl(k,iphy)=LARGE
                    END IF
                  ELSE
                    Gt_sl(k,iphy)=LARGE
                  END IF
!
!  Phosphorus limited growth rate.
!
                  IF ((HsPO4(iphy).lt.LARGE).and.                    &
     &                (Bio(k,iPhyP(iphy)).gt.0.0_r8)) THEN
                    FV1=Bio(k,iPhyC(iphy))/                           &
     &                  (Bio(k,iPhyP(iphy))+Bio_new(k,iPhyP(iphy)))
                    Gt_pl(k,iphy)=mu_bar_p(k,iphy)*                   &
     &                            (1.0_r8-ImaxC2pALG(iphy)*FV1)
                    Gt_pl(k,iphy)=MAX(0.0_r8,                           &
     &                                MIN(Gt_pl(k,iphy),                &
     &                                    GtALG(k,iphy)))
                  ELSE
                    Gt_pl(k,iphy)=LARGE
                  END IF
!
!  MFR - Only if iron is being calculated...

                  IF(IRON==1)THEN
!
!  Iron limited growth rate
!
                   IF ((HsFe(iphy).lt.LARGE).and.                     &
     &                 (Bio(k,iPhyF(iphy)).gt.0.0_r8)) THEN
                     FV1=Bio(k,iPhyC(iphy))/                           &
     &                   (Bio(k,iPhyF(iphy))+Bio_new(k,iPhyF(iphy)))
                     Gt_fl(k,iphy)=mu_bar_f(k,iphy)*                   &
     &                             (1.0_r8-ImaxC2FeALG(iphy)*FV1)
                     Gt_fl(k,iphy)=MAX(0.0_r8,                           &
     &                                 MIN(Gt_fl(k,iphy),                &
     &                                     GtALG(k,iphy)))
                   ELSE
                     Gt_fl(k,iphy)=LARGE
                   END IF

                  ELSE 
                    Gt_fl(k,iphy)=LARGE
                  END IF  !IRON

!
!  Realized growth rate is minimum of light or nutrient limited rate.
!
                  GtALG_r(k,iphy)=MIN(Gt_ll(k,iphy),Gt_nl(k,iphy),    &
     &                                  Gt_sl(k,iphy),Gt_pl(k,iphy),    &
     &                                  Gt_fl(k,iphy))
                  IF (GtALG_r(k,iphy).ge.LARGE) THEN
                    GtALG_r(k,iphy)=0.0_r8
                  END IF
!
!  Carbon growth calculations.
!
                  FV1=Bio(k,iPhyC(iphy))*GtALG_r(k,iphy)
                  Bio_new(k,iPhyC(iphy))=Bio_new(k,iPhyC(iphy))+    &
     &                                     FV1
                  Bio_new(k,iDIC_)=Bio_new(k,iDIC_)-                &
     &                               FV1
!       
!  Pigment growth calculations.
!
 		  DO ipig=1,Npig
		    IF (iPigs(iphy,ipig).gt.0) THEN
                      itrc=iPigs(iphy,ipig)
                        IF (Bio(k,iPhyC(iphy)).gt.0.0_r8) THEN
                          FV1=Bio(k,itrc)*GtALG_r(k,iphy)
                          Bio_new(k,itrc)=Bio_new(k,itrc)+FV1
                        END IF
		      END IF
		    END DO
                  END IF
                END DO
              END DO
!            END DO
!
!-----------------------------------------------------------------------
!  Bacterioplankton carbon growth terms.
!-----------------------------------------------------------------------
!
!          DO i=1,nea
!            if(idry_e(i)==1) cycle
            DO k=kbe(i)+1,nvrt
              Het_BAC=0.0_r8
              RelDOC1=0.0_r8
              RelDON1=0.0_r8
              RelDOP1=0.0_r8
              RelFe=0.0_r8
!
!  NOTE: Only DOC2/DON2 formation is in this section.
!        Take colored excretion off the top. 03/18/00
!        also, not excreting any DOP or Fe
!  REMEMBER, if excreting DOP and Fe, must address changes in growth if
!        tests. (see DON equations). 03/21/00.
!

              DO ibac=1,Nbac

!  MFR - Only if recalcitrant organic matter is being calculated...
                IF(Ndom==2)THEN

                 FV1=NupDOC_ba(k,ibac)*ExBAC_c*                    &
     &               (1.0_r8-cDOCfrac_c(irct))
                 FV2=NupDOC_ba(k,ibac)*ExBAC_c*                    &
     &               cDOCfrac_c(irct)
                 FV3=NupDON_ba(k,ibac)*ExBAC_n
                           
                 Bio_new(k,iDOMC(irct))=Bio_new(k,iDOMC(irct))+      &
     &                                    FV1
                 Bio_new(k,iCDMC(irct))=Bio_new(k,iCDMC(irct))+      &
     &                                    FV2
                 Bio_new(k,iDOMN(irct))=Bio_new(k,iDOMN(irct))+      &
     &                                    FV3
                
!
!  As we are taking it off the top, must remove from DOMN1 now. No other
!  organisms use DOMC1, so net term (totDOC_d) can be used in budgeting
!  below. This saves cycles, but makes code difficult to read. WPB
!
                 Bio_new(k,iDOMN(ilab))=Bio_new(k,iDOMN(ilab))-      &
     &                                   FV3
!
!  Remove from uptake.
!
                 NupDOC_ba(k,ibac)=NupDOC_ba(k,ibac)-                &
     &                              (FV1+FV2)
                 NupDON_ba(k,ibac)=NupDON_ba(k,ibac)-                &
     &                              FV3

                END IF  !Ndom
!
!  Determine growth limitation. Assuming 100% efficiency for N, P, Fe.
!  If DOMC=0, or DOMN=0, or DOMP=0, then NupDOC_ba = NupDON_ba =
!  NupDOP_ba = 0 and none of the divisions below are accessed. WPB
!
                Bac_G(1)=NupDOC_ba(k,ibac)*Bac_Ceff
                Bac_G(2)=(NupDON_ba(k,ibac)+                          &
     &                    NupNH4_ba(k,ibac))*                         &
     &                   C2nBAC
                Bac_G(3)=(NupDOP_ba(k,ibac)+                          &
     &                    NupPO4_ba(k,ibac))*                         &
     &                   C2pBAC

                IF(IRON==1)THEN
                  Bac_G(4)=NupFe_ba(k,ibac)*C2FeBAC
                ELSE 
                  Bac_G(4)= LARGE
                ENDIF

!  Energy limited case. All excess nutrients returned in inorganic form.
!
                IF ((Bac_G(1).le.Bac_G(2)).and.                         &
     &              (Bac_G(1).le.Bac_G(3)).and.                         &
     &              (Bac_G(1).le.Bac_G(4))) THEN
                  Het_BAC=Bac_G(1)
                  FV1=Bac_G(1)*N2cBAC
                  FV2=Bac_G(1)*P2cBAC
                  FV3=Bac_G(1)*Fe2cBAC
                  Bio_new(k,iBacN(ibac))=Bio_new(k,iBacN(ibac))+    &
     &                                     FV1
                  Bio_new(k,iBacP(ibac))=Bio_new(k,iBacP(ibac))+    &
     &                                     FV2
                  IF(IRON==1) Bio_new(k,iBacF(ibac))=                 &
     &                                     Bio_new(k,iBacF(ibac))+    &
     &                                     FV3
!
!  Uptake arrays should probably now be negative. If NH4 or PO4 is
!  positive, then there is some uptake of inorganic forms, but this
!  value will be less than the original Nup value because of IF test.
!
                  NupNH4_ba(k,ibac)=FV1-NupDON_ba(k,ibac)
                  NupPO4_ba(k,ibac)=FV2-NupDOP_ba(k,ibac)
!
!  Because Fe is considered to be all inorganic, only net uptake of Fe
!  is needed.
!
                  RelFe=NupFe_ba(k,ibac)-FV3
                  NupFe_ba(k,ibac)=FV3
!
!  Nitrogen limited case. Excess nutrients returned in organic form
!  first, inorganic second.
!
                ELSE IF ((Bac_G(2).le.Bac_G(3)).and.                    &
     &                   (Bac_G(2).le.Bac_G(4))) THEN
                  Het_BAC=Bac_G(2)
                  FV2=Bac_G(2)*P2cBAC
                  FV3=Bac_G(2)*Fe2cBAC
                  Bio_new(k,iBacN(ibac))=Bio_new(k,iBacN(ibac))+    &
     &                                     (NupDON_ba(k,ibac)+        &
     &                                      NupNH4_ba(k,ibac))
                  Bio_new(k,iBacP(ibac))=Bio_new(k,iBacP(ibac))+    &
     &                                     FV2
                  IF(IRON==1) Bio_new(k,iBacF(ibac))=                 &
     &                                     Bio_new(k,iBacF(ibac))+    &
     &                                     FV3
!
!  Uptake arrays will now reflect release of inorganic and organic
!  revision of uptake.
!
                  FV1=(Bac_G(1)-Bac_G(2))*I_Bac_Ceff
                  NupDOC_ba(k,ibac)=NupDOC_ba(k,ibac)-FV1
                  RelDOC1=FV1
!
!  To get accurate DOP from C2pDOC, must add back excreted DOC.
!
                  FV4=FV1*R_ExBAC_c*                                &
!!   &                DOC_frac(k)*                                    &
     &                Bio(k,iDOMP(ilab))/                             &
     &                Bio(k,iDOMC(ilab))
                  FV5=FV2-(NupDOP_ba(k,ibac)+                         &
                           NupPO4_ba(k,ibac)-FV4)
!
!  If FV5 is positive then released DOP is required for bacteria growth.
!
                  IF (FV5.lt.0.0_r8) THEN
                    RelDOP1=FV4
                    NupPO4_ba(k,ibac)=NupPO4_ba(k,ibac)+FV5
                  ELSE
                    RelDOP1=FV4-FV5
                  END IF
                  NupDOP_ba(k,ibac)=NupDOP_ba(k,ibac)-RelDOP1
!
!  Release Fe.
!
                  RelFe=NupFe_ba(k,ibac)-FV3
                  NupFe_ba(k,ibac)=FV3
!
!  Phosphorous limited case. Excess nutrients returned in organic form
!  first, inorganic second.
!
                ELSE IF (Bac_G(3).le.Bac_G(4)) THEN
                  Het_BAC=Bac_G(3)
                  FV2=Bac_G(3)*N2cBAC
                  FV3=Bac_G(3)*Fe2cBAC
                  Bio_new(k,iBacN(ibac))=Bio_new(k,iBacN(ibac))+    &
     &                                     FV2
                  Bio_new(k,iBacP(ibac))=Bio_new(k,iBacP(ibac))+    &
     &                                     (NupDOP_ba(k,ibac)+        &
     &                                      NupPO4_ba(k,ibac))
                  IF(IRON==1) Bio_new(k,iBacF(ibac))=                 &
     &                                     Bio_new(k,iBacF(ibac))+    &
     &                                     FV3
!
!  Uptake arrays will now reflect release of inorganic and organic
!  revision of uptake.
!
                  FV1=(Bac_G(1)-Bac_G(3))*I_Bac_Ceff
                  NupDOC_ba(k,ibac)=NupDOC_ba(k,ibac)-FV1
                  RelDOC1=FV1
!
!  To get accurate DON from C2nDOC, must add back excreted DOC.
!
                  FV4=FV1*R_ExBAC_c*                                &
!!   &                DOC_frac(k)*                                    &
     &                (Bio(k,iDOMN(ilab))/                            &
     &                 Bio(k,iDOMC(ilab)))*Frac_ExBAC_n
                  FV5=FV2-(NupDON_ba(k,ibac)+                         &
     &                     NupNH4_ba(k,ibac)-FV4)
!
!  If FV5 is positive then released DON is required for bacteria growth.
!
                  IF (FV5.lt.0.0_r8) THEN
                    RelDON1=FV4
                    NupNH4_ba(k,ibac)=NupNH4_ba(k,ibac)+FV5
                  ELSE
                    RelDON1=FV4-FV5
                  END IF
                  NupDON_ba(k,ibac)=NupDON_ba(k,ibac)-RelDON1
!
!  Release Fe.
!
                  RelFe=NupFe_ba(k,ibac)-FV3
                  NupFe_ba(k,ibac)=FV3
!
!  Fe limited case. Excess nutrients returned in organic form
!  first, inorganic second.
!
                ELSE
                  Het_BAC=Bac_G(4)
                  FV2=Bac_G(4)*N2cBAC
                  FV3=Bac_G(4)*P2cBAC
                  Bio_new(k,iBacN(ibac))=Bio_new(k,iBacN(ibac))+    &
     &                                     FV2
                  Bio_new(k,iBacP(ibac))=Bio_new(k,iBacP(ibac))+    &
     &                                     FV3
                  IF(IRON==1) Bio_new(k,iBacF(ibac))=                 &
     &                                     Bio_new(k,iBacF(ibac))+    &
     &                                     NupFe_ba(k,ibac)
!
!  Uptake arrays will now reflect release of inorganic and organic
!  revision of uptake.
!
                  FV1=(Bac_G(1)-Bac_G(4))*I_Bac_Ceff
                  NupDOC_ba(k,ibac)=NupDOC_ba(k,ibac)-FV1
                  RelDOC1=FV1
!
!  To get accurate DON from C2nDOC, must add back excreted DOC.
!
                  FV4=FV1*R_ExBAC_c*                                &
!!   &                DOC_frac(k)*                                    &
     &                Bio(k,iDOMN(ilab))/                             &
     &                Bio(k,iDOMC(ilab))*Frac_ExBAC_n
                  FV5=FV2-(NupDON_ba(k,ibac)+                         &
     &                     NupNH4_ba(k,ibac)-FV4)
!
!  If FV5 is positive then released DON is required for bacteria growth.
!
                  IF (FV5.lt.0.0_r8) THEN
                    RelDON1=FV4
                    NupNH4_ba(k,ibac)=NupNH4_ba(k,ibac)+FV5
                  ELSE
                    RelDON1=FV4-FV5
                  END IF
                  NupDON_ba(k,ibac)=NupDON_ba(k,ibac)-RelDON1
!
!  To get accurate DOP from C2pDOC, must add back excreted DOC.
!
                  FV4=FV1*R_ExBAC_c*                                &
!!   &               DOC_frac(k)*                                     &
     &                Bio(k,iDOMP(ilab))/                             &
     &                Bio(k,iDOMC(ilab))
                  FV5=FV2-(NupDOP_ba(k,ibac)+                         &
     &                     NupPO4_ba(k,ibac)-FV4)
!
!  If FV5 is positive then released DOP is required for bacteria growth.
!
                  IF (FV5.lt.0.0_r8) THEN
                    RelDOP1=FV4
                    NupPO4_ba(k,ibac)=NupPO4_ba(k,ibac)+FV5
                  ELSE
                    RelDOP1=FV4-FV5
                  END IF
                  NupDOP_ba(k,ibac)=NupDOP_ba(k,ibac)-RelDOP1
                END IF
!
!  Increment nutrient arrays.
!
                Bio_new(k,iBacC(ibac))=Bio_new(k,iBacC(ibac))+      &
     &                                   Het_BAC
                FV1=NupDOC_ba(k,ibac)-Het_BAC
                Bio_new(k,iDIC_)=Bio_new(k,iDIC_)+                  &
     &                             FV1
! 
!  NOTE: to be strictly accurate we should remove RelDOC1 from DOCNP1,
!       and then add it back, since NupDOC_ba is a net term. This should
!       wash out in the budgeting.
!
                Bio_new(k,iDOMC(ilab))=Bio_new(k,iDOMC(ilab))-      &
     &                                   (totDOC_d(k)-RelDOC1)
!!   &                                   (totDOC_d(k)-RelDOC1)*       &
!!   &                                   DOC_frac(k)
!!              IF(CDOC==1) Bio_new(k,iCDMC(ilab))=                   &
!!   &                                   Bio_new(k,iCDMC(ilab))-      &
!!   &                                   (totDOC_d(k)-RelDOC1)*       &
!!   &                                   (1.0_r8-DOC_frac(k))
!!
!  This is inclusive of RelDOX1, excretion of DON1 removed above.
!
                Bio_new(k,iDOMN(ilab))=Bio_new(k,iDOMN(ilab))-      &
     &                                   NupDON_ba(k,ibac)
                Bio_new(k,iDOMP(ilab))=Bio_new(k,iDOMP(ilab))-      &
     &                                   NupDOP_ba(k,ibac)
                Bio_new(k,iNH4_)=Bio_new(k,iNH4_)-                  &
     &                             NupNH4_ba(k,ibac)
                Bio_new(k,iPO4_)=Bio_new(k,iPO4_)-                  &
     &                             NupPO4_ba(k,ibac)
                IF(IRON==1) Bio_new(k,iFeO_)=Bio_new(k,iFeO_)-      &
     &                             NupFe_ba(k,ibac)
              END DO
            END DO
!          END DO


! Marta Rodrigues
!-----------------------------------------------------------------------
!  Zooplankton growth terms
!-----------------------------------------------------------------------
!
          DO izoo=1,Nzoo
!            DO i=1,nea
!              if(idry_e(i)==1) cycle
	      DO k=kbe(i)+1, nvrt	      
!
!  Carbon growth calculations
!                 	      
	        FV1=Bio(k,iZooC(izoo))*Gt_z(k,izoo)
		Bio_new(k,iZooC(izoo))=Bio_new(k,iZooC(izoo))+      &
		                         FV1

!
!  Nitrogen growth calculations
!
                IF(foodC(k,izoo)==0.d0)THEN
                   FV2=0.d0
                ELSE
          	   FV2=foodN(k,izoo)/foodC(k,izoo)
                ENDIF
		Bio_new(k,iZooN(izoo))=Bio_new(k,iZooN(izoo))+      &
		                         FV1*FV2

!
!  Phosphorus growth calculations
!
                IF(foodC(k,izoo)==0.d0)THEN
                   FV3=0.d0
                ELSE
		   FV3=foodP(k,izoo)/foodC(k,izoo)
                ENDIF 
		Bio_new(k,iZooP(izoo))=Bio_new(k,iZooP(izoo))+      &
		                         FV1*FV3
	      END DO
!	    END DO
	  END DO		 	 

!-----------------------------------------------------------------------

!
!-----------------------------------------------------------------------
!  Phytoplankton Losses.
!-----------------------------------------------------------------------
!

          DO iphy=1,Nphy
!            DO i=1,nea
!              if(idry_e(i)==1) cycle
              DO k=kbe(i)+1,nvrt
!
!  Excretion.
!
                IF ((C2nALG(k,iphy).ge.                               &
     &               C2nALGminABS(iphy)).and.                        &
     &              (C2pALG(k,iphy).ge.                               &
     &               C2pALGminABS(iphy)).and.                        &
     &              (HsSiO(iphy).gt.LARGE)) THEN
                  FV1=Bio(k,iPhyC(iphy))*ExALG(iphy)
                  Bio_new(k,iPhyC(iphy))=Bio_new(k,iPhyC(iphy))-    &
     &                                     FV1
!
!  No excretion of CDOC.
!
                  Bio_new(k,iDOMC(ilab))=Bio_new(k,iDOMC(ilab))+    &
     &                                     FV1
                  ELSE IF ((C2nALG(k,iphy).ge.                          &
     &                    C2nALGminABS(iphy)).and.                   &
     &                   (C2pALG(k,iphy).ge.                          &
     &                    C2pALGminABS(iphy)).and.                   &
     &                   (C2sALG(k,iphy).ge.                          &
     &                    C2SiALGminABS(iphy))) THEN
                  FV1=Bio(k,iPhyC(iphy))*ExALG(iphy)
                  Bio_new(k,iPhyC(iphy))=Bio_new(k,iPhyC(iphy))-    &
     &                                     FV1
!
!  No excretion of CDOC.
!
                  Bio_new(k,iDOMC(ilab))=Bio_new(k,iDOMC(ilab))+    &
     &                                     FV1
                END IF
!
!  Grazing.
!  And lises (Marta Rodrigues)

                IF (Bio(k,iPhyC(iphy)).gt.refuge(k,iphy)) THEN
!
!  Carbon calculations.
!
! Marta Rodrigues
!  Changed to account for zooplankton

! Calculate grazing of phytoplankton due to zooplankton (carbon)
	          DO izoo=1,Nzoo  
		     FV1=Gt_zoo_phy(k,izoo,iphy)*Bio(k,izooC(izoo))
		     !graz_act(k,i,iphy)=graz_act(k,i,iphy)*FV1
                     graz_act(k,iphy)=FV1
		  END DO 

!MFR          FV1=graz_act(k,iphy)*Bio(k,iPhyC(iphy))  
                      
		  FV1=HsGRZ(iphy)*Bio(k,iphyC(iphy))
                  FV2=HsGRZ(iphy)*Bio(k,iphyC(iphy))+graz_act(k,iphy)
               
                  Bio_new(k,iPhyC(iphy))=Bio_new(k,iPhyC(iphy))-    &
     &                                     FV2
                  Bio_new(k,iFecC(isfc))=Bio_new(k,iFecC(isfc))+    &
     &                                     FecPEL(iphy,isfc)*FV1
                  IF(Nfec==2) Bio_new(k,iFecC(iffc))=                 &
     &                                     Bio_new(k,iFecC(iffc))+    &
     &                                     FecPEL(iphy,iffc)*FV1
                  FV3=FecDOC(iphy)*FV1
                  Bio_new(k,iDOMC(ilab))=Bio_new(k,iDOMC(ilab))+    &
     &                                     (1.0_r8-cDOCfrac_c(ilab))*&
     &                                     FV3
                  IF(CDOC==1) Bio_new(k,iCDMC(ilab))=                 &
     &                                     Bio_new(k,iCDMC(ilab))+    &
     &                                     cDOCfrac_c(ilab)*FV3
                  Bio_new(k,iDIC_)=Bio_new(k,iDIC_)+                &
     &                               FecCYC(iphy)*FV1
!
!  Nitrogen calculations.
!
!MFR          FV2=graz_act(k,iphy)*Bio(k,iPhyN(iphy))

                  IF(C2nALG(k,iphy).gt.0.d0)THEN
                   FV3=HsGRZ(iphy)*Bio(k,iphyN(iphy))+              &
		       graz_act(k,iphy)/C2nALG(k,iphy)
                  ELSE
                   FV3=HsGRZ(iphy)*Bio(k,iphyN(iphy)) 
                  ENDIF

                  FV2=HsGRZ(iphy)*Bio(k,iphyN(iphy))

                  Bio_new(k,iPhyN(iphy))=Bio_new(k,iPhyN(iphy))-    &
     &                                     FV3
                  Bio_new(k,iFecN(isfc))=Bio_new(k,iFecN(isfc))+    &
     &                                     FecPEL(iphy,isfc)*FV2
                  IF(Nfec==2) Bio_new(k,iFecN(iffc))=                 &
     &                                     Bio_new(k,iFecN(iffc))+    &
     &                                     FecPEL(iphy,iffc)*FV2
                  Bio_new(k,iDOMN(ilab))=Bio_new(k,iDOMN(ilab))+    &
     &                                     FecDOC(iphy)*FV2
                  Bio_new(k,iNH4_)=Bio_new(k,iNH4_)+                &
     &                               FecCYC(iphy)*FV2
!
!  Silica calculations.
!
                  IF (iPhyS(iphy).gt.0) THEN
!MFR            FV2=graz_act(k,i,iphy)*Bio(k,i,iPhyS(iphy))

                    IF(C2sALG(k,iphy).gt.0.d0)THEN 
                     FV3=HsGRZ(iphy)*Bio(k,iphyS(iphy))+            &
		         graz_act(k,iphy)/C2sALG(k,iphy)
                    ELSE
                     FV3=HsGRZ(iphy)*Bio(k,iphyS(iphy))
                    ENDIF

                    FV2=HsGRZ(iphy)*Bio(k,iphyS(iphy)) 

                    Bio_new(k,iPhyS(iphy))=Bio_new(k,iPhyS(iphy))-  &
     &                                       FV3
!
!  Assuming that the fraction of material lost via sloppy feeding/cell
!  lysis also results in silica tests being put into FecS pool.
!
                    Bio_new(k,iFecS(isfc))=Bio_new(k,iFecS(isfc))+  &
     &                                       FecDOC(iphy)*FV2
                    IF(Nfec==2) Bio_new(k,iFecS(iffc))=               &
     &                                       Bio_new(k,iFecS(iffc))+  &
     &                                       (1.0_r8-FecDOC(iphy))*     &
     &                                       FV2
                  END IF
!
!  Phosphorus calculations.
!
!MFR          FV2=graz_act(k,iphy)*Bio(k,iPhyP(iphy))

                  IF(C2pALG(k,iphy).gt.0.d0)THEN
                    FV3=HsGRZ(iphy)*Bio(k,iphyP(iphy))+            &
		        graz_act(k,iphy)/C2pALG(k,iphy)  
		  ELSE
                    FV3=HsGRZ(iphy)*Bio(k,iphyP(iphy))
                  ENDIF

                  FV2=HsGRZ(iphy)*Bio(k,iphyP(iphy))                  

                  Bio_new(k,iPhyP(iphy))=Bio_new(k,iPhyP(iphy))-    &
     &                                     FV3
                  Bio_new(k,iFecP(isfc))=Bio_new(k,iFecP(isfc))+    &
     &                                     FecPEL(iphy,isfc)*FV2
                  IF(Nfec==2) Bio_new(k,iFecP(iffc))=                 &
     &                                     Bio_new(k,iFecP(iffc))+    &
     &                                     FecPEL(iphy,iffc)*FV2
                  Bio_new(k,iDOMP(ilab))=Bio_new(k,iDOMP(ilab))+    &
     &                                     FecDOC(iphy)*FV2
                  Bio_new(k,iPO4_)=Bio_new(k,iPO4_)+                &
     &                               FecCYC(iphy)*FV2
!
! MFR - Only if Fe is being calculated

                  IF(IRON==1)THEN
!
!  Iron calculations. Assuming no DOMF.
!
!MFR          FV2=graz_act(k,iphy)*Bio(k,iPhyF(iphy))

                   IF(C2fALG(k,iphy).gt.0.d0)THEN
                     FV3=HsGRZ(iphy)*Bio(k,iphyF(iphy))+            &
		         graz_act(k,iphy)/C2fALG(k,iphy)
                   ELSE
                     FV3=HsGRZ(iphy)*Bio(k,iphyF(iphy))
                   ENDIF 

                   FV2=HsGRZ(iphy)*Bio(k,iphyF(iphy))
		  
                   Bio_new(k,iPhyF(iphy))=Bio_new(k,iPhyF(iphy))-    &
     &                                      FV3
                   Bio_new(k,iFecF(isfc))=Bio_new(k,iFecF(isfc))+    &
     &                                      FecPEL(iphy,isfc)*FV2
                   IF(Nfec==2) Bio_new(k,iFecF(iffc))=                 &
     &                                      Bio_new(k,iFecF(iffc))+    &
     &                                      FecPEL(iphy,iffc)*FV2
                   Bio_new(k,iFeO_)=Bio_new(k,iFeO_)+                &
     &                                   (FecCYC(iphy)+                  &
     &                                   FecDOC(iphy))*FV2

                  END IF  !IRON
                END IF  !refuge
              END DO
!            END DO
          END DO
!
!  Pigment Grazing.  No fecal or dissolved terms for pigments.
!
          DO ipig=1,Npig
            DO iphy=1,Nphy
              IF (iPigs(iphy,ipig).gt.0) THEN
                itrc=iPigs(iphy,ipig)
!                DO i=1,nea
!                  if(idry_e(i)==1) cycle
                  DO k=kbe(i)+1,nvrt                  
                    IF (Bio(k,iPhyC(iphy)).gt.refuge(k,iphy)) THEN
                      FV1=HsGRZ(iphy)*Bio(k,itrc)+                 &
		          graz_act(k,iphy)*Bio(k,itrc)/Bio(k,iPhyC(iphy))  
		      Bio_new(k,itrc)=Bio_new(k,itrc) - FV1
                    END IF
                  END DO
 !               END DO
              END IF
            END DO
          END DO


! Marta Rodrigues
! Respiration
 
          DO iphy=1,Nphy
	    DO k=kbe(i)+1,nvrt
	    
! Carbon calculations	    
	    	
	      respPhy(iphy)=basalPhy(iphy)*QPhy(iphy)**((Bio(k,itemp)-10.d0)/10.d0)*&
	                    Bio(k,iPhyC(iphy))+gamaPhy(iphy)*(GtALG_r(k,iphy)        &
		            *Bio(k,iPhyC(iphy))-ExALG(iphy)*Bio(k,iPhyC(iphy)))
	      
              Bio_new(k,iPhyC(iphy))=Bio_new(k,iPhyC(iphy))-respPhy(iphy)
	      
	      Bio_new(k,iDIC_)=Bio_new(k,iDIC_)+respPhy(iphy)
	       
              DO ipig=1,Npig
		   IF (iPigs(iphy,ipig).gt.0) THEN
                      itrc=iPigs(iphy,ipig)
                      Bio_new(k,itrc)=Bio_new(k,itrc)-respPhy(iphy)
		   END IF
	      END DO
	      
! Oxygen calculations

              Bio_new(k,iDO_)=Bio_new(k,iDO_)+omegaO2C*&
                                     (GtALG_r(k,iphy)*Bio(k,iPhyC(iphy))- &
                                       respPhy(iphy))
	    END DO
	  END DO  	     


! Marta Rodrigues (! To check later...)
! Sinking
! MFR, Nov/2015 - Sinking is now computed with itr_met=3
!        IF((kbe(i)+1).ne.nvrt)THEN
!          DO iphy=1,Nphy
!            DO i=1,nea
!              if(idry_e(i)==1) cycle
!
!              DO k=kbe(i)+1,nvrt-1
!            
!                Bio_new(k,iPhyC(iphy))=Bio_new(k,iPhyC(iphy))- &
!                                      WS(iphy)*((Bio(k+1,iPhyC(iphy))- &
!                                      Bio(k-1,iPhyC(iphy)))/(ze(k+1,i)-   &
!                                      ze(k-1,i)))
!
!                Bio_new(k,iPhyN(iphy))=Bio_new(k,iPhyN(iphy))- &
!                                      WS(iphy)*((Bio(k+1,iPhyN(iphy))- &
!                                      Bio(k-1,iPhyN(iphy)))/(ze(k+1,i)-   &
!                                      ze(k-1,i)))
!
!                Bio_new(k,iPhyP(iphy))=Bio_new(k,iPhyP(iphy))- &
!                                      WS(iphy)*((Bio(k+1,iPhyP(iphy))- &
!                                      Bio(k-1,iPhyP(iphy)))/(ze(k+1,i)-   &
!                                      ze(k-1,i)))
!
!                IF(IRON==1)THEN
!                 Bio_new(k,iPhyF(iphy))=Bio_new(k,iPhyF(iphy))- &
!                                      WS(iphy)*((Bio(k+1,iPhyF(iphy))- &
!                                      Bio(k-1,iPhyF(iphy)))/(ze(k+1,i)-   &
!                                      ze(k-1,i)))
!                END IF  !IRON
!
!                IF(iphy<=2) THEN
!
!                  Bio_new(k,iPhyS(iphy))=Bio_new(k,iPhyS(iphy))- &
!                                        WS(iphy)*((Bio(k+1,iPhyS(iphy))- &
!                                        Bio(k-1,iPhyS(iphy)))/(ze(k+1,i)-   &
!                                        ze(k-1,i)))
!                END IF
!
!                DO ipig=1,Npig
!                  IF (iPigs(iphy,ipig).gt.0) THEN
!                    itrc=iPigs(iphy,ipig)
!
!                    Bio_new(k,itrc)=Bio_new(k,itrc)-- &
!                                    WS(iphy)*((Bio(k+1,itrc)- &
!                                    Bio(k-1,itrc))/(ze(k+1,i)-   &
!                                    ze(k-1,i)))
!                  END IF
!  
!                END DO
!              END DO  
!          END DO 
!          END DO ! Nphy
!        ENDIF !3D
!
!
!-----------------------------------------------------------------------
!  Bacterial losses.
!-----------------------------------------------------------------------
!
!  NOTE: Bacterial growth is completely reminerialized.
!
          DO ibac=1,Nbac
 !           DO i=1,nea
 !             if(idry_e(i)==1) cycle
              DO k=kbe(i)+1,nvrt             
!
!  Grazing calculation. (All fecal material to slow sinking pool.)
!
!! WPB - There appears to be some rounding errors that cause bacteria
!!       populations to drop just below initialization values.  Once
!!       they do, they never recover and the new lower values propagate
!!       through the model.  Only evident in Bac_P1 at the moment.
!!
!!                FV1=BacCYC(ng)*Bio_new(k,iBacC(ibac))
!!                FV2=BacPEL(ng)*Bio_new(k,iBacC(ibac))
!!                FV3=BacDOC(ng)*Bio_new(k,iBacC(ibac))
!!                FV4=FV1+FV2+FV3
!
!  Carbon calculations.
!
                Bio_new(k,iBacC(ibac))=Bio_new(k,iBacC(ibac))-      &
!!   &                                   FV4
     &                                   Bio_new(k,iBacC(ibac))
                Bio_new(k,iFecC(isfc))=Bio_new(k,iFecC(isfc))+      &
!!   &                                   FV2
     &                                   Bio_new(k,iBacC(ibac))*      &
     &                                   BacPEL
                Bio_new(k,iDOMC(ilab))=Bio_new(k,iDOMC(ilab))+      &
     &                                   (1.0_r8-cDOCfrac_c(ilab))*  &
!!   &                                   FV3
     &                                   Bio_new(k,iBacC(ibac))*      &
     &                                   BacDOC
                IF(CDOC==1) Bio_new(k,iCDMC(ilab))=                   &
     &                                   Bio_new(k,iCDMC(ilab))+      &
     &                                   cDOCfrac_c(ilab)*           &
!!   &                                   FV3
     &                                   Bio_new(k,iBacC(ibac))*      &
     &                                   BacDOC
                Bio_new(k,iDIC_)=Bio_new(k,iDIC_)+                  &
!!   &                             FV1
     &                             Bio_new(k,iBacC(ibac))*            &
     &                             BacCYC
!
!  Nitrogen calculations.
!
                Bio_new(k,iBacN(ibac))=Bio_new(k,iBacN(ibac))-      &
!!   &                                   N2cBAC(ng)*FV4
     &                                   Bio_new(k,iBacN(ibac))
                Bio_new(k,iFecN(isfc))=Bio_new(k,iFecN(isfc))+      &
!!   &                                   N2cBAC(ng)*FV2
     &                                   Bio_new(k,iBacN(ibac))*      &
     &                                   BacPEL
                Bio_new(k,iDOMN(ilab))=Bio_new(k,iDOMN(ilab))+      &
!!   &                                   N2cBAC(ng)*FV3
     &                                   Bio_new(k,iBacN(ibac))*      &
     &                                   BacDOC
                Bio_new(k,iNH4_)=Bio_new(k,iNH4_)+                  &
!!   &                             N2cBAC(ng)*FV1
     &                             Bio_new(k,iBacN(ibac))*            &
     &                             BacCYC
!
!  Phosphorous calculations.
!
                Bio_new(k,iBacP(ibac))=Bio_new(k,iBacP(ibac))-      &
!!   &                                   P2cBAC(ng)*FV4
     &                                   Bio_new(k,iBacP(ibac))
                Bio_new(k,iFecP(isfc))=Bio_new(k,iFecP(isfc))+      &
!!   &                                   P2cBAC(ng)*FV2
     &                                   Bio_new(k,iBacP(ibac))*      &
     &                                   BacPEL
                Bio_new(k,iDOMP(ilab))=Bio_new(k,iDOMP(ilab))+      &
!!   &                                   P2cBAC(ng)*FV3
     &                                   Bio_new(k,iBacP(ibac))*      &
     &                                   BacDOC
                Bio_new(k,iPO4_)=Bio_new(k,iPO4_)+                  &
!!   &                             P2cBAC(ng)*FV1
     &                             Bio_new(k,iBacP(ibac))*            &
     &                             BacCYC
!
! MFR - Only if Fe is being calculated...
                
                IF(IRON==1)THEN     

!  Iron calculations.
!
                 Bio_new(k,iBacF(ibac))=Bio_new(k,iBacF(ibac))-      &
!!   &                                    Fe2cBAC(ng)*FV4
     &                                    Bio_new(k,iBacF(ibac))
                 Bio_new(k,iFecF(isfc))=Bio_new(k,iFecF(isfc))+      &
!!   &                                    Fe2cBAC(ng)*FV2
     &                                    Bio_new(k,iBacF(ibac))*      &
     &                                    BacPEL
                 Bio_new(k,iFeO_)=Bio_new(k,iFeO_)+                  &
!!   &                              Fe2cBAC(ng)*(FV1+FV3)
     &                              Bio_new(k,iBacF(ibac))*            &
     &                              (BacDOC+BacCYC)

                END IF  !IRON                

              END DO
!            END DO
          END DO
	  
! Marta Rodrigues
!  Respiration
 
          DO ibac=1,Nbac
	    DO k=kbe(i)+1,nvrt
	  
!  Carbon calculations
	      	
	      fBac(ibac)=Bio(k,iDO_)**3/(Bio(k,iDO_)**3+HsBacO(ibac))
	    
	      respBac(ibac)=basalBac(ibac)*QBac(ibac)**((Bio(k,itemp)-10.d0)/10.d0)*&
	                    Bio(k,iBacC(ibac))+(1.d0-Bac_Ceff-GEE0C(ibac)*     &
			    (1.d0-fBac(ibac)))*Het_BAC

              Bio_new(k,iBacC(ibac))=Bio_new(k,iBacC(ibac))-respBac(ibac)
	      
	      Bio_new(k,iDIC_)=Bio_new(k,iDIC_)+respBac(ibac)
	      
!  Nitrogen calculations
!  Denitrification

              IF(DENIT_flag==1)THEN   ! Use Vichi et al., 2007
                Denit(k)=RtDenit*(1.d0*omegaO2C*(1.d0-fBac(ibac)*respBac(ibac)/MDenit))* &
		         Bio(k,iNO3_)
	      ELSE  ! Use Hang et al., 2008    
	        !Denit(k)=RtDenit*exp(1.d0-20.d0/Bio(k,itemp))*Bio(k,iNO3_)
		call parallel_abort('Only DENIT_flag=1 is available at the moment')
	      END IF
	      
	      Bio_new(k,iNO3_)=Bio_new(k,iNO3_)-Denit(k)
	      
!  Oxygen calculations

              Bio_new(k,iDO_)=Bio_new(k,iDO_)-omegaO2C*respBac(ibac)
	      
	      Bio_new(k,iCOD_)=Bio_new(k,iCOD_)+omegaS2O*omegaO2C*(1.d0*fBac(ibac))*    &
	                       respBac(ibac)-omegaS2O*omegaO2NDenit*Denit(k)

	    END DO
	  END DO	  


! Marta Rodrigues
!-----------------------------------------------------------------------
!  Zooplankton losses
!-----------------------------------------------------------------------

!  Excretion.
        
	DO izoo=1,Nzoo
!          DO i=1,nea
!            if(idry_e(i)==1) cycle
	    DO k=kbe(i)+1, nvrt
	    
!  Carbon calculations
	    
              FV1=Bio(k,iZooC(izoo))*ExZOO(izoo)
              Bio_new(k,iZooC(izoo))=Bio_new(k,iZooC(izoo))-    &
     &                                 FV1  

!  Excretion of DOC, CDOC and FecC.
!
              FV2=ZooDOC(izoo)*FV1    
	      Bio_new(k,iDOMC(ilab))=Bio_new(k,iDOMC(ilab))+    &
	                               (1.0_r8-cDOCfrac_c(ilab))*   &
     &                                 FV1
              IF(CDOC==1) Bio_new(k,iCDMC(ilab))=                 &
     &                                 Bio_new(k,iCDMC(ilab))+    &
     &                                 cDOCfrac_c(ilab)*FV1     
              Bio_new(k,iFecC(isfc))=Bio_new(k,iFecC(isfc))+    &
     &                                 ZooPEL(izoo,isfc)*FV1
              IF(Nfec==2) Bio_new(k,iFecC(iffc))=                 &
     &                                 Bio_new(k,iFecC(iffc))+    &
     &                                 ZooPEL(izoo,iffc)*FV1
     
!  Nitrogen calculations

              FV1=Bio(k,iZooN(izoo))*ExZOO(izoo)
              Bio_new(k,iZooN(izoo))=Bio_new(k,iZooN(izoo))-    &
     &                                 FV1  

!  Excretion of DON, FecN, NH4

              Bio_new(k,iDOMN(ilab))=Bio_new(k,iDOMN(ilab))+    &
     &                                 ZooDOC(izoo)*FV1
              Bio_new(k,iFecN(isfc))=Bio_new(k,iFecN(isfc))+    &
     &                                 ZooPEL(izoo,isfc)*FV1
              IF(Nfec==2) Bio_new(k,iFecN(iffc))=                 &
     &                                 Bio_new(k,iFecN(iffc))+    &
     &                                 ZooPEL(izoo,iffc)*FV1
              Bio_new(k,iNH4_)=Bio_new(k,iNH4_)+    &
     &                                 ZooCYC(izoo)*FV1          
     
!  Phosphorus calculations

              FV1=Bio(k,iZooP(izoo))*ExZOO(izoo)
              Bio_new(k,iZooP(izoo))=Bio_new(k,iZooP(izoo))-    &
     &                                 FV1  

!  Excretion of DOP, FecP, PO4

              Bio_new(k,iDOMP(ilab))=Bio_new(k,iDOMP(ilab))+    &
     &                                 ZooDOC(izoo)*FV1
              Bio_new(k,iFecP(isfc))=Bio_new(k,iFecP(isfc))+    &
     &                                 ZooPEL(izoo,isfc)*FV1
              IF(Nfec==2) Bio_new(k,iFecP(iffc))=                 &
     &                                 Bio_new(k,iFecP(iffc))+    &
     &                                 ZooPEL(izoo,iffc)*FV1
              Bio_new(k,iPO4_)=Bio_new(k,iPO4_)+    &
     &                                 ZooCYC(izoo)*FV1  

     
            END DO
!          END DO
	END DO
	                                       

!  Grazing.
 
        DO izoo=1,Nzoo
!         DO i=1,nea
!            if(idry_e(i)==1) cycle
	    DO k=kbe(i)+1, nvrt
	    
!  Carbon calculations
              
              FV1=Bio(k,iZooC(izoo))*GZ(izoo)
              Bio_new(k,iZooC(izoo))=Bio_new(k,iZooC(izoo))-    &
     &                                 FV1 
              Bio_new(k,iFecC(isfc))=Bio_new(k,iFecC(isfc))+    &
     &                                 ZooPEL(izoo,isfc)*FV1
              IF(Nfec==2) Bio_new(k,iFecC(iffc))=                 &
     &                                 Bio_new(k,iFecC(iffc))+    &
     &                                 ZooPEL(izoo,iffc)*FV1  
     
!  DIC

              Bio_new(k,iDIC_)=Bio_new(k,iDIC_)+    &
     &                           EfcPrd(izoo)*FV1         

!  Nitrogen calculations
              
	      FV1=Bio(k,iZooN(izoo))*GZ(izoo)
              Bio_new(k,iZooN(izoo))=Bio_new(k,iZooN(izoo))-    &
     &                                 FV1
      
              Bio_new(k,iFecN(isfc))=Bio_new(k,iFecN(isfc))+    &
     &                                 ZooPEL(izoo,isfc)*FV1
              IF(Nfec==2) Bio_new(k,iFecN(iffc))=                 &
     &                                 Bio_new(k,iFecN(iffc))+    &
     &                                 ZooPEL(izoo,iffc)*FV1
     
!  Phosphorus calculations

              FV1=Bio(k,iZooP(izoo))*GZ(izoo)
              Bio_new(k,iZooP(izoo))=Bio_new(k,iZooP(izoo))-    &
     &                                 FV1  
     
              Bio_new(k,iFecP(isfc))=Bio_new(k,iFecP(isfc))+    &
     &                                 ZooPEL(izoo,isfc)*FV1
              IF(Nfec==2) Bio_new(k,iFecP(iffc))=                 &
     &                                 Bio_new(k,iFecP(iffc))+    &
     &                                 ZooPEL(izoo,iffc)*FV1     
    
            END DO
!	  END DO
	END DO
	
! Respiration
 
          DO izoo=1,Nzoo
	    DO k=kbe(i)+1,nvrt   
	    
! Carbon calculations
	    	
	      respZoo(izoo)=basalZoo(izoo)*QZoo(izoo)**((Bio(k,itemp)-10.d0)/10.d0)*&
	                    Bio(k,iZooC(izoo))+(1.d0-(ExZOO(izoo)/Gt_z(k,izoo)))*   &
			    (1.d0-etaZoo(izoo))*Gt_z(k,izoo)*Bio(k,iZooC(izoo))
	      
              Bio_new(k,iZooC(izoo))=Bio_new(k,iZooC(izoo))-respZoo(izoo)
	      
	      Bio_new(k,iDIC_)=Bio_new(k,iDIC_)+respZoo(izoo)

	      
! Oxygen calculations

              Bio_new(k,iDO_)=Bio_new(k,iDO_)-omegaO2C*respZoo(izoo)
	      
	    END DO
	  END DO 	    
 
!-----------------------------------------------------------------------

!
!-----------------------------------------------------------------------
!  Fecal pellet remineralization.
!-----------------------------------------------------------------------
!
          DO ifec=1,Nfec
!            DO i=1,nea
!              if(idry_e(i)==1) cycle
              DO k=kbe(i)+1,nvrt

!
!  Carbon calculations.  All carbon goes to CO2.
!
                FV3=Regen_C(k,ifec)*Bio(k,iFecC(ifec))
                Bio_new(k,iFecC(ifec))=Bio_new(k,iFecC(ifec))-      &
     &                                   FV3
                Bio_new(k,iDIC_)=Bio_new(k,iDIC_)+                  &
     &                             FV3

!
!  Nitrogen calculations.  Nitrogen goes to NH4.
!
                FV2=Regen_N(k,ifec)*Bio(k,iFecN(ifec))
                Bio_new(k,iFecN(ifec))=Bio_new(k,iFecN(ifec))-      &
     &                                   FV2
                Bio_new(k,iNH4_)=Bio_new(k,iNH4_)+                  &
     &                             FV2
!
!  Silica calculations.
!
                FV2=Regen_S(k,ifec)*Bio(k,iFecS(ifec))
                Bio_new(k,iFecS(ifec))=Bio_new(k,iFecS(ifec))-      &
     &                                   FV2
                Bio_new(k,iSiO_)=Bio_new(k,iSiO_)+                  &
     &                             FV2
!
!  Phosphorous calculations.
!
                FV2=Regen_P(k,ifec)*Bio(k,iFecP(ifec))
                Bio_new(k,iFecP(ifec))=Bio_new(k,iFecP(ifec))-      &
     &                                   FV2
                Bio_new(k,iPO4_)=Bio_new(k,iPO4_)+                  &
     &                             FV2
!
! MFR - Only if Fe is being calculated...
!
                IF(IRON==1)THEN

!  Iron calculations.
!
                  FV2=Regen_F(k,ifec)*Bio(k,iFecF(ifec))
                  Bio_new(k,iFecF(ifec))=Bio_new(k,iFecF(ifec))-      &
       &                                   FV2
                  Bio_new(k,iFeO_)=Bio_new(k,iFeO_)+                  &
       &                             FV2

                 END IF  !IRON
              END DO
!            END DO
          END DO

! Marta Rodrigues (! To check later...)
! ----------------------------------------------------------------------
! Fecal matter sinking.
! ----------------------------------------------------------------------
! MFR, Nov/2015 - Sinking is now computed with itr_met=3
!          IF((kbe(i)+1).ne.nvrt)THEN
!            DO ifec=1,Nfec
!             DO i=1,nea
!                if(idry_e(i)==1) cycle
!                DO k=kbe(i)+1,nvrt-1
!
!                Bio_new(k,iFecC(ifec))=Bio_new(k,iFecC(ifec))- &
!                                      WF(ifec)*((Bio(k+1,iFecC(ifec))- &
!                                      Bio(k-1,iFecC(ifec)))/(ze(k+1,i)-   &
!                                      ze(k-1,i)))
!
!                Bio_new(k,iFecN(ifec))=Bio_new(k,iFecN(ifec))- &
!                                      WF(ifec)*((Bio(k+1,iFecN(ifec))- &
!                                      Bio(k-1,iFecN(ifec)))/(ze(k+1,i)-   &
!                                      ze(k-1,i)))
!
!                Bio_new(k,iFecP(ifec))=Bio_new(k,iFecP(ifec))- &
!                                      WF(ifec)*((Bio(k+1,iFecP(ifec))- &
!                                      Bio(k-1,iFecP(ifec)))/(ze(k+1,i)-   &
!                                      ze(k-1,i)))
!
!                IF(IRON==1)THEN
!                 Bio_new(k,iFecF(ifec))= Bio_new(k,iFecF(ifec))- &
!                                      WF(ifec)*((Bio(k+1,iFecF(ifec))- &
!                                      Bio(k-1,iFecF(ifec)))/(ze(k+1,i)-   &
!                                      ze(k-1,i)))
!                END IF
!
!                Bio_new(k,iFecS(ifec))=Bio_new(k,iFecS(ifec))- &
!                                      WF(ifec)*((Bio(k+1,iFecS(ifec))- &
!                                      Bio(k-1,iFecS(ifec)))/(ze(k+1,i)-   &
!                                      ze(k-1,i)))
!
!             END DO
!          END DO 
!          END DO ! Nfec
!        END IF
!
! ----------------------------------------------------------------------

!
!-----------------------------------------------------------------------
!  CDMC photolysis calculations.
!-----------------------------------------------------------------------
!
          IF (RtUVR_flag) THEN
!            DO i=1,nea
!              if(idry_e(i)==1) cycle
!
!  If Ed_nz(i,N(ng)) > zero, then there is sunlight. Standardizing rate
!  to 1500 umol quanta m-2 s-1.
!
              IF (Ed_nz(nvrt).ge.0.01) THEN
!
                FV1=RtUVR_DIC*Ed_nz(nvrt)/1500.0_r8
                FV2=RtUVR_DOC*Ed_nz(nvrt)/1500.0_r8

!
!  FV4 equals the CDMC1 absorption at 410 nm. 0.012 converts to g m-3.
!  FV5 equals the CDMC2 absorption at 410 nm.
!  Weighted average attenuation of UVB of water at 300 nm = 0.2 m-1.
!
                FV4=Bio(nvrt,iCDMC(ilab))*0.012_r8*aDOC410(ilab)
                FV5=Bio(nvrt,iCDMC(irct))*0.012_r8*aDOC410(irct)
                photo_decay=0.5_r8*Hz(nvrt)*                       &
     &                      (0.2_r8+(FV4+FV5)*aDOC300(ilab))
                FV3=EXP(-photo_decay)
                photo_decay=2.0_r8*photo_decay
!
!  Do not photolyze below the euphotic zone.
!
                DO k=nvrt,Keuphotic,-1 !número de verticais
                  IF (FV3.gt.0.01_r8) THEN
                    FV6=FV5+FV4
                    IF (FV6.gt.0.0_r8) THEN
                      FV7=FV4/FV6
                      photo_DIC=FV3*FV1*FV6
                      photo_DOC=FV3*FV2*FV6
                      total_photo=photo_DIC+photo_DOC
!
!  NOTE: not testing for excess photolysis (CDOC going negative).
!
                      FV4=(1.0_r8-FV7)*total_photo
                      Bio_new(k,iCDMC(irct))=Bio_new(k,iCDMC(irct))-&
     &                                         FV4
                      Bio_new(k,iDOMC(ilab))=Bio_new(k,iDOMC(ilab))+&
     &                                         photo_DOC
                      Bio_new(k,iCDMC(ilab))=Bio_new(k,iCDMC(ilab))-&
     &                                         FV7*total_photo
                      Bio_new(k,iDIC_)=Bio_new(k,iDIC_)+            &
     &                                         photo_DIC
                    END IF


!
!  FV4 equals the CDMC1 absorption at 410 nm. 0.012 converts to g m-3.
!  FV5 equals the CDMC2 absorption at 410 nm.
!  Weighted average attenuation of UVB of water at 300 nm = 0.2 m-1.
!
                    FV4=Bio(k,iCDMC(ilab))*0.012_r8*aDOC410(ilab)
                    FV5=Bio(k,iCDMC(irct))*0.012_r8*aDOC410(irct)
                    FV7=photo_decay+                                    &
     &                 0.5_r8*Hz(k)*(0.2_r8+(FV4+FV5)*aDOC300(ilab))
!
!  If k is greater than the bottom of the euphotic zone (and by
!  by extension the bottom boundary) or the decay constant is
!  greater than 4.61 (or < 1% photolysis zone) then exit do loop.
!
                    FV3=EXP(-FV7)
!
!  Store value for passage through entire Hz(i,j,k).
!
                    photo_decay=photo_decay+2.0_r8*FV7
                  END IF
                END DO
              END IF
!            END DO
          END IF
!
!-----------------------------------------------------------------------
!  Create optimal pigment ratios.
!-----------------------------------------------------------------------
!
!          DO i=1,nea
!            if(idry_e(i)==1) cycle
            IF (Keuphotic.le.nvrt) THEN
              DO iphy=1,Nphy
!
!  Carbon to chlorophyll a ratio
!  This statement says that nutrient limitation of C2CHL ratio overides
!  light adaptation. Minimum of two functions may be more ecologically
!  accurate?
!
                DO k=nvrt,Keuphotic,-1
                  IF (b_C2Cn(iphy).lt.0.0_r8+SMALL) THEN
                    C2CHL_w(k,iphy)=MIN((b_C2Cl(iphy)+               &
     &                                   mxC2Cl(iphy)*E0_nz(k)),   &
     &                                  C2CHL_max(iphy))
                  ELSE IF (C2nALG(k,iphy).gt.                         &
     &                     minC2nALG(iphy)+SMALL) THEN
                    C2CHL_w(k,iphy)=b_C2Cn(iphy)+                    &
     &                              mxC2Cn(iphy)*                    &
     &                              (C2nALG(k,iphy)-                  &
     &                               minC2nALG(iphy))
                  ELSE
                    C2CHL_w(k,iphy)=MIN((b_C2Cl(iphy)+               &
     &                                   mxC2Cl(iphy)*E0_nz(k)),   &
     &                                  C2CHL_max(iphy))
                  END IF
                END DO
!
!  Chlorophyll a concentation per species. form g CHL a / g C
!
                DO k=nvrt,Keuphotic,-1
                  Pigs_w(k,iphy,ichl)=1.0_r8/C2CHL_w(k,iphy)
                END DO
!
!  Chlorophyll b concentration per species. form g CHL b / g C
!
                IF (iPigs(iphy,2).gt.0) THEN
                  DO k=nvrt,Keuphotic,-1
                    Pigs_w(k,iphy,2)=b_ChlB(iphy)+                   &
     &                               mxChlB(iphy)*                   &
     &                               (C2CHL_w(k,iphy)-                  &
     &                                b_C2Cl(iphy))
                    Pigs_w(k,iphy,2)=Pigs_w(k,iphy,2)*                  &
     &                               Pigs_w(k,iphy,ichl)
                  END DO
                END IF
!
!  Chlorophyll c concentration per species. form g CHL c / g C
!
                IF (iPigs(iphy,3).gt.0) THEN
                  DO k=nvrt,Keuphotic,-1
                    Pigs_w(k,iphy,3)=b_ChlC(iphy)+                   &
     &                               mxChlC(iphy)*                   &
     &                               (C2CHL_w(k,iphy)-                  &
     &                                 b_C2Cl(iphy))
                    Pigs_w(k,iphy,3)=Pigs_w(k,iphy,3)*                  &
     &                               Pigs_w(k,iphy,ichl)
                  END DO
                END IF
!
!  Photosynthetic caroteniods per species. form g PSC / g C
!
                IF (iPigs(iphy,4).gt.0) THEN
                  DO k=nvrt,Keuphotic,-1
                    Pigs_w(k,iphy,4)=b_PSC(iphy)+                    &
     &                               mxPSC(iphy)*                    &
     &                               (C2CHL_w(k,iphy)-                  &
     &                                b_C2Cl(iphy))
                    Pigs_w(k,iphy,4)=Pigs_w(k,iphy,4)*                  &
     &                               Pigs_w(k,iphy,ichl)
                  END DO
                END IF
!
!  Photoprotective caroteniods per species. form g PPC / g C
!
                IF (iPigs(iphy,5).gt.0) THEN
                  DO k=nvrt,Keuphotic,-1
                    Pigs_w(k,iphy,5)=b_PPC(iphy)+                    &
     &                               mxPPC(iphy)*                    &
     &                               (C2CHL_w(k,iphy)-                  &
     &                                b_C2Cl(iphy))
                    Pigs_w(k,iphy,5)=Pigs_w(k,iphy,5)*                  &
     &                               Pigs_w(k,iphy,ichl)
                  END DO
                END IF
!
!  Low Urobilin Phycoeurythin concentration per species. g LPUB / g C
!
                IF (iPigs(iphy,6).gt.0) THEN
                  DO k=nvrt,Keuphotic,-1
                    Pigs_w(k,iphy,6)=b_LPUb(iphy)+                   &
     &                               mxLPUb(iphy)*                   &
     &                               (C2CHL_w(k,iphy)-                  &
     &                                b_C2Cl(iphy))
                    Pigs_w(k,iphy,6)=Pigs_w(k,iphy,6)*                  &
     &                               Pigs_w(k,iphy,ichl)
                  END DO
                END IF
!
!  High Urobilin Phycoeurythin concentration per species (g HPUB / g C).
!
                IF (iPigs(iphy,7).gt.0) THEN
                  DO k=nvrt,Keuphotic,-1
                    Pigs_w(k,iphy,7)=b_HPUb(iphy)+                   &
     &                               mxHPUb(iphy)*                   &
     &                               (C2CHL_w(k,iphy)-                  &
     &                                b_C2Cl(iphy))
                    Pigs_w(k,iphy,7)=Pigs_w(k,iphy,7)*                  &
     &                               Pigs_w(k,iphy,ichl)
                  END DO
                END IF
              END DO
!
!  Calculate pigment ratio changes.
!  NOTE: 12 factor to convert to ugrams (mg m-3)
!
              DO ipig=1,Npig
                DO iphy=1,Nphy
                  IF (iPigs(iphy,ipig).gt.0) THEN
                    itrc=iPigs(iphy,ipig)
                    DO k=nvrt,Keuphotic,-1
                      IF ((Bio(k,iPhyC(iphy)).gt.0.0_r8).and.         &
     &                    (Bio(k,itrc).gt.0.0_r8)) THEN
                        FV1=Bio(k,iPhyC(iphy))*12.0_r8
                        FV2=GtALG_r(k,iphy) !-respPhy(iphy)  ! MFR-Check respiration here
                        FV3=FV1/                                        &
     &                      (FV2/Pigs_w(k,iphy,ipig)+                   &
     &                       FV1*(1.0_r8-FV2)/                          &
     &                       Bio(k,itrc))
                        Bio_new(k,itrc)=Bio_new(k,itrc)+            &
     &                                    (FV3-Bio(k,itrc))
                      END IF
                    END DO
                  END IF
                END DO
              END DO
            END IF
!          END DO

! Marta Rodrigues
!-----------------------------------------------------------------------
! Oxygen calculations
!-----------------------------------------------------------------------
!
          DO k=kbe(i)+1,nvrt 
            reox=reox_COD*QCOD**((Bio(k,itemp)-10.d0)/10.d0)*      &
	         Bio(k,iDO_)*Bio(k,iCOD_)/(Bio(k,iDO_)+HsCOD) 
            Bio_new(k,iCOD_)=Bio_new(k,iCOD_)-reox

! Nitrification and reoxidation
            Bio_new(k,iDO_)=Bio_new(k,iDO_)-omegaO2N*NH4toNO3(k)-   &
                            (1/omegaS2O)*reox
          END DO

!-----------------------------------------------------------------------
!
!-----------------------------------------------------------------------
!  Update the tendency arrays
!-----------------------------------------------------------------------
!MFR,Apr/2017 - Removed BioIter loop
!          DO ibio=1,NBIT
!            itrc=idbio(ibio)
!            DO i=1,nea
!              if(idry_e(i)==1) cycle
!              DO k=kbe(i)+1,nvrt
!MFR             Bio(k,i,itrc)=Bio(k,i,itrc)+dtbio*Bio_new(k,i,itrc)
!                Bio(k,itrc)=Bio_old(k,itrc)+dtbio*Bio_new(k,itrc)
!              END DO
!            END DO
!          END DO
!        END DO ITER_LOOP
!
!-----------------------------------------------------------------------
!  Update global tracer variables.
!  MFR - Calculates body source for SELFE
!-----------------------------------------------------------------------
!
        DO ibio=1,NBIT
          itrc=idbio(ibio)
!          DO i=1,nea
!            if(idry_e(i)==1) cycle
            DO k=kbe(i)+1,nvrt

! MFR ... To check tracer concentration
!               Bio(k,itrc)=Bio_old(k,itrc)+Bio_new(k,itrc)
               Bio(k,itrc)=Bio_old(k,itrc)+Bio_new(k,itrc)*dt
               IF(Bio(k,itrc)<MinVal)THEN
                 bdy_frc(itrc+irange_tr(1,6)-1,k,i)=0.d0 !MFR - Check later...
               ELSE 
!                 bdy_frc(itrc+irange_tr(1,6)-1,k,i)=Bio_new(k,itrc)/dt ! MFR: removed BioIter loop
                  bdy_frc(itrc+irange_tr(1,6)-1,k,i)=Bio_new(k,itrc)
                 !IF(isnan(bdy_frc(itrc,k,i)))THEN
                 !ZYL: isnan is not recognized on some system
                 IF(bdy_frc(itrc+irange_tr(1,6)-1,k,i)/=bdy_frc(itrc+irange_tr(1,6)-1,k,i))THEN
                   WRITE(600,*) myrank,bdy_frc(itrc+irange_tr(1,6)-1,k,i),itrc,k,i
                 ENDIF
               END IF
! MFR               t(k,i,itrc)=MAX(MinVal,                           &
! MFR     &                                t(k,i,itrc)+               &
! MFR     &                                Hz(k,i)*                        &
! MFR     &                                (Bio(k,i,itrc)-Bio_old(k,i,itrc)))

            END DO
!          END DO
        END DO

      END DO  ! For each elemet (i=1,nea)

! Marta Rodrigues
! ---------------------------------------------------------------------
! Calculates bottom flux for SELFE
! ---------------------------------------------------------------------
!       flx_bt(ntracers,nea)
        flx_bt(irange_tr(1,6):irange_tr(2,6),:) = 0.d0
!
! ---------------------------------------------------------------------
! Calculates surface flux for SELFE
! ---------------------------------------------------------------------
!       flx_sf(ntracers,nea)
        flx_sf(irange_tr(1,6):irange_tr(2,6),:) = 0.d0
	
! Include surface flux from reaeration for DO and DIC
        wind_speed = 0.d0 
        DO i=1,nea
            ! Compute uwind and vwind @ elements
            u_wind1=sum(windx1(elnode(1:i34(i),i)))/i34(i)
            u_wind2=sum(windx2(elnode(1:i34(i),i)))/i34(i)
            v_wind1=sum(windy1(elnode(1:i34(i),i)))/i34(i)
            v_wind2=sum(windy2(elnode(1:i34(i),i)))/i34(i)

            uwind(i)=u_wind1+wtratio*(u_wind2-u_wind1)
            vwind(i)=v_wind1+wtratio*(v_wind2-v_wind1)
	    wind_speed(i)=SQRT(uwind(i)*uwind(i)+vwind(i)*vwind(i))
            ! To prevent unrealistic Kreaer
	    IF(wind_speed(i)>=10.0d0) wind_speed(i)=10.0d0         
	    temp(i)=tr_el(1,nvrt,i) !tsel(1,nvrt,i)
	    salt(i)=tr_el(2,nvrt,i) !tsel(2,nvrt,i)

	    IF(REAER_flag==1)then !Use Wanninkhof,1992
	      
	       ScO2=(1800.6-120.1*temp(i)+3.7818*temp(i)**2-0.047608*&
	            &temp(i)**3)*(1+0.00314*salt(i))
	       
	       ScCO2=(2073.1-125.62*temp(i)+3.6276*temp(i)**2-&
	             &0.043219*temp(i)**3)
	    
	    ELSE !Thomann and Muller, 1987; Antonopoulos and Gianniou, 2003
	    	    
	       IF(wind_speed(i)<=3.5)THEN
	          K20(i)=0.2*wind_speed(i)
	       ELSE
	          K20(i)=0.057*wind_speed(i)*wind_speed(i)
	       ENDIF

               Kreaer(i)=K20(i)*QWind**(temp(i)-20)
	       Kreaer(i)=Kreaer(i)/86400.d0

            ENDIF  !REAER_flag
	    
! DO	    
	    DOsat(i)=14.6244-0.367134*temp(i)+0.0044972*temp(i)**2-&
	             0.0966*salt(i)+0.00205*salt(i)*temp(i)+0.0002739*salt(i)**2
	    
            DOsat(i)=DOsat(i)*31.25        !mmolO2/m3		     
		     
	    OW(i)=psiWind*exp(alfaWind*wind_speed(i))
	    
	    temp_scale(i)=log((298.15-temp(i))/(273.15+temp(i)))
	    
	    solub(i)=2.00856+3.224*temp_scale(i)+3.99063*temp_scale(i)**2+    &
	             4.80299*temp_scale(i)**3+0.978188*temp_scale(i)**4+      &
		     1.71069*temp_scale(i)**5+salt(i)*(-0.00624067-0.00693498*&
		     temp_scale(i)-0.00690358*temp_scale(i)**2-0.00429155*    &
		     temp_scale(i)**3)-0.00000031168*salt(i)**2
	    
            solub(i)=exp(solub(i))
            
            solub(i)=solub(i)/22391.6
           
	    solub(i)=solub(i)*32000.d0     !mg/l
	    
	    DOW(i)=OW(i)*solub(i)*0.0143
	    DOW(i)=DOW(i)*31.25            !mmolO2/m3
	    
	    IF(REAER_flag==1) Kreaer(i)=(0.24d0/86400.d0)*0.31*wind_speed(i)**2*&
	                                (ScO2(i)/660.d0)**(-0.5d0)  
	    
	    
	    !reaerDO(i)=Kreaer(i)*(DOsat(i)+DOW(i)-tr_el(iDO_,nvrt,i))
	    reaerDO(i)=Kreaer(i)*(DOsat(i)+DOW(i)-tr_el(iDO_-1+irange_tr(1,6),nvrt,i))

	    flx_sf(iDO_+irange_tr(1,6)-1,i)=reaerDO(i)
	    
! DIC (CO2)

            IF(REAER_flag==1) Kreaer(i)=(0.24d0/86400.d0)*0.31*wind_speed(i)**2*&
	                                (ScCO2(i)/660.d0)**(-0.5d0)
					
            solubCO2(i)=exp(-60.2409+9345.17/(temp(i)+273.15)+23.3585*log((temp(i)+273.15)/&
	                100.d0)+salt(i)*(0.023517-0.00023656*(temp(i)+273.15)+&
			0.00000047036*(temp(i)+273.15)**2))

	    solubCO2(i)=solubCO2(i)*eqstate(7,temp(i),salt(i))/1000.d0
	    
	    kdiss1(i)=exp(-2307.1266/(temp(i)+273.15)+2.83655-1.5529413*log(temp(i)+273.15)+&
	              ((-4.0484/(temp(i)+273.15))-0.20760841)*salt(i)**0.5+0.08468345*salt(i)-&
		      0.00654208*salt(i)**(3.d0/2.d0)+log(1-0.001005*salt(i)))*&
                      eqstate(8,temp(i),salt(i))/1000.d0
			
	    kdiss2(i)=exp(-3351.6106/(temp(i)+273.15)-9.226508-0.2005743*log(temp(i)+273.15)+&
	              ((-23.9722/(temp(i)+273.15))-0.106901773)*salt(i)**0.5+0.1130822*salt(i)-&
		      0.00846934*salt(i)**(3.d0/2.d0)+log(1-0.001005*salt(i)))*&
                      eqstate(9,temp(i),salt(i))/1000.d0
	    
	    Hconc=10.d0**(-pH)
	    
	    CO2star(i)=1.d0/(Hconc*Hconc+kdiss1(i)*Hconc+kdiss1(i)*kdiss2(i))
            !CO2star(i)=CO2star(i)*Hconc**2.d0*tr_el(iDIC_,nvrt,i)*10.d0**(-6.d0) 	
            CO2star(i)=CO2star(i)*Hconc**2.d0*tr_el(iDIC_-1+irange_tr(1,6),nvrt,i)*10.d0**(-6.d0) 	
					
	    pCO2w(i)=CO2star(i)/solubCO2(i)				
										
            reaerCO2(i)=Kreaer(i)*solubCO2(i)*(pCO2w(i)-pCO2a)
	    
	    flx_sf(iDIC_+irange_tr(1,6)-1,i)=reaerCO2(i)
	    
	END DO     	    
!----------------------------------------------------------------------	   
	      
!
! ---------------------------------------------------------------------

      RETURN
      END SUBROUTINE
