!Note: most arrays in this file are from SCHISM directly (too account for
!quads)
#include "wwm_functions.h"
#ifdef SCHISM
!**********************************************************************
!*   This routine is for RADFLAG=LON (Longuet-Higgins)
!* March 2022, LRU team : update with the new roller model
!**********************************************************************
      SUBROUTINE RADIATION_STRESS_SCHISM
        USE DATAPOOL
        USE schism_glbl, only: errmsg,hmin_radstress,npa,nsa,idry_s,isidenode 
        USE schism_msgp 
        IMPLICIT NONE

        INTEGER      :: IP, IS, ID
        REAL(rkind)  :: COSE2, SINE2, COSI2, ELOC, ErLOC, WN, HTOT
        REAL(rkind)  :: ACLOC(MSC,MDC), RACLOC(MSC,MDC)
        REAL(rkind)  :: DSXX3D(2,NVRT,nsa), DSYY3D(2,NVRT,nsa), DSXY3D(2,NVRT,nsa)

! GD: imet_dry allows to choose between 2 different methods to compute the
! derivative at the sides between wet and dry elements:
!
! imet_dry=1 : only the values at the 2 nodes of the side are used to
! compute the derivative (this older method showed to provide inconsistent
! wave force at the wet/dry interface).
!
! imet_dry=2 : a 4-point stencil (the 3 wet nodes and an artificial
! node at the center of the side) is used to compute the derivative.
! This method is similar to using shape functions to compute the
! derivative at the center of the element and assigning this value to the
! the side center.

        ! Initialisation of the arrays
        IMET_DRY = 2
        HTOT = DMIN
        RSXX = ZERO; RSXY = ZERO; RSYY = ZERO
        SXX3D = ZERO; SYY3D = ZERO; SXY3D = ZERO
        DSXX3D = ZERO; DSYY3D = ZERO; DSXY3D = ZERO
        WWAVE_FORCE = ZERO

        DO IP = 1, MNP
          IF (idry(IP)==1) CYCLE 
          ACLOC = AC2(:,:,IP)

          ! Wave contribution [m^3/s^2]
          DO ID = 1, MDC
            COSE2 = COS(SPDIR(ID))**TWO
            SINE2 = SIN(SPDIR(ID))**TWO
            COSI2 = COS(SPDIR(ID)) * SIN(SPDIR(ID))
            DO IS = 2, MSC
              !Added grav in ELOC
              ELOC  = 0.5_rkind*G9*(SPSIG(IS)*ACLOC(IS,ID)+SPSIG(IS-1)*ACLOC(IS-1,ID))*DS_INCR(IS)*DDIR
              WN    = CG(IS,IP) / ( SPSIG(IS)/WK(IS,IP) )
              RSXX(IP) = RSXX(IP) + ( WN*COSE2 + WN - 0.5_rkind) * ELOC ![m^3/s^2]
              RSXY(IP) = RSXY(IP) + ( WN*COSI2          ) * ELOC
              RSYY(IP) = RSYY(IP) + ( WN*SINE2 + WN - 0.5_rkind) * ELOC
            END DO
          END DO

          ! Surface roller contribution [m^3/s^2] (NB: already contains G9, check in wwm_compute_roller.F90)
          ! Reference: e.g. Apotsos et al. (2007)
          IF (IROLLER == 1) THEN
            RSXX(IP) = RSXX(IP) + 2.0D0*EROL2(IP)*COS(DROLP(IP))**TWO
            RSXY(IP) = RSXY(IP) + 2.0D0*EROL2(IP)*COS(DROLP(IP))*SIN(DROLP(IP))
            RSYY(IP) = RSYY(IP) + 2.0D0*EROL2(IP)*SIN(DROLP(IP))**TWO
          END IF

          ! Storing depth-averaged radiation stresses [m^3/s^2]
          SXX3D(:,IP) = RSXX(IP)
          SXY3D(:,IP) = RSXY(IP)
          SYY3D(:,IP) = RSYY(IP)
        END DO

        ! Computing gradients of the depth-averaged radiation stress terms  [m^2/s^2]
        CALL hgrad_nodes(IMET_DRY,0,NVRT,MNP,nsa,SXX3D,DSXX3D)   ! (dSxx/dx , dSxx/dy )
        CALL hgrad_nodes(IMET_DRY,0,NVRT,MNP,nsa,SYY3D,DSYY3D)   ! (dSyy/dx , dSyy/dy )
        CALL hgrad_nodes(IMET_DRY,0,NVRT,MNP,nsa,SXY3D,DSXY3D)   ! (dSxy/dx , dSxy/dy )
        CALL exchange_s3d_2(DSXX3D)
        CALL exchange_s3d_2(DSYY3D)
        CALL exchange_s3d_2(DSXY3D)

        ! Computing the wave forces; these are noted Rsx, Rsy in Rolland et al. (2012), see Eq. (9)
        ! These are stored in wwave_force(:,1:nsa,1:2)  [m/s^2]
        DO IS = 1, nsa
          IF(idry_s(IS) == 1) CYCLE

          ! Total water depth at sides
          HTOT = MAX((DEP(isidenode(1,IS)) + DEP(isidenode(2,IS)))/2.0D0,hmin_radstress)

          ! Wave forces
          WWAVE_FORCE(1,:,IS) = WWAVE_FORCE(1,:,IS) - (DSXX3D(1,:,IS) + DSXY3D(2,:,IS)) / HTOT
          WWAVE_FORCE(2,:,IS) = WWAVE_FORCE(2,:,IS) - (DSXY3D(1,:,IS) + DSYY3D(2,:,IS)) / HTOT
        END DO !IS

      END SUBROUTINE RADIATION_STRESS_SCHISM

!**********************************************************************
!*  This routine is used with RADFLAG=VOR (3D vortex formulation, after Bennis et al., 2011)
!*  => Computation of the wave-induced pressure term at nodes (the gradient is computed directly 
!*  at sides when calculating the forces) and the Stokes drift velocities. The latter are 
!*  computed at all levels, at nodes and sides, and for both the wave and roller (kept separated).
!**********************************************************************
      SUBROUTINE STOKES_STRESS_INTEGRAL_SCHISM
        USE DATAPOOL
        USE schism_glbl, ONLY: errmsg, hmin_radstress, ns, kbs, kbe, nea, idry_e, &
                            &  isdel, indel, elnode, dldxy, zs, area, nsa, idry_s, &
                            &  isidenode, nne
        USE schism_msgp
        IMPLICIT NONE

        INTEGER     :: IP, k, ID, IS, IL
        REAL(rkind) :: D_loc, k_loc, kD_loc, z_loc, E_loc, Er_loc, JPress_loc
        REAL(rkind) :: Uint, Vint, Urint, Vrint
        REAL(rkind) :: USTOKES_loc(NVRT), VSTOKES_loc(NVRT), UrSTOKES_loc(NVRT), VrSTOKES_loc(NVRT)

        integer     :: IE, isd, j, l, n1, n2, n3, icount
        real(rkind) :: tmp0, tmp1, tmp2, ztmp, ubar, vbar, dhdx, dhdy
        real(rkind) :: STOKES_WVEL_ELEM(NVRT,MNE), ws_tmp1(NVRT,nsa),ws_tmp2(NVRT,nsa)
        real(rkind) :: dr_dxy_loc(2,NVRT,nsa)

!...    Computing Stokes drift horizontal velocities at nodes and pressure term
        DO IP = 1, MNP
          IF(idry(IP) == 1) CYCLE

          ! Total water depth at the node
          D_loc = MAX( DEP(IP) , hmin_radstress )

          ! Initialization of the local Stokes drift and J variables
          USTOKES_loc  = 0.D0; VSTOKES_loc  = 0.D0;
          UrSTOKES_loc = 0.D0; VrSTOKES_loc = 0.D0;
          JPress_loc   = 0.D0;

          ! Loop on the frequencies
          DO IS = 1, MSC
            Uint   = 0.D0; Vint   = 0.D0; Urint   = 0.D0; Vrint   = 0.D0
            k_loc  = MIN(KDMAX/DEP(IP),WK(IS,IP))
            kD_loc = MIN(KDMAX,WK(IS,IP)*D_loc)

            ! Loop on the directions
            DO ID = 1, MDC
              E_loc      = AC2(IS,ID,IP)*SPSIG(IS)*DDIR*DS_INCR(IS)
              JPress_loc = JPress_loc + G9*k_loc*E_loc/DSINH(2.D0*kD_loc)
              Uint       = Uint + SPSIG(IS)*k_loc*COSTH(ID)*E_loc
              Vint       = Vint + SPSIG(IS)*k_loc*SINTH(ID)*E_loc
            END DO !MDC

            ! Loop on the vertical nodes
            DO IL = KBP(IP), NVRT
              ! Here we need to compute z+h of Eq. C.1 of Bennis et al. (2011)
              ! In her framework, z varies from -h to eta, meaning that z+h corresponds to the distance to the bed
              ! -ZETA(KBP(IP),IP) corresponds to h, the depth at node IP (not the total water depth)
              ! Waves
              z_loc           = ZETA(IL,IP) - ZETA(KBP(IP),IP)
              USTOKES_loc(IL) = USTOKES_loc(IL) + Uint*DCOSH(2.D0*k_loc*z_loc)/DSINH(kD_loc)**2
              VSTOKES_loc(IL) = VSTOKES_loc(IL) + Vint*DCOSH(2.D0*k_loc*z_loc)/DSINH(kD_loc)**2
            END DO !NVRT
          END DO !MSC

          ! Surface roller contribution to horizontal Stokes drift velocities
          ! NB: we do not just add the contribution and keep separated arrays.
          ! This is motivated by the fact that we do not want this contribution to
          ! influence Wst, which is computed from the continuity equation for waves only
          
          IF (IROLLER == 1) THEN
            Urint = 2.D0*COS(DROLP(IP))*EROL2(IP)/(CROLP(IP)*D_loc)
            Vrint = 2.D0*SIN(DROLP(IP))*EROL2(IP)/(CROLP(IP)*D_loc)

            ! Homogeneous across depth
            UrSTOKES_loc = Urint
            VrSTOKES_loc = Vrint

            ! Making sure, the Stokes drift velocities do not blow up in very shallow water
            IF (D_loc < 2.D0*hmin_radstress) THEN
              UrSTOKES_loc = SIGN(MIN(0.1D0*SQRT(G9*D_loc),ABS(Urint)),Urint)
              VrSTOKES_loc = SIGN(MIN(0.1D0*SQRT(G9*D_loc),ABS(Vrint)),Vrint)
            END IF
          END IF

          ! Storing Stokes drift horizontal velocities
          STOKES_HVEL(1,:,IP) = USTOKES_loc
          STOKES_HVEL(2,:,IP) = VSTOKES_loc
          ! Surface rollers
          IF (IROLLER == 1) THEN
            ! Smoothing the roller contribution to the Stokes drift velocity near the shoreline
            ! With this profile, U_st < 10% of computed U_st at h < DMIN, and U_st > 95% of computed U_st at h > 2.25*DMIN
            IF (D_loc < 1.5D0*DMIN) THEN
              ROLLER_STOKES_HVEL(1,:,IP) = UrSTOKES_loc*(SINH(DEP(IP))/SINH(1.5D0))**2
              ROLLER_STOKES_HVEL(2,:,IP) = VrSTOKES_loc*(SINH(DEP(IP))/SINH(1.5D0))**2
            ELSE
              ROLLER_STOKES_HVEL(1,:,IP) = UrSTOKES_loc
              ROLLER_STOKES_HVEL(2,:,IP) = VrSTOKES_loc
            END IF
          END IF

          ! Storing pressure term
          JPRESS(IP) = JPress_loc
        END DO !MNP

!...    Computing Stokes drift horizontal velocities at sides (in pframe if ics=2)
        ! The average of the values from vertically adjacent nodes is taken
        STOKES_HVEL_SIDE = 0.D0; ROLLER_STOKES_HVEL_SIDE = 0.D0
        DO IS = 1,nsa
          IF(idry_s(IS) == 1) CYCLE

          ! Indexes of surrounding nodes
          n1 = isidenode(1,IS); n2 = isidenode(2,IS)
          DO k = kbs(IS),NVRT
            ! Waves
            STOKES_HVEL_SIDE(1,k,IS) = (STOKES_HVEL(1,k,n1) + STOKES_HVEL(1,k,n2))/2.D0
            STOKES_HVEL_SIDE(2,k,IS) = (STOKES_HVEL(2,k,n1) + STOKES_HVEL(2,k,n2))/2.D0

            ! Surface rollers
            IF (IROLLER == 1) THEN
              ROLLER_STOKES_HVEL_SIDE(1,k,IS) = (ROLLER_STOKES_HVEL(1,k,n1) + ROLLER_STOKES_HVEL(1,k,n2))/2.D0
              ROLLER_STOKES_HVEL_SIDE(2,k,IS) = (ROLLER_STOKES_HVEL(2,k,n1) + ROLLER_STOKES_HVEL(2,k,n2))/2.D0
            END IF
          END DO
        END DO !nsa

!...    Compute bottom Stokes drift z-vel. at elements
        STOKES_WVEL_ELEM = 0.D0
        DO IE = 1,nea
           IF(idry_e(IE) == 1) CYCLE

           ! Index of the surrounding nodes
           n1 = elnode(1,IE)
           n2 = elnode(2,IE)
           n3 = elnode(3,IE)
           IF(kbe(IE) == 0) THEN
             WRITE(errmsg,*)'Error: kbe(i) == 0'
             CALL parallel_abort(errmsg)
           END IF

           ubar = (STOKES_HVEL(1,max(kbp(n1),kbe(IE)),n1) + STOKES_HVEL(1,max(kbp(n2),kbe(IE)),n2) &
               & + STOKES_HVEL(1,max(kbp(n3),kbe(IE)),n3))/3.D0 !average bottom stokes-x-vel
           vbar = (STOKES_HVEL(2,max(kbp(n1),kbe(IE)),n1) + STOKES_HVEL(2,max(kbp(n2),kbe(IE)),n2) &
               & + STOKES_HVEL(2,max(kbp(n3),kbe(IE)),n3))/3.D0 !average bottom stokes-y-vel
           dhdx = DEP8(n1)*dldxy(1,1,IE) + DEP8(n2)*dldxy(2,1,IE) + DEP8(n3)*dldxy(3,1,IE) !eframe
           dhdy = DEP8(n1)*dldxy(1,2,IE) + DEP8(n2)*dldxy(2,2,IE) + DEP8(n3)*dldxy(3,2,IE)
           STOKES_WVEL_ELEM(kbe(IE),IE) = -dhdx*ubar - dhdy*vbar
        END DO !nea

!...    Compute bottom Stokes z-vel. at nodes
        STOKES_WVEL = 0.D0
        DO IP = 1,np !residents only
           IF(idry(IP) == 1) CYCLE

           !Bottom Stokes z-vel.
           tmp0 = 0.D0
           DO j = 1,nne(IP)
              ie = indel(j,IP)
              IF(idry_e(ie)==0) THEN
                STOKES_WVEL(kbp(IP),IP) = STOKES_WVEL(kbp(IP),IP) + STOKES_WVEL_ELEM(kbe(ie),ie)*area(ie)
              END IF
              tmp0 = tmp0 + area(ie)
           END DO !j
           STOKES_WVEL(kbp(IP),IP) = STOKES_WVEL(kbp(IP),IP)/tmp0
        END DO !np

!...    Compute horizontal gradient of Stokes x and y-vel. (to compute Stokes z-vel.)
        ws_tmp1 = 0.D0; ws_tmp2 = 0.D0
        CALL hgrad_nodes(2,0,NVRT,MNP,nsa,STOKES_HVEL(1,:,:),dr_dxy_loc)
        ws_tmp1(:,:) = dr_dxy_loc(1,:,:) !valid only in resident
        CALL hgrad_nodes(2,0,NVRT,MNP,nsa,STOKES_HVEL(2,:,:),dr_dxy_loc)
        ws_tmp2(:,:) = dr_dxy_loc(2,:,:)

!...    Compute Stokes z-vel. at all levels: STOKES_WVEL_SIDE(NVRT,nsa)
        STOKES_WVEL_SIDE = 0.D0
        DO IS = 1,ns !residents only
          IF(idry_s(IS) == 1) CYCLE
          n1 = isidenode(1,IS)
          n2 = isidenode(2,IS)

          !Bottom Stokes z-vel.
          STOKES_WVEL_SIDE(kbs(IS),IS) = (STOKES_WVEL(max(kbs(IS),kbp(n1)),n1) + STOKES_WVEL(max(kbs(IS),kbp(n2)),n2))/2.D0

          !Stokes z-vel. at all levels
          DO k = kbs(IS)+1, NVRT
            ztmp = zs(k,IS) - zs(k-1,IS)
            STOKES_WVEL_SIDE(k,IS) = STOKES_WVEL_SIDE(k-1,IS) &
                                   & -(ws_tmp1(k,IS)+ws_tmp1(k-1,IS))/2.D0*ztmp &
                                   & -(ws_tmp2(k,IS)+ws_tmp2(k-1,IS))/2.D0*ztmp
           END DO
        END DO !ns

      END SUBROUTINE STOKES_STRESS_INTEGRAL_SCHISM

!**********************************************************************
!*  This routine is used with RADFLAG=VOR (3D vortex formulation, after
!Bennis et al., 2011)
!*  => Computation of the conservative terms A1 and B1 from Eq. (11) and
!(12) respectively
!**********************************************************************
            SUBROUTINE COMPUTE_CONSERVATIVE_VF_TERMS_SCHISM
        USE DATAPOOL
        USE schism_glbl, ONLY: kbs, ns, idry_e, isdel, elnode, dldxy, cori, zs, su2, sv2,nsa,idry_s
        USE schism_msgp
        IMPLICIT NONE

        integer     :: IS, IE, k, l, icount
        real(rkind) :: dJ_dx_loc, dJ_dy_loc, du_loc, dv_loc, dz_loc, Ust_loc, Vst_loc, &
                       &VF_x_loc,VF_y_loc, STCOR_x_loc, STCOR_y_loc
        real(rkind) :: du_dxy(2,NVRT,nsa), dv_dxy(2,NVRT,nsa)

!...    Initialisation
        WWAVE_FORCE = 0.D0

!...    Computing the spatial derivative of horizontal velocities
        CALL hgrad_nodes(2,0,NVRT,MNP,nsa,uu2,du_dxy)
        CALL hgrad_nodes(2,0,NVRT,MNP,nsa,vv2,dv_dxy)

!...    Main loop over the sides
        DO IS = 1,ns !resident
          IF(idry_s(IS) == 1) CYCLE

          !------------------------
          ! Pressure term (grad(J))
          
          icount = 0; dJ_dx_loc = 0; dJ_dy_loc = 0

          IF (fwvor_gradpress == 1) THEN ! BM
            DO l = 1,2 !elements
              IE = isdel(l,IS)
              IF(ie /= 0 .AND. idry_e(IE) == 0) THEN
                icount = icount + 1
                dJ_dx_loc = dJ_dx_loc + dot_product(JPRESS(elnode(1:3,IE)),dldxy(1:3,1,IE)) !in eframe
                dJ_dy_loc = dJ_dy_loc + dot_product(JPRESS(elnode(1:3,IE)),dldxy(1:3,2,IE))
              END IF
            END DO
            ! Averaging the values from the two surrounding elements
            IF(icount > 2) CALL parallel_abort('Pressure term:icount>2')
            IF(icount == 2) THEN
              dJ_dx_loc = dJ_dx_loc/2.D0
              dJ_dy_loc = dJ_dy_loc/2.D0
            END IF
          END IF

          !---------------------------------------
          ! Depth-varying conservative wave forces 
          
          du_loc = 0; dv_loc = 0; dz_loc = 1
          
          DO k = kbs(IS),NVRT
          
            IF (fwvor_advz_stokes == 1) THEN ! BM
              ! du/dz and dv/dz terms
              IF (k == kbs(IS) .OR. k == kbs(IS)+1) THEN
                dz_loc = zs(kbs(IS)+2,IS) - zs(kbs(IS)+1,IS)
                du_loc = su2(kbs(IS)+2,IS) - su2(kbs(IS)+1,IS)
                dv_loc = sv2(kbs(IS)+2,IS) - sv2(kbs(IS)+1,IS)
              ELSE IF (k == NVRT) THEN
                dz_loc = zs(k,IS) - zs(k-1,IS)
                du_loc = su2(k,IS) - su2(k-1,IS)
                dv_loc = sv2(k,IS) - sv2(k-1,IS)
              ELSE
                dz_loc = zs(k+1,IS) - zs(k-1,IS)
                du_loc = su2(k+1,IS) - su2(k-1,IS)
                dv_loc = sv2(k+1,IS) - sv2(k-1,IS)
              END IF
            END IF

            ! Stokes drift velocity
            ! LRU team : switch off roller contribution, which is only accounted 
            ! for within continuity equation. This is motivated by the fact that VF
            ! arises from the irrotational part of the wave motion as opposed
            ! to surface rollers.
            Ust_loc = 0.D0; Vst_loc = 0.D0
            IF (fwvor_advxy_stokes == 1) THEN
              !IF (IROLLER == 1) THEN
              !  Ust_loc = STOKES_HVEL_SIDE(1,k,IS) + ROLLER_STOKES_HVEL_SIDE(1,k,IS)
              !  Vst_loc = STOKES_HVEL_SIDE(2,k,IS) + ROLLER_STOKES_HVEL_SIDE(2,k,IS)
              !ELSE
              !  Ust_loc = STOKES_HVEL_SIDE(1,k,IS)
              !  Vst_loc = STOKES_HVEL_SIDE(2,k,IS)
              !END IF
              Ust_loc = STOKES_HVEL_SIDE(1,k,IS)
              Vst_loc = STOKES_HVEL_SIDE(2,k,IS)
            END IF

            ! Vortex force 
            !  x axis : -du/dy*v_s + dv/dx*v_s - W_s*du/dz
            !  y axis : +du/dy*u_s - dv/dx*u_s - W_s*dv/dz
            VF_x_loc =  - du_dxy(2,k,IS)*Vst_loc + dv_dxy(1,k,IS)*Vst_loc - STOKES_WVEL_SIDE(k,IS)*du_loc/dz_loc
            VF_y_loc =  + du_dxy(2,k,IS)*Ust_loc - dv_dxy(1,k,IS)*Ust_loc - STOKES_WVEL_SIDE(k,IS)*dv_loc/dz_loc
            
            ! Stokes-Coriolis
            ! x axis : f*v_s
            ! y axis : -f*U_st
            STCOR_x_loc = cori(IS)*Vst_loc
            STCOR_y_loc = -cori(IS)*Ust_loc
            
            ! Saving wave forces
            WWAVE_FORCE(1,k,IS) = WWAVE_FORCE(1,k,IS) + VF_x_loc + STCOR_x_loc - dJ_dx_loc
            WWAVE_FORCE(2,k,IS) = WWAVE_FORCE(2,k,IS) + VF_y_loc + STCOR_y_loc - dJ_dy_loc
            
          END DO
        END DO !ns

        ! Exchange between ghost regions
        CALL exchange_s3d_2(WWAVE_FORCE)

      END SUBROUTINE COMPUTE_CONSERVATIVE_VF_TERMS_SCHISM


!**********************************************************************
!*  This routine is used with RADFLAG=VOR (3D vortex formulation, after Bennis et al., 2011)
!*  => Computation of the non-conservative terms due to depth-induced breaking (term Fb from Eq. (11) and (12))
!*  March 2022 : update LRU team
!    Accounts for depth-induced breaking, roller (if turned on) and whitecapping contribution
!**********************************************************************
      SUBROUTINE COMPUTE_BREAKING_VF_TERMS_SCHISM
        USE DATAPOOL
        USE schism_glbl, ONLY: hmin_radstress, kbs, ns, isbs, dps, h0, out_wwm, &
                               &zs,nsa,idry_s,isidenode
        USE schism_msgp 
        IMPLICIT NONE

        INTEGER     :: IS, isd, k, j, l, n1, n2, n3, icount
        REAL(rkind) :: eta_tmp, tmp0, htot, sum_2D, sum_3D
        REAL(rkind) :: Fdb_x_loc, Fdb_y_loc, Fds_x_loc, Fds_y_loc
        REAL(rkind) :: swild_2D(NVRT), swild_3D(NVRT)
        
        ! Compute sink of momentum due to wave breaking 
        DO IS = 1, ns
        
          ! Check IF dry segment or open bnd segment
          IF(idry_s(IS) == 1 .or. isbs(IS) > 0) CYCLE
          
          ! Water depth at side
          n1 = isidenode(1,IS); n2 = isidenode(2,IS)
          eta_tmp = (eta2(n1) + eta2(n2))/2.D0
          !htot = max(h0,dps(IS)+eta_tmp,hmin_radstress) ! KM
          htot = max(h0,dps(IS)+eta_tmp)

          IF(kbs(IS)+1 == NVRT) THEN !2D
          
            Fdb_x_loc = 0.D0 ; Fdb_y_loc = 0.D0
            Fds_x_loc = 0.D0 ; Fds_y_loc = 0.D0

            ! N.B. average between the two adjacent nodes
            
            ! Depth-induced breaking and roller contribution
            IF (IROLLER == 1) THEN
              Fdb_x_loc = -((1.D0-ALPROL)*(SBR(1,n1) + SBR(1,n2)) + SROL(1,n1) + SROL(1,n2))/2.D0/htot 
              Fdb_y_loc = -((1.D0-ALPROL)*(SBR(2,n1) + SBR(2,n2)) + SROL(2,n1) + SROL(2,n2))/2.D0/htot
            ELSE
              Fdb_x_loc = - (SBR(1,n1) + SBR(1,n2))/2.D0/htot
              Fdb_y_loc = - (SBR(2,n1) + SBR(2,n2))/2.D0/htot
            ENDIF
            
            ! Whitecapping contribution
            Fds_x_loc = -(SDS(1,n1) + SDS(1,n2))/2.D0/htot
            Fds_y_loc = -(SDS(2,n1) + SDS(2,n2))/2.D0/htot
            
            ! Save breaking wave force
            WWAVE_FORCE(1,:,IS) = WWAVE_FORCE(1,:,IS) + Fdb_x_loc + Fds_x_loc
            WWAVE_FORCE(2,:,IS) = WWAVE_FORCE(2,:,IS) + Fdb_y_loc + Fds_y_loc
            
          ELSE !3D

            ! Threshold on Hs
            tmp0 = (out_wwm(n1,1) + out_wwm(n2,1))/2.D0 !Hs
            IF(tmp0 <= 0.005D0) CYCLE

            ! Vertical distribution function of qdm (due to wave breaking)
            swild_3D = 0.D0
            DO k = kbs(IS), NVRT
              ! Homogeneous vertical distribution
              IF (ZPROF_BREAK == 1) swild_3D(k) = 1.D0 
              ! Hyperbolic distribution
              IF (ZPROF_BREAK == 2) swild_3D(k) = cosh((dps(IS)+zs(k,IS))/(0.2D0*tmp0))
              IF (ZPROF_BREAK == 3) swild_3D(k) = 1.D0 - dtanh(((eta_tmp-zs(k,IS))/(0.5D0*tmp0))**2.D0)
              IF (ZPROF_BREAK == 4) swild_3D(k) = 1.D0 - dtanh(((eta_tmp-zs(k,IS))/(0.5D0*tmp0))**4.D0)
              IF (ZPROF_BREAK == 5) swild_3D(k) = 1.D0 - dtanh(((eta_tmp-zs(k,IS))/(0.5D0*tmp0))**8.D0)
              ! All in the two surface layers
              IF (ZPROF_BREAK == 6 .AND. k .GE. NVRT-1) swild_3D(k)=1.D0
            END DO !k

            ! In shallow depths, we make the vertical profile tend to a vertical-uniform one
            ! Objectives: 1 - vertical mixing; 2 - numerical stability
            !IF (htot .LT. 5.D0*DMIN_SCHISM) swild_3D = 1.D0 + (swild_3D - 1.D0)*tanh((0.2D0*htot/DMIN_SCHISM)**8.D0)
            !IF (htot .LT. 2.0D0) swild_3D = 1.D0 + (swild_3D - 1.D0)*tanh((htot/2.0D0)**8.D0)
 
            ! Integral of the vertical distribution function
            sum_3D = 0.0D0
            DO k = kbs(IS), NVRT-1
              sum_3D = sum_3D + (swild_3D(k+1) + swild_3D(k))/2.D0*(zs(k+1,IS) - zs(k,IS))
            END DO !NVRT-1

            DO k = kbs(IS), NVRT
            
              Fdb_x_loc = 0.D0 ; Fdb_y_loc = 0.D0
              Fds_x_loc = 0.D0 ; Fds_y_loc = 0.D0
              
              ! Depth-induced breaking and roller contribution
              IF (IROLLER == 1) THEN
                Fdb_x_loc = -swild_3D(k)*((1.D0-ALPROL)*(SBR(1,n1) + SBR(1,n2)) + SROL(1,n1) + SROL(1,n2))/2.D0/sum_3D 
                Fdb_y_loc = -swild_3D(k)*((1.D0-ALPROL)*(SBR(2,n1) + SBR(2,n2)) + SROL(2,n1) + SROL(2,n2))/2.D0/sum_3D
              ELSE
                Fdb_x_loc = -swild_3D(k)*(SBR(1,n1) + SBR(1,n2))/2.D0/sum_3D
                Fdb_y_loc = -swild_3D(k)*(SBR(2,n1) + SBR(2,n2))/2.D0/sum_3D
              ENDIF
              ! Whitecapping contribution
              Fds_x_loc = -swild_3D(k)*(SDS(1,n1) + SDS(1,n2))/2.D0/sum_3D
              Fds_y_loc = -swild_3D(k)*(SDS(2,n1) + SDS(2,n2))/2.D0/sum_3D
              ! Save breaking wave force
              WWAVE_FORCE(1,k,IS) = WWAVE_FORCE(1,k,IS) + Fdb_x_loc + Fds_x_loc
              WWAVE_FORCE(2,k,IS) = WWAVE_FORCE(2,k,IS) + Fdb_y_loc + Fds_y_loc
            END DO
          END IF !2D/3D
          
          !! Smoothing wave forces near the shoreline
          !! With this profile, F < 10% of computed F at h < DMIN, and F > 95% of computed F at h > 2.25*DMIN
          !!IF (htot < 8.*DMIN) WWAVE_FORCE(:,:,IS) = WWAVE_FORCE(:,:,IS)*tanh((0.5D0*htot/DMIN)**8.D0)
          !IF (htot < 1.5D0) WWAVE_FORCE(1,:,IS) = WWAVE_FORCE(1,:,IS)*(SINH(htot)/SINH(1.5D0))**2
          !IF (htot < 0.8D0) WWAVE_FORCE(2,:,IS) = WWAVE_FORCE(2,:,IS)*(SINH(htot)/SINH(0.8D0))**2

        END DO !nsa

        ! Exchange between ghost regions
        CALL exchange_s3d_2(WWAVE_FORCE)

      END SUBROUTINE COMPUTE_BREAKING_VF_TERMS_SCHISM

!**********************************************************************
!*  This routine is used with RADFLAG=VOR (3D vortex formulation, after Bennis et al., 2011)
!*  => Computation of the non-conservative terms due to bottom friction (see Uchiyama et al., 2010)
!*  TO DO : pass the vertical distribution in option similar to breaking wave force
!**********************************************************************
      SUBROUTINE COMPUTE_STREAMING_VF_TERMS_SCHISM
        ! MP
        USE DATAPOOL
        USE schism_glbl, ONLY: hmin_radstress, kbs, ns, isbs, dps, h0, out_wwm,&
                               &zs, idry_s, isidenode, nchi, rough_p, iwbl, delta_wbl
        USE schism_msgp 
        IMPLICIT NONE

        INTEGER     :: IS, isd, k, j, l, n1, n2, n3, icount
        REAL(rkind) :: eta_tmp, tmp0, tmp1, tmp2, htot, sum_2D, sum_3D
        REAL(rkind) :: Fws_x_loc,Fws_y_loc
        REAL(rkind) :: swild_2D(NVRT), swild_3D(NVRT)

        ! Compute sink of momentum due to wave breaking 
        DO IS = 1, ns
          ! Check IF dry segment or open bnd segment
          IF(idry_s(IS) == 1 .or. isbs(IS) > 0) CYCLE
          
          ! Water depth at side
          n1 = isidenode(1,IS); n2 = isidenode(2,IS)
          eta_tmp = (eta2(n1) + eta2(n2))/2.D0
          !htot = max(h0,dps(IS)+eta_tmp,hmin_radstress) ! KM
          htot = max(h0,dps(IS)+eta_tmp)
   
          IF(kbs(IS)+1 == NVRT) THEN !2D
          
            ! N.B. average between the two adjacent nodes
            Fws_x_loc = - (SBF(1,n1) + SBF(1,n2))/2.D0/htot
            Fws_y_loc = - (SBF(2,n1) + SBF(2,n2))/2.D0/htot
            ! Saving wave streaming
            WWAVE_FORCE(1,:,IS) = WWAVE_FORCE(1,:,IS) + Fws_x_loc 
            WWAVE_FORCE(2,:,IS) = WWAVE_FORCE(2,:,IS) + Fws_y_loc

          ELSE !3D
          
            ! Threshold on WBBL
            ! 1/kwd = awd * delta_wbl
            ! we take awd = 1 but literature suggests awd>1
            ! we note 1/kwd = tmp0
            tmp0 = (delta_wbl(n1) + delta_wbl(n2))/2.D0
            IF(tmp0 .LT. SMALL) CYCLE
            
            ! Vertical distribution function of qdm
            swild_3D = 0.D0
            DO k = kbs(IS), NVRT
              ! Homogeneous vertical distribution
              swild_3D(k) = 1.D0
              ! Hyperbolic distribution - Type of profile 1
              !swild_3D(k) = cosh((eta_tmp-zs(k,IS))/tmp0)
              ! Hyperbolic distribution - Type of profile 2
              swild_3D(k) = 1.D0 - dtanh(((dps(IS)+zs(k,IS))/tmp0)**2.D0)
              !swild_3D(k) = 1.D0 - dtanh(((dps(IS)+zs(k,IS))/tmp0)**4.D0)
              !swild_3D(k) = 1.D0 - dtanh(((dps(IS)+zs(k,IS))/tmp0)**8.D0)
            END DO !k
            
            ! In shallow depths, we make the vertical profile tend to a vertical-uniform one
            ! Objectives: 1 - vertical mixing; 2 - numerical stability
            !IF (htot .LT. 2.0D0) swild_3D = 1.D0 + (swild_3D - 1.D0)*tanh((htot/2.0D0)**8.D0)

            ! Integral of the vertical distribution function
            sum_3D = 0.0D0
            DO k = kbs(IS), NVRT-1
              sum_3D = sum_3D + (swild_3D(k+1) + swild_3D(k))/2.D0*(zs(k+1,IS) - zs(k,IS))
            END DO !NVRT-1

            DO k = kbs(IS), NVRT
              Fws_x_loc = - swild_3D(k)*(SBF(1,n1) + SBF(1,n2))/2.D0/sum_3D
              Fws_y_loc = - swild_3D(k)*(SBF(2,n1) + SBF(2,n2))/2.D0/sum_3D
              ! Saving wave streaming
              WWAVE_FORCE(1,k,IS) = WWAVE_FORCE(1,k,IS) + Fws_x_loc 
              WWAVE_FORCE(2,k,IS) = WWAVE_FORCE(2,k,IS) + Fws_y_loc              
            END DO
          END IF !2D/3D

        END DO !MNS

        ! Exchange between ghost regions
        CALL exchange_s3d_2(WWAVE_FORCE)

      END SUBROUTINE COMPUTE_STREAMING_VF_TERMS_SCHISM

!**********************************************************************
!*   This routine fixes the wave forces to the barotropic gradient at the numerical shoreline (boundary between dry and wet elements)			
!**********************************************************************
      SUBROUTINE SHORELINE_WAVE_FORCES
        USE DATAPOOL
        USE schism_glbl, only: errmsg,hmin_radstress,idry_e,idry_s,isidenode, & 
                             & isdel,elnode,i34,dldxy,grav,thetai,ns 
        USE schism_msgp 
        IMPLICIT NONE

        INTEGER      :: IP,IS,INODE_1,INODE_2,IELEM
        REAL(rkind)  :: TMP_X,TMP_Y,BPGR(2)

!... 	Loop on the resident sides
        do IS = 1,ns
          if(idry_s(IS) == 1) cycle

          ! Adjacent nodes index
          INODE_1 = isidenode(1,IS) ! Side node #1
          INODE_2 = isidenode(2,IS) ! Side node #2
          if(idry(INODE_1) == 1 .OR. idry(INODE_2) == 1) cycle  

          ! Sides we are not interested in
          if(isdel(1,IS) == 0 .OR. isdel(2,IS) == 0) cycle ! Boundaries
          if(idry_e(isdel(1,IS)) == 0 .AND. idry_e(isdel(2,IS)) == 0) cycle ! Case where both adjacent elements are wet    
          if(idry_e(isdel(1,IS)) == 1 .AND. idry_e(isdel(2,IS)) == 1) cycle ! Case where both adjacent elements are dry (should never occur anyway)
          !if(isbnd(1,INODE_1) == 1 .OR. isbnd(1,INODE_2) == 1) cycle ! Case where the side touches open boundaries
          if(isbnd(1,INODE_1)>0.OR.isbnd(1,INODE_2)>0) cycle ! Case where the side touches open boundaries

          ! Reinitializing the wave force
          WWAVE_FORCE(:,:,IS) = 0

          ! We are left with sides that belong to one dry element and one wet element
          ! We store the elements indexes for future use
          if(idry_e(isdel(1,IS)) == 0 .AND. idry_e(isdel(2,IS)) == 1) then 
            IELEM = isdel(1,IS)
          elseif(idry_e(isdel(2,IS)) == 0 .AND. idry_e(isdel(1,IS)) == 1) then 
            IELEM = isdel(2,IS)
          else
            cycle
          endif

          ! We compute the barotropic gradient at the shoreline (only the wet element is used)
          BPGR = 0
          do IP = 1,i34(IELEM)
            ! x and y-components of grad(eta)
            TMP_X = (1-thetai)*eta1(elnode(IP,IELEM))*dldxy(IP,1,IELEM) + thetai*eta2(elnode(IP,IELEM))*dldxy(IP,1,IELEM)
            TMP_Y = (1-thetai)*eta1(elnode(IP,IELEM))*dldxy(IP,2,IELEM) + thetai*eta2(elnode(IP,IELEM))*dldxy(IP,2,IELEM)
            ! Barotropic gradient = g*grad(eta)
            BPGR(1) = BPGR(1) + grav*TMP_X
            BPGR(2) = BPGR(2) + grav*TMP_Y
          enddo !IP

          ! Overwriting wave forces to balance out pressure gradient
          WWAVE_FORCE(1,:,IS) = BPGR(1)
          WWAVE_FORCE(2,:,IS) = BPGR(2)
        enddo !IS

        call exchange_s3d_2(WWAVE_FORCE)

      END SUBROUTINE SHORELINE_WAVE_FORCES


!**********************************************************************
!* This routine, called in main, applies a ramp to wave fores starting
!  from the open boundary.
!  This ramp is defined in the input file wafo_ramp.gr3, and is read 
!  in wwm_initio  (subroutine READ_WAFO_OPBND_RAMP).
!  If wafo_obcramp==1, APPLY_WAFO_OPBND_RAMP is called in main at
!  each time step.
!  Authors: X. Bertin & B. Mengual (05/2020)
!**********************************************************************
      SUBROUTINE APPLY_WAFO_OPBND_RAMP

        USE DATAPOOL
        USE schism_glbl, only : ns,kbs,idry_s,nsa,wafo_opbnd_ramp
        USE schism_msgp

        IMPLICIT NONE

        INTEGER      :: IS,k

        DO IS = 1,nsa 
          IF(idry_s(IS) == 1) CYCLE
          DO k = kbs(IS),NVRT
            WWAVE_FORCE(1,k,IS) = WWAVE_FORCE(1,k,IS)*wafo_opbnd_ramp(IS)
            WWAVE_FORCE(2,k,IS) = WWAVE_FORCE(2,k,IS)*wafo_opbnd_ramp(IS)
          END DO ! NVRT
        END DO ! nsa

      END SUBROUTINE APPLY_WAFO_OPBND_RAMP


#endif /*SCHISM*/
!**********************************************************************
!*                                                                    *
!**********************************************************************
