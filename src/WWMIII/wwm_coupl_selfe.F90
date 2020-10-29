!Note: most arrays in this file are from SCHISM directly (too account for
!quads)
#include "wwm_functions.h"
#ifdef SCHISM
!**********************************************************************
!*   This routine is for RADFLAG=LON (Longuet-Higgins)
!**********************************************************************
      SUBROUTINE RADIATION_STRESS_SCHISM
        USE DATAPOOL
        use schism_glbl, only: errmsg,hmin_radstress,npa,nsa,idry_s,isidenode 
        USE schism_msgp 
        IMPLICIT NONE

        INTEGER :: IP,IS,ID
        REAL(rkind)  :: ACLOC(MSC,MDC),RACLOC(MSC,MDC)
        REAL(rkind)  :: COSE2, SINE2, COSI2
        REAL(rkind)  :: EWK(MNP),EWS(MNP),EWN(MNP),ETOT(MNP),MDIR(MNP)
        REAL(rkind)  :: m0, m0d, tmp, EHFR, ELOC, ErLOC, EFTAIL, DVEC2RAD
        REAL(rkind)  :: DS, ETOTS, ETOTC, EWSIG
        REAL(rkind)  :: WNTMP,WKTMP,WCGTMP,WCTMP,WN,WKDEPTMP
        REAL(rkind)  :: WSTMP, DEPLOC

        REAL(rkind)  :: HTOT
        REAL(rkind)  :: DSXX3D(2,NVRT,MNS), DSYY3D(2,NVRT,MNS),DSXY3D(2,NVRT,MNS)


!GD: imet_dry allows to choose between 2 different methods to compute the
!derivative at the sides between wet and dry elements:
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

        EFTAIL = ONE / (PTAIL(1) - ONE)

        ETOT = ZERO
        MDIR = ZERO

        IF (LETOT) THEN
!AR: Estimate zeroth moment m0, mean wave direction, dominant wave number, dominant sigma ...
          DO IP = 1, MNP
!            IF (ABS(IOBP(IP)) .GT. 0) CYCLE
            !IF (DEP(IP) .LT. DMIN) CYCLE
            IF (idry(IP)==1) CYCLE

            DEPLOC = DEP(IP)
            ACLOC = AC2(:,:,IP)
            m0    = ZERO
            EWSIG  = ZERO
            ETOTS  = ZERO
            ETOTC  = ZERO
            IF (MSC .GE. 2) THEN
              DO ID = 1, MDC
                m0d = ZERO
                DO IS = 2, MSC
                  tmp = 0.5_rkind*(SPSIG(IS)*ACLOC(IS,ID)+SPSIG(IS-1)*ACLOC(IS-1,ID))*DS_INCR(IS)*DDIR
                  m0 = m0 + tmp
                  EWSIG  = EWSIG  + SPSIG(IS) * tmp
                  m0d = m0d + tmp
                END DO
                IF (MSC > 3) THEN
                  EHFR = ACLOC(MSC,ID) * SPSIG(MSC)
                  m0 = m0 + DDIR * EHFR * SPSIG(MSC) * EFTAIL
                endif
                ETOTC  = ETOTC + m0d * COS(SPDIR(ID))
                ETOTS  = ETOTS + m0d * SIN(SPDIR(ID))
              END DO
            ELSE
              DS = SGHIGH - SGLOW
              DO ID = 1, MDC
                m0d = ACLOC(1,ID) * DS * DDIR
                m0 = m0 + m0d
              END DO
            END IF
            ETOT(IP) = m0
            IF (m0 .GT. small .and. .not. dep(ip) .lt. dmin) then
              EWS(IP) = EWSIG/m0
              WSTMP = EWSIG/m0
              CALL ALL_FROM_TABLE(WSTMP,DEPLOC,WKTMP,WCGTMP,WKDEPTMP,WNTMP,WCTMP)
              EWN(IP) = WNTMP
              EWK(IP) = WKTMP
              MDIR(IP) = DVEC2RAD (ETOTC, ETOTS)
            ELSE
              EWS(IP)  = ZERO 
              EWN(IP)  = ZERO 
              EWK(IP)  = 10. 
              MDIR(IP) = ZERO 
            END IF 
          END DO !IP
        END IF !LETOT

!AR: Here comes the whole story ... 
! Etot = 1/16 * Hs² = 1/8 * Hmono² => Hs² = 2 * Hmono² => Hs = sqrt(2) * Hmono => Hmono = Hs / SQRT(2) 
! Etot = 1/16 * Hs² = 1/16 * (4 * sqrt(m0))² = m0 
! Etot = 1/8 * Hmono² ... so the problem for the analytical solution evolved because we treat the Etot from Hs and Hmono there is a factor of 2 between this!
! Or in other words for the analytical solution we impose a Hs = X[m], we integrate m0 out of it and get Etot, since this Etot is a function of Hs and not Hmono^X^O
! it needs the factor of 2 between it! This should make now things clear forever. So the question is not how we calculate the total energy the question is 
! what is defined on the boundary that means we should always recalculate the boundary in terms of Hs =  SQRT(2) * Hmono !!!
! Or saying it again in other words our boundary conditions is wrong if we impose Hmono in wwminput.nml !!!


        IF (RADFLAG .EQ. 'LON') THEN
          DO IP = 1, MNP
            IF (idry(IP)==1) CYCLE

            IF (.NOT. LETOT) THEN
              ACLOC = AC2(:,:,IP)
              IF (IROLLER == 1) RACLOC = RAC2(:,:,IP)
              DO ID = 1, MDC
                COSE2 = COS(SPDIR(ID))**TWO
                SINE2 = SIN(SPDIR(ID))**TWO
                COSI2 = COS(SPDIR(ID)) * SIN(SPDIR(ID))
                DO IS = 2, MSC
                  ELOC  = 0.5_rkind * (SPSIG(IS)*ACLOC(IS,ID)+SPSIG(IS-1)*ACLOC(IS-1,ID))*DS_INCR(IS)*DDIR
                  WN    = CG(IS,IP) / ( SPSIG(IS)/WK(IS,IP) )
                  ! Wave contribution
                  RSXX(IP) = RSXX(IP) + ( WN * COSE2 + WN - 0.5_rkind) * ELOC   ! Units = [ 1/s + 1/s - 1/s ] * m²s = m²
                  RSXY(IP) = RSXY(IP) + ( WN * COSI2          ) * ELOC
                  RSYY(IP) = RSYY(IP) + ( WN * SINE2 + WN - 0.5_rkind) * ELOC
                  ! Roller contribution
                  IF (IROLLER == 1) THEN
                    ErLOC  = 0.5_rkind *(SPSIG(IS)*RACLOC(IS,ID)+SPSIG(IS-1)*RACLOC(IS-1,ID))*DS_INCR(IS)*DDIR
                    RSXX(IP) = RSXX(IP) + 2 * COSE2 * ErLOC  ! Units = [ 1/s+ 1/s - 1/s ] * m²s = m²
                    RSXY(IP) = RSXY(IP) + 2 * COSI2 * ErLOC  ! For the '2' factor, check Stive and de Vriend (1994) and the Appendix by Rolf Deigaard
                    RSYY(IP) = RSYY(IP) + 2 * SINE2 * ErLOC
                  END IF !IROLLER
                END DO !IS
              END DO !ID

            ELSE IF (LETOT) THEN
              RSXX(IP) =  ETOT(IP) * (EWN(IP)*((EWK(IP)*SIN(MDIR(IP)))**TWO/EWK(IP)**TWO+ONE)-0.5_rkind)
              RSXY(IP) =  ETOT(IP) *  EWN(IP)* EWK(IP)*SIN(MDIR(IP))*EWK(IP)*COS(MDIR(IP))* ONE/EWK(IP)
              RSYY(IP) =  ETOT(IP) * (EWN(IP)*((EWK(IP)*COS(MDIR(IP)))**TWO/EWK(IP)**TWO+ONE)-0.5_rkind)
            END IF !LETOT

          END DO

        ELSE
          call parallel_abort('R.S.: unknown R.S. model') 
        END IF !RADFLAG 

        ! Transforming into depth-averaged radiation stress terms,
        ! varying in the vertical (unit: m^2/s/s)
        DO IP = 1, MNP
          IF (idry(IP) == 1) CYCLE
          SXX3D(:,IP) = RSXX(IP) * G9
          SXY3D(:,IP) = RSXY(IP) * G9
          SYY3D(:,IP) = RSYY(IP) * G9
        END DO

        ! Computing gradients of the depth-averaged radiation stress
        ! terms (unit: m^2/s/s)
        CALL hgrad_nodes(IMET_DRY,0,NVRT,MNP,MNS,SXX3D,DSXX3D)   ! (dSxx/dx , dSxx/dy )
        CALL hgrad_nodes(IMET_DRY,0,NVRT,MNP,MNS,SYY3D,DSYY3D)   ! (dSyy/dx , dSyy/dy )
        CALL hgrad_nodes(IMET_DRY,0,NVRT,MNP,MNS,SXY3D,DSXY3D)   ! (dSxy/dx , dSxy/dy )
        CALL exchange_s3d_2(DSXX3D)
        CALL exchange_s3d_2(DSYY3D)
        CALL exchange_s3d_2(DSXY3D)

        ! Computing the wave forces; these are noted Rsx, Rsy in Rolland
        ! et al. (2012), see Eq. (9)
        ! These are stored in wwave_force(:,1:MNS,1:2) (unit: m/s/s)
        DO IS = 1, MNS
          IF(idry_s(IS) == 1) CYCLE

          ! Total water depth at sides
          HTOT = MAX((DEP(isidenode(1,IS)) + DEP(isidenode(2,IS)))/2.0D0,hmin_radstress)

          ! Wave forces
          WWAVE_FORCE(1,:,IS) = WWAVE_FORCE(1,:,IS) - (DSXX3D(1,:,IS) + DSXY3D(2,:,IS)) / HTOT
          WWAVE_FORCE(2,:,IS) = WWAVE_FORCE(2,:,IS) - (DSXY3D(1,:,IS) + DSYY3D(2,:,IS)) / HTOT
        END DO !IS


      END SUBROUTINE RADIATION_STRESS_SCHISM

!**********************************************************************
!*  This routine is used with RADFLAG=VOR (3D vortex formulation, after
!Bennis et al., 2011)
!*  => Computation of the wave-induced pressure term at nodes (the
!gradient is computed directly
!*  at sides when calculating the forces) and the Stokes drift
!velocities. The latter are
!*  computed at all levels, at nodes and sides, and for both the wave
!and roller (kept separated).
!**********************************************************************
      SUBROUTINE STOKES_STRESS_INTEGRAL_SCHISM
 
        USE DATAPOOL
        USE schism_glbl, ONLY: errmsg, hmin_radstress, ns, kbs, kbe, nea, idry_e, &
                            &  isdel, indel, elnode, dldxy, zs, area
        USE schism_msgp
        IMPLICIT NONE

        INTEGER     :: IP, k, ID, IS, IL
        REAL(rkind) :: D_loc, k_loc, kD_loc, z_loc, E_loc, Er_loc, JPress_loc
        REAL(rkind) :: Uint, Vint, Urint, Vrint
        REAL(rkind) :: USTOKES_loc(NVRT), VSTOKES_loc(NVRT),UrSTOKES_loc(NVRT), VrSTOKES_loc(NVRT)

        integer     :: IE, isd, j, l, n1, n2, n3, icount
        real(rkind) :: tmp0, tmp1, tmp2, ztmp, ubar, vbar, dhdx, dhdy
        real(rkind) :: STOKES_WVEL_ELEM(NVRT,MNE), ws_tmp1(NVRT,MNS),ws_tmp2(NVRT,MNS)

        real(rkind) :: dr_dxy_loc(2,NVRT,MNS)


!...    Computing Stokes drift horizontal velocities at nodes and
!pressure term
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
              IF (IROLLER == 1) THEN
                Er_loc   = RAC2(IS,ID,IP)*SPSIG(IS)*DDIR*DS_INCR(IS)
                Urint    = Urint + SPSIG(IS)*k_loc*COSTH(ID)*Er_loc
                Vrint    = Vrint + SPSIG(IS)*k_loc*SINTH(ID)*Er_loc
              END IF
            END DO !MDC

            ! Loop on the vertical nodes
            DO IL = KBP(IP), NVRT
              ! Here we need to compute z+h of Eq. C.1 of Bennis et al.
              ! (2011)
              ! In her framework, z varies from -h to eta, meaning that
              ! z+h corresponds to the distance to the bed
              ! -ZETA(KBP(IP),IP) corresponds to h, the depth at node IP
              ! (not the total water depth)
              ! Waves
              z_loc            = ZETA(IL,IP) - ZETA(KBP(IP),IP)
              USTOKES_loc(IL)  = USTOKES_loc(IL)  + Uint*DCOSH(2.D0*k_loc*z_loc)/DSINH(kD_loc)**2
              VSTOKES_loc(IL)  = VSTOKES_loc(IL)  + Vint*DCOSH(2.D0*k_loc*z_loc)/DSINH(kD_loc)**2
              ! Surface rollers
              IF (IROLLER == 1) THEN
                UrSTOKES_loc(IL) = UrSTOKES_loc(IL) + Urint*DCOSH(2.D0*k_loc*z_loc)/DSINH(kD_loc)**2
                VrSTOKES_loc(IL) = VrSTOKES_loc(IL) + Vrint*DCOSH(2.D0*k_loc*z_loc)/DSINH(kD_loc)**2
              END IF
            END DO !NVRT
          END DO !MSC

          ! Notes: by default, the roller contribution to the Stokes
          ! drift velocity has the same vertical
          ! variation as the classic Stokes drift velocities. This is
          ! unrealistic since the roller particles,
          ! and hence the transport, are concentrated around the mean
          ! water level.
          ! Kumar et al. (2012) found that it had a tiny influence on
          ! the results, we hence leave it as it is for now

          ! Storing Stokes drift horizontal velocities and J term
          ! variables
          ! Waves
          STOKES_HVEL(1,:,IP) = USTOKES_loc
          STOKES_HVEL(2,:,IP) = VSTOKES_loc
          ! Surface rollers
          IF (IROLLER == 1) THEN
            ! Smoothing the roller contribution to the Stokes drift
            ! velocity near the shoreline
            ! With this profile, U_st < 10% of computed U_st at h <
            ! DMIN, and U_st > 95% of computed U_st at h > 2.25*DMIN
            IF (D_loc < 3*DMIN) THEN
              ROLLER_STOKES_HVEL(1,:,IP) = UrSTOKES_loc*tanh((0.5D0*D_loc/DMIN)**8.D0)
              ROLLER_STOKES_HVEL(2,:,IP) = VrSTOKES_loc*tanh((0.5D0*D_loc/DMIN)**8.D0)
            ELSE
              ROLLER_STOKES_HVEL(1,:,IP) = UrSTOKES_loc
              ROLLER_STOKES_HVEL(2,:,IP) = VrSTOKES_loc
            END IF
          END IF
          ! Pressure term
          JPRESS(IP) = JPress_loc
        END DO !MNP

!...    Computing Stokes drift horizontal velocities at sides (in pframe
!if ics=2)
        ! The average of the values from vertically adjacent nodes is
        ! taken
        STOKES_HVEL_SIDE = 0.D0; ROLLER_STOKES_HVEL_SIDE = 0.D0
        DO IS = 1,MNS
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
        END DO !MNS

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

!...    Compute horizontal gradient of Stokes x and y-vel. (to compute
!Stokes z-vel.)
        ws_tmp1 = 0.D0; ws_tmp2 = 0.D0
        CALL hgrad_nodes(2,0,NVRT,MNP,MNS,STOKES_HVEL(1,:,:),dr_dxy_loc)
        ws_tmp1(:,:) = dr_dxy_loc(1,:,:) !valid only in resident
        CALL hgrad_nodes(2,0,NVRT,MNP,MNS,STOKES_HVEL(2,:,:),dr_dxy_loc)
        ws_tmp2(:,:) = dr_dxy_loc(2,:,:)

!...    Compute Stokes z-vel. at all levels: STOKES_WVEL_SIDE(NVRT,MNS)
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
        USE schism_glbl, ONLY: kbs, ns, idry_e, isdel, elnode, dldxy, cori, zs, su2, sv2
        USE schism_msgp
        IMPLICIT NONE

        integer     :: IS, IE, k, l, icount
        real(rkind) :: dJ_dx_loc, dJ_dy_loc, du_loc, dv_loc, dz_loc, Ust_loc, Vst_loc
        real(rkind) :: du_dxy(2,NVRT,MNS), dv_dxy(2,NVRT,MNS)

!...    Initialisation
        WWAVE_FORCE = 0.D0

!...    Computing the spatial derivative of horizontal velocities
        CALL hgrad_nodes(2,0,NVRT,MNP,MNS,uu2,du_dxy)
        CALL hgrad_nodes(2,0,NVRT,MNP,MNS,vv2,dv_dxy)

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
            END DO !l

            ! Averaging the values from the two surrounding elements
            IF(icount > 2) CALL parallel_abort('Pressure term:icount>2')
            IF(icount == 2) THEN
              dJ_dx_loc = dJ_dx_loc/2.D0
              dJ_dy_loc = dJ_dy_loc/2.D0
            END IF

          END IF

          !------------------------
          ! Saving wave forces: loop over the vertical
          ! Description of the terms:
          ! 1 - Terms with Coriolis force and the spatial derivative of
          ! horizontal velocities
          ! 2 - Term of the spatial variation de the wave-induced
          ! pressure (J)
          ! 3 - Term -w_s*(du/dz,dv/dz)
          du_loc = 0; dv_loc = 0; dz_loc = 1
          DO k = kbs(IS),NVRT

            IF (fwvor_advz_stokes == 1) THEN ! BM

              ! du/dz and dv/dz terms
              IF (k == kbs(IS)) THEN
                dz_loc = zs(k+1,IS) - zs(k,IS)
                du_loc = su2(k+1,IS) - su2(k,IS)
                dv_loc = sv2(k+1,IS) - sv2(k,IS)
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

            ! Are surface rollers modelled?
            Ust_loc = 0.D0; Vst_loc = 0.D0

            IF (fwvor_advxy_stokes == 1) THEN ! BM

              IF (IROLLER == 1) THEN
                Ust_loc = STOKES_HVEL_SIDE(1,k,IS) + ROLLER_STOKES_HVEL_SIDE(1,k,IS)
                Vst_loc = STOKES_HVEL_SIDE(2,k,IS) + ROLLER_STOKES_HVEL_SIDE(2,k,IS)
              ELSE
                Ust_loc = STOKES_HVEL_SIDE(1,k,IS)
                Vst_loc = STOKES_HVEL_SIDE(2,k,IS)
              END IF

            END IF

            ! Saving wave forces
            WWAVE_FORCE(1,k,IS) = WWAVE_FORCE(1,k,IS) + (cori(IS) + dv_dxy(1,k,IS) - du_dxy(2,k,IS))*Vst_loc   &
                                &                     - dJ_dx_loc &
                                &                     - STOKES_WVEL_SIDE(k,IS)*du_loc/dz_loc
            WWAVE_FORCE(2,k,IS) = WWAVE_FORCE(2,k,IS) - (cori(IS) + dv_dxy(1,k,IS) - du_dxy(2,k,IS))*Ust_loc   &
                                &                     - dJ_dy_loc &
                                &                     - STOKES_WVEL_SIDE(k,IS)*dv_loc/dz_loc
          END DO
        END DO !ns

        ! Exchange between ghost regions
        CALL exchange_s3d_2(WWAVE_FORCE)

      END SUBROUTINE COMPUTE_CONSERVATIVE_VF_TERMS_SCHISM


!**********************************************************************
!*  This routine is used with RADFLAG=VOR (3D vortex formulation, after Bennis et al., 2011)
!*  => Computation of the non-conservative terms due to depth-induced breaking (term Fb from Eq. (11) and (12))
!**********************************************************************
      SUBROUTINE COMPUTE_BREAKING_VF_TERMS_SCHISM
        USE DATAPOOL
        USE schism_glbl, ONLY: hmin_radstress, kbs, ns, isbs, dps, h0, out_wwm, zs
        USE schism_msgp 
        IMPLICIT NONE

        INTEGER     :: IS, isd, k, j, l, n1, n2, n3, icount
        REAL(rkind) :: eta_tmp, tmp0, htot, sum_2D, sum_3D
        REAL(rkind) :: swild_2D(NVRT), swild_3D(NVRT)


        ! Apply lpp_filter
        IF (LPP_FILT_FLAG) CALL LPP_FILT(SBR(1,:))  
        IF (LPP_FILT_FLAG) CALL LPP_FILT(SBR(2,:))

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

            ! Breaking acceleration: average between the two adjacent nodes
            IF (IROLLER == 1) THEN
              WWAVE_FORCE(1,:,IS) = WWAVE_FORCE(1,:,IS) - (SBR(1,n1) + SBR(1,n2) + SROL(1,n1) + SROL(1,n2))/2.D0/htot
              WWAVE_FORCE(2,:,IS) = WWAVE_FORCE(2,:,IS) - (SBR(2,n1) + SBR(2,n2) + SROL(2,n1) + SROL(2,n2))/2.D0/htot
            ELSE
              WWAVE_FORCE(1,:,IS) = WWAVE_FORCE(1,:,IS) - (SBR(1,n1) + SBR(1,n2))/2.D0/htot
              WWAVE_FORCE(2,:,IS) = WWAVE_FORCE(2,:,IS) - (SBR(2,n1) + SBR(2,n2))/2.D0/htot
            END IF

          ELSE !3D

            ! Threshold on Hs
            tmp0 = (out_wwm(n1,1) + out_wwm(n2,1))/2.D0 !Hs

            IF(tmp0 <= 0.05D0) CYCLE
            IF(tmp0/htot < 0.1D0) CYCLE

            ! Vertical distribution function of qdm (due to wave breaking)
            !swild_2D = 0.D0; 
            swild_3D = 0.D0
            DO k = kbs(IS), NVRT

              ! swild_2D(k) = 1.D0
              ! Homogeneous vertical distribution
              IF (ZPROF_BREAK == 1) swild_3D(k) = 1.D0 

              IF (ZPROF_BREAK == 2) swild_3D(k) = cosh((dps(IS)+zs(k,IS))/(0.2D0*tmp0))

              IF (ZPROF_BREAK == 3) swild_3D(k) = 1.D0 - dtanh(((eta_tmp-zs(k,IS))/(0.5D0*tmp0))**2.D0)

              IF (ZPROF_BREAK == 4) swild_3D(k) = 1.D0 - dtanh(((eta_tmp-zs(k,IS))/(0.5D0*tmp0))**4.D0)

              IF (ZPROF_BREAK == 5) swild_3D(k) = 1.D0 - dtanh(((eta_tmp-zs(k,IS))/(0.5D0*tmp0))**8.D0)

              ! All in the two surface layers
              IF (ZPROF_BREAK == 6 .AND. k .GE. NVRT-1) swild_3D(k)=1.D0

            END DO !k
 
            ! Integral of the vertical distribution function
            !sum_2D = 0.0D0
            sum_3D = 0.0D0
            DO k = kbs(IS), NVRT-1
              !sum_2D = sum_2D + (swild_2D(k+1) + swild_2D(k))/2.D0*(zs(k+1,IS) - zs(k,IS))
              sum_3D = sum_3D + (swild_3D(k+1) + swild_3D(k))/2.D0*(zs(k+1,IS) - zs(k,IS))
            END DO !NVRT-1
            !IF(sum_2D*sum_3D == 0) CALL parallel_abort('Vertical profile in wave breaking-induced force: integral=0')
            IF(sum_3D == 0) CALL parallel_abort('Vertical profile in wave breaking-induced force: integral=0')
!'

            DO k = kbs(IS), NVRT
              ! Breaking acceleration
              IF (IROLLER == 1) THEN
                WWAVE_FORCE(1,k,IS) = WWAVE_FORCE(1,k,IS) - swild_3D(k)*(SBR(1,n1) + SBR(1,n2))/2.D0/sum_3D    &
                                                        & - swild_3D(k)*(SROL(1,n1) + SROL(1,n2))/2.D0/sum_3D
                WWAVE_FORCE(2,k,IS) = WWAVE_FORCE(2,k,IS) - swild_3D(k)*(SBR(2,n1) + SBR(2,n2))/2.D0/sum_3D    &
                                                        & - swild_3D(k)*(SROL(2,n1) + SROL(2,n2))/2.D0/sum_3D
              ELSE
                WWAVE_FORCE(1,k,IS) = WWAVE_FORCE(1,k,IS) - swild_3D(k)*(SBR(1,n1) + SBR(1,n2))/2.D0/sum_3D
                WWAVE_FORCE(2,k,IS) = WWAVE_FORCE(2,k,IS) - swild_3D(k)*(SBR(2,n1) + SBR(2,n2))/2.D0/sum_3D
              END IF
            END DO
          END IF !2D/3D

          ! Smoothing wave forces near the shoreline
          ! With this profile, F < 10% of computed F at h < DMIN, and F > 95% of computed F at h > 2.25*DMIN
          IF (htot < 3*DMIN) WWAVE_FORCE(:,:,IS) = WWAVE_FORCE(:,:,IS)*tanh((0.5D0*htot/DMIN)**8.D0)

        END DO !MNS

        ! Exchange between ghost regions
        CALL exchange_s3d_2(WWAVE_FORCE)

      END SUBROUTINE COMPUTE_BREAKING_VF_TERMS_SCHISM

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
!*                                                                    *
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
        USE schism_glbl, only : ns,kbs
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

!        CALL exchange_s3d_2(WWAVE_FORCE)

      END SUBROUTINE APPLY_WAFO_OPBND_RAMP


#endif /*SCHISM*/
!**********************************************************************
!*                                                                    *
!**********************************************************************
