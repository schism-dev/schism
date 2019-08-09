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

MODULE sed2d_morpho
!--------------------------------------------------------------------
! This module contains several subroutines to compute the bottom
! evolution
!
! - base_scheme: Node-centered method where sediment fluxes are
!                considered constant within each element, during
!                a morphological time step
! - weno_scheme: Node-centered WENO-based method
!                (subroutines lin_polyn and det_3x3 are related to
!                this weno_scheme subroutine)
!
! Authors: Thomas Guerin   (thomas.guerin@univ-lr.fr)
! Date: 11/2016
!--------------------------------------------------------------------

  IMPLICIT NONE

CONTAINS

  SUBROUTINE base_scheme(it)
!--------------------------------------------------------------------
! This subroutine computes bottom evolution based on Exner equation
!
! Adapted from sediment_v8.F90 (MORSELFE, L. Pinto) 
!
! Authors: Guillaume Dodet (guillaume.dodet@univ-brest.fr)
!          Thomas Guerin   (thomas.guerin@univ-lr.fr)
!
! History:
! 02/2013 - G.Dodet: - Removed arguments shared in schism_glbl; 
!                    - Removed useless ghost exchange; 
! 03/2013 - G.Dodet: - Corrected argument dimension for solve_jcg;
! 04/2013 - G.Dodet: - Added element-centered method;
!                    - Added information for debugging;
! 07/2013 - G.Dodet: - Flux time integration is now done in 
!                      sed2d_main based on dtsed2d
! 03/2014 - T.Guerin: - Adapted the routine for mutli-class multi-
!                       layer mode: computation of new bathymetry is
!                       done in sed2d_main_mcml for this mode
! 10/2016 - T.Guerin: Modifications related to the merging of single-
!                     class and multi-class routines
!--------------------------------------------------------------------

  USE schism_glbl, ONLY : area,dt,dp,eta2,idry,idry_e,indel,iplg,isdel,     &
                        isidenode,elside,nsa,moitn0,mxitn0,nea,i34,elnode,nne,mnei_p,np,&
                        npa,rkind,rtol0,xctr,xnd,yctr,ynd
  USE schism_msgp, ONLY : exchange_e2d,exchange_p2d,parallel_abort
  USE sed2d_mod, ONLY : bc_flag,bc_val,bed_delta,bed_del_e,h0_sed,   &
                        imeth,mcoefd,poro,qdt_e,nb_class

  IMPLICIT NONE

!- Arguments --------------------------------------------------------
  INTEGER, INTENT(IN) :: it          ! time step
!- Local variables --------------------------------------------------
  INTEGER :: i,iside,j  
  REAL(rkind) :: area_tot,flux,htot,x21,xp,y21,yp,tmp2 
  REAL(rkind), DIMENSION(nea,2) :: qtotdt_e
  REAL(rkind), DIMENSION(nsa,2) :: qtotdt_s
  REAL(rkind), DIMENSION(npa) :: qsan2d
  REAL(rkind), DIMENSION(np) :: tmp
!--------------------------------------------------------------------

!- Compute bed change ----------------------------------------------
!-
!- Integrate sediment transport (m^2) along node-centered control volume 
!- and solve equation system with jcg solver to obtain bed delta (m).
!-
!- Here the sediment transport is considered constant inside each element.
!--------------------------------------------------------------------

  qsan2d = 0.d0
  DO i = 1,nea
     htot = 0.d0
     DO j = 1,i34(i)
        htot = htot+(eta2(elnode(j,i))+dp(elnode(j,i)))/i34(i)
     ENDDO !j

     IF(idry_e(i) == 1.OR.htot.LE.h0_sed) CYCLE

!Error: YJZ - need to account for quads
     xp  = (xnd(elnode(2,i))+xnd(elnode(3,i)))/2.d0
     yp  = (ynd(elnode(2,i))+ynd(elnode(3,i)))/2.d0
     flux= qdt_e(i,1)*(yp-yctr(i))-qdt_e(i,2)*(xp-xctr(i))
     qsan2d(elnode(2,i)) = qsan2d(elnode(2,i))-flux
     qsan2d(elnode(3,i)) = qsan2d(elnode(3,i))+flux

     xp  = (xnd(elnode(3,i))+xnd(elnode(1,i)))/2.d0
     yp  = (ynd(elnode(3,i))+ynd(elnode(1,i)))/2.d0
     flux= qdt_e(i,1)*(yp-yctr(i))-qdt_e(i,2)*(xp-xctr(i))
     qsan2d(elnode(3,i)) = qsan2d(elnode(3,i))-flux
     qsan2d(elnode(1,i)) = qsan2d(elnode(1,i))+flux

     xp  = (xnd(elnode(1,i))+xnd(elnode(2,i)))/2.d0
     yp  = (ynd(elnode(1,i))+ynd(elnode(2,i)))/2.d0
     flux= qdt_e(i,1)*(yp-yctr(i))-qdt_e(i,2)*(xp-xctr(i))
     qsan2d(elnode(1,i)) = qsan2d(elnode(1,i))-flux
     qsan2d(elnode(2,i)) = qsan2d(elnode(2,i))+flux
  ENDDO !nea
  CALL exchange_p2d(qsan2d)

  tmp = 0.d0
  DO i = 1,np
     tmp(i) = qsan2d(i)/(1.d0-poro)
  ENDDO !np

!- Solve Exner equation ---------------------------------------------
  bed_delta = 0.d0 
  CALL solve_jcg(mnei_p,np,npa,it,moitn0,mxitn0,rtol0,mcoefd,bed_delta,tmp,bc_val,bc_flag)
  CALL exchange_p2d(bed_delta)

  END SUBROUTINE base_scheme

  SUBROUTINE weno_scheme(it)
!--------------------------------------------------------------------
! This subroutine computes bottom evolution using a WENO-base scheme
! to solve the Exner equation
!
! Authors: Thomas Guerin   (thomas.guerin@univ-lr.fr)
! Date: 11/2016
!--------------------------------------------------------------------

  USE schism_glbl, ONLY : area,dp,dt,elnode,eta2,idry,idry_e,indel,mnei_p,&
                         &moitn0,mxitn0,nea,nne,np,npa,rkind,rtol0,xctr,xnd,yctr,ynd
  USE schism_msgp, ONLY : exchange_p2d
  USE sed2d_mod, ONLY : bc_flag,bc_val,bed_delta,h0_sed,mcoefd,morfac,poro,qdt_e,vc_area

  IMPLICIT NONE

!- Arguments --------------------------------------------------------
  INTEGER, INTENT(IN) :: it
!- Parameters -------------------------------------------------------
  INTEGER, PARAMETER :: nb_comb = 8
  REAL(rkind), PARAMETER :: ONETHIRD = 1.d0/3.d0
  REAL(rkind), PARAMETER :: c = 0.5d0 + DSQRT(3.d0)/6.d0
  REAL(rkind), PARAMETER :: epsi = 1.e-10
  REAL(rkind), PARAMETER :: rp = 1.d0
  REAL(rkind), PARAMETER :: beta_comp = 2.d0
!- Local variables --------------------------------------------------
  INTEGER :: el1,el2,el3,i,j,k,m,nd1,nd2,nd3,p
  REAL(rkind) :: alpha,flux,h1,h2,h3,htot,OItmp,OIx,OIxtmp,OIy,OIytmp,          &
                 P1,P1G1,P1G2,P1nd,P2,P2G1,P2G2,P2nd,Px1G1,Px1G2,Px1nd,Px2G1,Px2G2,Px2nd,  &
                 Py1G1,Py1G2,Py1nd,Py2G1,Py2G2,Py2nd,r_flux,w_sum,xG1,xG2,xi,xnd1,xnd2, &
                 xnd3,xp,yG1,yG2,yi,ynd1,ynd2,ynd3,yp
  REAL(rkind), DIMENSION(npa) :: fluxtot
  REAL(rkind), DIMENSION(np) :: fluxtmp
  REAL(rkind), DIMENSION(npa,56,3) :: Cx,Cy
  REAL(rkind), DIMENSION(3) :: Cxtmp,Cytmp
  REAL(rkind), DIMENSION(56) :: OI
  REAL(rkind), DIMENSION(npa,nb_comb) :: w
!--------------------------------------------------------------------


!- Initialization
  Cx(:,:,:) = 0.d0; Cy(:,:,:) = 0.d0
  w(:,:) = 0.d0
  DO i=1,np
     IF(idry(i)==1) CYCLE

     IF(nne(i)<3) CYCLE

     !Compute linear polynom for all stencils
     m=1
     OI(:) = 0.d0

     DO j=1,nne(i)
        el1 = indel(j,i)
        IF(el1<1) CYCLE

        IF (j<=nne(i)-2) THEN
          el2 = indel(j+1,i)
          el3 = indel(j+2,i)

        ELSEIF (j==nne(i)-1) THEN
          el2 = indel(j+1,i)
          el3 = indel(1,i)

        ELSEIF (j==nne(i)) THEN
          el2 = indel(1,i)
          el3 = indel(2,i)
        ENDIF

        IF((el2<1).or.(el3<1)) CYCLE

        CALL lin_polyn(xctr(el1),yctr(el1),xctr(el2),yctr(el2),&
            &xctr(el3),yctr(el3),qdt_e(el1,:),qdt_e(el2,:),  &
            &qdt_e(el3,:),Cx(i,m,:),Cy(i,m,:))

        !Compute oscillation indicator
        OIx = DSQRT(vc_area(i)*(Cx(i,m,1)**2.d0+Cx(i,m,2)**2.d0)&
                   &/((area(el1)+area(el2)+area(el3))/3.d0))

        OIy = DSQRT(vc_area(i)*(Cy(i,m,1)**2.d0+Cy(i,m,2)**2.d0)&
                   &/((area(el1)+area(el2)+area(el3))/3.d0))

        OI(m) = OIx+OIy

        !Sorting oscillation indicator
        IF (m>1) THEN
          DO p=1,m-1
             IF (OI(m)<OI(p)) THEN
               OItmp = OI(p)
               Cxtmp(:) = Cx(i,p,:); Cytmp(:) = Cy(i,p,:)

               OI(p) = OI(m)
               Cx(i,p,:) = Cx(i,m,:); Cy(i,p,:) = Cy(i,m,:)

               OI(m) = OItmp
               Cx(i,m,:) = Cxtmp(:); Cy(i,m,:) = Cytmp(:)
             ENDIF
          ENDDO
        ENDIF

        m=m+1

     ENDDO!j

     w_sum = 0.d0
     IF (nne(i)<=nb_comb) THEN
       DO j=1,nne(i)
          w_sum = w_sum + 1.d0/((epsi+OI(j))**rp)
       ENDDO
     ELSE
       DO j=1,nb_comb
          w_sum = w_sum + 1.d0/((epsi+OI(j))**rp)
       ENDDO
     ENDIF

     IF (nne(i)<=nb_comb) THEN
       DO j=1,nne(i)
          w(i,j) = 1.d0/((epsi+OI(j))**rp)/w_sum
       ENDDO
     ELSE
       DO j=1,nb_comb
          w(i,j) = 1.d0/((epsi+OI(j))**rp)/w_sum
       ENDDO
     ENDIF

  ENDDO !i=1,np

  DO j=1,nb_comb
     CALL exchange_p2d(w(:,j))
     DO k=1,3
        CALL exchange_p2d(Cx(:,j,k)); CALL exchange_p2d(Cy(:,j,k))
     ENDDO
  ENDDO

  fluxtot = 0.d0
  DO i = 1,nea
     xi = xctr(i); yi = yctr(i)
     nd1 = elnode(1,i); nd2 = elnode(2,i); nd3 = elnode(3,i)
     xnd1 = xnd(nd1); xnd2 = xnd(nd2); xnd3 = xnd(nd3)
     ynd1 = ynd(nd1); ynd2 = ynd(nd2); ynd3 = ynd(nd3)

     !Compute total water depth at the 3 nodes
     h1 = dp(nd1)+eta2(nd1)
     h2 = dp(nd2)+eta2(nd2)
     h3 = dp(nd3)+eta2(nd3)
     htot = (h1+h2+h3)*ONETHIRD

     IF((idry_e(i)==1).OR.(htot.LE.h0_sed)) CYCLE

     !1st segment
     xp = 0.5d0*(xnd1+xnd2)
     yp = 0.5d0*(ynd1+ynd2)
     xG1 = c*xp+(1.d0-c)*xi; yG1 = c*yp+(1.d0-c)*yi
     xG2 = c*xi+(1.d0-c)*xp; yG2 = c*yi+(1.d0-c)*yp

     Px1G1 = 0.d0; Py1G1 = 0.d0; Px2G1 = 0.d0; Py2G1 = 0.d0
     Px1G2 = 0.d0; Py1G2 = 0.d0; Px2G2 = 0.d0; Py2G2 = 0.d0
     Px1nd = 0.d0; Py1nd = 0.d0; Px2nd = 0.d0; Py2nd = 0.d0
     IF (nne(nd1)>3) THEN
       DO j=1,nb_comb
          Px1G1 = Px1G1 + w(nd1,j)*(Cx(nd1,j,1)*xG1+Cx(nd1,j,2)*yG1+  &
                                    Cx(nd1,j,3))
          Py1G1 = Py1G1 + w(nd1,j)*(Cy(nd1,j,1)*xG1+Cy(nd1,j,2)*yG1+  &
                                    Cy(nd1,j,3))
          Px1G2 = Px1G2 + w(nd1,j)*(Cx(nd1,j,1)*xG2+Cx(nd1,j,2)*yG2+  &
                                    Cx(nd1,j,3))
          Py1G2 = Py1G2 + w(nd1,j)*(Cy(nd1,j,1)*xG2+Cy(nd1,j,2)*yG2+  &
                                    Cy(nd1,j,3))
          Px1nd = Px1nd + w(nd1,j)*(Cx(nd1,j,1)*xnd1+Cx(nd1,j,2)*ynd1+&
                                    Cx(nd1,j,3))
          Py1nd = Py1nd + w(nd1,j)*(Cy(nd1,j,1)*xnd1+Cy(nd1,j,2)*ynd1+&
                                    Cy(nd1,j,3))
       ENDDO
     ELSEIF (nne(nd1)==3) THEN
       Px1G1 = Cx(nd1,1,1)*xG1+Cx(nd1,1,2)*yG1+Cx(nd1,1,3)
       Py1G1 = Cy(nd1,1,1)*xG1+Cy(nd1,1,2)*yG1+Cy(nd1,1,3)
       Px1G2 = Cx(nd1,1,1)*xG2+Cx(nd1,1,2)*yG2+Cx(nd1,1,3)
       Py1G2 = Cy(nd1,1,1)*xG2+Cy(nd1,1,2)*yG2+Cy(nd1,1,3)
       Px1nd = Cx(nd1,1,1)*xnd1+Cx(nd1,1,2)*ynd1+Cx(nd1,1,3)
       Py1nd = Cy(nd1,1,1)*xnd1+Cy(nd1,1,2)*ynd1+Cy(nd1,1,3)
     ELSEIF (nne(nd1)<3) THEN
       Px1G1 = qdt_e(i,1); Py1G1 = qdt_e(i,2)
       Px1G2 = qdt_e(i,1); Py1G2 = qdt_e(i,2)
       Px1nd = qdt_e(i,1); Py1nd = qdt_e(i,2)
     ENDIF

     IF (nne(nd2)>3) THEN
       DO j=1,nb_comb
          Px2G1 = Px2G1 + w(nd2,j)*(Cx(nd2,j,1)*xG1+Cx(nd2,j,2)*yG1+  &
                                    Cx(nd2,j,3))
          Py2G1 = Py2G1 + w(nd2,j)*(Cy(nd2,j,1)*xG1+Cy(nd2,j,2)*yG1+  &
                                    Cy(nd2,j,3))
          Px2G2 = Px2G2 + w(nd2,j)*(Cx(nd2,j,1)*xG2+Cx(nd2,j,2)*yG2+  &
                                    Cx(nd2,j,3))
          Py2G2 = Py2G2 + w(nd2,j)*(Cy(nd2,j,1)*xG2+Cy(nd2,j,2)*yG2+  &
                                    Cy(nd2,j,3))
          Px2nd = Px2nd + w(nd2,j)*(Cx(nd2,j,1)*xnd2+Cx(nd2,j,2)*ynd2+&
                                    Cx(nd2,j,3))
          Py2nd = Py2nd + w(nd2,j)*(Cy(nd2,j,1)*xnd2+Cy(nd2,j,2)*ynd2+&
                                    Cy(nd2,j,3))
       ENDDO
     ELSEIF (nne(nd2)==3) THEN
       Px2G1 = Cx(nd2,1,1)*xG1+Cx(nd2,1,2)*yG1+Cx(nd2,1,3)
       Py2G1 = Cy(nd2,1,1)*xG1+Cy(nd2,1,2)*yG1+Cy(nd2,1,3)
       Px2G2 = Cx(nd2,1,1)*xG2+Cx(nd2,1,2)*yG2+Cx(nd2,1,3)
       Py2G2 = Cy(nd2,1,1)*xG2+Cy(nd2,1,2)*yG2+Cy(nd2,1,3)
       Px2nd = Cx(nd2,1,1)*xnd2+Cx(nd2,1,2)*ynd2+Cx(nd2,1,3)
       Py2nd = Cy(nd2,1,1)*xnd2+Cy(nd2,1,2)*ynd2+Cy(nd2,1,3)
     ELSEIF (nne(nd2)<3) THEN
       Px2G1 = qdt_e(i,1); Py2G1 = qdt_e(i,2)
       Px2G2 = qdt_e(i,1); Py2G2 = qdt_e(i,2)
       Px2nd = qdt_e(i,1); Py2nd = qdt_e(i,2)
     ENDIF

     P1G1 = Px1G1*(yi-yp)-Py1G1*(xi-xp)
     P1G2 = Px1G2*(yi-yp)-Py1G2*(xi-xp)
     P1 = 0.5d0*(P1G1+P1G2)
     P2G1 = Px2G1*(yi-yp)-Py2G1*(xi-xp)
     P2G2 = Px2G2*(yi-yp)-Py2G2*(xi-xp)
     P2 = 0.5d0*(P2G1+P2G2)
     P1nd = Px1nd*(yi-yp)-Py1nd*(xi-xp)
     P2nd = Px2nd*(yi-yp)-Py2nd*(xi-xp)

     !compute the r_flux parameter (quantifying the local upwinding applied)
     IF (0.5d0*(h1+h2)>0.d0) THEN
       r_flux = abs(h2-h1)/(0.5d0*(h1+h2))
     ELSE
       r_flux = 0.d0
     ENDIF

     !compute weight parameter for numerical flux
     alpha = MAX(0.d0,MIN(r_flux,beta_comp))

     !compute the numerical flux
     IF (h2.LE.h1) THEN
       IF (P1nd<=P2nd) THEN
         flux = P1+0.5d0*alpha*(P1nd-P1)
       ELSE
         flux = P2+0.5d0*alpha*(P2nd-P2)
       ENDIF
     ELSE
       IF (P1nd>=P2nd) THEN
         flux = P1+0.5d0*alpha*(P1nd-P1)
       ELSE
         flux = P2+0.5d0*alpha*(P2nd-P2)
       ENDIF
     ENDIF

     fluxtot(nd1) = fluxtot(nd1) + flux
     fluxtot(nd2) = fluxtot(nd2) - flux

     !2nd segment
     xp = 0.5d0*(xnd(nd2)+xnd(nd3))
     yp = 0.5d0*(ynd(nd2)+ynd(nd3))
     xG1 = c*xp+(1.d0-c)*xi; yG1 = c*yp+(1.d0-c)*yi
     xG2 = c*xi+(1.d0-c)*xp; yG2 = c*yi+(1.d0-c)*yp

     Px1G1 = 0.d0; Py1G1 = 0.d0; Px2G1 = 0.d0; Py2G1 = 0.d0
     Px1G2 = 0.d0; Py1G2 = 0.d0; Px2G2 = 0.d0; Py2G2 = 0.d0
     Px1nd = 0.d0; Py1nd = 0.d0; Px2nd = 0.d0; Py2nd = 0.d0
     IF (nne(nd2)>3) THEN
       DO j=1,nb_comb
          Px1G1 = Px1G1 + w(nd2,j)*(Cx(nd2,j,1)*xG1+Cx(nd2,j,2)*yG1+  &
                                    Cx(nd2,j,3))
          Py1G1 = Py1G1 + w(nd2,j)*(Cy(nd2,j,1)*xG1+Cy(nd2,j,2)*yG1+  &
                                    Cy(nd2,j,3))
          Px1G2 = Px1G2 + w(nd2,j)*(Cx(nd2,j,1)*xG2+Cx(nd2,j,2)*yG2+  &
                                    Cx(nd2,j,3))
          Py1G2 = Py1G2 + w(nd2,j)*(Cy(nd2,j,1)*xG2+Cy(nd2,j,2)*yG2+  &
                                    Cy(nd2,j,3))
          Px1nd = Px1nd + w(nd2,j)*(Cx(nd2,j,1)*xnd2+Cx(nd2,j,2)*ynd2+&
                                    Cx(nd2,j,3))
          Py1nd = Py1nd + w(nd2,j)*(Cy(nd2,j,1)*xnd2+Cy(nd2,j,2)*ynd2+&
                                    Cy(nd2,j,3))
       ENDDO
     ELSEIF (nne(nd2)==3) THEN
       Px1G1 = Cx(nd2,1,1)*xG1+Cx(nd2,1,2)*yG1+Cx(nd2,1,3)
       Py1G1 = Cy(nd2,1,1)*xG1+Cy(nd2,1,2)*yG1+Cy(nd2,1,3)
       Px1G2 = Cx(nd2,1,1)*xG2+Cx(nd2,1,2)*yG2+Cx(nd2,1,3)
       Py1G2 = Cy(nd2,1,1)*xG2+Cy(nd2,1,2)*yG2+Cy(nd2,1,3)
       Px1nd = Cx(nd2,1,1)*xnd2+Cx(nd2,1,2)*ynd2+Cx(nd2,1,3)
       Py1nd = Cy(nd2,1,1)*xnd2+Cy(nd2,1,2)*ynd2+Cy(nd2,1,3)
     ELSEIF (nne(nd2)<3) THEN
       Px1G1 = qdt_e(i,1); Py1G1 = qdt_e(i,2)
       Px1G2 = qdt_e(i,1); Py1G2 = qdt_e(i,2)
       Px1nd = qdt_e(i,1); Py1nd = qdt_e(i,2)
     ENDIF

     IF (nne(nd3)>3) THEN
       DO j=1,nb_comb
          Px2G1 = Px2G1 + w(nd3,j)*(Cx(nd3,j,1)*xG1+Cx(nd3,j,2)*yG1+  &
                                     Cx(nd3,j,3))
          Py2G1 = Py2G1 + w(nd3,j)*(Cy(nd3,j,1)*xG1+Cy(nd3,j,2)*yG1+  &
                                    Cy(nd3,j,3))
          Px2G2 = Px2G2 + w(nd3,j)*(Cx(nd3,j,1)*xG2+Cx(nd3,j,2)*yG2+  &
                                    Cx(nd3,j,3))
          Py2G2 = Py2G2 + w(nd3,j)*(Cy(nd3,j,1)*xG2+Cy(nd3,j,2)*yG2+  &
                                    Cy(nd3,j,3))
          Px2nd = Px2nd + w(nd3,j)*(Cx(nd3,j,1)*xnd3+Cx(nd3,j,2)*ynd3+&
                                    Cx(nd3,j,3))
          Py2nd = Py2nd + w(nd3,j)*(Cy(nd3,j,1)*xnd3+Cy(nd3,j,2)*ynd3+&
                                    Cy(nd3,j,3))
       ENDDO
     ELSEIF (nne(nd3)==3) THEN
       Px2G1 = Cx(nd3,1,1)*xG1+Cx(nd3,1,2)*yG1+Cx(nd3,1,3)
       Py2G1 = Cy(nd3,1,1)*xG1+Cy(nd3,1,2)*yG1+Cy(nd3,1,3)
       Px2G2 = Cx(nd3,1,1)*xG2+Cx(nd3,1,2)*yG2+Cx(nd3,1,3)
       Py2G2 = Cy(nd3,1,1)*xG2+Cy(nd3,1,2)*yG2+Cy(nd3,1,3)
       Px2nd = Cx(nd3,1,1)*xnd3+Cx(nd3,1,2)*ynd3+Cx(nd3,1,3)
       Py2nd = Cy(nd3,1,1)*xnd3+Cy(nd3,1,2)*ynd3+Cy(nd3,1,3)
     ELSEIF (nne(nd3)<3) THEN
       Px2G1 = qdt_e(i,1); Py2G1 = qdt_e(i,2)
       Px2G2 = qdt_e(i,1); Py2G2 = qdt_e(i,2)
       Px2nd = qdt_e(i,1); Py2nd = qdt_e(i,2)
     ENDIF

     P1G1 = Px1G1*(yi-yp)-Py1G1*(xi-xp)
     P1G2 = Px1G2*(yi-yp)-Py1G2*(xi-xp)
     P1 = 0.5d0*(P1G1+P1G2)
     P2G1 = Px2G1*(yi-yp)-Py2G1*(xi-xp)
     P2G2 = Px2G2*(yi-yp)-Py2G2*(xi-xp)
     P2 = 0.5d0*(P2G1+P2G2)
     P1nd = Px1nd*(yi-yp)-Py1nd*(xi-xp)
     P2nd = Px2nd*(yi-yp)-Py2nd*(xi-xp)

     !compute the r_flux parameter (quantifying the local upwinding applied)
     IF (0.5d0*(h2+h3)>0.d0) THEN
       r_flux = abs(h3-h2)/(0.5d0*(h2+h3))
     ELSE
       r_flux = 0.d0
     ENDIF

     !compute weight parameter for numerical flux
     alpha = MAX(0.d0,MIN(r_flux,beta_comp))

     !compute the numerical flux
     IF (h3.LE.h2) THEN
       IF (P1nd<=P2nd) THEN
         flux = P1+0.5d0*alpha*(P1nd-P1)
       ELSE
         flux = P2+0.5d0*alpha*(P2nd-P2)
       ENDIF
     ELSE
       IF (P1nd>=P2nd) THEN
         flux = P1+0.5d0*alpha*(P1nd-P1)
       ELSE
         flux = P2+0.5d0*alpha*(P2nd-P2)
       ENDIF
     ENDIF

     fluxtot(nd2) = fluxtot(nd2) + flux
     fluxtot(nd3) = fluxtot(nd3) - flux

     !3rd segment
     xp = 0.5d0*(xnd(nd3)+xnd(nd1))
     yp = 0.5d0*(ynd(nd3)+ynd(nd1))
     xG1 = c*xp+(1.d0-c)*xi; yG1 = c*yp+(1.d0-c)*yi
     xG2 = c*xi+(1.d0-c)*xp; yG2 = c*yi+(1.d0-c)*yp

     Px1G1 = 0.d0; Py1G1 = 0.d0; Px2G1 = 0.d0; Py2G1 = 0.d0
     Px1G2 = 0.d0; Py1G2 = 0.d0; Px2G2 = 0.d0; Py2G2 = 0.d0
     Px1nd = 0.d0; Py1nd = 0.d0; Px2nd = 0.d0; Py2nd = 0.d0
     IF (nne(nd3)>3) THEN
       DO j=1,nb_comb
          Px1G1 = Px1G1 + w(nd3,j)*(Cx(nd3,j,1)*xG1+Cx(nd3,j,2)*yG1+  &
                                    Cx(nd3,j,3))
          Py1G1 = Py1G1 + w(nd3,j)*(Cy(nd3,j,1)*xG1+Cy(nd3,j,2)*yG1+  &
                                    Cy(nd3,j,3))
          Px1G2 = Px1G2 + w(nd3,j)*(Cx(nd3,j,1)*xG2+Cx(nd3,j,2)*yG2+  &
                                    Cx(nd3,j,3))
          Py1G2 = Py1G2 + w(nd3,j)*(Cy(nd3,j,1)*xG2+Cy(nd3,j,2)*yG2+  &
                                    Cy(nd3,j,3))
          Px1nd = Px1nd + w(nd3,j)*(Cx(nd3,j,1)*xnd3+Cx(nd3,j,2)*ynd3+&
                                    Cx(nd3,j,3))
          Py1nd = Py1nd + w(nd3,j)*(Cy(nd3,j,1)*xnd3+Cy(nd3,j,2)*ynd3+&
                                    Cy(nd3,j,3))
       ENDDO
     ELSEIF (nne(nd3)==3) THEN
       Px1G1 = Cx(nd3,1,1)*xG1+Cx(nd3,1,2)*yG1+Cx(nd3,1,3)
       Py1G1 = Cy(nd3,1,1)*xG1+Cy(nd3,1,2)*yG1+Cy(nd3,1,3)
       Px1G2 = Cx(nd3,1,1)*xG2+Cx(nd3,1,2)*yG2+Cx(nd3,1,3)
       Py1G2 = Cy(nd3,1,1)*xG2+Cy(nd3,1,2)*yG2+Cy(nd3,1,3)
       Px1nd = Cx(nd3,1,1)*xnd3+Cx(nd3,1,2)*ynd3+Cx(nd3,1,3)
       Py1nd = Cy(nd3,1,1)*xnd3+Cy(nd3,1,2)*ynd3+Cy(nd3,1,3)
     ELSEIF (nne(nd3)<3) THEN
       Px1G1 = qdt_e(i,1); Py1G1 = qdt_e(i,2)
       Px1G2 = qdt_e(i,1); Py1G2 = qdt_e(i,2)
       Px1nd = qdt_e(i,1); Py1nd = qdt_e(i,2)
     ENDIF

     IF (nne(nd1)>3) THEN
       DO j=1,nb_comb
          Px2G1 = Px2G1 + w(nd1,j)*(Cx(nd1,j,1)*xG1+Cx(nd1,j,2)*yG1+  &
                                    Cx(nd1,j,3))
          Py2G1 = Py2G1 + w(nd1,j)*(Cy(nd1,j,1)*xG1+Cy(nd1,j,2)*yG1+  &
                                    Cy(nd1,j,3))
          Px2G2 = Px2G2 + w(nd1,j)*(Cx(nd1,j,1)*xG2+Cx(nd1,j,2)*yG2+  &
                                    Cx(nd1,j,3))
          Py2G2 = Py2G2 + w(nd1,j)*(Cy(nd1,j,1)*xG2+Cy(nd1,j,2)*yG2+  &
                                    Cy(nd1,j,3))
          Px2nd = Px2nd + w(nd1,j)*(Cx(nd1,j,1)*xnd1+Cx(nd1,j,2)*ynd1+&
                                    Cx(nd1,j,3))
          Py2nd = Py2nd + w(nd1,j)*(Cy(nd1,j,1)*xnd1+Cy(nd1,j,2)*ynd1+&
                                    Cy(nd1,j,3))
       ENDDO
     ELSEIF (nne(nd1)==3) THEN
       Px2G1 = Cx(nd1,1,1)*xG1+Cx(nd1,1,2)*yG1+Cx(nd1,1,3)
       Py2G1 = Cy(nd1,1,1)*xG1+Cy(nd1,1,2)*yG1+Cy(nd1,1,3)
       Px2G2 = Cx(nd1,1,1)*xG2+Cx(nd1,1,2)*yG2+Cx(nd1,1,3)
       Py2G2 = Cy(nd1,1,1)*xG2+Cy(nd1,1,2)*yG2+Cy(nd1,1,3)
       Px2nd = Cx(nd1,1,1)*xnd1+Cx(nd1,1,2)*ynd1+Cx(nd1,1,3)
       Py2nd = Cy(nd1,1,1)*xnd1+Cy(nd1,1,2)*ynd1+Cy(nd1,1,3)
     ELSEIF (nne(nd1)<3) THEN
       Px2G1 = qdt_e(i,1); Py2G1 = qdt_e(i,2)
       Px2G2 = qdt_e(i,1); Py2G2 = qdt_e(i,2)
       Px2nd = qdt_e(i,1); Py2nd = qdt_e(i,2)
     ENDIF

     P1G1 = Px1G1*(yi-yp)-Py1G1*(xi-xp)
     P1G2 = Px1G2*(yi-yp)-Py1G2*(xi-xp)
     P1 = 0.5d0*(P1G1+P1G2)
     P2G1 = Px2G1*(yi-yp)-Py2G1*(xi-xp)
     P2G2 = Px2G2*(yi-yp)-Py2G2*(xi-xp)
     P2 = 0.5d0*(P2G1+P2G2)
     P1nd = Px1nd*(yi-yp)-Py1nd*(xi-xp)
     P2nd = Px2nd*(yi-yp)-Py2nd*(xi-xp)

     !compute the r_flux parameter (quantifying the local upwinding applied)
     IF (0.5d0*(h3+h1)>0.d0) THEN
       r_flux = abs(h1-h3)/(0.5d0*(h3+h1))
     ELSE
       r_flux = 0.d0
     ENDIF

     !compute weight parameter for numerical flux
     alpha = MAX(0.d0,MIN(r_flux,beta_comp))

     !compute the numerical flux
     IF (h1.LE.h3) THEN
       IF (P1nd<=P2nd) THEN
         flux = P1+0.5d0*alpha*(P1nd-P1)
       ELSE
         flux = P2+0.5d0*alpha*(P2nd-P2)
       ENDIF
     ELSE
       IF (P1nd>=P2nd) THEN
         flux = P1+0.5d0*alpha*(P1nd-P1)
       ELSE
         flux = P2+0.5d0*alpha*(P2nd-P2)
       ENDIF
     ENDIF

     fluxtot(nd3) = fluxtot(nd3) + flux
     fluxtot(nd1) = fluxtot(nd1) - flux

  ENDDO !nea

  CALL exchange_p2d(fluxtot)

  fluxtmp = 0.d0
  DO i=1,np
     fluxtmp(i) = fluxtot(i)/(1.d0-poro)
  ENDDO

!- Compute bed change at nodes
  bed_delta = 0.d0
  CALL solve_jcg(mnei_p,np,npa,it,moitn0,mxitn0,rtol0,mcoefd,bed_delta,fluxtmp,bc_val,bc_flag)
  CALL exchange_p2d(bed_delta)

  END SUBROUTINE weno_scheme

  SUBROUTINE lin_polyn(x1,y1,x2,y2,x3,y3,q1,q2,q3,Cx,Cy)
!--------------------------------------------------------------------
! Compute coefficients for bilinear interpolation of vectors
!--------------------------------------------------------------------
  USE schism_glbl, only : rkind
  USE schism_msgp, ONLY : parallel_abort

  IMPLICIT NONE

!- Arguments --------------------------------------------------------
  REAL(rkind), INTENT(IN) :: x1,y1,x2,y2,x3,y3
  REAL(rkind), DIMENSION(2), INTENT(IN) :: q1,q2,q3
  REAL(rkind), DIMENSION(3), INTENT(OUT) :: Cx,Cy
!- Local variables --------------------------------------------------
  REAL(rkind) :: detA
  REAL(rkind), DIMENSION(2) :: det1,det2,det3
!--------------------------------------------------------------------

  CALL det_3x3(x1,y1,1.d0,x2,y2,1.d0,x3,y3,1.d0,detA)

  CALL det_3x3(q1(1),y1,1.d0,q2(1),y2,1.d0,q3(1),y3,1.d0,det1(1))
  CALL det_3x3(q1(2),y1,1.d0,q2(2),y2,1.d0,q3(2),y3,1.d0,det1(2))

  CALL det_3x3(x1,q1(1),1.d0,x2,q2(1),1.d0,x3,q3(1),1.d0,det2(1))
  CALL det_3x3(x1,q1(2),1.d0,x2,q2(2),1.d0,x3,q3(2),1.d0,det2(2))

  CALL det_3x3(x1,y1,q1(1),x2,y2,q2(1),x3,y3,q3(1),det3(1))
  CALL det_3x3(x1,y1,q1(2),x2,y2,q2(2),x3,y3,q3(2),det3(2))

  IF (detA /= 0.d0) THEN
    Cx(1) = det1(1)/detA
    Cy(1) = det1(2)/detA
    Cx(2) = det2(1)/detA
    Cy(2) = det2(2)/detA
    Cx(3) = det3(1)/detA
    Cy(3) = det3(2)/detA
  ELSE
    CALL parallel_abort('detA=0 in lin_polyn subroutine')
  ENDIF

  END SUBROUTINE lin_polyn

  SUBROUTINE det_3x3(a,b,c,d,e,f,g,h,i,det)
!--------------------------------------------------------------------
! Compute determinant of a 3x3 matrice
!--------------------------------------------------------------------
  USE schism_glbl, only : rkind

  IMPLICIT NONE

!- Arguments --------------------------------------------------------
  REAL(rkind), INTENT(IN) :: a,b,c,d,e,f,g,h,i
  REAL(rkind), INTENT(OUT) :: det
!--------------------------------------------------------------------

  det = a*(e*i-f*h)-b*(d*i-f*g)+c*(d*h-e*g)

  END SUBROUTINE det_3x3
!--------------------------------------------------------------------
END MODULE sed2d_morpho
