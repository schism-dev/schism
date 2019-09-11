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

MODULE sed2d_filter
!--------------------------------------------------------------------
! This module contains several subroutines to filter the bathymetry:
! - sed2d_filter_extrema
! - sed2d_filter_slope
! - sed2d_filter_diffu
!
! Author: guillaume dodet (gdodet01@univ-lr.fr, gdodet@lnec.pt)
! Date: 18/01/2013 
!
! History:
! 01/2013 - G.Dodet: added sed2d_filter_main
! 02/2013 - G.Dodet: removed sed2d_filter_main (useless) 
! 03/2013 - G.Dodet: solved a bug in sed2d_filter_extrema
!--------------------------------------------------------------------

  IMPLICIT NONE

CONTAINS

!--------------------------------------------------------------------
  SUBROUTINE sed2d_filter_extrema(dp1,dp2)
!--------------------------------------------------------------------
! This subroutine filters iteratively nodes whose depth presents a 
! local extremum (with respect to the surrounding nodes).
!
! Adapted from filternl.f (SAND2D, A. Fortunato)
!
! Two options are available (str_filt flag): 
! - the strong filter will compute the average with node with the 
!   largest amplitude (this is default)
! - the weak filter will compute the average with node with the
!   smallest amplitude (set str_filt to 2 in this routine) 
!
!                     i                    i is a local maximum 
!                     +                           
!                    / \                   strong will use i+1
!              +----+   \     +-----+
!                  i-1   \   /             weak will use i-1
!                         \ /
!                          +               to compute an average 
!                          i+1             
!                                          value
!
!                                                                  
! Author: guillaume dodet (gdodet01@univ-lr.fr, gdodet@lnec.pt)    
! Date:   20/11/2012
!
! History: 
! 01/2013 - G.Dodet: Included this subroutine in module sed2d_filter 
! 03/2013 - G.Dodet: Solved a bug in the computation of the node-
!                    associated area (formerly area of element i was
!                    used as the area of node i... nonsense!)
!                    Removed STR_FILT from input file and included it
!                    as a local variable
! 04/2013 - G.Dodet: - Solved a bug in the computation of area (can
!                      be negative if outside augmented domain). The 
!                      volume control area is now computed in 
!                      sed2d_init;
!                    - Exchange information between processors after
!                      each iteration;
!                    - Synchronized processors and mpi_reduced some
!                      variables for consistency;
!                    - Added debug options; 
!                    - Removed SAVE statement (useless?);
!--------------------------------------------------------------------

  USE schism_glbl, ONLY : area,iplg,indel,indnd,nne,nnp,np,npa,rkind
  USE schism_msgp, ONLY : comm,exchange_p2d,ierr,itype,myrank,nproc
  USE schism_msgp, ONLY : parallel_barrier 
  USE sed2d_mod, ONLY : bc_flag,idsed2d,vc_area

  IMPLICIT NONE

  INCLUDE 'mpif.h'

!- Arguments --------------------------------------------------------
  REAL(rkind), DIMENSION(npa), INTENT(IN)  :: dp1
  REAL(rkind), DIMENSION(npa), INTENT(OUT) :: dp2 
!- Local variables --------------------------------------------------
  INTEGER :: i,iflag,iflag_gb,imax,imin,iter,j,nij,nnode,nnode_gb
  REAL(rkind) :: dpij,dpmax,dpmed,dpmin
!- User-defined parameter -------------------------------------------
  INTEGER, PARAMETER :: str_filt = 1 !1=strong filter / 2=weak filter
  REAL(rkind), PARAMETER :: epsi = 0.10d0 !Threshold difference  
!--------------------------------------------------------------------
  iter = 0 
  iflag = 1
  iflag_gb = 1
  dp2 = dp1
!- Look for local extrema -------------------------------------------
  DO WHILE (iflag_gb .EQ. 1)
     iflag = 0
     iflag_gb = 0
     nnode = 0
     iter  = iter+1
     DO i = 1, np
        imin = indnd(1,i)
        imax = imin
        dpmin = dp2(imin)
        dpmax = dpmin
        DO j = 2, nnp(i)
           nij  = indnd(j,i)
           dpij = dp2(nij)
           IF (dpij .LT. dpmin) THEN
              imin = nij
              dpmin = dpij
           ELSEIF (dpij .GT. dpmax) THEN
              imax = nij
              dpmax = dpij
           ENDIF
        ENDDO !nnp(i)

!- Apply the strong filter ------------------------------------------
        IF ((str_filt.EQ.1).AND.((.NOT.bc_flag(i)).AND.              &
           (.NOT.bc_flag(imax)).AND.(.NOT.bc_flag(imin)))) THEN
           IF (dp2(i)-dpmin .LT. -epsi) THEN !i is a local minimum
              iflag = 1
              dpmed = (dp2(i)*vc_area(i)+dpmax*vc_area(imax))/       &
                      (vc_area(i)+vc_area(imax))
              dp2(imax) = dpmed
              dp2(i) = dpmed
              nnode = nnode+1
#ifdef USE_DEBUG
              WRITE(12,*)'Strong min. found at node:',iplg(i)
              !IF(isnan(dp2(i)) .EQV. .TRUE.) THEN
              IF(dp2(i)/=dp2(i)) THEN
                WRITE(12,*)'NaN found (NL filter 1)',dp2(i),dpmax
              ENDIF
#endif
           ENDIF
           IF (dp2(i)-dpmax .GT. epsi) THEN !i is a local maximum
              iflag = 1
              dpmed = (dp2(i)*vc_area(i)+dpmin*vc_area(imin))/       &
                      (vc_area(i)+vc_area(imin))           
              dp2(imin) = dpmed
              dp2(i) = dpmed
              nnode = nnode+1
#ifdef USE_DEBUG
              WRITE(12,*)'Strong max. found at node:',iplg(i)
              !IF(isnan(dp2(i)) .EQV. .TRUE.) THEN
              IF(dp2(i)/=dp2(i)) THEN
                WRITE(12,*)'NaN found (NL filter 2)',dp2(i),dpmin
              ENDIF
#endif
           ENDIF

!- Apply the weak filter --------------------------------------------
!- If the node is on a boundary the filter is always set to weak
!--------------------------------------------------------------------
        ELSE
           IF (dp2(i)-dpmin .LT. -epsi) THEN !i is a local minimum
              iflag = 1
              dpmed = (dp2(i)*vc_area(i)+dpmin*vc_area(imin))/       &
                      (vc_area(i)+vc_area(imin))
              dp2(imin) = dpmed
              dp2(i) = dpmed
              nnode = nnode+1
#ifdef USE_DEBUG
              WRITE(12,*)'Weak min. found at node:',iplg(i)
              !IF(isnan(dp2(i)) .EQV. .TRUE.) THEN
              IF(dp2(i)/=dp2(i)) THEN
                WRITE(12,*)'NaN found (NL filter 3)',dp2(i),dpmin
              ENDIF
#endif
           ENDIF
           IF (dp2(i)-dpmax .GT. epsi) THEN ! i is a local maximum
              iflag = 1
              dpmed = (dp2(i)*vc_area(i)+dpmax*vc_area(imax))/       &
                      (vc_area(i)+vc_area(imax))
              dp2(imax) = dpmed
              dp2(i) = dpmed
              nnode = nnode+1
#ifdef USE_DEBUG
              WRITE(12,*)'Weak max. found at node:',iplg(i)
              !IF(isnan(dp2(i)) .EQV. .TRUE.) THEN
              IF(dp2(i)/=dp2(i)) THEN
                WRITE(12,*)'NaN found (NL filter 4)',dp2(i),dpmax
              ENDIF
#endif
           ENDIF
        ENDIF
     ENDDO !np
     CALL exchange_p2d(dp2)
     CALL mpi_reduce(nnode,nnode_gb,1,itype,MPI_SUM,0,comm,ierr)
     CALL mpi_allreduce(iflag,iflag_gb,1,itype,MPI_MAX,comm,ierr) 
     IF(myrank==0 .AND. nnode_gb>0) THEN
       WRITE(idsed2d,*)'Iter #:',iter,'Nb of extrema:',nnode_gb
     ENDIF
  ENDDO !iflag_gb

!  IF(myrank==0 .AND. iter>1) THEN
!    WRITE(idsed2d,*)'Total number of iterations in NL filter:',iter
!  ENDIF

  END SUBROUTINE sed2d_filter_extrema
!--------------------------------------------------------------------

!--------------------------------------------------------------------
  SUBROUTINE sed2d_filter_slope(dp1,dp2)
!--------------------------------------------------------------------
! This routine computes a new bathymetry to account for avalanching. 
! The bed slopes are computed at each element. When the slope 
! exceeds a critical value (specified in slope_cr.gr3 or default values below), the  
! bathymetry of each element node is modified to obtain a slope lower
! than the critical threshold. This method conserves the volume.     
!
! Adapted from filter.f (SAND2D, A. Fortunato)                      
! 
! Author: florian ganthy (fganthy@lnec.pt ; florian.ganthy@gmail.com)
! Date:   2013/01/16                                                 
!                                                                    
! History:                                                          
! 01/2013 - G.Dodet: Adapted this subroutine to sed2d environment and
!                    included it in module sed2d_filter.
!                    Critical values are read in slope_cr.gr3 (if no
!                    file, default valut is 1.
! 03/2013 - G.Dodet: Added user-defined threshold value for maximum
!                    slope in wet and dry area, in case they are not
!                    provided in slope_cr.gr3
! 05/2013 - G.Dodet: - Synchronized processors and mpi_reduced some
!                      variables for consistency; 
!                    - Removed computation of volume control area and
!                      slope_cr.gr3 reading (redundant). This is now 
!                      done in sed2d_init
! 04/2017 - T.Guerin: Added user-defined dryslope and wetslope
!--------------------------------------------------------------------

  USE schism_glbl, ONLY: dldxy,idry,nea,elnode,npa,rkind,xnd,ynd, &
 &in_dir,out_dir,len_in_dir,len_out_dir
  USE schism_msgp, ONLY: comm,exchange_p2d,ierr,itype,myrank,nproc,    &
                       parallel_abort

  USE sed2d_mod, ONLY: idsed2d,slope_cr,vc_area,dryslope,wetslope

  IMPLICIT NONE

  INCLUDE 'mpif.h'

!- Arguments --------------------------------------------------------
  REAL(rkind), DIMENSION(npa), INTENT(IN) :: dp1
  REAL(rkind), DIMENSION(npa), INTENT(OUT) :: dp2 
!- Local variables --------------------------------------------------
  INTEGER :: i,iflag,iflag_gb,iter,j,n1,n2,n3
  REAL(rkind) :: det,h1,h2,h3,h1p,h2p,h3p,m11,m12,m13,m21,m22,m23,   &
                 m31,m32,m33,mas1,mas2,mas3,slope,slopemax,tmp,vec1, &
                 vec2,vec3,xtmp,ytmp
  REAL(rkind), DIMENSION(nea) :: slope_cr_e
  REAL(rkind), DIMENSION(nea,2) :: dpdxy_el
  LOGICAL :: lexist
!- User-defined parameters ------------------------------------------
  INTEGER, PARAMETER :: maxiter = 10 !Max. number of iterations
  REAL(rkind), PARAMETER :: epsi = 0.01d0
!--------------------------------------------------------------------

!- Check if user file is provided -----------------------------------
  INQUIRE(FILE=in_dir(1:len_in_dir)//'slope_cr.gr3',EXIST=lexist)

!Error: YJZ - not working with quads
!- Start iterative procedure ----------------------------------------
  dp2 = dp1
  iter  = 0
  iflag = 1
  iflag_gb = 1
  DO WHILE (iflag_gb.EQ.1)
     iflag = 0
     iter  = iter+1
     dpdxy_el = 0.d0
     slope_cr_e = 0.d0
     DO i = 1,nea

!- Compute bed slope at element center ------------------------------
        DO j=1,3
           dpdxy_el(i,1) = dpdxy_el(i,1)+dp2(elnode(j,i))*dldxy(j,1,i)
           dpdxy_el(i,2) = dpdxy_el(i,2)+dp2(elnode(j,i))*dldxy(j,2,i)
           slope_cr_e(i) = slope_cr_e(i)+slope_cr(elnode(j,i))/3.d0
        ENDDO
        slope = DSQRT(dpdxy_el(i,1)*dpdxy_el(i,1)+dpdxy_el(i,2)*     &
                      dpdxy_el(i,2))

!- Check if the 3 nodes of the element are dry ----------------------
        IF(lexist) THEN
          slopemax = slope_cr_e(i)
        ELSE
          IF(idry(elnode(1,i))+idry(elnode(2,i))+idry(elnode(3,i)) == 3)THEN
            slopemax = dryslope
          ELSE
            slopemax = wetslope
          ENDIF
        ENDIF

        IF(slope-slopemax.LE.epsi) CYCLE
!- Preparation of equation system -----------------------------------
        n1 = elnode(1,i)
        n2 = elnode(2,i)
        n3 = elnode(3,i)
        h1 = dp2(n1)
        h2 = dp2(n2)
        h3 = dp2(n3)
        m11 = vc_area(n1)
        m12 = vc_area(n2)
        m13 = vc_area(n3)
        m21 = ynd(n2)-ynd(n3)
        m22 = ynd(n3)-ynd(n1)
        m23 = ynd(n1)-ynd(n2)
        m31 = xnd(n3)-xnd(n2)
        m32 = xnd(n1)-xnd(n3)
        m33 = xnd(n2)-xnd(n1)
        vec1 = h1*m11+h2*m12+h3*m13
        vec2 = slopemax/slope*(h1*m21+h2*m22+h3*m23)
        vec3 = slopemax/slope*(h1*m31+h2*m32+h3*m33)

!- Solving the system following Cramer's rule -----------------------
        det = m11*(m22*m33-m32*m23)-                                 &
              m12*(m21*m33-m31*m23)+                                 &
              m13*(m21*m32-m31*m22)

        IF(det.EQ.0.0d0) CALL parallel_abort('sed2d_filter: det=0.0')
!- Compute new depth at nodes ---------------------------------------
        dp2(n1) = (vec1*(m22*m33-m32*m23)-                           &
                  m12*(vec2*m33-vec3*m23)+                           &
                  m13*(vec2*m32-vec3*m22))/det
        dp2(n2) = (m11*(vec2*m33-vec3*m23)-                          &
                  vec1*(m21*m33-m31*m23)+                            &
                  m13*(m21*vec3-m31*vec2))/det
        dp2(n3) = (m11*(m22*vec3-m32*vec2)-                          &
                  m12*(m21*vec3-m31*vec2)+                           &
                  vec1*(m21*m32-m31*m22))/det
        iflag = 1
      ENDDO !nea
      CALL exchange_p2d(dp2)
      CALL mpi_allreduce(iflag,iflag_gb,1,itype,MPI_MAX,comm,ierr)
      IF (iter.GE.maxiter) THEN
         iflag_gb = 0
         IF(myrank.EQ.0)WRITE(idsed2d,*)'Warning: maximum iterations &
     &number reached in slope filter. Increase maxiter in sed2d_filter_slope'
      ENDIF
  ENDDO !iflag
  IF(iter > 1)THEN
    IF(myrank.EQ.0)WRITE(idsed2d,*)'# of iter. in slope filter:',iter
  ENDIF
  END SUBROUTINE sed2d_filter_slope
!--------------------------------------------------------------------

!--------------------------------------------------------------------
  SUBROUTINE sed2d_filter_diffu(var1,var2,ndim)
!--------------------------------------------------------------------
! This subroutine applies a diffusive filter on any variable by
! interpolating from node (resp. side, resp. element) to element 
! (resp. node, resp. node) and vice-versa.
!
! Author: guillaume dodet (gdodet01@univ-lr.fr, gdodet@lnec.pt)
! Date:   23/01/2013
!
! History:
! 04/2013 - G.Dodet: Extended the filter to all kind of variables 
!                    (eg. defined on node, side, elment).
!--------------------------------------------------------------------

  USE schism_glbl, ONLY : indel,isidenode,nea,nne,i34,elnode,np,npa,nsa,rkind
  USE schism_msgp, ONLY : exchange_e2d,exchange_p2d,exchange_s2d,      &
                        parallel_abort

  IMPLICIT NONE

!- Arguments --------------------------------------------------------
  INTEGER, INTENT(IN) :: ndim
  REAL(rkind), DIMENSION(ndim), INTENT(IN)  :: var1
  REAL(rkind), DIMENSION(ndim), INTENT(OUT) :: var2
!- Local variables --------------------------------------------------
  INTEGER :: i,iel,inode,j
  REAL(rkind), DIMENSION(nea) :: tmp_e
  REAL(rkind), DIMENSION(nsa) :: tmp_s
  REAL(rkind), DIMENSION(npa) :: neigh,tmp_n
!--------------------------------------------------------------------
  var2 = 0.D0
  IF(ndim == npa) THEN !Node-element-node interpolation
    tmp_e = 0.D0
    DO i = 1,nea
       DO j = 1,i34(i)
          inode = elnode(j,i)
          tmp_e(i) = tmp_e(i)+var1(inode)/i34(i)
       ENDDO
    ENDDO
    CALL exchange_e2d(tmp_e)

    DO i = 1,np
       DO j = 1,nne(i)
          iel = indel(j,i)
          var2(i) = var2(i)+tmp_e(iel)/nne(i)
       ENDDO
    ENDDO
    CALL exchange_p2d(var2)

  ELSEIF(ndim == nea) THEN !Element-node-element interpolation
    tmp_n = 0.d0
    DO i = 1,np
       DO j = 1,nne(i)
          iel = indel(j,i)
          tmp_n(i) = tmp_n(i)+var1(iel)/nne(i)
       ENDDO
    ENDDO
    CALL exchange_p2d(tmp_n)

    DO i = 1,nea
       DO j = 1,i34(i)
          inode = elnode(j,i)
          var2(i) = var2(i)+tmp_n(inode)/i34(i)
       ENDDO
    ENDDO
    CALL exchange_e2d(var2)
           
  ELSEIF(ndim == nsa) THEN !Side-node-side interpolation
    tmp_n = 0.d0
    neigh = 0.d0
    DO i = 1,nsa
       DO j = 1,2
          inode = isidenode(j,i) 
          tmp_n(inode) = tmp_n(inode)+var1(i)
          neigh(inode) = neigh(inode)+1  
       ENDDO
    ENDDO
    DO i = 1,npa
       tmp_n(i) = tmp_n(i)/neigh(i)
    ENDDO
    CALL exchange_p2d(tmp_n) 

    DO i = 1,nsa
       DO j = 1,2
          inode = isidenode(j,i)
          var2(i) = var2(i)+tmp_n(inode)/2.d0
       ENDDO
    ENDDO
    CALL exchange_s2d(var2)
  ELSE
    CALL parallel_abort('Wrong dimension in sed2d_filter_diffu')
  ENDIF
  
  END SUBROUTINE sed2d_filter_diffu
!--------------------------------------------------------------------
END MODULE sed2d_filter
