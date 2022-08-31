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

!=====================================================================
!=====================================================================
! MORSELFE FILTER SUBROUTINES
!
! subroutine sed_avalanching
! subroutine sed2d_filter_diffu
!
!=====================================================================
!=====================================================================

      SUBROUTINE sed_avalanching(dhnd)
!--------------------------------------------------------------------!
! This routine updates the bathymetry to account for avalanching     !
! The bed slopes are computed at each node. When a critical slope for!
! element is exceeded, bathymetry of each element node is modified   !
! in order to obtain a slope=critical threshold. The      
! method used here conserves the volume.                             !
! Adapted from filter.f (SAND2D, A. Fortunato)                       !
!                                                                    !
! Currently it does not account for modification to bed sediment
! characteristics related to modification of bathymetry              !
!                                                                    !
! Author: florian ganthy (fganthy@lnec.pt ; florian.ganthy@gmail.com)!
! Date:   2013/01/16                                                 !
!                                                                    !
! History:                                                           !
! 2013/01 - F.Ganthy : Modification of wet/dry element consideration !
!                      to be more physical.                          !
! 2013/05 - F.Ganthy : Synchronized processors and mpi_reduced some  !
!                      variables for consistency.                    !
!                      Removed computation of volume control, now    ! 
!                      done in sed_init                              !
!                                                                    !
!--------------------------------------------------------------------!
! Details of algorithm/eqs.
! (0) split quads into pair of tri's
! (1) conservation of volume in an elem.
! \sum(Ai*xi)=\sum(Ai*hi) : Ai area of dual graph, hi is depth at node of elem,
!    and xi is the new hi;
! (2) Force slope to become slope_cr, and preserve the original direction
! \sum(xi*dldxy(,1:2,))=slope_cr*slope_[xy]/slope : slope_[xy] are \nabla h,
! and slope=|\nabla h|

      USE schism_glbl, ONLY: dldxy,idry,nea,i34,elnode,npa,rkind,area,    &
     &                     errmsg,dp,np,xel,yel,nxq
      USE schism_msgp, ONLY: comm,exchange_p2d,ierr,itype,myrank,      &
     &                     nproc,parallel_abort
      USE sed_mod,   ONLY: dry_slope_cr,wet_slope_cr,vc_area

      IMPLICIT NONE

      INCLUDE 'mpif.h'

      REAL(rkind) :: signa

!- Arguments --------------------------------------------------------!
      REAL(rkind), DIMENSION(npa), INTENT(inout) :: dhnd
!- Local variables --------------------------------------------------!
      INTEGER :: i,j,iter,iflag,iflag_gb,n1,n2,n3,ndry,m,jj,nwild(3),nwild2(3)
      REAL(rkind) :: slope_cr,slope,h1,h2,h3,vec1,vec2,vec3,         &
                     m11,m12,m13,m21,m22,m23,m31,m32,m33,det,        &
                     h1p,h2p,h3p,ar2,dldxy2(3,2)
      REAL(rkind), DIMENSION(:)   :: dp0(npa),dp1(npa),area2(npa),   &
                                     dph(npa)
      REAL(rkind), DIMENSION(:,:) :: dpdxy_el(nea,2)
!- User-defined parameters ------------------------------------------      
      INTEGER,PARAMETER :: maxiter = 30 !Maximum number of iterations
      !REAL(rkind), PARAMETER :: epsi  = 0.01d0 !too large??
      REAL(rkind), PARAMETER :: epsi  = 1.e-4

!- Start Statement --------------------------------------------------!

!YJZ: Error: this routine will violate bare rock limit
!--------------------------------------------------------------------!
! * Compute bed changes due to sediment transport
!--------------------------------------------------------------------!
      DO i=1,npa
        dp0(i) = dp(i)+dhnd(i)
      ENDDO
      dp1=dp0

!--------------------------------------------------------------------!
! * Start iterative procedure
!--------------------------------------------------------------------!
      iter     = 0
      iflag    = 1
      iflag_gb = 1
      
      DO WHILE (iflag_gb.EQ.1)
        iflag = 0
        iter  = iter+1
        dpdxy_el = 0.0d0
        DO i=1,nea
!--------------------------------------------------------------------!
! * Compute bed slope at element center
!--------------------------------------------------------------------!
          vec1=dot_product(dp1(elnode(1:i34(i),i)),dldxy(1:i34(i),1,i)) !/dx
          vec2=dot_product(dp1(elnode(1:i34(i),i)),dldxy(1:i34(i),2,i)) !/dy
          vec3=sqrt(vec1*vec1+vec2*vec2) !slope

!--------------------------------------------------------------------!
! * Testing for critical slope value (dry or wet)
!   Here we consider that element is dry only if all nodes are 
!   dry. In the former version, the test was done over idry_e (=1), 
!   but this potentially overestimated the slumping at the wet-dry
!   limit
!--------------------------------------------------------------------!
          ndry = sum(idry(elnode(1:i34(i),i)))
          IF(ndry==i34(i)) THEN
            slope_cr = dry_slope_cr
          ELSE
            slope_cr = wet_slope_cr
          ENDIF

          !IF(slope<=slope_cr+epsi) CYCLE
          IF(vec3<=slope_cr+epsi) CYCLE

!         Adjust depths
          iflag = 1
          do m=1,i34(i)-2 !split quads into 2 tri's
            !Calc derivatives of shape function if quads
            if(i34(i)==3) then
              nwild(1:3)=(/1,2,3/) !local indices
              nwild2(1:3)=elnode(nwild(1:3),i) !3 nodes
              dldxy2(1:3,1:2)=dldxy(1:3,1:2,i)
            else !quad
              if(m==1) then !node 1,2,3
                nwild(1:3)=(/1,2,3/) !local indices
              else !1,3,4
                nwild(1:3)=(/1,3,4/)
              endif !m

              nwild2(1:3)=elnode(nwild(1:3),i) !3 nodes of the tri
!              ar2=signa(xnd(nwild2(1)),xnd(nwild2(2)),xnd(nwild2(3)),ynd(nwild2(1)),ynd(nwild2(2)),ynd(nwild2(3)))
              ar2=signa(xel(nwild(1),i),xel(nwild(2),i),xel(nwild(3),i),yel(nwild(1),i),yel(nwild(2),i),yel(nwild(3),i))
              if(ar2<=0) call parallel_abort('SED_FILTER: ar2<=0')
              do jj=1,3
                !Elem. type is 3 not 4!!
                dldxy2(jj,1)=(yel(nwild(nxq(1,jj,3)),i)-yel(nwild(nxq(2,jj,3)),i))/2/ar2 !dL/dx
                dldxy2(jj,2)=(xel(nwild(nxq(2,jj,3)),i)-xel(nwild(nxq(1,jj,3)),i))/2/ar2 !dL/dy
              enddo !jj
            endif !i34

            dpdxy_el(i,1)=dot_product(dp1(nwild2(1:3)),dldxy2(1:3,1)) !dh/dx
            dpdxy_el(i,2)=dot_product(dp1(nwild2(1:3)),dldxy2(1:3,2))
!            DO j=1,i34(i)
!              dpdxy_el(i,1) = dpdxy_el(i,1)+dp1(elnode(j,i))*dldxy(j,1,i)
!              dpdxy_el(i,2) = dpdxy_el(i,2)+dp1(elnode(j,i))*dldxy(j,2,i)
!            ENDDO
            slope=sqrt(dpdxy_el(i,1)*dpdxy_el(i,1)+dpdxy_el(i,2)*dpdxy_el(i,2))

!--------------------------------------------------------------------!
! * Preparation of system equation
!--------------------------------------------------------------------!
            n1 = nwild2(1)
            n2 = nwild2(2)
            n3 = nwild2(3)
            h1 = dp1(n1)
            h2 = dp1(n2)
            h3 = dp1(n3)
            !Matrix
            m11 = vc_area(n1)
            m12 = vc_area(n2)
            m13 = vc_area(n3)
            m21 = dldxy2(1,1) !dL_1/dx
            m22 = dldxy2(2,1) 
            m23 = dldxy2(3,1) 
            m31 = dldxy2(1,2) 
            m32 = dldxy2(2,2) 
            m33 = dldxy2(3,2) 
            vec1 = h1*m11 + h2*m12 + h3*m13 !RHS
            !slope checked
            !vec2 = slope_cr/slope * (h1*m21 + h2*m22 + h3*m23)
            !vec3 = slope_cr/slope * (h1*m31 + h2*m32 + h3*m33)
            vec2 = slope_cr*dpdxy_el(i,1)/slope
            vec3 = slope_cr*dpdxy_el(i,2)/slope
!--------------------------------------------------------------------!
! * Solving the system by Cramer's rule
!--------------------------------------------------------------------!
            det=m11*(m22*m33-m32*m23)-m12*(m21*m33-m31*m23)+m13*(m21*m32-m31*m22)

            IF(det==0.0d0) THEN
              WRITE(errmsg,*)'SED_AVALANCHING: det=0.0'
              CALL parallel_abort(errmsg)
            ENDIF
!--------------------------------------------------------------------!
! * Compute new depth at nodes
!--------------------------------------------------------------------!
            dp1(n1)=(vec1*(m22*m33-m32*m23)-m12*(vec2*m33-vec3*m23)+m13*(vec2*m32-vec3*m22))/det
            dp1(n2)=(m11*(vec2*m33-vec3*m23)-vec1*(m21*m33-m31*m23)+m13*(m21*vec3-m31*vec2))/det
            dp1(n3)=(m11*(m22*vec3-m32*vec2)-m12*(m21*vec3-m31*vec2)+vec1*(m21*m32-m31*m22))/det
          enddo !m=1,i34(i)-2
        ENDDO !i=1,nea
        CALL exchange_p2d(dp1)
        CALL mpi_allreduce(iflag,iflag_gb,1,itype,MPI_MAX,comm,ierr)

        !No CPU dependency as _allreduce is a barrier
        IF(iter>=maxiter) THEN
          iflag_gb = 0 !reset flag for exit
          IF(myrank.EQ.0) WRITE(16,*)'Warning: max iterations     &
     &    number reached in sed_avalanching.'
        ENDIF
        
      ENDDO ! End do while on iflag

!--------------------------------------------------------------------!
! * Write number of iteration within mirror.out
!--------------------------------------------------------------------!
      IF(iter > 1.and.myrank==0) THEN
        WRITE(16,*)'# of iter. in sed_avalanching:',iter
      ENDIF
!--------------------------------------------------------------------!
! * Apply depth changes
!--------------------------------------------------------------------!
      DO i=1,npa
        dhnd(i) = dhnd(i)+dp1(i)-dp0(i)
        IF(dhnd(i)/=dhnd(i)) THEN
          WRITE(errmsg,*) 'Avalanching: dhnd is NaN',myrank,dhnd(i), &
                          dp1(i),dp0(i)
          CALL parallel_abort(errmsg)
        ENDIF
      ENDDO
!      CALL exchange_p2d(dhnd)

      END SUBROUTINE sed_avalanching

!=====================================================================

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
! 02/2020 - B.Mengual: implementation in SED3D
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


!=====================================================================

