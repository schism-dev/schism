!     Module to compute wave force using Longuet-Higgins Stewart formulation, 
!     to be used for ESMF coupler.
      module compute_wave_force
      implicit none

      contains

      subroutine compute_wave_force_lon(RSXX,RSXY,RSYY)
      use schism_glbl, only : rkind,nsa,npa,nvrt,idry,idry_s,dps,hmin_radstress, &
     &WWAVE_FORCE
      use schism_msgp
      REAL(rkind), intent(inout) :: RSXX(npa),RSXY(npa),RSYY(npa) !dimension: m*m/s/s
     
      REAL(rkind), allocatable :: DSXX3D(:,:,:),DSXY3D(:,:,:),DSYY3D(:,:,:)
      integer :: IS
      REAL(rkind) :: HTOT
    
      allocate(DSXX3D(2,NVRT,nsa), DSYY3D(2,NVRT,nsa),DSXY3D(2,NVRT,nsa))

      where(idry==1) 
        RSXX=0.d0
        RSXY=0.d0
        RSYY=0.d0
      end where

      ! Computing gradients of the depth-averaged radiation stress
      ! terms (unit: m^2/s/s)
      CALL hgrad_nodes(2,0,NVRT,npa,nsa,RSXX,DSXX3D)   !C(dSxx/dx , dSxx/dy )
      CALL hgrad_nodes(2,0,NVRT,npa,nsa,RSYY,DSYY3D)   !C(dSyy/dx , dSyy/dy )
      CALL hgrad_nodes(2,0,NVRT,npa,nsa,RSXY,DSXY3D)   !C(dSxy/dx , dSxy/dy )
      CALL exchange_s3d_2(DSXX3D)
      CALL exchange_s3d_2(DSYY3D)
      CALL exchange_s3d_2(DSXY3D)

      ! Computing the wave forces; these are noted Rsx, Rsy in Rolland
      ! et al. (2012), see Eq. (9)
      ! These are stored in wwave_force(:,1:nsa,1:2) (unit: m/s/s)
      WWAVE_FORCE=0.d0 !m/s/s
      DO IS=1,nsa
        IF(idry_s(IS)==1) CYCLE

        ! Total water depth at sides
        HTOT=MAX(dps(IS),hmin_radstress)

        ! Wave forces
        WWAVE_FORCE(1,:,IS)=WWAVE_FORCE(1,:,IS)-(DSXX3D(1,:,IS)+DSXY3D(2,:,IS))/HTOT
        WWAVE_FORCE(2,:,IS)=WWAVE_FORCE(2,:,IS)-(DSXY3D(1,:,IS)+DSYY3D(2,:,IS))/HTOT
      ENDDO !IS

      deallocate(DSXX3D,DSYY3D,DSXY3D)
      end subroutine compute_wave_force_lon

      end module compute_wave_force
