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

!**********************************************************************
!*  Author: Ivica Janeković [ivica.jan@gmail.com]                                                                   *
!**********************************************************************
! This is set of subroutines for handling wind and mslp from netcdf WRF
! to get Sschism forcing fields.
! Before running model you have to compute interpolation weights using 
! matlab or any other program. To compute weights you need information 
! of atmo model grid (lon, lat) and Schism grid (lon, lat)
! see example Utility/Pre-Processing/make_schism_interp_coefs.m
! *********************************************************************
! list of SUBROUTINES:

! READ_REC_ATMO_FD	   read fields and interpolate them on the Schism grid
! READ_REC_ATMO_FEM     read fileds already interpolated on the Schism grid
! GENERIC_NETCDF_ERROR  for error handling
!
!**********************************************************************
!*                                                                    *
!**********************************************************************      
      SUBROUTINE READ_REC_ATMO_FD(rec, windx, windy, pr)
!      input is record in netcdf that you want to read, outputs are U10,V10,MSLP     
      USE NETCDF
      USE schism_glbl, only : npa, rkind, cf_x2, cf_x1, cf_y2, cf_y1, cf_i, cf_j, cf_denom, &
     &in_dir,out_dir,len_in_dir,len_out_dir
!      those are already computed weights externaly and read in the schism_init 
      USE schism_msgp, only : parallel_abort, myrank

      IMPLICIT NONE
      INTEGER, INTENT(in)          :: rec
      REAL(rkind), INTENT(out)     :: windx(npa), windy(npa), pr(npa)
      INTEGER                      :: i,istat, fid, varid, ndims, dimids(3), dim1, dim2, dim3
      REAL(rkind)			     :: cf_add_offset_uwind, cf_scale_factor_uwind
      REAL(rkind)			     :: cf_add_offset_vwind, cf_scale_factor_vwind
      REAL(rkind)			     :: cf_add_offset_pr,    cf_scale_factor_pr
      REAL(rkind), ALLOCATABLE     :: uwind_fd(:,:), vwind_fd(:,:), mslp_fd(:,:)
      CHARACTER(len = *), parameter       :: callfct="READ_REC_ATMO_FD"
      CHARACTER(len = 100)         :: chrerr
      LOGICAL                      :: got_uwind, got_vwind, got_mslp
      
      got_uwind = .false.
      got_vwind = .false.
      got_mslp  = .false.
        
! Open NC file
! hardcoded WIND MSLP file name to UVP.nc !
      istat = NF90_OPEN(in_dir(1:len_in_dir)//'UVP.nc', nf90_nowrite, fid)
	IF (istat /= 0) THEN
            CALL PARALLEL_ABORT('READ_REC_ATMO_FD: missing UVP.nc file as defined nws = 5 in param.in')
	ENDIF
! Reading Uwind attributes
      istat = NF90_INQ_VARID(fid, "Uwind", varid)
	IF (istat /= 0) THEN
	      IF(myrank == 0)  WRITE(16,*) 'There is no Uwind variable defined in UVP.nc, setting all to zero!'
            windx(:) = 0.0_rkind
	ELSE
! scale_factor
            istat = NF90_GET_ATT(fid, varid, "scale_factor", cf_scale_factor_uwind)
	      IF (istat /= 0) cf_scale_factor_uwind = 1.0_rkind
	      IF(myrank==0) WRITE(16,'(A,F6.3)') '++ netCDF scale factor for Uwind = ', cf_scale_factor_uwind
! add_offset
            istat = NF90_GET_ATT(fid, varid, "add_offset", cf_add_offset_uwind)
	      IF (istat /= 0)  cf_add_offset_uwind = 0.0_rkind
	      IF(myrank==0) WRITE(16,'(A,F6.3)') '++ netCDF offset for Uwind       = ', cf_add_offset_uwind
! read Uwind
            istat = NF90_INQUIRE_VARIABLE (fid, varid, ndims=ndims, dimids=dimids)
            IF (ndims /= 3) CALL PARALLEL_ABORT('Uwind must be 3D variable in format (time,nlat,nlon)')
            istat = NF90_INQUIRE_DIMENSION(fid, dimids(1), len = dim1)
            istat = NF90_INQUIRE_DIMENSION(fid, dimids(2), len = dim2)
            istat = NF90_INQUIRE_DIMENSION(fid, dimids(3), len = dim3)
            IF(myrank==0)  WRITE(16,'(A,I0,A,I0,A,I0,A)') '++ Uwind dimensions: (', dim1, ',', dim2, ',', dim3, ')'
            ALLOCATE(uwind_fd(dim1,dim2), stat=istat)
            istat = NF90_GET_VAR(fid, varid, uwind_fd, start = (/ 1, 1, rec /), count = (/ dim1, dim2, 1 /))
            uwind_fd(:,:)= cf_add_offset_uwind + cf_scale_factor_uwind * uwind_fd(:,:)
	      IF(myrank==0)  WRITE(16,'(A,I0)') '++ Done with reading Uwind from netCDF file UVP.nc and record ', rec
            got_uwind = .true.
      ENDIF
! reading Vwind attributes
      istat = NF90_INQ_VARID(fid, "Vwind", varid)
	IF (istat /= 0) THEN
	      IF(myrank==0)  WRITE(16,*) 'There is no Vwind variable defined in UVP.nc, setting all to zero!'
            windy(:) = 0.0_rkind
	ELSE
! scale_factor 
            istat = NF90_GET_ATT(fid, varid, "scale_factor", cf_scale_factor_vwind)
            IF (istat /= 0) cf_scale_factor_vwind = 1.0_rkind
            IF(myrank==0) WRITE(16,'(A,F6.3)') '++ netCDF scale factor for Vwind = ', cf_scale_factor_vwind
! add_offset
            istat = NF90_GET_ATT(fid, varid, "add_offset", cf_add_offset_vwind)
            IF (istat /= 0)  cf_add_offset_vwind = 0.0_rkind
            IF(myrank==0) WRITE(16,'(A,F6.3)') '++ netCDF offset for Vwind       = ', cf_add_offset_vwind
! read Vwind
            istat = NF90_INQUIRE_VARIABLE (fid, varid, ndims=ndims, dimids=dimids)
            IF (ndims /= 3) CALL PARALLEL_ABORT('Vwind must be 3D variable in format (time,nlat,nlon)')
            istat = NF90_INQUIRE_DIMENSION(fid, dimids(1), len = dim1)
            istat = NF90_INQUIRE_DIMENSION(fid, dimids(2), len = dim2)
            istat = NF90_INQUIRE_DIMENSION(fid, dimids(3), len = dim3)
            IF(myrank==0)  WRITE(16,'(A,I0,A,I0,A,I0,A)') '++ Vwind dimensions: (', dim1, ',', dim2, ',', dim3, ')'
            ALLOCATE(vwind_fd(dim1,dim2), stat=istat)
            istat = NF90_GET_VAR(fid, varid, vwind_fd, start = (/ 1, 1, rec /), count = (/ dim1, dim2, 1 /))
            vwind_fd(:,:)= cf_add_offset_vwind + cf_scale_factor_vwind * vwind_fd(:,:)
            IF(myrank==0)  WRITE(16,'(A,I0)') '++ Done with reading Vwind from netCDF file UVP.nc and record ', rec
            got_vwind = .true.
      ENDIF
! reading Pair attributes
      istat = NF90_INQ_VARID(fid, "Pair", varid)
      IF (istat /= 0) THEN
            IF(myrank==0)  WRITE(16,*) 'There is no Pair variable defined in UVP.nc, setting all to zero!'
            pr(:) = 101325.0_rkind          
      ELSE
! scale_factor 
            istat = NF90_GET_ATT(fid, varid, "scale_factor", cf_scale_factor_pr)
            IF (istat /= 0) cf_scale_factor_pr = 1.0_rkind
            IF(myrank==0) WRITE(16,'(A,F6.3)') '++ netCDF scale factor for Pair = ', cf_scale_factor_pr
! add_offset
            istat = NF90_GET_ATT(fid, varid, "add_offset", cf_add_offset_pr)
            IF (istat /= 0)  cf_add_offset_pr = 0.0_rkind
            IF(myrank==0) WRITE(16,'(A,F12.3)') '++ netCDF offset for Pair      = ', cf_add_offset_pr
! read Pair
            istat = NF90_INQUIRE_VARIABLE (fid, varid, ndims=ndims, dimids=dimids)
            IF (ndims /= 3) CALL PARALLEL_ABORT('Pair must be 3D variable in format (time,nlat,nlon)')
            istat = NF90_INQUIRE_DIMENSION(fid, dimids(1), len = dim1)
            istat = NF90_INQUIRE_DIMENSION(fid, dimids(2), len = dim2)
            istat = NF90_INQUIRE_DIMENSION(fid, dimids(3), len = dim3)
            IF(myrank==0)  WRITE(16,'(A,I0,A,I0,A,I0,A)') '++ Pair dimensions: (', dim1, ',', dim2, ',', dim3, ')'
            ALLOCATE(mslp_fd(dim1,dim2), stat=istat)
            istat = NF90_GET_VAR(fid, varid, mslp_fd, start = (/ 1, 1, rec /), count = (/ dim1, dim2, 1 /))
            mslp_fd(:,:)= 100.0_rkind * (cf_add_offset_pr + cf_scale_factor_pr * mslp_fd(:,:))  !! original data is in mb !!
            IF(myrank==0)  THEN
                  WRITE(16,'(A,I0)') '++ Done with reading Pair from netCDF file UVP.nc and record ', rec
                  FLUSH(16)
            ENDIF  
            got_mslp = .true.
        ENDIF
! close UVP.nc 
      istat = NF90_CLOSE(fid)

! test for grids and interpolation
      IF ( (maxval(cf_i(:)) > dim1) .or. (maxval(cf_j(:)) > dim2) ) THEN
            IF(myrank==0) THEN 
            WRITE(16,*) 'Max(cf_i)= ',maxval(cf_i(:)),' dim1= ',dim1,' Max(cf_j)= ',maxval(cf_j(:)),' dim2= ',dim2 
            FLUSH(16)
            ENDIF
            CALL PARALLEL_ABORT('READ_REC_ATMO_FD: max(cf_i)>dim1 or max(cf_j)>dim2')
	ENDIF        
! construct interpolated fields for mapping on the FEM grid, uwind(nx,ny) -> windx(i)
! uwind_fd, vwind_fd and mslp_fd MUST HAVE first dimension X, then Y ! Transpose them if not!

      DO I = 1, npa
            windx(I)=(  uwind_fd( cf_i(I)  , cf_j(I)  )*cf_x2(I)*cf_y2(I)+        &
            &           uwind_fd( cf_i(I)+1, cf_j(I)  )*cf_x1(I)*cf_y2(I)+        &
            &           uwind_fd( cf_i(I)  , cf_j(I)+1)*cf_x2(I)*cf_y1(I)+        &
            &           uwind_fd( cf_i(I)+1, cf_j(I)+1)*cf_x1(I)*cf_y1(I) )/cf_denom(I)
            windy(I)=(  vwind_fd( cf_i(I)  , cf_j(I)  )*cf_x2(I)*cf_y2(I)+        &
            &           vwind_fd( cf_i(I)+1, cf_j(I)  )*cf_x1(I)*cf_y2(I)+        &
            &           vwind_fd( cf_i(I)  , cf_j(I)+1)*cf_x2(I)*cf_y1(I)+        &
            &           vwind_fd( cf_i(I)+1, cf_j(I)+1)*cf_x1(I)*cf_y1(I) )/cf_denom(I)
            pr(I)   =(  mslp_fd(  cf_i(I)  , cf_j(I)  )*cf_x2(I)*cf_y2(I)+        &
            &           mslp_fd(  cf_i(I)+1, cf_j(I)  )*cf_x1(I)*cf_y2(I)+        &
            &           mslp_fd(  cf_i(I)  , cf_j(I)+1)*cf_x2(I)*cf_y1(I)+        &
            &           mslp_fd(  cf_i(I)+1, cf_j(I)+1)*cf_x1(I)*cf_y1(I) )/cf_denom(I)
      END DO
      IF(myrank==0) THEN
            WRITE(16,'(A,I0)') '++ Done with reading/interpolating atmo forcing from UVP.nc for record ',rec
            WRITE(16,'(A,F6.1,2x,F6.1)') '++ Uwind @ atmo model, min/max = ', minval(uwind_fd), maxval(uwind_fd)
            WRITE(16,'(A,F6.1,2x,F6.1)') '++ Vwind @ atmo model, min/max = ', minval(vwind_fd), maxval(vwind_fd)
            WRITE(16,'(A,F8.1,2x,F8.1)') '++ Pair  @ atmo model, min/max = ', minval(mslp_fd),  maxval(mslp_fd)
            WRITE(16,'(A,F6.1,2x,F6.1)') '++ U10   @ model grid, min/max = ', minval(windx(:)), maxval(windx(:))
            WRITE(16,'(A,F6.1,2x,F6.1)') '++ V10   @ model grid, min/max = ', minval(windy(:)), maxval(windy(:))
            WRITE(16,'(A,F8.1,2x,F8.1)') '++ Pair  @ model grid, min/max = ', minval(pr(:)), maxval(pr(:))
            FLUSH(16)
      ENDIF
      DEALLOCATE(uwind_fd, vwind_fd, mslp_fd)
      END SUBROUTINE READ_REC_ATMO_FD

!**********************************************************************
!*                                                                    *
!**********************************************************************      

      SUBROUTINE READ_REC_ATMO_FEM(rec, windx, windy, pr)
! This subroute reads atmo model forcing from UVP_direct.nc file with variables already interpolated onto
! Schism grid. All dimensions are size of number of nodes in the grid. Input is record that you wish to read
! and outputs are windx, windy, pressure for the record at nodes for the tile.
      USE NETCDF
      USE schism_glbl, only : npa, np_global, ipgl,rkind,in_dir,out_dir,len_in_dir,len_out_dir
      USE schism_msgp, only : myrank, parallel_abort
       
      IMPLICIT NONE
      INTEGER, INTENT(in)           :: rec
      REAL(rkind), INTENT(out)      :: windx(npa), windy(npa), pr(npa)
      INTEGER                       :: i, itmp, istat, fid, varid, ndims, dimids(3), dim1, dim2, dim3
      REAL(rkind)                   :: cf_add_offset_uwind, cf_scale_factor_uwind
      REAL(rkind)                   :: cf_add_offset_vwind, cf_scale_factor_vwind
      REAL(rkind)                   :: cf_add_offset_pr, cf_scale_factor_pr
      REAL(rkind)                   :: pin(1:np_global), uin(1:np_global), vin(1:np_global)

! open NC file, the name is hardcoded
      istat = NF90_OPEN(in_dir(1:len_in_dir)//'UVP_direct.nc', NF90_NOWRITE, fid)
	IF (istat /= 0) THEN
            CALL PARALLEL_ABORT('READ_REC_ATMO_FD: missing UVP_direct.nc file')
	ENDIF
! Reading Uwind attributes
      istat = NF90_INQ_VARID(fid, "Uwind", varid)
	IF (istat /= 0) THEN
	     IF(myrank == 0)  WRITE(16,*) 'There is no Uwind variable in UVP_direct.nc file, setting variable to zero!'
            uin(:) = 0.0_rkind
	ELSE
! read scale_factor
            istat = NF90_GET_ATT(fid, varid, "scale_factor", cf_scale_factor_uwind)
            IF (istat /= 0) cf_scale_factor_uwind = 1.0_rkind
            IF(myrank==0) WRITE(16,'(A,F6.3)') '++ netCDF scale factor for Uwind = ', cf_scale_factor_uwind
! add_offset
            istat = NF90_GET_ATT(fid, varid, "add_offset", cf_add_offset_uwind)
            IF (istat /= 0)  cf_add_offset_uwind = 0.0_rkind
            IF(myrank==0) WRITE(16,'(A,F6.3)') '++ netCDF offset for Uwind       = ', cf_add_offset_uwind
! read Uwind
            istat = NF90_GET_VAR(fid, varid, uin, start = (/ 1, rec /), count = (/ np_global, 1 /))
      ENDIF
! Reading Vwind attributes
      istat = nf90_inq_varid(fid, "Vwind", varid)
      IF (istat /= 0) THEN
	IF(myrank==0)  WRITE(16,*) 'There is no Vwind variable in UVP_direct.nc file, setting variable to zero!'
            vin(:) = 0.0_rkind
	ELSE
! read scale_factor 
            istat = NF90_GET_ATT(fid, varid, "scale_factor", cf_scale_factor_vwind)
            IF (istat /= 0) cf_scale_factor_vwind = 1.0_rkind
            IF(myrank==0) WRITE(16,'(A,F6.3)') '++ netCDF scale factor for Vwind = ', cf_scale_factor_vwind
! add_offset
            istat = NF90_GET_ATT(fid, varid, "add_offset", cf_add_offset_vwind)
            IF (istat /= 0)  cf_add_offset_vwind = 0.0_rkind
            IF(myrank==0) WRITE(16,'(A,F6.3)') '++ netCDF offset for Vwind       = ', cf_add_offset_vwind
! read Vwind
            istat = NF90_GET_VAR(fid, varid, vin, start = (/ 1, rec /), count = (/ np_global, 1 /))
      ENDIF
! Reading Pair attributes
      istat = NF90_INQ_VARID(fid, "Pair", varid)
	IF (istat /= 0) THEN
            IF(myrank==0)  WRITE(16,*) 'There is no Pair variable in UVP_direct.nc file, setting variable to 1013.25!'
            pin(:) = 1013.25_rkind          
	ELSE
! read scale_factor 
            istat = NF90_GET_ATT(fid, varid, "scale_factor", cf_scale_factor_pr)
            IF (istat /= 0) cf_scale_factor_pr = 1.0_rkind
            IF(myrank==0) WRITE(16,'(A,F6.3)') '++ netCDF scale factor for Pair = ', cf_scale_factor_pr
! add_offset
            istat = NF90_GET_ATT(fid, varid, "add_offset", cf_add_offset_pr)
            IF (istat /= 0)  cf_add_offset_pr = 0.0_rkind
            IF(myrank==0) WRITE(16,'(A,F12.3)') '++ netCDF offset for Pair      = ', cf_add_offset_pr
! Read Pair
            istat = NF90_GET_VAR(fid, varid, pin, start = (/ 1, rec /), count = (/ np_global, 1 /))
      ENDIF
      istat = NF90_CLOSE(FID)
 
! MSLP in the WRF/NCEP/ROMS file is in hPa but have to use Pa -> multiply with 100.0
! Get local nodes for the tile
      DO I=1, np_global
            IF(ipgl(I)%rank==myrank) THEN
                  itmp=ipgl(I)%id
                  windx(itmp) =  cf_add_offset_uwind + cf_scale_factor_uwind * uin(I)
                  windy(itmp) =  cf_add_offset_vwind + cf_scale_factor_vwind * vin(I)
                  pr(itmp) = 100.0_rkind * (cf_add_offset_pr + cf_scale_factor_pr * pin(I))
            ENDIF
      END DO
      IF(myrank==0) THEN
            WRITE(16,'(A,F6.1,2x,F6.1)') '++ Uwind @ model min/max = ', cf_add_offset_vwind + cf_scale_factor_vwind * minval(uin), &
      &                                                                 cf_add_offset_vwind + cf_scale_factor_vwind * maxval(uin)
            WRITE(16,'(A,F6.1,2x,F6.1)') '++ Vwind @ model min/max = ', cf_add_offset_vwind + cf_scale_factor_vwind * minval(vin), &
      &                                                                 cf_add_offset_vwind + cf_scale_factor_vwind * maxval(vin)
            WRITE(16,'(A,F9.3,2x,F9.3)') '++ Pair  @ model min/max = ', cf_add_offset_pr + cf_scale_factor_pr * minval(pin), &
      &                                                                 cf_add_offset_pr + cf_scale_factor_pr * maxval(pin)
            WRITE(16,'(A, I0)') '++ Done with READ_REC_ATMO_FEM from UVP_direct.nc and record = ', rec
            FLUSH(16)
      ENDIF
      END SUBROUTINE READ_REC_ATMO_FEM
      
!**********************************************************************
!*                                                                    *
!**********************************************************************      
      SUBROUTINE GENERIC_NETCDF_ERROR(CallFct, idx, iret)
      USE NETCDF
! Generic netcdf error handling subroute
      implicit none
      integer, intent(in)           :: iret, idx
      character(*), intent(in)      :: CallFct
      character(len=500)            :: CHRERR
      IF (iret .NE. nf90_noerr) THEN
            CHRERR = nf90_strerror(iret)
            WRITE(16,*) TRIM(CallFct), ' -', idx, '-', CHRERR
      ENDIF
      END SUBROUTINE GENERIC_NETCDF_ERROR
