!     Prototye for converting sflux data into .nc; can be driven by gen_atmos.pl
!     Output: sflux_[air,rad,prc]_1.???.nc

!     Compile with: 
!     ifort -assume byterecl -O2 -o gen_atmos gen_atmos.f90 -I$NETCDF_FORTRAN/include -L$NETCDF_FORTRAN/lib -lnetcdff
      program fgennc
      implicit none
      include 'netcdf.inc'

! error status return
      integer  iret
! netCDF id
      integer  ncid
! dimension ids
      integer  nx_grid_dim
      integer  ny_grid_dim
      integer  time_dim
! dimension lengths
!      integer  nx_grid_len
!      integer  ny_grid_len
      integer  time_len
!      parameter (nx_grid_len = 349)
!      parameter (ny_grid_len = 277)
      parameter (time_len = NF_UNLIMITED)
! variable ids
      integer  time_id
      integer  lon_id
      integer  lat_id
      integer  uwind_id
      integer  vwind_id
      integer  prmsl_id
      integer  stmp_id
      integer  spfh_id
      integer  dlwrf_id
      integer  dswrf_id
      integer  prate_id

! rank (number of dimensions) for each variable
      integer  time_rank
      integer  lon_rank
      integer  lat_rank
      integer  uwind_rank
      integer  vwind_rank
      integer  prmsl_rank
      integer  stmp_rank
      integer  spfh_rank
      integer  dlwrf_rank
      integer  dswrf_rank
      integer  prate_rank
      parameter (time_rank = 1)
      parameter (lon_rank = 2)
      parameter (lat_rank = 2)
      parameter (uwind_rank = 3)
      parameter (vwind_rank = 3)
      parameter (prmsl_rank = 3)
      parameter (stmp_rank = 3)
      parameter (spfh_rank = 3)
      parameter (dlwrf_rank = 3)
      parameter (dswrf_rank = 3)
      parameter (prate_rank = 3)

! variable shapes
      integer  time_dims(time_rank)
      integer  lon_dims(lon_rank)
      integer  lat_dims(lat_rank)
      integer  uwind_dims(uwind_rank)
      integer  vwind_dims(vwind_rank)
      integer  prmsl_dims(prmsl_rank)
      integer  stmp_dims(stmp_rank)
      integer  spfh_dims(spfh_rank)
      integer  dlwrf_dims(dlwrf_rank)
      integer  dswrf_dims(dswrf_rank)
      integer  prate_dims(prate_rank)

! data variables
      real, allocatable ::  lon(:,:),lat(:,:),times(:),data1(:,:)
! attribute vectors
      integer  intval(4)

!     My own variables
      character(len=14) fname
      character(len=21) date_in
      character(len=4) year_char
      character(len=2) mon_char
      character(len=2) day_char
      integer num_times !# of time records in each file
      integer ni,nj,i_time,data_start(3), data_count(3),i,j,k,ifiletype
!     CWB RC NFS
      integer nday,iday,iii,jjj,iiday
      character*34 fn1,fn2,fn3
      real, allocatable :: pre(:,:),uu(:,:),vv(:,:),sst(:,:),t2k(:,:)
      real, allocatable :: sph(:,:),rh(:,:),dsw(:,:),dlw(:,:),prc(:,:)
      real, allocatable :: uout(:,:),vout(:,:)
      real es,ea,mr
      character*27 outfn
      integer mday(12)
      character*6 infn
    

!     Inputs
      open(10,file='gen_atmos.in',status='old')
      read(10,*) ifiletype !1: air; 2: rad; 3: prc
      read(10,'(a)')fname !output name; e.g., '../sflux/sflux_air_1.002.nc'
      read(10,*) intval(1:3) !year (e.g. 1990); month (no leading 0s); day (no leading 0s)
      read(10,*) nday  !how many days to convert
      close(10)

      if(ifiletype<1.or.ifiletype>3) stop "Unknown ifiletype"
      intval(4)=0 !starting hour offset from UTC
      write(year_char,'(i4)')intval(1)
      write(mon_char,'(i2.2)')intval(2)
      write(day_char,'(i2.2)')intval(3)
      date_in='days since '//year_char//'-'//mon_char//'-'//day_char !in UTC

!     num_times=24 !# of time records in each file
      ni=791 !dimension in longitude
      nj=403 !dimension in lattitude

!     define mday
      mday(1)=31
      mday(2)=28
      if (mod(intval(1),4).eq.0) mday(2)=29
      mday(3)=31
      mday(4)=30
      mday(5)=31
      mday(6)=30
      mday(7)=31
      mday(8)=31
      mday(9)=30
      mday(10)=31
      mday(11)=30
      mday(12)=31

!     num_times=4*mday(intval(2))!+1!# of time records in each file
      num_times=24 !5160
      write(*,*) num_times
!     if (ifiletype.eq.2) num_times=24*(mday(intval(2))+5)+1

!     Allocate arrays
      allocate(lon(ni,nj),lat(ni,nj),data1(ni,nj),times(num_times)) 
      allocate(pre(ni,nj),uu(ni,nj),vv(ni,nj),sst(ni,nj),t2k(ni,nj))
      allocate(sph(ni,nj),rh(ni,nj),dsw(ni,nj),dlw(ni,nj),prc(ni,nj)) 
      allocate(uout(ni,nj),vout(ni,nj))

!     Write times (>=0 but can exceed 1)
!     times=(/0., 0.0416667, 0.5, 0.75/) !in days; UTC
      do i=1,num_times
         times(i)=(i-1)*1./24.+00./24 !1800 run,change this for different run
      end do

!     Write lon/lat (degrees_east and degrees_north)
      open(9,file='./lonlat.dat')
      do j=1,nj
        do i=1,ni
          read(9,*) lon(i,j),lat(i,j)
        enddo !j
      enddo !i
      close(9)

!  open filelist
!     write(infn,'(i4,i2.2)') intval(1),intval(2)
      infn='EC2018'
      open(20,file=infn//'_prmsl.str',form="unformatted")
      open(21,file=infn//'_uwind.str',form="unformatted")
      open(22,file=infn//'_vwind.str',form="unformatted")
!      open(23,file='')
      open(24,file=infn//'_stmp.str',form="unformatted")
      open(25,file=infn//'_spfh.str',form="unformatted")
      open(26,file=infn//'_dswrf.str',form="unformatted")
      open(27,file=infn//'_dlwrf.str',form="unformatted")
      open(28,file=infn//'_prc.str',form="unformatted")

!  do loop for different input

      do iday=1,nday

      write(outfn,'(a14,i4.4,a3)') fname,iday,'.nc'

      do iiday=1,1
      if (iday/=1) then
        intval(3)=intval(3)+1
        if (mday(intval(2)).lt.intval(3)) then
         intval(3)=1
         intval(2)=intval(2)+1
        end if
        if (intval(2).eq.13) then
         intval(2)=1
         intval(1)=intval(1)+1
        end if
      write(year_char,'(i4)')intval(1)
      write(mon_char,'(i2.2)')intval(2)
      write(day_char,'(i2.2)')intval(3)
      date_in='days since '//year_char//'-'//mon_char//'-'//day_char !in UTC
      endif
      end do

! enter define mode
      iret = nf_create(trim(outfn), NF_CLOBBER, ncid)
      call check_err(iret)
! define dimensions
      iret = nf_def_dim(ncid, 'nx_grid', ni, nx_grid_dim)
      call check_err(iret)
      iret = nf_def_dim(ncid, 'ny_grid', nj, ny_grid_dim)
      call check_err(iret)
      iret = nf_def_dim(ncid, 'time', NF_UNLIMITED, time_dim)
      call check_err(iret)
! define variables
      time_dims(1) = time_dim
      iret = nf_def_var(ncid, 'time', NF_REAL, time_rank, time_dims, time_id)
      call check_err(iret)
      lon_dims(2) = ny_grid_dim
      lon_dims(1) = nx_grid_dim
      iret = nf_def_var(ncid, 'lon', NF_REAL, lon_rank, lon_dims, lon_id)
      call check_err(iret)
      lat_dims(2) = ny_grid_dim
      lat_dims(1) = nx_grid_dim
      iret = nf_def_var(ncid, 'lat', NF_REAL, lat_rank, lat_dims, lat_id)
      call check_err(iret)
      uwind_dims(3) = time_dim
      uwind_dims(2) = ny_grid_dim
      uwind_dims(1) = nx_grid_dim

      if(ifiletype==1) then
!       air
        iret = nf_def_var(ncid, 'uwind', NF_REAL, uwind_rank, uwind_dims,uwind_id)
        call check_err(iret)
        vwind_dims(3) = time_dim
        vwind_dims(2) = ny_grid_dim
        vwind_dims(1) = nx_grid_dim
        iret = nf_def_var(ncid, 'vwind', NF_REAL, vwind_rank, vwind_dims,vwind_id)
        call check_err(iret)
        prmsl_dims(3) = time_dim
        prmsl_dims(2) = ny_grid_dim
        prmsl_dims(1) = nx_grid_dim
        iret = nf_def_var(ncid, 'prmsl', NF_REAL, prmsl_rank, prmsl_dims,prmsl_id)
        call check_err(iret)
        stmp_dims(3) = time_dim
        stmp_dims(2) = ny_grid_dim
        stmp_dims(1) = nx_grid_dim
        iret = nf_def_var(ncid, 'stmp', NF_REAL, stmp_rank, stmp_dims, stmp_id)
        call check_err(iret)
        spfh_dims(3) = time_dim
        spfh_dims(2) = ny_grid_dim
        spfh_dims(1) = nx_grid_dim
        iret = nf_def_var(ncid, 'spfh', NF_REAL, spfh_rank, spfh_dims, spfh_id)
        call check_err(iret)
      else if(ifiletype==2) then
!       rad
        dlwrf_dims(3) = time_dim
        dlwrf_dims(2) = ny_grid_dim
        dlwrf_dims(1) = nx_grid_dim
        iret = nf_def_var(ncid, 'dlwrf', NF_REAL, dlwrf_rank, dlwrf_dims,dlwrf_id)
        call check_err(iret)
        dswrf_dims(3) = time_dim
        dswrf_dims(2) = ny_grid_dim
        dswrf_dims(1) = nx_grid_dim
        iret = nf_def_var(ncid, 'dswrf', NF_REAL, dswrf_rank, dswrf_dims,dswrf_id)
        call check_err(iret)
      else
!       prcp
        prate_dims(3) = time_dim
        prate_dims(2) = ny_grid_dim
        prate_dims(1) = nx_grid_dim
        iret = nf_def_var(ncid, 'prate', NF_REAL, prate_rank, prate_dims,prate_id)
        call check_err(iret)
      endif

! assign attributes
      iret = nf_put_att_text(ncid, time_id, 'long_name', 4, 'Time')
      call check_err(iret)
      iret = nf_put_att_text(ncid, time_id, 'standard_name', 4, 'time')
      call check_err(iret)
      iret = nf_put_att_text(ncid, time_id, 'units', 21, date_in)
!'
      call check_err(iret)
      iret = nf_put_att_int(ncid, time_id, 'base_date', NF_INT, 4, intval)
      call check_err(iret)
      iret = nf_put_att_text(ncid, lon_id, 'long_name', 9, 'Longitude')
      call check_err(iret)
      iret = nf_put_att_text(ncid, lon_id, 'standard_name', 9, 'longitude')
      call check_err(iret)
      iret = nf_put_att_text(ncid, lon_id, 'units', 12, 'degrees_east')
      call check_err(iret)
      iret = nf_put_att_text(ncid, lat_id, 'long_name', 8, 'Latitude')
      call check_err(iret)
      iret = nf_put_att_text(ncid, lat_id, 'standard_name', 8, 'latitude')
      call check_err(iret)
      iret = nf_put_att_text(ncid, lat_id, 'units', 13, 'degrees_north')
      call check_err(iret)

      if(ifiletype==1) then !air
        iret = nf_put_att_text(ncid, uwind_id, 'long_name', 39, 'Surface Eastward Air Velocity (10m AGL)')
        call check_err(iret)
        iret = nf_put_att_text(ncid, uwind_id, 'standard_name', 13, 'eastward_wind')
        call check_err(iret)
        iret = nf_put_att_text(ncid, uwind_id, 'units', 3, 'm/s')
        call check_err(iret)
        iret = nf_put_att_text(ncid, vwind_id, 'long_name', 40, 'Surface Northward Air Velocity (10m AGL)')
        call check_err(iret)
        iret = nf_put_att_text(ncid, vwind_id, 'standard_name', 14, 'northward_wind')
        call check_err(iret)
        iret = nf_put_att_text(ncid, vwind_id, 'units', 3, 'm/s')
        call check_err(iret)
        iret = nf_put_att_text(ncid, prmsl_id, 'long_name', 23, 'Pressure reduced to MSL')
        call check_err(iret)
        iret = nf_put_att_text(ncid, prmsl_id, 'standard_name', 25, 'air_pressure_at_sea_level')
        call check_err(iret)
        iret = nf_put_att_text(ncid, prmsl_id, 'units', 2, 'Pa')
        call check_err(iret)
        iret = nf_put_att_text(ncid, stmp_id, 'long_name', 32, 'Surface Air Temperature (2m AGL)')
        call check_err(iret)
        iret = nf_put_att_text(ncid, stmp_id, 'standard_name', 15, 'air_temperature')
        call check_err(iret)
        iret = nf_put_att_text(ncid, stmp_id, 'units', 1, 'K')
        call check_err(iret)
        iret = nf_put_att_text(ncid, spfh_id, 'long_name', 34, 'Surface Specific Humidity (2m AGL)')
        call check_err(iret)
        iret = nf_put_att_text(ncid, spfh_id, 'standard_name', 17, 'specific_humidity')
        call check_err(iret)
        iret = nf_put_att_text(ncid, spfh_id, 'units', 1, '1')
        call check_err(iret)

      else if(ifiletype==2) then !rad

        iret = nf_put_att_text(ncid, dlwrf_id, 'long_name', 33, 'Downward Long Wave Radiation Flux')
        call check_err(iret)
        iret = nf_put_att_text(ncid, dlwrf_id, 'standard_name', 40, 'surface_downwelling_longwave_flux_in_air')
        call check_err(iret)
        iret = nf_put_att_text(ncid, dlwrf_id, 'units', 5, 'W/m^2')
        call check_err(iret)
        iret = nf_put_att_text(ncid, dswrf_id, 'long_name', 34, 'Downward Short Wave Radiation Flux')
        call check_err(iret)
        iret = nf_put_att_text(ncid, dswrf_id, 'standard_name', 41, 'surface_downwelling_shortwave_flux_in_air')
        call check_err(iret)
        iret = nf_put_att_text(ncid, dswrf_id, 'units', 5, 'W/m^2')
        call check_err(iret)

      else !prcp
        iret = nf_put_att_text(ncid, prate_id, 'long_name', 26, 'Surface Precipitation Rate')
        call check_err(iret)
        iret = nf_put_att_text(ncid, prate_id, 'standard_name', 18, 'precipitation_flux')
        call check_err(iret)
        iret = nf_put_att_text(ncid, prate_id, 'units', 8, 'kg/m^2/s')
        call check_err(iret)

      endif

      iret = nf_put_att_text(ncid, NF_GLOBAL, 'Conventions', 6, 'CF-1.0')
      call check_err(iret)
!     leave define mode
      iret = nf_enddef(ncid)
      call check_err(iret)
       
!     Write times
      iret = nf_put_vara_real(ncid, time_id, 1, num_times, times)
      call check_err(iret)

!     Write lon/lat
      data_start(1:2)=1
      data_count(1)=ni
      data_count(2)=nj
      iret = nf_put_vara_real(ncid, lon_id, data_start, data_count, lon)
      iret = nf_put_vara_real(ncid, lat_id, data_start, data_count, lat)
      call check_err(iret)

!     Write record variables for time loop
      do i_time=1,num_times
        data_start(1) = 1
        data_start(2) = 1
        data_start(3) = i_time
        data_count(1) = ni
        data_count(2) = nj
        data_count(3) = 1

        if(ifiletype==1) then !air
!------------------------------------------------------------------------------
!       Calculate uwind: Surface Eastward Air Velocity (10m AGL) in m/s
!       store values in data1(1:ni,1:nj) (at lon(1:ni,1:nj) and lat(1:ni,1:nj)).
!       for invalid value use 9.999e+20

        
        read(21) uout
        read(22) vout
!       data1=uu
        data1=uout

        iret = nf_put_vara_real(ncid,uwind_id,data_start,data_count,data1)

!       Calculate vwind: Surface Northward Air Velocity (10m AGL) in m/s
!       store values in data1(1:ni,1:nj) (at lon(1:ni,1:nj) and lat(1:ni,1:nj)).
!       for invalid value use 9.999e+20
!       read(22,'(a)') fn3
!       open(32,file='./wrf02data/'//fn3)
!       read(32,'(6f12.4)') ((vv(i,j),i=1,ni),j=1,nj)
!       close(32)
!       data1=vv
        data1=vout

        iret = nf_put_vara_real(ncid,vwind_id,data_start,data_count,data1)

!       Calculate prmsl: Pressure reduced to MSL in Pa
!       store values in data1(1:ni,1:nj) (at lon(1:ni,1:nj) and lat(1:ni,1:nj)).
!       for invalid value use 9.999e+20
        read(20) pre
        pre(:,nj)=pre(:,nj-1) !to avoid ncl output error
        data1=pre

        iret = nf_put_vara_real(ncid,prmsl_id,data_start,data_count,data1)

!       Calculate stmp: Surface Air Temperature (2m AGL) in K
!       store values in data1(1:ni,1:nj) (at lon(1:ni,1:nj) and lat(1:ni,1:nj)).
!       for invalid value use 9.999e+20
!       do i=1,ni
!         do j=1,nj
!           data1(i,j)=298
!         enddo !j
!       enddo !i
        read(24) t2k
        data1=t2k

        iret = nf_put_vara_real(ncid,stmp_id,data_start,data_count,data1)

!       Calculate spfh: Surface Specific Humidity (2m AGL) (dimensionless)
!       store values in data1(1:ni,1:nj) (at lon(1:ni,1:nj) and lat(1:ni,1:nj)).
!       for invalid value use 9.999e+20
!       do i=1,ni
!         do j=1,nj
!           data1(i,j)=1.e-2
!         enddo !j
!       enddo !i
        read(25) sph

        data1=sph

        iret = nf_put_vara_real(ncid,spfh_id,data_start,data_count,data1)

!------------------------------------------------------------------------------
        else if(ifiletype==2) then !rad
!------------------------------------------------------------------------------
!       Calculate dlwrf: Downward Long Wave Radiation Flux in W/m^2
!       store values in data1(1:ni,1:nj) (at lon(1:ni,1:nj) and lat(1:ni,1:nj)).
!       for invalid value use 9.999e+20
!       do i=1,ni
!         do j=1,nj
!           data1(i,j)=abs(lon(i,j)+lat(i,j))/10 !400
!         enddo !j
!       enddo !i
        read(27) dlw
        data1=dlw

        iret = nf_put_vara_real(ncid,dlwrf_id,data_start,data_count,data1)

!       Calculate dswrf: Downward Short Wave Radiation Flux in W/m^2
!       store values in data1(1:ni,1:nj) (at lon(1:ni,1:nj) and lat(1:ni,1:nj)).
!       for invalid value use 9.999e+20
!       do i=1,ni
!         do j=1,nj
!           data1(i,j)=abs(lon(i,j)-lat(i,j))/10 !40
!         enddo !j
!       enddo !i
        read(26) dsw
        data1=dsw

        iret = nf_put_vara_real(ncid,dswrf_id,data_start,data_count,data1)

!------------------------------------------------------------------------------
        else  !prcp
!------------------------------------------------------------------------------
!       Calculate prate: Surface Precipitation Rate in kg/m^2/s
!       store values in data1(1:ni,1:nj) (at lon(1:ni,1:nj) and lat(1:ni,1:nj)).
!       for invalid value use 9.999e+20
!       do i=1,ni
!         do j=1,nj
!           data1(i,j)=abs(lon(i,j)+lat(i,j))/10 !1.e-6
!         enddo !j
!       enddo !i
        read(28) prc
        data1=prc

        iret = nf_put_vara_real(ncid,prate_id,data_start,data_count,data1)

!------------------------------------------------------------------------------
        endif
      enddo !i_time=1,num_times
       
      iret = nf_close(ncid)
      call check_err(iret)

      end do  !iday loop

      stop
      end
       
      subroutine check_err(iret)
      integer iret
      include 'netcdf.inc'
      if (iret .ne. NF_NOERR) then
      print *, nf_strerror(iret)
      stop
      endif
      end

!     Program getqs.f
!     Calculates the saturation specific humidity for water
!     vapor in air for a given temperature in deg. C. and
!     air pressure in mbars (hPa).
!     Uses sixth-order polynomial fit from Flatau et al. (1992).
!     Note, this code gives a value of es = 40.42 hPa at T = 29 deg. C 
!     and 1 atmosphere.
      subroutine getqs(temp,pres,es)
      implicit none
      integer*2 i
      real*4 ac(7),temp,rm,qs
      real*4 es,ws,pres
!     Initialize coefficients:
      data ac(1) /6.1117675/
      data ac(2) /0.44398606/
      data ac(3) /1.430533e-02/
      data ac(4) /2.650272e-04/
      data ac(5) /3.02247e-06/
      data ac(6) /2.038863e-08/
      data ac(7) /6.3878097e-10/
!
!     Mass ratio of H2O to AIR :
      rm = 18.015/28.965
!     Saturation vapor pressure in hPa = es
      es = ac(1)
      do i = 2,7
        es = es + ac(i)*temp**(i-1)
      enddo
!     Correction for sea water :
!     es = 0.980*es
!     Water vapor mixing ratio (saturation) :
!     ws = rm*es/(pres-es)
!     Specific humidity (saturation):
!     qs = rm*es/(pres-(1.0-rm)*es)
      return
      end

        subroutine uv2uv(uin,vin,uout,vout,lon,m,n,it)

!       it=0
!         convert (uin,vin) at model grid to (uout,vout) at lat, lon
!       it=1
!         convert (uin,vin) at lat/lon to (uout,vout) at model grid

        real lon(m,n),a,b,stdlta,stdltb,stdlon,pi,d2r,radius,gcon
        real uin(m,n),vin(m,n),uout(m,n),vout(m,n)

        stdlta=10.
        stdltb=40.
        stdlon=120.

        pi  = 4.0*atan(1.0)
        d2r = pi/180.0

        radius = 6371229.0
            gcon = (log(sin((90.0-stdlta)*d2r)) &
                 -log(sin((90.0-stdltb)*d2r))) &
                  /(log(tan((90.0-stdlta)*0.5*d2r)) &
                   -log(tan((90.0-stdltb)*0.5*d2r)))

        if(it.eq.0)then
!         convert (uin,vin) at model grid to (uout,vout) at lat, lon
          do i=1,m
          do j=1,n
            a=sin(gcon*(lon(i,j)-stdlon)*d2r)
            b=cos(gcon*(lon(i,j)-stdlon)*d2r)
            Uout(i,j)=(b*uin(i,j)+a*vin(i,j))/(a**2+b**2)
            Vout(i,j)=(-1.*a*uin(i,j)+b*vin(i,j))/(a**2+b**2)
          enddo
          enddo
        else if(it.eq.1)then
!         convert (uin,vin) at lat/lon to (uout,vout) at model grid
          do i=1,m
          do j=1,n
            a=sin(gcon*(lon(i,j)-stdlon)*d2r)
            b=cos(gcon*(lon(i,j)-stdlon)*d2r)
            Uout(i,j)=b*uin(i,j)-a*vin(i,j)
            Vout(i,j)=a*uin(i,j)+b*vin(i,j)
          enddo
          enddo
        else
          print*,' Please input the it'
          stop
        endif

        return
        end

