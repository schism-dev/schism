!$Id: bio_save.F90,v 1.3 2004-07-30 09:22:20 hb Exp $
#include"cppdefs.h"
!-----------------------------------------------------------------------
!BOP
! !ROUTINE: Storing the results
!
! !INTERFACE:
   subroutine bio_save(nlev,h,totn)
!
! !DESCRIPTION:
!  Here, storing of the sediment profiles to an ascii or a
!  netCDF file is managed.
!
! !USES:
   use bio_var
   use output, only: out_fmt,ts
#ifdef NETCDF_FMT
   use ncdfout, only: ncid
   use ncdfout, only: lon_dim,lat_dim,z_dim,time_dim,dims
   use ncdfout, only: define_mode,new_nc_variable,set_attributes,store_data
#endif
   IMPLICIT NONE
#ifdef NETCDF_FMT
#include "netcdf.inc"
#endif
!
! !INPUT PARAMETERS:
   integer, intent(in)                 :: nlev
   REALTYPE, intent(in)                :: h(0:nlev)
   REALTYPE, intent(in)                :: totn
!
! !REVISION HISTORY:
!  Original author(s): Hans Burchard & Karsten Bolding
!
! !LOCAL VARIABLES:
   logical, save             :: first=.true.
   integer, save             :: nn
   integer, save             :: totn_id
   integer                   :: iret
   integer                   :: out_unit=67,rc
   REALTYPE                  :: zz
   integer                   :: i,j,n
!EOP
!-----------------------------------------------------------------------
!BOC
   select case (out_fmt)
      case (ASCII)
         if(first) then
            open(out_unit,file='bio.out',status='unknown')
            nn = ubound(cc(1,:),1)
            first = .false.
         end if
         write(out_unit,*)
         write(out_unit,*) trim(ts)
         zz = _ZERO_
         do i=nn,1,-1
            zz=zz+0.5*h(i)
            write(out_unit,115) zz,(cc(j,i) , j=1,numc)
            zz=zz+0.5*h(i)
         end do
115 format(F10.4,100(1x,E10.4E2))

      case (NETCDF)
#ifdef NETCDF_FMT
         if(first) then
            first = .false.
            dims(1) = lon_dim
            dims(2) = lat_dim
            dims(3) = z_dim
            dims(4) = time_dim

            iret = define_mode(ncid,.true.)

            do n=1,numc
               iret = new_nc_variable(ncid,var_names(n),NF_REAL, &
                                      4,dims,var_ids(n))
               iret = set_attributes(ncid,var_ids(n),       &
                                     units=var_units(n),    &
                                     long_name=var_long(n))
            end do

            nn = ubound(cc(1,:),1)

            dims(1) = time_dim
            iret = new_nc_variable(ncid,'totn',NF_REAL,1,dims,totn_id)
            iret = set_attributes(ncid,totn_id,units='mmol/m**2',    &
                   long_name='total N')

            iret = define_mode(ncid,.false.)
         end if

         do n=1,numc
            iret = store_data(ncid,var_ids(n),XYZT_SHAPE,nn,array=cc(n,:))
         end do
!KBK         iret = store_data(ncid,phy_id,XYZT_SHAPE,nn,array=cc(2,:)+P0)
!KBK         iret = store_data(ncid,zoo_id,XYZT_SHAPE,nn,array=cc(3,:)+Z0)

         iret = store_data(ncid,totn_id,T_SHAPE,1,scalar=totn)
#endif
      case default
         FATAL 'A non valid output format has been chosen'
         stop 'bio_save'
   end select

   return
   end subroutine bio_save
!EOC

!-----------------------------------------------------------------------
! Copyright by the GOTM-team under the GNU Public License - www.gnu.org
!-----------------------------------------------------------------------
