!$Id: seagrass.F90,v 1.5 2005-06-27 13:44:07 kbk Exp $
#include"cppdefs.h"
!-----------------------------------------------------------------------
!BOP
!
! !MODULE: seagrass --- sea grass dynamics \label{sec:seagrass}
!
! !INTERFACE:
   module seagrass
!
! !DESCRIPTION:
! In this module, seagrass canopies are treated as Lagrangian tracers,
! which either advect passively with the horizontal current speed or
! rest at their excursion limits and thus exert friction on the mean flow,
! see \cite{VerduinBackhaus2000}.
! Turbulence generation due to seagrass friction is possible, see
! namelist file {\tt seagrass.inp}. The extra production term
! in the balance of TKE, \eq{tkeA}, is included as described in
! \sect{sec:production}.

!
! !USES:
   use meanflow, only:     u,v,h,drag,xP
   use output,   only:     out_fmt,write_results,ts

!  default: all is private.
   private

!
! !PUBLIC MEMBER FUNCTIONS:
   public init_seagrass, calc_seagrass, end_seagrass

!
! !REVISION HISTORY:!
!  Original author(s): Hans Burchard & Karsten Bolding
!  $Log: seagrass.F90,v $
!  Revision 1.5  2005-06-27 13:44:07  kbk
!  modified + removed traling blanks
!
!  Revision 1.4  2003/03/28 09:20:34  kbk
!  added new copyright to files
!
!  Revision 1.3  2003/03/28 08:28:36  kbk
!  removed tabs
!
!  Revision 1.2  2003/03/10 09:13:09  gotm
!  Improved documentation
!
!  Revision 1.1.1.1  2001/02/12 15:55:57  gotm
!  initial import into CVS
!
!
!EOP
!-----------------------------------------------------------------------
!
!  private data members
   REALTYPE, dimension(:), allocatable :: xx,yy
   REALTYPE, dimension(:), allocatable :: exc,vfric,grassz

!  from a namelist
   logical                   :: grass_calc=.false.
   character(len=PATH_MAX)   :: grassfile='seagrass.dat'
   REALTYPE                  :: XP_rat
   integer                   :: grassind
   integer                   :: grassn
   integer                   :: out_unit

!-----------------------------------------------------------------------

   contains

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: Initialise the sea grass module
!
! !INTERFACE:
   subroutine init_seagrass(namlst,fname,unit,nlev,h)
!
! !DESCRIPTION:
! Here, the seagrass namelist {\tt seagrass.inp} is read
! and memory is allocated
! for some relevant vectors. Afterwards, excursion limits and friction
! coefficients are read from a file. The uppermost grid related index
! for the seagrass canopy is then calculated.

!
! !USES:
   IMPLICIT NONE

!
! !INPUT PARAMETERS:
   integer,          intent(in)   :: namlst
   character(len=*), intent(in)   :: fname
   integer,          intent(in)   :: unit
   integer,          intent(in)   :: nlev
   REALTYPE,         intent(in)   :: h(0:nlev)

!
! !REVISION HISTORY:
!  Original author(s): Karsten Bolding & Hans Burchard
!
!EOP
!
! !LOCAL VARIABLES:
   integer                   :: i,rc
   REALTYPE                  :: z
   namelist /canopy/  grass_calc,grassfile,XP_rat
!-----------------------------------------------------------------------
!BOC


   LEVEL1 'init_seagrass'

!  Open and read the namelist
   open(namlst,file=fname,action='read',status='old',err=98)
   read(namlst,nml=canopy,err=99)
   close(namlst)

   if (grass_calc) then
      out_unit=unit

      allocate(xx(0:nlev),stat=rc)
      if (rc /= 0) STOP 'init_seagrass: Error allocating (xx)'

      allocate(yy(0:nlev),stat=rc)
      if (rc /= 0) STOP 'init_seagrass: Error allocating (yy)'

      allocate(exc(0:nlev),stat=rc)
      if (rc /= 0) STOP 'init_seagrass: Error allocating (exc)'

      allocate(vfric(0:nlev),stat=rc)
      if (rc /= 0) STOP 'init_seagrass: Error allocating (vfric)'

      allocate(grassz(0:nlev),stat=rc)
      if (rc /= 0) STOP 'init_seagrass: Error allocating (grassz)'

      xx = -999.0
      yy = -999.0

      open(unit,status='unknown',file=grassfile)

      read(unit,*) grassn

      do i=1,grassn
         read(unit,*) grassz(i),exc(i),vfric(i)
      end do

      z=0.5*h(1)
      do i=2,nlev
         z=z+0.5*(h(i-1)+h(i))
         if (grassz(grassn).gt.z) grassind=i
      end do

      xx = _ZERO_
      yy = _ZERO_

      close(unit)

      LEVEL1 'Seagrass initialised ...'

   end if
   return

98 LEVEL2 'I could not open seagrass.inp'
   LEVEL2 'Ill continue but set grass_calc to false.'
   LEVEL2 'If thats not what you want you have to supply seagrass.inp'
   LEVEL2 'See the Seagrass example on www.gotm.net for a working seagrass.inp'
   grass_calc = .false.
   return
99 FATAL 'I could not read seagrass.inp'
   stop 'init_seagrass'
   end subroutine init_seagrass
!EOC

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: Update the sea grass model
!
! !INTERFACE:
   subroutine calc_seagrass(nlev,dt)
!
! !DESCRIPTION:
!
!  Here the time depending seagrass equation suggested by
!  \cite{VerduinBackhaus2000} is calculated. In order to
!  explain the basic principle, an idealised example is examined here
!  with a simplified momentum equation,
!
!  \begin{equation}
!  \partial_t u - \partial_z(\nu_t \partial_z u) = -g\partial_x\zeta-C_fu|u|
!  \comma
!  \end{equation}
!  and the Lagrangian tracer equation for seagrass,
!  \begin{equation}
!  \partial_t X =
!  \left\{
!  \begin{array}{ll}
!  u & \mbox{ for } |X|< X_{\max} \mbox{ or } X \cdot u <0,\\
!  0 & \mbox{ else}
!  \comma
!  \end{array}
!  \right.
!  \end{equation}
!  where $X$ is the Langrangian horizontal excursion of the seagrass.
!  The seagrass friction coefficient, $C_f$, is only non--zero at heights
!  where seagrass tracers are at their excursion limits:
!
!  \begin{equation}
!  C_f =
!  \left\{
!  \begin{array}{ll}
!  C_f^{\max} & \mbox{ for } |X|=X_{\max} \comma \\
!  0 & \mbox{ else} \point
!  \end{array}
!  \right.
!  \end{equation}
!
!  The maximum excursion limits $X_{\max}$ and the friction coefficients
!  $C_f^{\max}$ are read from a file.
!
!  The production of turbulence is calculated here as the sum of shear
!  production and friction loss at the seagrass leaves,
!  \begin{equation}
!   \label{sgProduction}
!    X_P = \alpha_{sg} C_f |u|^3
!   \comma
!  \end{equation}
!  which is added to the usual shear--production term as illustrated in
!  \eq{computeP}. The efficiency coefficient of turbulence production
!  by sea--grass friction, $\alpha_{sg}$, is denoted as {\tt xP\_rat}
!  in the code. It has to be read--in from the {\tt canopy} namelist.
!  For details and example calculations, see \cite{BurchardBolding2000}.

!
! !USES:
   IMPLICIT NONE

!
! !INPUT PARAMETERS:
   integer,  intent(in)     :: nlev
   REALTYPE, intent(in)     :: dt

!
! !REVISION HISTORY:
!  Original author(s): Karsten Bolding & Hans Burchard
!
!BOC
! !LOCAL VARIABLES:
   integer                   :: i
   REALTYPE                  :: dist
   REALTYPE                  :: grassfric(0:nlev)
   REALTYPE                  :: excur(0:nlev)
   REALTYPE                  :: z(0:nlev)
   REALTYPE                  :: xxP(0:nlev)
!EOP
!-----------------------------------------------------------------------

   if (grass_calc) then

      z(1)=0.5*h(1)
      do i=2,nlev
         z(i)=z(i-1)+0.5*(h(i-1)+h(i))
      end do

! Interpolate excursion limits and frcition to actual grid.

      call gridinterpol(grassn,1,grassz,exc,nlev,z,excur)
      call gridinterpol(grassn,1,grassz,vfric,nlev,z,grassfric)

      do i=1,nlev
         xx(i)=xx(i)+dt*u(i)                ! Motion of seagrass elements with
         yy(i)=yy(i)+dt*v(i)                ! mean flow.
         dist=sqrt(xx(i)*xx(i)+yy(i)*yy(i))
         if (dist .gt. excur(i)) then       ! Excursion limit reached
            xx(i)= excur(i)/dist * xx(i)
            yy(i)= excur(i)/dist * yy(i)

   ! Increased drag by seagrass friction
            drag(i)=drag(i)+grassfric(i)

   ! Extra turbulence production by seagrass friction
            xxP(i)=xP_rat*grassfric(i)*(sqrt(u(i)**2+v(i)**2))**3
         else
            xxP(i)=0.
         end if
      end do

   ! Interpolate onto turbulence grid points
      do i=1,nlev-1
         xP(i)=0.5*(xxP(i)+xxP(i+1))
      end do

      if (write_results) then
         call save_seagrass()
      end if

   end if

   return
   end subroutine calc_seagrass
!EOC
!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: Finish the sea grass calculations
!
! !INTERFACE:
   subroutine end_seagrass
!
! !DESCRIPTION:
!  Nothing done yet --- supplied for completeness.

!
! !USES:
   IMPLICIT NONE

!
! !REVISION HISTORY:
!  Original author(s): Karsten Bolding & Hans Burchard
!
!EOP
!-----------------------------------------------------------------------
!BOC

   return
   end subroutine end_seagrass
!EOC
!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: Storing the results
!
! !INTERFACE:
   subroutine save_seagrass
!
! !DESCRIPTION:
!  Here, storing of the sediment profiles to an ascii or a
!  netCDF file is managed.

!
! !USES:
   use output, only: out_fmt
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
! !REVISION HISTORY:
!  Original author(s): Karsten Bolding & Hans Burchard
!
!EOP
!
! !LOCAL VARIABLES:
   logical, save             :: first=.true.
   integer, save             :: x_excur_id,y_excur_id,n
   integer                   :: i,iret
   REALTYPE                  :: zz
   REALTYPE                  :: miss_val
!-----------------------------------------------------------------------
!BOC

   select case (out_fmt)
      case (ASCII)
         if(first) then
            open(out_unit,file='seagrass.out',status='unknown')
            first = .false.
         end if
         write(out_unit,*)
         write(out_unit,*) trim(ts)
         zz = _ZERO_
         do i=1,grassind
            zz=zz+0.5*h(i)
            write(out_unit,*) zz,xx(i),yy(i)
            zz=zz+0.5*h(i)
         end do
      case (NETCDF)
#ifdef NETCDF_FMT
         if(first) then
            dims(1) = lon_dim
            dims(2) = lat_dim
            dims(3) = z_dim
            dims(4) = time_dim
            miss_val = -999.0
            iret = define_mode(ncid,.true.)
            iret = new_nc_variable(ncid,'x-excur',NF_REAL,4,dims,x_excur_id)
            iret = set_attributes(ncid,x_excur_id,units='m',    &
                   long_name='seagrass excursion(x)',missing_value=miss_val)
            iret = new_nc_variable(ncid,'y-excur',NF_REAL,4,dims,y_excur_id)
            iret = set_attributes(ncid,y_excur_id,units='m',    &
                   long_name='seagrass excursion(y)',missing_value=miss_val)
            iret = define_mode(ncid,.false.)
            n = ubound(xx,1)
            first = .false.
         end if
         iret = store_data(ncid,x_excur_id,XYZT_SHAPE,n,array=xx)
         iret = store_data(ncid,y_excur_id,XYZT_SHAPE,n,array=yy)
#endif
      case default
         FATAL 'A non valid output format has been chosen'
         stop 'save_seagrass'
   end select
   return
   end subroutine save_seagrass
!EOC

!-----------------------------------------------------------------------

   end module seagrass

!-----------------------------------------------------------------------
! Copyright by the GOTM-team under the GNU Public License - www.gnu.org
!-----------------------------------------------------------------------
