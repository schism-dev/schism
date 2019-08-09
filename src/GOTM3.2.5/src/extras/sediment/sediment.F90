!$Id: sediment.F90,v 1.10 2005-06-27 13:44:07 kbk Exp $
#include"cppdefs.h"
!-----------------------------------------------------------------------
!BOP
!
! !MODULE: sediment --- suspended sediment dynamics\label{sec:sediment}
!
! !INTERFACE:
   module sediment
!
! !DESCRIPTION:
!  This subroutine computes the transport of sediment, given
!  by its concentration, $C$. The transport equation for this
!  quantity can be written as
!  \begin{equation}
!   \label{CEq}
!    \partder{C}{t} +  \partder{w_s C}{z}
!    = {\cal D}_C
!    \comma
!  \end{equation}
!  where $w_s$ is the constant sinking speed of the sediment, and
!  ${\cal D}_C$ denotes the sum of the turbulent and viscous transport
!   terms modelled according to
!  \begin{equation}
!   \label{DC}
!    {\cal D}_C
!    = \frstder{z}
!     \left(
!        \nu_t \partder{C}{z}
!      \right)
!    \point
!  \end{equation}
!  For simplicity, we set the turbulent diffusivity of sediment
!  equal to $\nu_t$, the turbulent diffusivity of momentum, see
!   \sect{sec:turbulenceIntro}. Surface fluxes and inner sources or
!  sinks are not considered.
!
!  The sinking speed, $w_s$, is negative by definition, and may depend
!  on the grain diameter, molecular viscosity of water
!  and sediment density as discussed in \sect{sec:initSediment}.
!  Diffusion is discretized implicitly as discussed in \sect{sec:diffusionMean},
!  advection (settling) is treated explicitely,
!  see \sect{sec:advectionMean}.
!
!  There is an interesting stationary solution of \eq{CEq} for the case of
!  advection by vertical settling and mixing by diffusion to balance exactly. If
!  one considers only the region of a turbulent flow near to a rigid wall,
!  where the law-of-the-wall
!  relation $\nu_t = \kappa u_* (z+z_0)$ holds, it is easy to show that
!  \begin{equation}
!   \label{rouseProfile}
!     \dfrac{C}{C_0} = \left( \dfrac{z+z_0}{z_0} \right)^R
!  \end{equation}
!  is a solution of \eq{CEq}. Here, $C_0$ is the reference sediment
!  concentration at $z=0$ and $R=w_s/(\kappa u_*)$ is the so-called
!  Rouse number. \eq{rouseProfile} can be used to derive boundary
!  conditions for \eq{CEq}.
!
! !USES:
   use util
   use meanflow,     only:    depth,u_taub,gravity,rho_0,z0b,za
   use meanflow,     only:    h,avh,NN,buoy
   use meanflow,     only:    ho,w_grid,grid_method
   use turbulence,   only:    kappa,num,nuh
   use observations, only:    w_adv_discr,w_adv_method
   use output,       only:    out_fmt,write_results,ts
!
   IMPLICIT NONE
!
!  default: all is private
   private
!
!  !PUBLIC MEMBER FUNCTIONS:
   public init_sediment, do_sediment, end_sediment
!
!  !DEFINED PARAMETERS:
!  how to compute bottom sediment flux or concentration
   integer, parameter      :: NoFlux          = 1
   integer, parameter      :: SmithMcLean     = 2

! !REVISION HISTORY:
!  Original author(s): Hans Burchard & Karsten Bolding
!
!  $Log: sediment.F90,v $
!  Revision 1.10  2005-06-27 13:44:07  kbk
!  modified + removed traling blanks
!
!  Revision 1.9  2004/06/29 12:56:36  lars
!  removed tabs
!
!  Revision 1.8  2004/03/22 10:14:25  kbk
!  cleaned, store old index -> much faster, fixed conc. calc.
!
!  Revision 1.7  2004/03/04 09:34:54  kbk
!  selection between eularian and lagrangian solver
!
!  Revision 1.6  2004/01/13 10:20:21  lars
!  removed small bug in namelist
!
!  Revision 1.5  2004/01/13 10:00:52  lars
!  partical re-write using new adv.-diff. routines
!
!  Revision 1.4  2003/03/28 09:20:34  kbk
!  added new copyright to files
!
!  Revision 1.3  2003/03/28 08:24:01  kbk
!  removed tabs
!
!  Revision 1.2  2003/03/10 09:13:24  gotm
!  Improved documentation
!
!  Revision 1.1.1.1  2001/02/12 15:55:57  gotm
!  initial import into CVS
!
!  private data members from here on
!
!  the 'sedi' namelist
   logical         :: sedi_calc=.false.
   logical         :: sedi_eulerian=.true.
   logical         :: sedi_dens=.false.
   REALTYPE        :: rho_sed=2650.
   REALTYPE        :: size=62.5e-6
   REALTYPE        :: init_conc=0.0001
   integer         :: adv_method=1
   REALTYPE        :: cnpar=0.5
   integer         :: sedi_method
   integer         :: z0bMethod
   integer         :: sedi_npar=10000
   logical         :: take_mean=.true.
!
!  sediment concentration
   REALTYPE, dimension(:), allocatable :: C
!
!  sinking speed, observed cocentration, source term
   REALTYPE, dimension(:), allocatable :: wc,Cobs,Qsour
!
!  particle postions for lagragian simulation
   REALTYPE, dimension(:), allocatable :: zp
   integer,  dimension(:), allocatable :: zi
!
!  boundary condition types for diffusion and advection
   integer         :: DiffBcup
   integer         :: DiffBcdw
   integer         :: AdvBcup
   integer         :: AdvBcdw
!
!  boundary values for diffusion and advection
   REALTYPE        :: DiffCup
   REALTYPE        :: DiffCdw
   REALTYPE        :: AdvCup
   REALTYPE        :: AdvCdw

!  some parameter of the sediment model
   REALTYPE        :: ustarc,gs
!
!  output unit
   integer         :: out_unit
!EOP
!----------------------------------------------------------------------

   contains

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: Initialise the sediment module\label{sec:initSediment}
!
! !INTERFACE:
   subroutine init_sediment(namlst,fname,unit,nlev,g,rho_0)
!
! !DESCRIPTION:
!  This routine reads the sediment namelist from {\tt sediment.inp}
!  and allocates memory for the sediment-related vectors.
!  Further, depending on the sediment model, the settling velocity, $w_s$,
! and the critical friction velocity, $u_*^c$, are calculated here.
!  The settling velocity is based on a formula proposed
!  by \cite{Zanke77} which is valid for spherical particles.
!  The critical friction velocity is a function of the settling
!  velocity and the particle size (see \cite{SmithMcLean77}).
!
! !USES:
   IMPLICIT NONE
!
! !INPUT PARAMETERS:
   integer,intent(in)                  :: namlst
   character(len=*), intent(in)        :: fname
   integer,intent(in)                  :: unit
   integer,intent(in)                  :: nlev
   REALTYPE,intent(in)                 :: g,rho_0
!
! !REVISION HISTORY:
!  Original author(s): Karsten Bolding & Hans Burchard
!
! !LOCAL VARIABLES:
   integer          :: rc
   integer          :: i,n
   REALTYPE         :: x,Dsize,avmolu=1.3e-6
   namelist /sedi/ sedi_calc,sedi_eulerian,sedi_dens, &
                   rho_sed,size,init_conc, &
                   adv_method,cnpar,sedi_method,z0bMethod, &
                   sedi_npar,take_mean
   REALTYPE        :: zlev(0:nlev)
!EOP
!-----------------------------------------------------------------------
!BOC
   LEVEL1 'init_sediment'

!  Open and read the namelist
   open(namlst,file=fname,action='read',status='old',err=98)
   read(namlst,nml=sedi,err=99)
   close(namlst)

   if (sedi_calc) then
      out_unit = unit

      allocate(C(0:nlev),stat=rc)
      if (rc /= 0) STOP 'InitSediment: Error allocating (C)'
      C = _ZERO_

      allocate(wc(0:nlev),stat=rc)
      if (rc /= 0) STOP 'InitSediment: Error allocating (wc)'
      wc = _ZERO_

      allocate(Cobs(0:nlev),stat=rc)
      if (rc /= 0) STOP 'InitSediment: Error allocating (Cobs)'
      Cobs = _ZERO_

      allocate(Qsour(0:nlev),stat=rc)
      if (rc /= 0) STOP 'InitSediment: Error allocating (Qsour)'

      Qsour = _ZERO_
      C     = init_conc

      ! reduced gravity
      gs=g*(rho_sed-rho_0)/rho_0

      ! Zanke formula for fall velocity of sediment
      x = -10.0*avmolu/size*(sqrt(1+(0.01*gs*size**3)/avmolu/avmolu)-1.0)
      wc= x

      if (sedi_eulerian) then
         LEVEL2 'Using eulerian approach'
         select case(sedi_method)
            case(NoFlux)
               LEVEL2 'Assuming no net flux across the lowest interface'
            case(SmithMcLean)
               LEVEL2 'Computing sediment concentration in lowest box'
               LEVEL2 'according to Smith and McLean (1977)'
               Dsize=size*(gs/avmolu/avmolu)**0.3333333
               if (Dsize.gt.10.0) then
                  ustarc=-0.4*wc(1)
               else
                  ustarc=-4.0/Dsize*wc(1)
               endif
            case default
               FATAL "unkown method to compute sediment flux"
               stop "init sediment"
         end select
      else
         LEVEL2 'Using lagragian approach..'
         LEVEL3 sedi_npar,' particles'
         allocate(zp(sedi_npar),stat=rc)
         if (rc /= 0) stop 'init_sediment: Error allocating (zp)'
         allocate(zi(sedi_npar),stat=rc)
         if (rc /= 0) stop 'init_sediment: Error allocating (zi)'

         do n=1,sedi_npar
!           Equidist. particle distribution
            zp(n)=-depth+n/float(sedi_npar+1)*depth
         end do
         zlev(0)=-depth
         do n=1,nlev
            zlev(n)=zlev(n-1)+h(n)
         end do
         do n=1,sedi_npar
            do i=1,nlev
               if (zlev(i) .gt. zp(n)) EXIT
            end do
            zi(n)=i
         end do
      end if
   end if

   return
98 LEVEL2 'I could not open sediment.inp'
   LEVEL2 'Ill continue but set sedi_calc to false.'
   LEVEL2 'If thats not what you want you have to supply sediment.inp'
   LEVEL2 'See the Rouse example on www.gotm.net for a working sediment.inp'
   sedi_calc = .false.
   return
99 FATAL 'I could not read sediment.inp'
   stop 'init_sediment'
   end subroutine init_sediment
!EOC

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: do sediment calculations
!
! !INTERFACE:
   subroutine do_sediment(nlev,dt)
!
! !DESCRIPTION:
!
! !USES:
   IMPLICIT NONE
!
! !INPUT PARAMETERS:
   integer,  intent(in)                :: nlev
   REALTYPE, intent(in)                :: dt
!
! !REVISION HISTORY:
!  Original author(s): Karsten Bolding & Hans Burchard
!
! !LOCAL VARIABLES:
   integer         :: n
   REALTYPE        :: dcdz,rho_mean
   REALTYPE        :: rho(0:nlev)
!EOP
!-----------------------------------------------------------------------
!BOC
   if (sedi_calc) then
      if( sedi_eulerian ) then
         call sediment_eulerian(nlev,dt)
      else
         call sediment_lagrangian(nlev,dt)
      end if

      if (sedi_dens) then   ! Update buoyancy and NN.
         do n=1,nlev
            rho(n)  = rho_0*(1.0-buoy(n)/gravity)
            buoy(n) = -gravity*((1.0-C(n))*rho(n) &
                      +C(n)*rho_sed - rho_0)/rho_0
         end do
         do n=1,nlev-1
            rho_mean = 0.5*(rho(n+1)+rho(n))
            dcdz     = (C(n+1)-C(n))/(0.5*(h(n+1)+h(n)))
            NN(n)=(1-C(n))*NN(n) &
                 -gravity/rho_0*(rho_sed-rho_mean)*dcdz
         end do
      end if

      if (write_results) then
         call save_sediment()
      end if

   end if
   return
   end subroutine do_sediment
!EOC

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: Eularian sediment transport calculation\label{sec:calcSediment}
!
! !INTERFACE:
   subroutine sediment_eulerian(nlev,dt)
!
! !DESCRIPTION:
!  In this subroutine, the dynamical equation for sediment \eq{CEq}
!  including turbulent diffusion and settling of suspended matter  is updated.
!
!  The models to compute the boundary  conditions at the lowest grid box are
!  set by the parameter {\tt sedi\_method} in {\tt sediment.inp}.
!  Currently, there are the following models available in GOTM:
!  \begin{itemize}
!    \item zero-flux at the lowest interface ({\tt sedi\_method = 1}).
!    With this option, no sediment can
!    leave the domain. The net zero-flux boundary condition is implemented as zero-flux
!    for both, advection and diffusion. Note, that this boundary condition
!    results in the Rouse-profile \eq{rouseProfile}, however, with unkown
!    reference concentration $C_0$.
!    \item Dirichlet boundary condition suggested by \cite{SmithMcLean77}
!    ({\tt sedi\_method = 2}). These authors set the concentration
!    at the lowest interface to zero, unless
!    the friction velocity is larger than a threshold, $u_*^c$. Then, they
!    assume that the concentration at the lowest interface, $C_0$, is a
!    quadratic function of the ratio $u_*/u_*^c$. In GOTM, the Rouse profile
!    \eq{rouseProfile} is assumed to interpolate the value of $C_0$ (located
!   at the lowest interface) to the center of the lowest grid cell.
!  \end{itemize}
!
!  The sediment induced bottom roughness $z_a$ (also see \sect{sec:friction}) can
!  be either ignored ({\tt z0bMethod = 1}) or updated from empirical formulae
!  as suggested for example by \cite{SmithMcLean77} ({\tt z0bMethod = 2}).
!  If {\tt sedi\_dens = .true.}, the effect of sediment on the density
!  stratifiction (and hence on turbulence) is considered in the code.
!
! !USES:
   IMPLICIT NONE
!
! !INPUT PARAMETERS:
   integer,  intent(in)                :: nlev
   REALTYPE, intent(in)                :: dt

! !DEFINED PARAMETERS:
   REALTYPE, parameter                 :: g1 =  1.56E-3
   REALTYPE, parameter                 :: a0 = 26.3
!
! !REVISION HISTORY:
!  Original author(s): Karsten Bolding & Hans Burchard
!
! !LOCAL VARIABLES:
   REALTYPE                  :: CBott,Cbalg
   REALTYPE                  :: y,z,ya
   REALTYPE                  :: dcdz,rho_mean
   REALTYPE                  :: rho(0:nlev)
   REALTYPE,save             :: Cb
   REALTYPE                  :: RelaxT(0:nlev)
   integer                   :: i,flag
   REALTYPE                  :: Cint
   REALTYPE                  :: rouse
!EOP
!-----------------------------------------------------------------------
!BOC
   select case(sedi_method)
      case(NoFlux)

         za       = _ZERO_      ! no model to update za

         DiffBcup = Neumann
         DiffBcdw = Neumann
         DiffCup  = _ZERO_
         DiffCdw  = _ZERO_

         AdvBcup  = flux
         AdvBcdw  = flux
         AdvCup   = _ZERO_
         AdvCdw   = _ZERO_

      case(SmithMcLean)

         if (u_taub .ge. ustarc) then

            ! compute reference level
            if (z0bMethod.ne.1) then
               za = a0/gs*(u_taub**2-ustarc**2)
            else
               za = _ZERO_
            end if

            ! Compute reference concentration at the interface
            Cbott = g1*((u_taub/ustarc)**2-1)

            ! compute Rouse number
            rouse = wc(1)/kappa/u_taub
         else
            Cbott = _ZERO_
            za    = _ZERO_
            rouse = _ZERO_
         end if

         DiffBcup = Neumann
         DiffBcdw = Dirichlet
         DiffCup  = _ZERO_
         DiffCdw  = Cbott*((0.5*h(1)+z0b)/z0b)**rouse

         AdvBcup  = flux
         AdvBcdw  = oneSided     ! allow for sinking sediment to leave domain
         AdvCup   = _ZERO_
         AdvCdw   = _ZERO_       ! not used

      case default
         FATAL 'Invalid method for sedi_calc'
         stop  'sediment.F90'
   end select

!  Does not work for prescribed vertical current velocity and
!  adaptive grids!!

   if (w_adv_method .ne. 0) then
      FATAL 'w_adv_method=0 in obs.inp is required for sediment'
      stop 'sediment.F90'
   endif

   if (grid_method .eq. 3) then
      FATAL 'adaptive grids do not yet work with sediment.'
      stop 'sediment.F90'
   endif

   call advection_mean(nlev,dt,h,h,wc,AdvBcup,AdvBcdw,AdvCup,AdvCdw,adv_method,conserv,C)

   RelaxT = 1.e15  ! no relaxation to observed value

   call diffusion_mean(nlev,dt,cnpar,h,DiffBcup,DiffBcdw,DiffCup,DiffCdw,        &
        num,Qsour,RelaxT,Cobs,C)

   return
 end subroutine sediment_eulerian
 !EOC

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: Lagrangian sediment transport calculation\label{sec:calcSediment}
!
! !INTERFACE:
   subroutine sediment_lagrangian(nlev,dt)
!
! !DESCRIPTION:
!  In this subroutine, the dynamical equation for sediment \eq{CEq}
!  including turbulent diffusion and settling of suspended matter is
!  this routine uses a lagrangian approach.
!  Should contain more info (KBK).
!
! !USES:
   IMPLICIT NONE
!
! !INPUT PARAMETERS:
   integer,  intent(in)                :: nlev
   REALTYPE, intent(in)                :: dt

! !REVISION HISTORY:
!  Original author(s): Karsten Bolding & Hans Burchard
!
! !LOCAL VARIABLES:
   integer         :: n,np
   integer, save   :: count=0
   logical, save   :: set_C_zero=.true.
   logical         :: active(sedi_npar)
   REALTYPE        :: zlev(0:nlev)
!EOP
!-----------------------------------------------------------------------
!BOC
   zlev(0)=-depth
   do n=1,nlev
      zlev(n)=zlev(n-1)+h(n)
   end do

   active=.true.
   call lagrange(nlev,dt,zlev,nuh,wc(1),sedi_npar,active,zi,zp)
   active=.true.
   if (write_results .or. take_mean) then
      if (set_C_zero) then
         C=_ZERO_
         set_C_zero=.false.
      end if

      do np=1,sedi_npar
         if (active(np)) then
            n=zi(np)
            C(n) = C(n)+_ONE_
            active(np)=.false.
         end if
      end do
      if (take_mean) then
         count=count+1
      else
         count=1
      end if
      if (write_results) then
         do n=1,nlev
            C(n) = init_conc*C(n)/sedi_npar*depth/h(n)/count
         end do
         count=0
         set_C_zero=.true.
      end if
   end if

   return
   end subroutine sediment_lagrangian
!EOC

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: Finish  sediment calculations
!
! !INTERFACE:
   subroutine end_sediment
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
   end subroutine end_sediment
!EOC

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: Storing the results
!
! !INTERFACE:
   subroutine save_sediment
!
! !DESCRIPTION:
!  Here, the storing of the sediment profiles to an ascii or a
!  netCDF file is managed.
!
! !USES:
   use output,  only:    out_fmt
#ifdef NETCDF_FMT
   use ncdfout, only:    ncid
   use ncdfout, only:    lon_dim,lat_dim,z_dim,time_dim,dims
   use ncdfout, only:    define_mode,new_nc_variable,set_attributes,store_data
#endif
   IMPLICIT NONE
#ifdef NETCDF_FMT
   include "netcdf.inc"
#endif
!
! !REVISION HISTORY:
!  Original author(s): Karsten Bolding & Hans Burchard
!
!EOP
!
! !LOCAL VARIABLES:
   logical, save             :: first=.true.
   integer, save             :: sedi_id,n
   integer                   :: i,iret
   REALTYPE                  :: z
!
!-----------------------------------------------------------------------
!BOC
   select case (out_fmt)
      case (ASCII)
         if(first) then
            open(out_unit,file='sediment.out',status='unknown')
            n = ubound(C,1)
            first = .false.
         end if
         write(out_unit,*) trim(ts)
         z = _ZERO_
         do i=n,1,-1
            z=z-0.5*h(i)
            write(out_unit,*) z,C(i)
            z=z-0.5*h(i)
         end do
      case (NETCDF)
#ifdef NETCDF_FMT
         if(first) then
            dims(1) = lon_dim
            dims(2) = lat_dim
            dims(3) = z_dim
            dims(4) = time_dim
            iret = define_mode(ncid,.true.)
            iret = new_nc_variable(ncid,'sediment',NF_REAL,4,dims,sedi_id)
            iret = set_attributes(ncid,sedi_id,units='-',  &
                                  long_name='sediment concentration')
            iret = define_mode(ncid,.false.)
            n = ubound(C,1)
            first = .false.
         end if
         iret = store_data(ncid,sedi_id,XYZT_SHAPE,n,array=C)
#endif
      case default
         FATAL 'A non valid output format has been chosen'
         stop 'save_sediment'
   end select
   return
   end subroutine save_sediment
!EOC

!-----------------------------------------------------------------------

   end module sediment

!-----------------------------------------------------------------------
! Copyright by the GOTM-team under the GNU Public License - www.gnu.org
!-----------------------------------------------------------------------
