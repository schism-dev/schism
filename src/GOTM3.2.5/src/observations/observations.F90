!$Id: observations.F90,v 1.10 2005-08-15 11:54:01 hb Exp $
#include"cppdefs.h"
!-----------------------------------------------------------------------
!BOP
!
! !MODULE: observations --- the 'real' world \label{sec:observations}
!
! !INTERFACE:
   module observations
!
! !DESCRIPTION:
!  This module provides the necessary subroutines for communicating
!  `observations' to GOTM.
!  The module operates according to the general philosophy used in GOTM,
!  i.e.\ it provides {\tt init\_observ\-ations()} to be called in the overall
!  initialisation routine and {\tt get\_all\_obs()} to be called in the time
!  loop to actually obtain the `observations'.
!  In addition to these subroutines the module also provides two routines
!  for reading scalar-type observations and profile-type observations.
!  Each observation has a date stamp with the format {\tt yyyy-mm-dd hh:dd:mm}.
!  The module uses the {\tt time} module (see \sect{sec:time})
!  to convert the time string to the
!  internal time representation of GOTM.
!  Profiles are interpolated to the actual GOTM model grid.
!  Free format is used for reading-in the actual data.
!
! !USES:
   IMPLICIT NONE

!  default: all is private.
   private
!
! !PUBLIC MEMBER FUNCTIONS:
   public init_observations,get_all_obs,read_obs,read_profiles
!
! !PUBLIC DATA MEMBERS:
!
!  'observed' salinity profile
   REALTYPE, public, dimension(:), allocatable   :: sprof

!  'observed' temperature profile
   REALTYPE, public, dimension(:), allocatable   :: tprof

!  'observed' horizontal salinity  gradients
   REALTYPE, public, dimension(:), allocatable   :: dsdx,dsdy

!  'observed' horizontal temperarure  gradients
   REALTYPE, public, dimension(:), allocatable   :: dtdx,dtdy

!  internal horizontal pressure gradients
   REALTYPE, public, dimension(:), allocatable   :: idpdx,idpdy

!  horizontal velocity profiles
   REALTYPE, public, dimension(:), allocatable   :: uprof,vprof

!  observed profile of turbulent dissipation rates
   REALTYPE, public, dimension(:), allocatable   :: epsprof

!  ralaxation times for salinity and temperature
   REALTYPE, public, dimension(:), allocatable   :: SRelaxTau,TRelaxTau

!  sea surface elevation, sea surface gradients and height of velocity obs.
   REALTYPE, public          :: zeta=0.,dpdx=0.,dpdy=0.,h_press=0

!  vertical advection velocity
   REALTYPE, public          :: w_adv=0.,w_height

!  Parameters for water classification - default Jerlov type I
   REALTYPE, public          :: A=0.58,g1=0.35,g2=23.0

!------------------------------------------------------------------------------
!
! the following data are not all public,
! but have been included for clarity here
!
!------------------------------------------------------------------------------

!  Salinity profile(s)
   integer, public           :: s_prof_method=0
   integer, public           :: s_analyt_method=1
   character(LEN=PATH_MAX)   :: s_prof_file='sprof.dat'
   REALTYPE                  :: z_s1,s_1,z_s2,s_2
   REALTYPE                  :: s_obs_NN
   REALTYPE                  :: SRelaxTauM=0.
   REALTYPE                  :: SRelaxTauS=0.
   REALTYPE                  :: SRelaxTauB=0.
   REALTYPE                  :: SRelaxSurf=0.
   REALTYPE                  :: SRelaxBott=0.

!  Temperature profile(s)
   integer, public           :: t_prof_method=0
   integer, public           :: t_analyt_method=1
   character(LEN=PATH_MAX)   :: t_prof_file='tprof.dat'
   REALTYPE                  :: z_t1,t_1,z_t2,t_2
   REALTYPE                  :: t_obs_NN
   REALTYPE                  :: TRelaxTauM=0.
   REALTYPE                  :: TRelaxTauS=0.
   REALTYPE                  :: TRelaxTauB=0.
   REALTYPE                  :: TRelaxSurf=0.
   REALTYPE                  :: TRelaxBott=0.

!  External pressure - 'press' namelist
   integer, public           :: ext_press_method=0,PressMethod=0
   character(LEN=PATH_MAX)   :: ext_press_file=''
   REALTYPE, public          :: PressConstU=0.
   REALTYPE, public          :: PressConstV=0.
   REALTYPE, public          :: PressHeight=0.
   REALTYPE, public          :: PeriodM=44714.
   REALTYPE, public          :: AmpMu=0.
   REALTYPE, public          :: AmpMv=0.
   REALTYPE, public          :: PhaseMu=0.
   REALTYPE, public          :: PhaseMv=0.
   REALTYPE, public          :: PeriodS=43200.
   REALTYPE, public          :: AmpSu=0.
   REALTYPE, public          :: AmpSv=0.
   REALTYPE, public          :: PhaseSu=0.
   REALTYPE, public          :: PhaseSv=0.

!  Internal pressure - 'internal_pressure' namelist
   integer, public           :: int_press_method=0
   character(LEN=PATH_MAX)   :: int_press_file=''
   REALTYPE, public          :: const_dsdx=0.
   REALTYPE, public          :: const_dsdy=0.
   REALTYPE, public          :: const_dtdx=0.
   REALTYPE, public          :: const_dtdy=0.
   logical, public           :: s_adv=.false.
   logical, public           :: t_adv=.false.

!  Light extinction - the 'extinct' namelist
   integer                   :: extinct_method=1
   character(LEN=PATH_MAX)   :: extinct_file='extinction.dat'

!  Vertical advection velocity - 'w_advspec' namelist
   integer, public           :: w_adv_method=0
   REALTYPE, public          :: w_adv0=0.
   REALTYPE, public          :: w_adv_height0=0.
   character(LEN=PATH_MAX)   :: w_adv_file='w_adv.dat'
   integer, public           :: w_adv_discr=1

!  Sea surface elevations - 'zetaspec' namelist
   integer,public            :: zeta_method=0
   character(LEN=PATH_MAX)   :: zeta_file='zeta.dat'
   REALTYPE, public          :: zeta_0=0.
   REALTYPE, public          :: period_1=44714.
   REALTYPE, public          :: amp_1=0.
   REALTYPE, public          :: phase_1=0.
   REALTYPE, public          :: period_2=43200.
   REALTYPE, public          :: amp_2=0.
   REALTYPE, public          :: phase_2=0.

!  Observed velocity profile profiles - typically from ADCP
   integer                   :: vel_prof_method=0
   CHARACTER(LEN=PATH_MAX)   :: vel_prof_file='velprof.dat'
   REALTYPE, public          :: vel_relax_tau=3600.
   REALTYPE, public          :: vel_relax_ramp=86400.

!  Observed dissipation profiles
   integer                   :: e_prof_method=0
   REALTYPE                  :: e_obs_const=1.e-12
   CHARACTER(LEN=PATH_MAX)   :: e_prof_file='eprof.dat'

!  Buoyancy - 'bprofile' namelist
   REALTYPE, public          :: b_obs_surf=0.,b_obs_NN=0.
   REALTYPE, public          :: b_obs_sbf=0.

   REALTYPE,public, parameter:: pi=3.141592654

! !DEFINED PARAMETERS:

!  Unit numbers for reading observations/data.
   integer, parameter        :: s_prof_unit=30
   integer, parameter        :: t_prof_unit=31
   integer, parameter        :: ext_press_unit=32
   integer, parameter        :: int_press_unit=33
   integer, parameter        :: extinct_unit=34
   integer, parameter        :: w_adv_unit=35
   integer, parameter        :: zeta_unit=36
   integer, parameter        :: vel_prof_unit=37
   integer, parameter        :: e_prof_unit=38

!  pre-defined parameters
   integer, parameter        :: READ_SUCCESS=1
   integer, parameter        :: END_OF_FILE=-1
   integer, parameter        :: READ_ERROR=-2
   integer, parameter        :: NOTHING=0
   integer, parameter        :: ANALYTICAL=1
   integer, parameter        :: CONSTANT=1
   integer, parameter        :: FROMFILE=2
   integer, parameter        :: CONST_PROF=1
   integer, parameter        :: TWO_LAYERS=2
   integer, parameter        :: CONST_NN=3
!
! !REVISION HISTORY:
!  Original author(s): Karsten Bolding & Hans Burchard
!
!  $Log: observations.F90,v $
!  Revision 1.10  2005-08-15 11:54:01  hb
!  sequence of reading w_adv and w_height changed, w_adv_height0 introduced, documentation extended
!
!  Revision 1.9  2005/07/06 16:20:14  kbk
!  updated documentation - added const_NNT and const_NNS
!
!  Revision 1.8  2004/07/30 09:26:01  hb
!  Simple exponential light absorption added --> Wilfried Kuehn
!
!  Revision 1.7  2003/03/28 09:20:35  kbk
!  added new copyright to files
!
!  Revision 1.6  2003/03/28 08:08:21  kbk
!  removed tabs
!
!  Revision 1.5  2003/03/10 13:51:08  lars
!  changed intent(out) to intent(inout) for lines in read_profiles()
!
!  Revision 1.4  2003/03/10 08:51:58  gotm
!  Improved documentation and cleaned up code
!
!  Revision 1.3  2001/11/27 15:35:55  gotm
!  zeta_method now public - used by updategrid()
!
!  Revision 1.1.1.1  2001/02/12 15:55:58  gotm
!  initial import into CVS
!
!EOP
!
! !LOCAL VARIABLES:
   character(len=72)         :: cbuf
!
!-----------------------------------------------------------------------

   contains

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: Initialise the observation module
!
! !INTERFACE:
   subroutine init_observations(namlst,fn,julday,secs,                 &
                                depth,nlev,z,h,gravity,rho_0)
!
! !DESCRIPTION:
!  The {\tt init\_observations()} subroutine basically reads the {\tt obs.inp}
!  file with a number of different namelists and takes actions according
!  to the specifications in the different namelists.
!  In this routine also memory is allocated to hold the 'observations'.
!  Finally, all variables are initialised to sane values, either by
!  reading from files, by prescribing constant values, or by using analytical
!  expressions.
!
! !USES:
   IMPLICIT NONE
!
! !INPUT PARAMETERS:
   integer, intent(in)                 :: namlst
   character(len=*), intent(in)        :: fn
   integer, intent(in)                 :: julday,secs
   REALTYPE, intent(in)                :: depth
   integer, intent(in)                 :: nlev
   REALTYPE, intent(in)                :: z(0:nlev),h(0:nlev)
   REALTYPE, intent(in)                :: gravity,rho_0
!
!
! !REVISION HISTORY:
!  Original author(s): Karsten Bolding & Hans Burchard
!
!EOP
!
! !LOCAL VARIABLES:
   namelist /sprofile/                                          &
            s_prof_method,s_analyt_method,                      &
            z_s1,s_1,z_s2,s_2,s_prof_file,s_obs_NN,             &
            SRelaxTauM,SRelaxTauB,SRelaxTauS,                   &
            SRelaxBott,SRelaxSurf

   namelist /tprofile/                                          &
            t_prof_method,t_analyt_method,                      &
            z_t1,t_1,z_t2,t_2,t_prof_file,t_obs_NN,             &
            TRelaxTauM,TRelaxTauB,TRelaxTauS,                   &
            TRelaxBott,TRelaxSurf

   namelist /ext_pressure/                                      &
            ext_press_method,PressMethod,ext_press_file,        &
            PressConstU,PressConstV,PressHeight,                &
            PeriodM,AmpMu,AmpMv,PhaseMu,PhaseMv,                &
            PeriodS,AmpSu,AmpSv,PhaseSu,PhaseSv

   namelist /int_pressure/                                      &
            int_press_method,int_press_file,                    &
            const_dsdx,const_dsdy,const_dtdx,const_dtdy,        &
            s_adv,t_adv

   namelist /extinct/ extinct_method,extinct_file

   namelist /w_advspec/                                         &
            w_adv_method,w_adv_file,w_adv_height0,w_adv0,w_adv_discr

   namelist /zetaspec/                                          &
            zeta_method,zeta_file,zeta_0,                       &
            period_1,amp_1,phase_1,period_2,amp_2,phase_2

   namelist /velprofile/ vel_prof_method,vel_prof_file,         &
            vel_relax_tau,vel_relax_ramp

   namelist /eprofile/ e_prof_method,e_obs_const,e_prof_file

   namelist /bprofile/ b_obs_surf,b_obs_NN,b_obs_sbf

   integer                   :: rc,i
   REALTYPE                  :: ds,db
!
!-----------------------------------------------------------------------
!BOC
   LEVEL1 'init_observations'

   open(namlst,file=fn,status='old',action='read',err=80)
   read(namlst,nml=sprofile,err=81)
   read(namlst,nml=tprofile,err=82)
   read(namlst,nml=ext_pressure,err=83)
   read(namlst,nml=int_pressure,err=84)
   read(namlst,nml=extinct,err=85)
   read(namlst,nml=w_advspec,err=86)
   read(namlst,nml=zetaspec,err=87)
   read(namlst,nml=velprofile,err=88)
   read(namlst,nml=eprofile,err=89)
   read(namlst,nml=bprofile,err=90)
   close(namlst)

   allocate(sprof(0:nlev),stat=rc)
   if (rc /= 0) STOP 'init_observations: Error allocating (sprof)'
   sprof = 0.

   allocate(tprof(0:nlev),stat=rc)
   if (rc /= 0) STOP 'init_observations: Error allocating (tprof)'
   tprof = 0.

   allocate(dsdx(0:nlev),stat=rc)
   if (rc /= 0) STOP 'init_observations: Error allocating (dsdx)'
   dsdx = 0.

   allocate(dsdy(0:nlev),stat=rc)
   if (rc /= 0) STOP 'init_observations: Error allocating (dsdy)'
   dsdy = 0.

   allocate(dtdx(0:nlev),stat=rc)
   if (rc /= 0) STOP 'init_observations: Error allocating (dtdx)'
   dtdx = 0.

   allocate(dtdy(0:nlev),stat=rc)
   if (rc /= 0) STOP 'init_observations: Error allocating (dtdy)'
   dsdy = 0.

   allocate(idpdx(0:nlev),stat=rc)
   if (rc /= 0) STOP 'init_observations: Error allocating (idpdx)'
   idpdx = 0.

   allocate(idpdy(0:nlev),stat=rc)
   if (rc /= 0) STOP 'init_observations: Error allocating (idpdy)'
   idpdy = 0.

   allocate(SRelaxTau(0:nlev),stat=rc)
   if (rc /= 0) STOP 'init_observations: Error allocating (SRelaxTau)'
   SRelaxTau = 0.

   allocate(TRelaxTau(0:nlev),stat=rc)
   if (rc /= 0) STOP 'init_observations: Error allocating (TRelaxTau)'
   TRelaxTau = 0.

   allocate(uprof(0:nlev),stat=rc)
   if (rc /= 0) stop 'init_observations: Error allocating (uprof)'
   uprof = 0.

   allocate(vprof(0:nlev),stat=rc)
   if (rc /= 0) stop 'init_observations: Error allocating (vprof)'
   vprof = 0.

   allocate(epsprof(0:nlev),stat=rc)
   if (rc /= 0) stop 'init_observations: Error allocating (epsprof)'
   epsprof = _ZERO_

   db=_ZERO_
   ds=depth
   SRelaxTau(0)=SRelaxTauB
   TRelaxTau(0)=TRelaxTauB
   do i=1,nlev
      TRelaxTau(i)=TRelaxTauM
      SRelaxTau(i)=SRelaxTauM
      db=db+0.5*h(i)
      ds=ds-0.5*h(i)
      if (db.le.SRelaxBott) SRelaxTau(i)=SRelaxTauB
      if (ds.le.SRelaxSurf) SRelaxTau(i)=SRelaxTauS
      if (db.le.TRelaxBott) TRelaxTau(i)=TRelaxTauB
      if (ds.le.TRelaxSurf) TRelaxTau(i)=TRelaxTauS
      db=db+0.5*h(i)
      ds=ds-0.5*h(i)
      if ((s_prof_method.ne.0).and.(SRelaxTau(i).le.0.)) then
         LEVEL2 ''
         LEVEL2 '***************************************************'
         LEVEL2 'SRelaxTau at i=',i,' is not a positive value.'
         LEVEL2 'Please correct namelist.inp and rerun.'
         LEVEL2 'Program aborted.'
         LEVEL2 '***************************************************'
         stop 'init_observations'
      end if
      if ((t_prof_method.ne.0).and.(TRelaxTau(i).le.0.)) then
         LEVEL2 ''
         LEVEL2 '***************************************************'
         LEVEL2 'TRelaxTau at i=',i,' is not a positive value.'
         LEVEL2 'Please correct namelist.inp and rerun.'
         LEVEL2 'Program aborted.'
         LEVEL2 '***************************************************'
         stop 'init_observations'
      end if
   end do

!  The salinity profile
   select case (s_prof_method)
      case (NOTHING)
         sprof = _ZERO_
      case (ANALYTICAL)

         ! different ways to prescribe profiles analytically
         select case (s_analyt_method)
            case (CONST_PROF)
               sprof = s_1
            case (TWO_LAYERS)
               call analytical_profile(nlev,z,z_s1,s_1,z_s2,s_2,sprof)
            case (CONST_NN)

               if (.not.((t_prof_method       .eq. ANALYTICAL) .and.      &
                         (t_analyt_method .eq. CONST_PROF))   )  then
                  LEVEL2 ''
                  LEVEL2 '***************************************************'
                  LEVEL2 'For salinity profiles with NN=const. you have to   '
                  LEVEL2 'prescribe constant temperature.                    '
                  LEVEL2 'Please correct obs.inp and re-run.                 '
                  LEVEL2 'Program aborted.                                   '
                  LEVEL2 '***************************************************'
                  stop 'init_observations'
               endif

               call const_NNS(nlev,z,s_1,t_1,s_obs_NN,gravity,rho_0,sprof)
            case default
         end select

      case (FROMFILE)
         open(s_prof_unit,file=s_prof_file,status='unknown',err=101)
         LEVEL2 'Reading salinity profiles from:'
         LEVEL3 trim(s_prof_file)
         call get_s_profile(s_prof_unit,julday,secs,nlev,z)
      case default
   end select

!  The temperature profile
   select case (t_prof_method)
      case (NOTHING)
         tprof = 0.
      case (ANALYTICAL)

        ! different ways to prescribe profiles analytically
         select case (t_analyt_method)
         case (CONST_PROF)
               tprof = t_1
            case (TWO_LAYERS)
               call analytical_profile(nlev,z,z_t1,t_1,z_t2,t_2,tprof)
            case (CONST_NN)

               if (.not.((s_prof_method       .eq. ANALYTICAL) .and.      &
                         (s_analyt_method .eq. CONST_PROF))   )  then

                  LEVEL2 ''
                  LEVEL2 '***************************************************'
                  LEVEL2 'For temperature profiles with NN=const you have to '
                  LEVEL2 'prescribe constant salinity.                       '
                  LEVEL2 'Please correct obs.inp and re-run.                 '
                  LEVEL2 'Program aborted.                                   '
                  LEVEL2 '***************************************************'
                  stop 'init_observations'
               endif

               call const_NNT(nlev,z,t_1,s_1,t_obs_NN,gravity,rho_0,tprof)
            case default
         end select
      case (FROMFILE)
         open(t_prof_unit,file=t_prof_file,status='unknown',err=102)
         LEVEL2 'Reading temperature profiles from:'
         LEVEL3 trim(t_prof_file)
         call get_t_profile(t_prof_unit,julday,secs,nlev,z)
      case default
   end select

!  The external pressure
   select case (ext_press_method)
      case (FROMFILE)
         open(ext_press_unit,file=ext_press_file,status='unknown',err=103)
         LEVEL2 'Reading external pressure from:'
         LEVEL3 trim(ext_press_file)
      case default
   end select
   call get_ext_pressure(ext_press_method,ext_press_unit,julday,secs)

!  The internal pressure
   select case (int_press_method)
      case (CONSTANT)
         dsdx=const_dsdx
         dsdy=const_dsdy
         dtdx=const_dtdx
         dtdy=const_dtdy
      case (FROMFILE)
         open(int_press_unit,file=int_press_file,status='unknown',err=104)
         LEVEL2 'Reading internal pressure from:'
         LEVEL3 trim(int_press_file)
      case default
   end select
   call get_int_pressure(int_press_method,int_press_unit,julday,secs,nlev,z)

!  The light extinction profiles
   select case (extinct_method)
      case (0)
         open(extinct_unit,file=extinct_file,status='unknown',err=105)
         LEVEL2 'Reading extinction data from:'
         LEVEL3 trim(extinct_file)
         call read_extinction(extinct_unit,julday,secs)
      case (1)
         A=0.58;g1=0.35;g2=23.0
      case (2)
         A=0.68;g1=1.20;g2=28.0
      case (3)
         A=0.62;g1=0.60;g2=20.0
      case (4)
         A=0.67;g1=1.00;g2=17.0
      case (5)
         A=0.77;g1=1.50;g2=14.0
      case (6)
         A=0.78;g1=1.40;g2=7.9
      case (7)
         A=0.7;g1=0.40;g2=8.0 ! Adolf Stips - Lago Maggiore
      case default
   end select


!  The vertical advection velocity
   select case (w_adv_method)
      case (FROMFILE)
         open(w_adv_unit,file=w_adv_file,status='unknown',err=106)
         LEVEL2 'Reading vertical velocity observations from:'
         LEVEL3 trim(w_adv_file)
      case default
   end select
   call get_w_adv(w_adv_method,w_adv_unit,julday,secs)

!  The sea surface elevation
   select case (zeta_method)
      case (FROMFILE)
         open(zeta_unit,file=zeta_file,status='unknown',err=107)
         LEVEL2 'Reading sea surface elevations from:'
         LEVEL3 trim(zeta_file)
      case default
   end select
   call get_zeta(zeta_method,zeta_unit,julday,secs)

!  The observed velocity profile
   select case (vel_prof_method)
      case (0)
         uprof = 0.
         vprof = 0.
      case (2)
         open(vel_prof_unit,file=vel_prof_file,status='UNKNOWN',err=108)
         LEVEL2 'Reading velocity profiles from:'
         LEVEL3 trim(vel_prof_file)
         call get_vel_profile(vel_prof_unit,julday,secs,nlev,z)
      case default
   end select

!  The observed dissipation profile
   select case (e_prof_method)
      case (0)
         epsprof = 0.
      case (2)
         open(e_prof_unit,file=e_prof_file,status='UNKNOWN',err=109)
         LEVEL2 'Reading dissipation profiles from:'
         LEVEL3 trim(e_prof_file)
         call get_eps_profile(e_prof_unit,julday,secs,nlev,z)
      case default
   end select

   return

80 FATAL 'Unable to open "',trim(fn),'" for reading'
   stop 'init_observations'
81 FATAL 'I could not read "sprofile" namelist'
   stop 'init_observations'
82 FATAL 'I could not read "tprofile" namelist'
   stop 'init_observations'
83 FATAL 'I could not read "ext_pressure" namelist'
   stop 'init_observations'
84 FATAL 'I could not read "int_pressure" namelist'
   stop 'init_observations'
85 FATAL 'I could not read "extinct" namelist'
   stop 'init_observations'
86 FATAL 'I could not read "w_advspec" namelist'
   stop 'init_observations'
87 FATAL 'I could not read "zetaspec" namelist'
   stop 'init_observations'
88 FATAL 'I could not read "velprofile" namelist'
   stop 'init_observations'
89 FATAL 'I could not read "eprofile" namelist'
   stop 'init_observations'
90 FATAL 'I could not read "bprofile" namelist'
   stop 'init_observations'

101 FATAL 'Unable to open "',trim(s_prof_file),'" for reading'
   stop 'init_observations'
102 FATAL 'Unable to open "',trim(t_prof_file),'" for reading'
   stop 'init_observations'
103 FATAL 'Unable to open "',trim(ext_press_file),'" for reading'
   stop 'init_observations'
104 FATAL 'Unable to open "',trim(int_press_file),'" for reading'
   stop 'init_observations'
105 FATAL 'Unable to open "',trim(extinct_file),'" for reading'
   stop 'init_observations'
106 FATAL 'Unable to open "',trim(w_adv_file),'" for reading'
   stop 'init_observations'
107 FATAL 'Unable to open "',trim(zeta_file),'" for reading'
   stop 'init_observations'
108 FATAL 'Unable to open "',trim(vel_prof_file),'" for reading'
   stop 'init_observations'
109 FATAL 'Unable to open "',trim(e_prof_file),'" for reading'
   stop 'init_observations'

   return
   end subroutine init_observations
!EOC

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: get_all_obs
!
! !INTERFACE:
   subroutine get_all_obs(julday,secs,nlev,z)
!
! !DESCRIPTION:
!  During the time integration this subroutine is called each time step
!  to update all 'observation'. The routine is basically a wrapper
!  routine which calls the variable specific routines.
!  The only input to this routine is the time (in internal GOTM
!  representation) and the vertical grid. It is up to each of the individual
!  routines to use this information and to provide updated 'observations'.
!
! !USES:
   IMPLICIT NONE
!
! !INPUT PARAMETERS:
   integer, intent(in)                 :: julday,secs
   integer, intent(in)                 :: nlev
   REALTYPE, intent(in)                :: z(:)
!
! !REVISION HISTORY:
!  Original author(s): Karsten Bolding & Hans Burchard
!
!  See observation module
!
!EOP
!-----------------------------------------------------------------------
!BOC
   if(s_prof_method .eq. 2) then
      call get_s_profile(s_prof_unit,julday,secs,nlev,z)
   end if

   if(t_prof_method .eq. 2) then
      call get_t_profile(t_prof_unit,julday,secs,nlev,z)
   end if

   call get_ext_pressure(ext_press_method,ext_press_unit,julday,secs)

   call get_int_pressure(int_press_method,int_press_unit,julday,secs,nlev,z)

   if(extinct_method .eq. 0) then
      call read_extinction(extinct_unit,julday,secs)
   end if

   call get_w_adv(w_adv_method,w_adv_unit,julday,secs)

   call get_zeta(zeta_method,zeta_unit,julday,secs)

   if(vel_prof_method .eq. 2) then
      call get_vel_profile(vel_prof_unit,julday,secs,nlev,z)
   end if

   if(e_prof_method .eq. 2) then
      call get_eps_profile(e_prof_unit,julday,secs,nlev,z)
   end if

   return
   end subroutine get_all_obs
!EOC

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: read_obs
!
! !INTERFACE:
   subroutine read_obs(unit,yy,mm,dd,hh,min,ss,N,obs,ierr)
!
! !DESCRIPTION:
!  This routine will read all non-profile observations.
!  The routine allows for reading more than one scalar variable at a time.
!  The number of data to be read is specified by {\tt N}.
!  Data read-in are returned
!  in the 'obs' array. It is up to the calling routine to assign
!  meaning full variables to the individual elements in {\tt obs}.
!
! !USES:
   IMPLICIT NONE
!
! !INPUT PARAMETERS:
   integer, intent(in)                 :: unit
   integer, intent(in)                 :: N
!
! !OUTPUT PARAMETERS:
   integer, intent(out)                :: yy,mm,dd,hh,min,ss
   REALTYPE,intent(out)                :: obs(:)
   integer, intent(out)                :: ierr
!
! !REVISION HISTORY:
!  Original author(s): Karsten Bolding & Hans Burchard
!
!  See observation module
!
!EOP
!
! !LOCAL VARIABLES:
   character                 :: c1,c2,c3,c4
   integer                   :: i
!
!-----------------------------------------------------------------------
!BOC
   ierr=0
   read(unit,'(A72)',ERR=100,END=110) cbuf
   read(cbuf,900,ERR=100,END=110) yy,c1,mm,c2,dd,hh,c3,min,c4,ss
   read(cbuf(20:),*,ERR=100,END=110) (obs(i),i=1,N)

   return
100 ierr=READ_ERROR
   return
110 ierr=END_OF_FILE
   return
900 format(i4,a1,i2,a1,i2,1x,i2,a1,i2,a1,i2)
   end subroutine read_obs
!EOC

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: read_profiles
!
! !INTERFACE:
   subroutine read_profiles(unit,nlev,cols,yy,mm,dd,hh,min,ss,z, &
                            profiles,lines,ierr)
!
! !DESCRIPTION:
!  Similar to {\tt read\_obs()} but used for reading profiles instead of
!  scalar data.
!  The data will be interpolated on the grid specified by nlev and z.
!  The data can be read 'from the top' or 'from the bottom' depending on
!  a flag in the actual file.
!
! !USES:
   IMPLICIT NONE
!
! !INPUT PARAMETERS:
   integer, intent(in)                 :: unit
   integer, intent(in)                 :: nlev,cols
   REALTYPE, intent(in)                :: z(:)
!
! !INPUT/OUTPUT PARAMETERS:
   integer, intent(inout)              :: lines
!
! !OUTPUT PARAMETERS:
   integer, intent(out)                :: yy,mm,dd,hh,min,ss
   REALTYPE, intent(out)               :: profiles(:,:)
   integer, intent(out)                :: ierr
!
! !REVISION HISTORY:
!  Original author(s): Karsten Bolding & Hans Burchard
!
!  See observation module
!
!EOP
!
! !LOCAL VARIABLES:
   character                 :: c1,c2,c3,c4
   integer                   :: i,j,rc
   integer                   :: N,up_down
   REALTYPE,allocatable,dimension(:)   :: tmp_depth
   REALTYPE,allocatable,dimension(:,:) :: tmp_profs
!
!-----------------------------------------------------------------------
!BOC
   ierr=0
   read(unit,'(A72)',ERR=100,END=110) cbuf
   read(cbuf,900,ERR=100,END=110) yy,c1,mm,c2,dd,hh,c3,min,c4,ss
   read(cbuf(20:),*,ERR=100,END=110) N,up_down

   lines = lines+1

   allocate(tmp_depth(0:N),stat=rc)
   if (rc /= 0) stop 'read_profiles: Error allocating memory (tmp_depth)'
   allocate(tmp_profs(0:N,cols),stat=rc)
   if (rc /= 0) stop 'read_profiles: Error allocating memory (tmp_profs)'

   if(up_down .eq. 1) then
      do i=1,N
         lines = lines+1
         read(unit,*,ERR=100,END=110) tmp_depth(i),(tmp_profs(i,j),j=1,cols)
      end do
   else
      do i=N,1,-1
         lines = lines+1
         read(unit,*,ERR=100,END=110) tmp_depth(i),(tmp_profs(i,j),j=1,cols)
      end do
   end if

   call gridinterpol(N,cols,tmp_depth,tmp_profs,nlev,z,profiles)

   deallocate(tmp_depth)
   deallocate(tmp_profs)

   return
100 ierr=READ_ERROR
   return
110 ierr=END_OF_FILE
   return
900 format(i4,a1,i2,a1,i2,1x,i2,a1,i2,a1,i2)

   end subroutine read_profiles
!EOC

!-----------------------------------------------------------------------

   end module observations

!-----------------------------------------------------------------------
! Copyright by the GOTM-team under the GNU Public License - www.gnu.org
!-----------------------------------------------------------------------
