!====================================================================
!  Routines to calculate various types of z-coord. in SCHISM (single or
!  double precision)
!  zcor_SZ_[single,double]
!  get_vgrid_[single,double]

!  ifort -Bstatic -O3 -c compute_zcor.f90
!  pgf90 -O2 -mcmodel=medium  -Bstatic -c compute_zcor.f90
!====================================================================
!====================================================================
  module compute_zcor

!#ifdef USE_DOUBLE
!    integer,parameter,private :: RKIND=8
!#else
!    integer,parameter,private :: RKIND=4
!#endif

    public 

    contains

      subroutine zcor_SZ_double(dp,eta,h0,h_s,h_c,theta_b,theta_f,kz,nvrt,ztot,sigma,zcor,idry,kbp)
!     Routine to compute z coordinates for SCHISM's SZ vertical system (double
!     precision)
!     Inputs:
!             dp: depth;
!             eta: elevation;
!             h0: min. depth;
!             h_s: transition depth between S and Z layers;
!             h_c: transition depth between S and sigma
!             theta_b, theta_f: spacing const. in S coordinate system;
!             nvrt: total # of vertical levels (S+Z);
!             kz: # of Z levels (1 if pure S);
!             ztot(1:kz):  z coordinates for Z levels; note that ztot(kz)=-h_s;
!             sigma(1:nvrt-kz+1): sigma coordinates for S (or sigma) levels; note that sigma(1)=-1, sigma(nvrt-kz+1)=0;
!     Outputs:
!             idry: wet (0) or dry (1) flag;
!             kbp: bottom index (0 if dry);
!             zcor(kbp:nvrt): z coordinates (junk if dry);    
      implicit real(8)(a-h,o-z)

      integer, intent(in) :: kz,nvrt
      real(8), intent(in) :: dp,eta,h0,h_s,h_c,theta_b,theta_f,ztot(nvrt),sigma(nvrt)
      integer, intent(out) :: idry,kbp
      real(8), intent(out) :: zcor(nvrt)

      !Local
      real(8) :: s_con1,cs(nvrt)

      include 'zcor_SZ.txt'

      end subroutine zcor_SZ_double

      subroutine zcor_SZ_single(dp,eta,h0,h_s,h_c,theta_b,theta_f,kz,nvrt,ztot,sigma,zcor,idry,kbp)
!     Routine to compute z coordinates for SCHISM's SZ vertical system (single
!     precision)
!     Inputs:
!             dp: depth;
!             eta: elevation;
!             h0: min. depth;
!             h_s: transition depth between S and Z layers;
!             h_c: transition depth between S and sigma
!             theta_b, theta_f: spacing const. in S coordinate system;
!             nvrt: total # of vertical levels (S+Z);
!             kz: # of Z levels (1 if pure S);
!             ztot(1:kz):  z coordinates for Z levels; note that ztot(kz)=-h_s;
!             sigma(1:nvrt-kz+1): sigma coordinates for S (or sigma) levels; note that sigma(1)=-1, sigma(nvrt-kz+1)=0;
!     Outputs:
!             idry: wet (0) or dry (1) flag;
!             kbp: bottom index (0 if dry);
!             zcor(kbp:nvrt): z coordinates (junk if dry);    
      implicit real(4)(a-h,o-z)

      integer, intent(in) :: kz,nvrt
      real(4), intent(in) :: dp,eta,h0,h_s,h_c,theta_b,theta_f,ztot(nvrt),sigma(nvrt)
      integer, intent(out) :: idry,kbp
      real(4), intent(out) :: zcor(nvrt)

      !Local
      real(4) :: s_con1,cs(nvrt)

      include 'zcor_SZ.txt'

      end subroutine zcor_SZ_single

!====================================================================
      subroutine get_vgrid_double(vgrid,np,nvrt1,ivcor,kz,h_s,h_c,theta_b,theta_f,ztot,sigma,sigma_lcl,kbp)
!     Routine to read in vgrid.in (either in the current dir or ../) (double
!     precision)
!     Not all outputs have meaningful values, depending on ivcor.
      implicit real(8)(a-h,o-z)

      character(len=*), intent(in) :: vgrid
      integer, intent(in) :: np,nvrt1
      integer, intent(out) :: ivcor,kz,kbp(np)
      real(8), intent(out) :: h_s,h_c,theta_b,theta_f,ztot(nvrt1),sigma(nvrt1),sigma_lcl(nvrt1,np)

      !local
      logical :: lexist
      
      include 'get_vgrid.txt'

      end subroutine get_vgrid_double

      subroutine get_vgrid_single(vgrid,np,nvrt1,ivcor,kz,h_s,h_c,theta_b,theta_f,ztot,sigma,sigma_lcl,kbp)
!     Routine to read in vgrid.in (either in the current dir or ../) (single
!     precision)
!     Not all outputs have meaningful values, depending on ivcor.
      implicit real(4)(a-h,o-z)

      character(len=*), intent(in) :: vgrid
      integer, intent(in) :: np,nvrt1
      integer, intent(out) :: ivcor,kz,kbp(np)
      real(4), intent(out) :: h_s,h_c,theta_b,theta_f,ztot(nvrt1),sigma(nvrt1),sigma_lcl(nvrt1,np)

      !local
      logical :: lexist
      
      include 'get_vgrid.txt'

      end subroutine get_vgrid_single

!====================================================================
  end module compute_zcor
