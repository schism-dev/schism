!$Id: bio_var.F90,v 1.11 2007-03-14 12:46:07 kbk Exp $
#include"cppdefs.h"
!-----------------------------------------------------------------------
!BOP
!
! !MODULE: bio_var --- declaration of biological variables
!
! !INTERFACE:
   module bio_var
!
! !DESCRIPTION:
!  Here all variables necessary for the biogeochemical models are
!  declared, mostly as allocatable variables.
!
! !USES:
!  default: all is public.
   public
!
! !PUBLIC DATA MEMBERS:
   integer                               :: bio_model
   integer                               :: numc,numcc
   REALTYPE, dimension(:), allocatable   :: zlev
   REALTYPE, dimension(:), allocatable   :: par
   REALTYPE, dimension(:,:), allocatable :: cc,ws
   integer                               :: surface_flux_method=-1
   integer                               :: n_surface_fluxes=-1
   REALTYPE, dimension(:), allocatable   :: sfl_read
   REALTYPE, dimension(:), allocatable   :: sfl,bfl
   integer, dimension(:), allocatable    :: posconc
   logical, dimension(:), allocatable    :: mussels_inhale
   logical, dimension(:,:), allocatable  :: particle_active
   integer, dimension(:,:), allocatable  :: particle_indx
   REALTYPE, dimension(:,:), allocatable :: particle_pos

   integer, dimension(:), allocatable    :: var_ids
   character(len=64), dimension(:), allocatable :: var_names
   character(len=64), dimension(:), allocatable :: var_units
   character(len=64), dimension(:), allocatable :: var_long

   REALTYPE, parameter                   :: secs_pr_day=86400.

!  external variables - i.e. provided by the calling program but
!  made available via this module to the different biological models
!  the variables are copied via set_env_spm() in bio.F90
   REALTYPE, dimension(:), allocatable      :: h
   REALTYPE, dimension(:), allocatable      :: t
   REALTYPE, dimension(:), allocatable      :: s
   REALTYPE, dimension(:), allocatable      :: rho
   REALTYPE, dimension(:), allocatable      :: nuh
   REALTYPE, dimension(:), allocatable      :: w
   REALTYPE, dimension(:), allocatable      :: rad
   REALTYPE                                 :: wind
   REALTYPE                                 :: I_0
   integer                                  :: w_adv_ctr=0

!  external variables updated by the biological models
!  the variables are copied back to the calling program using
!  get_bio_updates()
   REALTYPE, dimension(:), allocatable      :: bioshade_
   REALTYPE, dimension(:), allocatable      :: abioshade_

   logical                                  :: init_saved_vars=.true.
!
! !REVISION HISTORY:
!  Original author(s): Hans Burchard & Karsten Bolding
!
!  $Log: bio_var.F90,v $
!  Revision 1.11  2007-03-14 12:46:07  kbk
!  proper cleaning after simulation
!
!  Revision 1.10  2006-11-17 07:13:17  kbk
!  rho amd wind-speed available via bio_var
!
!  Revision 1.9  2006-11-12 19:42:44  hb
!  vertical advection due to physical vertical velocities enabled for the bio module
!
!  Revision 1.8  2006-10-26 13:12:46  kbk
!  updated bio models to new ode_solver
!
!  Revision 1.7  2005-12-02 20:57:27  hb
!  Documentation updated and some bugs fixed
!
!  Revision 1.6  2005-11-17 09:58:18  hb
!  explicit argument for positive definite variables in diff_center()
!
!  Revision 1.5  2004/07/30 09:22:20  hb
!  use bio_var in specific bio models - simpliefied internal interface
!
!  Revision 1.4  2004/03/30 11:32:48  kbk
!  select between eulerian or lagrangian solver
!
!  Revision 1.3  2003/10/16 15:42:16  kbk
!  simple mussesl model implemented - filter only
!
!  Revision 1.2  2003/09/16 12:11:24  hb
!  added new biological model - bio_iow
!
!  Revision 1.1  2003/07/23 12:27:31  hb
!  more generic support for different bio models
!
!
!EOP
!-----------------------------------------------------------------------

   end module bio_var

!-----------------------------------------------------------------------
! Copyright by the GOTM-team under the GNU Public License - www.gnu.org
!----------------------------------------------------------------------- 
