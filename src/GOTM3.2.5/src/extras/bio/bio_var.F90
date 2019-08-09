!$Id: bio_var.F90,v 1.5 2004-07-30 09:22:20 hb Exp $
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
!  Remember this Hans
!
! !USES:
!  default: all is public.
   public
!
! !PUBLIC DATA MEMBERS:
   integer                               :: bio_model
   integer                               :: numc,numcc
   REALTYPE                              :: I_0
   REALTYPE, dimension(:), allocatable   :: zlev
   REALTYPE, dimension(:), allocatable   :: par
   REALTYPE, dimension(:,:), allocatable :: cc,ws
   integer                               :: surface_flux_method=-1
   integer                               :: n_surface_fluxes=-1
   REALTYPE, dimension(:), allocatable   :: sfl_read
   REALTYPE, dimension(:), allocatable   :: sfl,bfl
   logical, dimension(:), allocatable    :: mussels_inhale
   logical, dimension(:,:), allocatable  :: particle_active
   integer, dimension(:,:), allocatable  :: particle_indx
   REALTYPE, dimension(:,:), allocatable :: particle_pos

   integer, dimension(:), allocatable    :: var_ids
   character(len=64), dimension(:), allocatable :: var_names
   character(len=64), dimension(:), allocatable :: var_units
   character(len=64), dimension(:), allocatable :: var_long

   REALTYPE, parameter                   :: secs_pr_day=86400.

!
! !REVISION HISTORY:
!  Original author(s): Hans Burchard & Karsten Bolding
!
!  $Log: bio_var.F90,v $
!  Revision 1.5  2004-07-30 09:22:20  hb
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
