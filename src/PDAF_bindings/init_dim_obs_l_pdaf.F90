!$Id: init_dim_obs_l_pdaf.F90 82 2019-02-26 15:12:22Z lnerger $
!BOP
!
! !ROUTINE: init_dim_obs_l_pdaf --- Set dimension of local observation vector
!
! !INTERFACE:
SUBROUTINE init_dim_obs_l_pdaf(domain_p, step, dim_obs_f, dim_obs_l)

! !DESCRIPTION:
! User-supplied routine for PDAF.
! Used in the filters: LSEIK/LETKF/LESTKF
!
! The routine is called during the loop over
! all local analysis domains. It has to set 
! the dimension of the local observation vector 
! for the current local analysis domain.
!
! The routine is called by each filter process.
!
! !REVISION HISTORY:
! 2005-09 - Lars Nerger - Initial code
! Later revisions - see svn log
!
! !USES:
!   USE mod_assimilation, &
!        ONLY: local_range, coords_obs_f, obs_index_l, distance_l
  IMPLICIT NONE

! !ARGUMENTS:
  INTEGER, INTENT(in)  :: domain_p   ! Current local analysis domain
  INTEGER, INTENT(in)  :: step       ! Current time step
  INTEGER, INTENT(in)  :: dim_obs_f  ! Full dimension of observation vector
  INTEGER, INTENT(out) :: dim_obs_l  ! Local dimension of observation vector

! !CALLING SEQUENCE:
! Called by: PDAF_lseik_update   (as U_init_dim_obs_l)
! Called by: PDAF_lestkf_update  (as U_init_dim_obs_l)
! Called by: PDAF_letkf_update   (as U_init_dim_obs_l)
! Called by: PDAF_lnetf_update   (as U_init_dim_l)
!EOP


! **********************************************
! *** Initialize local observation dimension ***
! **********************************************

  ! Template reminder - delete when implementing functionality
  WRITE (*,*) 'TEMPLATE init_dim_obs_l_pdaf.F90: Set local observation dimension here!'

!   dim_obs_l = ??

! Also initialize array of distances here
!   distance_l = ?

END SUBROUTINE init_dim_obs_l_pdaf

