!$Id: prepoststep_ens.F90 82 2019-02-26 15:12:22Z lnerger $
!BOP
!
! !ROUTINE: prepoststep_ens --- Used-defined Pre/Poststep routine for PDAF
!
! !INTERFACE:
SUBROUTINE prepoststep_ens(step, dim_p, dim_ens, dim_ens_p, dim_obs_p, &
     state_p, Uinv, ens_p, flag)

! !DESCRIPTION:
! User-supplied routine for PDAF.
! Used in the filters: SEIK/LSEIK/ETKF/LETKF/ESTKF/LESTKF
! 
! The routine is called for global filters (e.g. SEIK)
! before the analysis and after the ensemble transformation.
! For local filters (e.g. LSEIK) the routine is called
! before and after the loop over all local analysis
! domains. Also it is called once at the initial time
! before any forecasts are computed.
! The routine provides full access to the state 
! estimate and the state ensemble to the user.
! Thus, user-controlled pre- and poststep 
! operations can be performed here. For example 
! the forecast and the analysis states and ensemble
! covariance matrix can be analized, e.g. by 
! computing the estimated variances. In addition, 
! the estimates can be written to disk. If a user 
! considers to perform adjustments to the 
! estimates (e.g. for balances), this routine is 
! the right place for it.
!
! The routine is called by all filter processes.
!
! !REVISION HISTORY:
! 2004-10 - Lars Nerger - Initial code
! Later revisions - see svn log
!
! !USES:
  IMPLICIT NONE

! !ARGUMENTS:
  INTEGER, INTENT(in) :: step        ! Current time step
     ! (When the routine is called before the analysis -step is provided.)
  INTEGER, INTENT(in) :: dim_p       ! PE-local state dimension
  INTEGER, INTENT(in) :: dim_ens     ! Size of state ensemble
  INTEGER, INTENT(in) :: dim_ens_p   ! PE-local size of ensemble
  INTEGER, INTENT(in) :: dim_obs_p   ! PE-local dimension of observation vector
  REAL, INTENT(inout) :: state_p(dim_p) ! PE-local forecast/analysis state
  ! The array 'state_p' is not generally not initialized in the case of SEIK.
  ! It can be used freely here.
  REAL, INTENT(inout) :: Uinv(dim_ens-1, dim_ens-1) ! Inverse of matrix U
  REAL, INTENT(inout) :: ens_p(dim_p, dim_ens)      ! PE-local state ensemble
  INTEGER, INTENT(in) :: flag        ! PDAF status flag

! !CALLING SEQUENCE:
! Called by: PDAF_get_state       (as U_prepoststep)
! Called by: PDAF_seik_update     (as U_prepoststep)
! Called by: PDAF_enkf_analysis   (as U_prepoststep)
! Called by: PDAF_lseik_update    (as U_prepoststep)
! Called by: PDAF_estkf_analysis  (as U_prepoststep)
! Called by: PDAF_lestkf_analysis (as U_prepoststep)
! Called by: PDAF_etkf_analysis   (as U_prepoststep)
! Called by: PDAF_letkf_analysis  (as U_prepoststep)
! Called by: PDAF_lenkf_analysis  (as U_prepoststep)
! Called by: PDAF_netf_analysis   (as U_prepoststep)
! Called by: PDAF_lnetf_analysis  (as U_prepoststep)
!EOP


! ****************************
! *** Perform pre/poststep ***
! ****************************

  ! Template reminder - delete when implementing functionality
  WRITE (*,*) 'TEMPLATE prepoststep_ens.F90: Implement prepoststep here!'


END SUBROUTINE prepoststep_ens
