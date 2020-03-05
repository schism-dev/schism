   !!==============================================================================
  SUBROUTINE sed_comp_poro_noncoh(frac_sed, poro_gravsan)

   !--------------------------------------------------------------------------
   !                 ***  SUBROUTINE sed_comp_poro_noncoh  ***
   !
   ! ** Purpose : Compute the porosity representative of a non-cohesive matrix
   !
   ! ** Called by : sed_init and sediment
   !
   ! ** References :
   !    - Wooster et al. (2008) : Sediment supply and relative size
   !                              distribution effects on fine sediment
   !                              infiltration into immobile gravels
   ! ** History :
   !       !  2020-02  B. Mengual
   !
   !--------------------------------------------------------------------------

    USE sed_mod
    USE schism_glbl, ONLY : rkind

    IMPLICIT NONE
    include 'mpif.h'

    !! * Arguments
    REAL(rkind), INTENT(IN)    :: frac_sed(ntr_l)
    REAL(rkind), INTENT(OUT)   :: poro_gravsan

    !! * Local declarations
    INTEGER         ::  ised
    REAL(rkind), PARAMETER :: eps = 1.0d-4
    REAL(rkind)   ::    psim,siggeo2,siggeo,stdgeo, &
                        frac_gravsan,diam_gravsan
    REAL(rkind)   ::    frac_noncoh(ntr_l), psi_sed(ntr_l)

   !!---------------------------------------------------------------------------
   !! * Executable part

    frac_gravsan=0.0d0
    DO ised=1,ntr_l
      IF (iSedtype(ised) > 0) THEN
        frac_gravsan=frac_gravsan+frac_sed(ised)
      END IF
    END DO

    IF (frac_gravsan .GT. eps) THEN

      frac_noncoh(:)=0.0d0
      DO ised=1,ntr_l
        IF (iSedtype(ised) > 0) THEN
          frac_noncoh(ised)=frac_sed(ised)/frac_gravsan
        END IF
      END DO

      psim=0.0d0
      psi_sed(:)=0.0d0
      DO ised=1,ntr_l
        IF (iSedtype(ised) > 0) THEN
          psi_sed(ised)=LOG(Sd50(ised)*1000.0d0)/LOG(2.0d0)
          psim=psim+(psi_sed(ised)*frac_noncoh(ised))
        END IF
      END DO

      siggeo2=0.0d0
      DO ised=1,ntr_l
        IF (iSedtype(ised) > 0) THEN
          siggeo2=siggeo2+( frac_noncoh(ised)*(psi_sed(ised)-psim)**2.0d0 )
        END IF
      END DO

      siggeo  = siggeo2**(0.5d0)
      stdgeo  = MAX(1.0d0,MIN(2.0d0**(siggeo),4.0d0)) ! to remain in the range of values described by Wooster et al. (2008)

      diam_gravsan = (2.0d0**(psim))/1000.0d0 ! psi_sed is computed with diam_sed in mm --> /1e3
      poro_gravsan = Awooster*(stdgeo**(Bwooster))

    ELSE

      poro_gravsan = Awooster

    END IF

  END SUBROUTINE sed_comp_poro_noncoh
   !!==============================================================================

