#ifdef CCSMCOUPLED
#define CESMCOUPLED
#endif
!=======================================================================

! Defines variable precision for all common data types
! Code originally based on kinds_mod.F in POP
!
! author: Elizabeth C. Hunke and William H. Lipscomb, LANL
! 2006: ECH converted to free source form (F90)
! 2020: TC added support for NO_I8 and NO_R16

      module icepack_kinds

!=======================================================================

      implicit none
      public

      integer, parameter :: char_len  = 80, &
                            char_len_long  = 256, &
                            log_kind  = kind(.true.), &
                            int_kind  = selected_int_kind(6), &
#if defined (NO_I8)
                            int8_kind = selected_int_kind(6), &
#else
                            int8_kind = selected_int_kind(13), &
#endif
                            real_kind = selected_real_kind(6), &
                            dbl_kind  = selected_real_kind(13), &
#if defined (NO_R16)
                            r16_kind  = selected_real_kind(13)
#else
                            r16_kind  = selected_real_kind(33,4931)
#endif

!=======================================================================

      end module icepack_kinds

!=======================================================================
