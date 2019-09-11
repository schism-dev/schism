!   Copyright 2014 College of William and Mary
!
!   Licensed under the Apache License, Version 2.0 (the "License");
!   you may not use this file except in compliance with the License.
!   You may obtain a copy of the License at
!
!     http://www.apache.org/licenses/LICENSE-2.0
!
!   Unless required by applicable law or agreed to in writing, software
!   distributed under the License is distributed on an "AS IS" BASIS,
!   WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or
!   implied.
!   See the License for the specific language governing permissions and
!   limitations under the License.

      MODULE fib_param

!
!!= Marta Rodrigues=====================================================
!!= National Laboratory for Civil Enginnering                          !
!!======================================================================

!!======================================================================
!! April, 2010 - Original code
!! June, 2010 - Sink term introduced                                   ! 
!!======================================================Marta Rodrigues=
!!                                                                     !
!                                                                      !
!  Fecal Indicator Bacteria (FIB) model parameters                     !
!                                                                      !
!  flag_fib       Flag that the defines how the decay rate is          !
!                 computed:                                            !
!                   1 - user-defined constant (in time) decay rate     !
!                      (can be horizontally varying)                   !
!                   2 - Canteras et al., 1995                          !
!                   3 - Servais et al., 2007
!  kk_fib         Decay rate associated with mortality (/s)            !
!  sink_fib       Particle settling velocity (m/s)                     !
!  fraction_fib   Fraction of microrganisms in the water column        !
!                   associated with suspended sediments                !        
!                                                                      !
!=======================================================================
!
        
      IMPLICIT NONE
     
!        integer :: flag_fib
        integer, parameter :: r8 = 8
        real(r8), allocatable  :: kk_fib(:,:)
        real(r8),allocatable :: sink_fib(:), fraction_fib(:)
      
      END MODULE fib_param
