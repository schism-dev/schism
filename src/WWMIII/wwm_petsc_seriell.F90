#include "wwm_functions.h"
#ifdef PETSC
      MODULE PETSC_SERIELL

#include "finclude/petscsys.h"
#include "finclude/petscvec.h"
#include "finclude/petscvec.h90"
#include "finclude/petscmat.h"
#include "finclude/petscksp.h"
#include "finclude/petscao.h"
#include "finclude/petscis.h"
#include "finclude/petscdraw.h"
#include "finclude/petscviewer.h"
#include "finclude/petscviewer.h90"
#include "finclude/petscviewerdef.h"


      KSPConvergedReason reason;
      PetscInt iterationen;

      contains

      SUBROUTINE PETSC_INIT_SERIELL()

        USE DATAPOOL
        use petscpool
        IMPLICIT NONE

! create sparse matrix
!        call MatCreateSeqAIJWithArrays(PETSC_COMM_SELF, MNP, MNP, IA_P, JA_P, ASPAR, matrix, petscErr);CHKERRQ(petscErr)
        call MatCreate(PETSC_COMM_SELF, matrix, petscErr);CHKERRQ(petscErr)
        call MatSetType(matrix, MATSEQAIJ , petscErr);CHKERRQ(petscErr)
        call MatSetSizes(matrix, PETSC_DECIDE, PETSC_DECIDE, MNP, MNP, petscErr);CHKERRQ(petscErr)

!create x vector
         call VecCreate(PETSC_COMM_SELF, myX, petscErr);CHKERRQ(petscErr)
         call VecSetSizes(myX, PETSC_DECIDE, MNP, petscErr);CHKERRQ(petscErr)
         call VecSetType(myX,VECSEQ, petscErr);CHKERRQ(petscErr)

! create vec myB
         call VecCreate(PETSC_COMM_SELF, myB, petscErr);CHKERRQ(petscErr)
         call VecSetSizes(myB, PETSC_DECIDE, MNP, petscErr);CHKERRQ(petscErr)
         call VecSetType(myB,VECSEQ, petscErr);CHKERRQ(petscErr)

! create solver
         call KSPCreate(PETSC_COMM_SELF, Solver, petscErr);CHKERRQ(petscErr)
!          call KSPSetOperators(Solver, matrix, matrix, SAME_NONZERO_PATTERN, petscErr); CHKERRQ(petscErr)
         call KSPSetOperators(Solver, matrix, matrix, 0, petscErr); CHKERRQ(petscErr)

        call KSPSetType(Solver,KSPBCGS, petscErr);CHKERRQ(petscErr)
!          call KSPSetType(Solver, KSPIBCGS, petscErr);CHKERRQ(petscErr) ! this solver create a segfault

! Create preconditioner
         call KSPGetPC(Solver, Prec, petscErr);CHKERRQ(petscErr)
         call PCSetType(Prec, PCILU, petscErr);CHKERRQ(petscErr)


      END SUBROUTINE

!**********************************************************************
!*                                                                    *
!**********************************************************************
      SUBROUTINE  EIMPS_PETSC_SERIELL(ISS, IDD)

         USE DATAPOOL
         use petscpool
         IMPLICIT NONE

         INTEGER, INTENT(IN) :: ISS,IDD

         INTEGER :: I, J
!          INTEGER :: KK

         INTEGER :: IP, IPGL1, IE, POS

         INTEGER :: I1, I2, I3, IPrel

         real(rkind) :: DTK, TMP3

         real(rkind) :: LAMBDA(2)
         real(rkind) :: FL11, FL12, FL21, FL22, FL31, FL32
         real(rkind) :: CRFS(3), K1, KM(3), K(3), TRIA03

         real(rkind) :: DELTAL(3,MNE)
         real(rkind) :: KP(3,MNE), NM(MNE)
         real(rkind) :: U(MNP), C(2,MNP)

         real(rkind) :: X(MNP)
         real(rkind) :: B(MNP)

         real(rkind) ::  ASPAR(NNZ)

#ifdef TIMINGS
         REAL    :: TIME1, TIME2, TIME3, TIME4, TIME5
#endif

         INTEGER :: POS_TRICK(3,2)


         ! petsc stuff
         integer :: counter
         integer :: nnznew
         integer :: temp
         PetscInt :: ncols
         PetscInt :: eCol, eRow
         PetscScalar :: eValue

#ifdef TIMINGS
         CALL WAV_MY_WTIME(TIME1)
#endif

         POS_TRICK(1,1) = 2
         POS_TRICK(1,2) = 3
         POS_TRICK(2,1) = 3
         POS_TRICK(2,2) = 1
         POS_TRICK(3,1) = 1
         POS_TRICK(3,2) = 2

         CALL CADVXY(ISS,IDD,C)
!
!        Calculate countour integral quantities ...
!
         DO IE = 1, MNE
           I1 = INE(1,IE)
           I2 = INE(2,IE)
           I3 = INE(3,IE)
           LAMBDA(1) = ONESIXTH * (C(1,I1)+C(1,I2)+C(1,I3))
           LAMBDA(2) = ONESIXTH * (C(2,I1)+C(2,I2)+C(2,I3))
           K(1)  = LAMBDA(1) * IEN(1,IE) + LAMBDA(2) * IEN(2,IE)
           K(2)  = LAMBDA(1) * IEN(3,IE) + LAMBDA(2) * IEN(4,IE)
           K(3)  = LAMBDA(1) * IEN(5,IE) + LAMBDA(2) * IEN(6,IE)
           KP(1,IE) = MAX(ZERO,K(1))
           KP(2,IE) = MAX(ZERO,K(2))
           KP(3,IE) = MAX(ZERO,K(3))
           KM(1) = MIN(ZERO,K(1))
           KM(2) = MIN(ZERO,K(2))
           KM(3) = MIN(ZERO,K(3))
           FL11 = C(1,I2)*IEN(1,IE)+C(2,I2)*IEN(2,IE)
           FL12 = C(1,I3)*IEN(1,IE)+C(2,I3)*IEN(2,IE)
           FL21 = C(1,I3)*IEN(3,IE)+C(2,I3)*IEN(4,IE)
           FL22 = C(1,I1)*IEN(3,IE)+C(2,I1)*IEN(4,IE)
           FL31 = C(1,I1)*IEN(5,IE)+C(2,I1)*IEN(6,IE)
           FL32 = C(1,I2)*IEN(5,IE)+C(2,I2)*IEN(6,IE)
           CRFS(1) =  - ONESIXTH *  (TWO *FL31 + FL32 + FL21 + TWO * FL22 )
           CRFS(2) =  - ONESIXTH *  (TWO *FL32 + TWO * FL11 + FL12 + FL31 )
           CRFS(3) =  - ONESIXTH *  (TWO *FL12 + TWO * FL21 + FL22 + FL11 )
           DELTAL(:,IE) = CRFS(:)- KP(:,IE)
           NM(IE)       = 1.0_rkind/MIN(-THR,SUM(KM(:)))
         END DO

#ifdef TIMINGS
         CALL WAV_MY_WTIME(TIME2)
#endif

         U(:) = AC2(ISS,IDD,:)

         J     = 0    ! Counter ...
         ASPAR = ZERO ! Mass matrix ...
         B     = ZERO ! Right hand side ...
!
! ... assembling the linear equation system ....
!
         DO IP = 1, NP_RES
           IF (IOBPD(IDD,IP) .EQ. 1 .AND. IOBWB(IP) .EQ. 1 .AND. DEP(IP) .GT. DMIN) THEN
             DO I = 1, CCON(IP)
               J = J + 1
               IE    =  IE_CELL(J)
               POS   =  POS_CELL(J)
               K1    =  KP(POS,IE) ! Flux Jacobian
               TRIA03 = ONETHIRD * TRIA(IE)
               DTK   =  K1 * DT4A * IOBPD(IDD,IP)
               TMP3  =  DTK * NM(IE)
               I1    =  POSI(1,J) ! Position of the recent entry in the ASPAR matrix ... ASPAR is shown in fig. 42, p.122
               I2    =  POSI(2,J)
               I3    =  POSI(3,J)
               ASPAR(I1) =  TRIA03 + DTK - TMP3 * DELTAL(POS             ,IE) + ASPAR(I1)  ! Diagonal entry
               ASPAR(I2) =               - TMP3 * DELTAL(POS_TRICK(POS,1),IE) + ASPAR(I2)  ! off diagonal entries ...
               ASPAR(I3) =               - TMP3 * DELTAL(POS_TRICK(POS,2),IE) + ASPAR(I3)
               B(IP)     =  B(IP) + TRIA03 * U(IP)
             END DO !I: loop over connected elements ...
           ELSE
             DO I = 1, CCON(IP)
               J = J + 1
               IE    =  IE_CELL(J)
               TRIA03 = ONETHIRD * TRIA(IE)
               I1    =  POSI(1,J) ! Position of the recent entry in the ASPAR matrix ... ASPAR is shown in fig. 42, p.122
               ASPAR(I1) =  TRIA03 + ASPAR(I1)  ! Diagonal entry
               B(IP)     =  0.!B(IP)  + TRIA03 * 0.
             END DO !I: loop over connected elements ...
           END IF
         END DO !IP

         IF (LBCWA .OR. LBCSP) THEN
           DO IP = 1, IWBMNP
             IF (LINHOM) THEN
               IPrel=IP
             ELSE
               IPrel=1
             ENDIF
             IPGL1 = IWBNDLC(IP)
             ASPAR(I_DIAG(IPGL1)) = SI(IPGL1) ! Set boundary on the diagonal
             B(IPGL1)             = SI(IPGL1) *  WBAC(ISS,IDD,IPrel)
           END DO
         END IF

         IF (ICOMP .GE. 2 .AND. SMETHOD .GT. 0) THEN
           DO IP = 1, NP_RES
             IF (IOBWB(IP) .EQ. 1) THEN
               ASPAR(I_DIAG(IP)) = ASPAR(I_DIAG(IP)) + IMATDAA(ISS,IDD,IP) * DT4A * SI(IP) ! Add source term to the diagonal
               B(IP)             = B(IP) + IMATRAA(ISS,IDD,IP) * DT4A * SI(IP) ! Add source term to the right hand side
             ENDIF
           END DO
         ENDIF

#ifdef TIMINGS
         CALL WAV_MY_WTIME(TIME3)
#endif

!> \todo for efficiency reason this will be moved out of this part when the code is working, AR.
! fill the matrix
         call MatZeroEntries(matrix, petscErr);CHKERRQ(petscErr)
         counter = 1
         nnznew=0

         do i = 1, MNP
           ncols = IA_P(i+1) - IA_P(i)
           eRow = i - 1

           ! insert col by col into matrix
           do j = 1, ncols
             ! the value we want to insert
             eValue=ASPAR(counter)
             ! get the col index in old counting
             temp = JA_P(counter)
             counter = counter + 1
             nnznew=nnznew + 1

             eCol = temp
             call MatSetValue(matrix, eRow, eCol, eValue, ADD_VALUES, petscErr);CHKERRQ(petscErr)
           end do
         end do ! np
         call MatAssemblyBegin(matrix,MAT_FINAL_ASSEMBLY, petscErr);CHKERRQ(petscErr);
         call MatAssemblyEnd(matrix,MAT_FINAL_ASSEMBLY, petscErr);CHKERRQ(petscErr);


!fill RHS vector
         call VecGetArrayF90(myB, myBtemp, petscErr);CHKERRQ(petscErr)
         myBtemp = B
         call VecRestoreArrayF90(myB, myBtemp, petscErr);CHKERRQ(petscErr)
         call VecAssemblyBegin(myB, petscErr);CHKERRQ(petscErr);
         call VecAssemblyEnd(myB, petscErr);CHKERRQ(petscErr);

#ifdef TIMINGS
         CALL WAV_MY_WTIME(TIME4)
#endif
! Solve
         call KSPSolve(Solver, myB, myX, petscErr);CHKERRQ(petscErr)

#ifdef TIMINGS
         CALL WAV_MY_WTIME(TIME5)
#endif

         !WRITE(*,*) TIME2-TIME1, TIME3-TIME2, TIME4-TIME3, TIME5-TIME4

         call KSPGetConvergedReason(Solver, reason, petscErr);CHKERRQ(petscErr)
         if (reason .LT. 0) then
           write(*,*) "Failure to converge\n"
          stop 'wwm_petsc_seriell l.260'
         else
         call KSPGetIterationNumber(Solver, iterationen, petscErr);CHKERRQ(petscErr)
            if(iterationen /= 0)  write(*,*) "Number of iterations", iss,idd,iterationen
!            write(*,*) "Number of iterations", iterationen
         endif

         ! write the soluten from vec into fortran array
         X    = ZERO
         call VecGetArrayF90(myX, myXtemp, petscErr);CHKERRQ(petscErr)
         X = myXtemp
         call VecRestoreArrayF90(myX, myXtemp, petscErr);CHKERRQ(petscErr)

         IF (SUM(X) .NE. SUM(X)) CALL WWM_ABORT('NaN in X')
!          AC2(ISS, IDD,:) = MAX(ZERO,X)

         AC2(ISS,IDD,:) = X(:)

#ifdef TIMINGS
!          CALL WAV_MY_WTIME(TIME4)
#endif
        END SUBROUTINE

      !> cleanup memory. You never need to call this function by hand. It will automaticly called by PETSC_FINALIZE()
      subroutine PETSC_FINALIZE_SERIELL()
        implicit none

      end subroutine
end module
#endif
