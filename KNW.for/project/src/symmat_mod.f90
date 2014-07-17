MODULE symmat_mod

CONTAINS

  SUBROUTINE calculate_eigs_symmetric_matrix(size,matrix,eigs)

    USE types
    IMPLICIT NONE

    INTEGER, intent(in) :: size
    REAL(kind=r_kind), intent(in) :: matrix(size,size)
    REAL(kind=r_kind), intent(out) :: eigs(size)
    REAL(kind=r_kind) :: copy_matrix(size,size)
    INTEGER :: lwork, info, II, JJ
    REAL(kind=r_kind) :: work1(1)
    REAL(kind=r_kind), ALLOCATABLE :: work2(:)



    !eigenvalues of real symmetrix matrix, used to check that input args are of correct type
    INTERFACE 
       SUBROUTINE dsyev(jobz, uplo, n, a, lda, w, work, lwork, info)
         CHARACTER(len=1) ::  jobz, uplo
         INTEGER          ::  info, lda, lwork, n
         REAL(KIND=8)     ::  w(*)
         REAL(KIND=8)   ::  a(lda, *), work(*)
       END SUBROUTINE dsyev
    END INTERFACE

    !copy_matrix = matrix
    !copy_matrix(1:size,1:size) = matrix(1:size,1:size)

    DO II = 1,size
       DO JJ = 1,size
          copy_matrix(II,JJ) = matrix(II,JJ)
       END DO
    END DO

!!$  WRITE(*,*) 'copy of matrix'
!!$  CALL print_matrix(size,copy_matrix)
!!$
!!$  WRITE(*,*) 'the matrix'
!!$  CALL print_matrix(size,matrix)

    lwork = -1
    CALL dsyev('N','U',size,copy_matrix,size,eigs,work1,lwork,info)

    lwork = work1(1)
    ALLOCATE(work2(lwork))
    CALL dsyev('N','U',size,copy_matrix,size,eigs,work2,lwork,info)

    !WRITE(*,*) '..copy of matrix'
    !CALL print_matrix(size,copy_matrix)

    !WRITE(*,*) '..the matrix'
    !CALL print_matrix(size,matrix)

!!$  WRITE(*,*) 'copy of matrix'
!!$  CALL print_matrix(size,copy_matrix)

    DEALLOCATE(work2)


  END SUBROUTINE calculate_eigs_symmetric_matrix



  SUBROUTINE generate_symmetrix_matrix(size,matrix)

    USE types

    IMPLICIT NONE

    INTEGER, intent(in) :: size
    REAL(kind=r_kind), intent(out) :: matrix(size,size)
    REAL(kind=r_kind) :: val

    INTEGER :: II, JJ


    DO II = 1,size
       DO JJ =II,size

          val = II*JJ
          matrix(II,JJ) = val
          matrix(JJ,II) = val
       END DO
    END DO



  END SUBROUTINE generate_symmetrix_matrix


  SUBROUTINE print_matrix(size,matrix)

    USE types

    IMPLICIT NONE

    INTEGER, intent(in) :: size
    REAL(kind=r_kind) :: matrix(size,size)

    INTEGER :: II, JJ

    DO II = 1,size
       DO JJ = 1,size
          WRITE(*,'(E18.8)',advance = "no") matrix(II,JJ)
       END DO
       WRITE(*,*)
    END DO



  END SUBROUTINE PRINT_MATRIX


END MODULE symmat_mod
