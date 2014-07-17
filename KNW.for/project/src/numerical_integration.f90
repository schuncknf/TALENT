MODULE numerical_integration

USE types

LOGICAL, protected :: is_init_grid_GH = .false.
INTEGER, protected :: grid_size_GH
REAL(KIND=r_kind), ALLOCATABLE, protected :: grid_points_GH(:)
REAL(KIND=r_kind), ALLOCATABLE, protected :: grid_weights_GH(:), grid_weights_GL(:)

INTEGER, parameter :: factorial_table_max = 300
REAL(KIND=r_kind), protected :: factorial_table(0:factorial_table_max)


CONTAINS

! initializes grid points and wheigts for Gauss-Hermite integration, n is the number of grid points between r=0,infty. The order of the Hermite polynomial for the quadrature is 2*n
SUBROUTINE init_grid_GH(n)

  IMPLICIT NONE

  INTEGER, intent(in) :: n

  REAL(kind=r_kind), parameter :: one_r = 1.0_r_kind
  INTEGER :: ndfact, ndherm
  REAL(kind=r_kind) :: eps
  REAL(kind=r_kind), allocatable :: gpwrk(:),gwwrk(:)

  INTEGER :: II

  grid_size_GH = n

  IF(allocated(grid_points_GH)) DEALLOCATE(grid_points_GH)
  IF(allocated(grid_weights_GH)) DEALLOCATE(grid_weights_GH)
  ALLOCATE(grid_points_GH(grid_size_GH),grid_weights_GH(grid_size_GH))
  
  

  CALL init_factorial_table

  eps = 1e-14_r_kind
  ndfact = factorial_table_max
  ndherm = ndfact
  
  IF(2*grid_size_GH-1 > ndherm) THEN
       WRITE(*,*) 'ERROR in subrotuine numerical_integration_ra, in module numerical_integration:'
       WRITE(*,*) 'grid_size_gh > ',(ndherm+1)/2,' is not supported, stopping.'
       STOP
    END IF


  IF(allocated(gpwrk)) DEALLOCATE(gpwrk)
  IF(allocated(gwwrk)) DEALLOCATE(gwwrk)
  ALLOCATE(gpwrk(ndherm),gwwrk(ndherm))
  gpwrk = 0.0_r_kind
  gwwrk = 0.0_r_kind

  CALL HERMIT(2*grid_size_GH, gpwrk, gwwrk, eps, factorial_table, ndherm, ndfact)
  
  grid_points_GH = gpwrk(1:grid_size_GH)
  grid_weights_GH = gwwrk(1:grid_size_GH)  

  DEALLOCATE(gpwrk,gwwrk)

!!$  !for test, order 2,4 OK
!!$  WRITE(*,*) 'order GH = ', 2*grid_size_GH
!!$  WRITE(*,*) 'grid pts    grid weights'
!!$  DO II=1,grid_size_GH
!!$     WRITE(*,'(2f16.10)') grid_points_GH(II), grid_weights_GH(II)
!!$  END DO
!!$  !end for test

  is_init_grid_GH = .true.

END SUBROUTINE init_grid_GH

SUBROUTINE init_factorial_table
  IMPLICIT NONE
  
  INTEGER :: II
  
  factorial_table(0) = 1.0_r_kind

  DO II=1,factorial_table_max
     factorial_table(II) = II*factorial_table(II-1)
  END DO


END SUBROUTINE init_factorial_table



  SUBROUTINE HERMIT(NN,X,A,EPS,FACT,NDHERM,NDFACT)
    USE types
    IMPLICIT REAL(KIND=r_kind) (A-H,O-Z) !<- Added to ensure that everything is double
    REAL(KIND=r_kind) PI,XTT,DDPN,PPN1
    INTEGER FN,NN,NDHERM,NDFACT
    REAL(KIND=r_kind) X,A,FACT,EPS
    DIMENSION X(NDHERM),A(NDHERM)
    DIMENSION FACT(0:NDFACT)
    !
    !=======================================================================
    !
    !         HERMIT CALCULATES ZEROES X(I) OF THE N-TH ORDER HERMITE
    !         POLYNOMIAL.   THE  LARGEST  ZERO  IS  STORED  IN  X(1).
    !         IT CALCULATES  ALSO THE CORRESPONDING COEFFICIENTS A(I)
    !         OF THE  NN-TH ORDER  GAUSS-HERMITE  QUADRATURE  FORMULA
    !         OF THE DEGREE 2*NN-1 .
    !
    !         (TRIVIAL MODIFICATIONS OF THE ROUTINE FROM STROUD BOOK)
    !
    !=======================================================================
    !
    PI=4.0_r_kind*ATAN(1.0_r_kind)
    !
    IF (NDHERM.GT.NDFACT) then
       write(*,*) 'module grid: ND/HERMI :',' NDHERM = ',NDHERM,' NDFACT = ',NDFACT
       STOP
    end IF
    IF (NN-1.GT.NDFACT) then 
       write(*,*) 'module grid: ND/HERMI :',' NN-1 = ',NN-1,' NDFACT = ',NDFACT
       STOP
    end IF
    
 !   IF (NDHERM.GT.NDFACT) stop 'module grid: ND/HERMI'
 !   IF (NN-1.GT.NDFACT) stop 'module grid: ND/HERMI'
    !
    EEPS=EPS
    !
    FN=NN
    N1=NN-1
    N2=(NN+1)/2
    !
    CC=SQRT(PI)*FACT(FN-1)/(2.0_r_kind**N1)
    S=(2.0_r_kind*FN+1)**.16667_r_kind
    !
    DO I=1,N2
       !
       IF (I-1) 10,1,2
       !
       !=======================================================================
       !            T H E   L A R G E S T    Z E R O
       !=======================================================================
       !
1      XT=S**3-1.85575_r_kind/S
       !
       GO TO 9
       !
2      IF (I-2) 10,3,4
       !
       !=======================================================================
       !            T H E     S E C O N D    Z E R O
       !=======================================================================
       !
3      XT=XT-1.14_r_kind*FN**.426_r_kind/XT
       !
       GO TO 9
       !
4      IF (I-3) 10,5,6
       !
       !=======================================================================
       !            T H E     T H I R D     Z E R O
       !=======================================================================
       !
5      XT=1.86_r_kind*XT-.86_r_kind*X(1)
       !
       GO TO 9
       !
6      IF (I-4) 10,7,8
       !
       !=======================================================================
       !            T H E     F O U R T H   Z E R O
       !=======================================================================
       !
7      XT=1.91_r_kind*XT-.91_r_kind*X(2)
       !
       GO TO 9
       !
       !=======================================================================
       !            A L L     O T H E R   Z E R O S
       !=======================================================================
       !
8      XT=2.0_r_kind*XT-X(I-2)
       !
9      XTT=XT
       !
       !=======================================================================
       CALL H_ROOT(XTT,NN,DDPN,PPN1,EEPS)
       !=======================================================================
       !
       DPN=DDPN
       PN1=PPN1
       XT =XTT
       !
       X(I)=XT
       A(I)=CC/DPN/PN1
       NI=NN-I+1
       X(NI)=-XT
       A(NI)=A(I)
       !
    END DO
    !
10  CONTINUE
    !
    !=======================================================================
    !
    RETURN
  END SUBROUTINE HERMIT

  SUBROUTINE H_ROOT(X,NN,DPN,PN1,EPS)
    USE types
    IMPLICIT REAL(KIND=r_kind) (A-H,O-Z) ! (A-D,F-H,O-Z)
    REAL(KIND=r_kind) X,DPN,PN1,EPS
    INTEGER NN,ITERMX
    DATA ITERMX / 25 /
    !
    !=======================================================================
    !        IMPROVES THE APPROXIMATED ROOT X. IN ADDITION WE ALSO OBTAIN
    !                                          PN1 = VALUE OF H(N-1) AT X
    !        (FROM THE MONOGRAPH BY STROUD)
    !=======================================================================
    !
    ITER=0
1   ITER=ITER+1
    !
    CALL HRECUR(P,DP,PN1,X,NN)
    !
    D=P/DP
    X=X-D
    !
    IF (ABS(D)-EPS) 3,3,2
    !
2   IF (ITER-ITERMX) 1,4,4
    !
3   DPN=DP
    !
    RETURN
    !
4   CONTINUE
    !
    DPN=DP
    !
    PRINT 5,ITER,X,DP
5   FORMAT(///,1X,'ITERATIONS IN H_ROOT, ITER = ',I2,'  X =',F20.10,&
         '  DP = ',F20.10,///)
    !
    !=======================================================================
    !
    RETURN
  END SUBROUTINE H_ROOT

    SUBROUTINE HRECUR(PN,DPN,PN1,X,NN)
    USE types
    IMPLICIT REAL(KIND=r_kind) (A-H,O-Z)
    INTEGER NN
    REAL(KIND=r_kind) PN,DPN,PN1,X
    !
    !=======================================================================
    !         AUXILIARY RECURRENCE ROUTINE (FROM THE MONOGRAPH BY STROUD)
    !=======================================================================
    !
    P1=1.0_r_kind
    P =X
    !
    DP1=0.0_r_kind
    DP =1.0_r_kind
    !
    DO J=2,NN
       !
       FJ =J
       FJ2=(FJ-1.0_r_kind)/2.0_r_kind
       !
       Q =X*P-FJ2*P1
       DQ=X*DP+P-FJ2*DP1
       !
       P1=P
       P =Q
       !
       DP1=DP
       DP =DQ
       !
    END DO
    !
    PN=P
    DPN=DP
    PN1=P1
    !
    !=======================================================================
    !
    RETURN
  END SUBROUTINE HRECUR




END MODULE NUMERICAL_INTEGRATION
