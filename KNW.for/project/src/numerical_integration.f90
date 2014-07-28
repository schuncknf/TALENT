MODULE numerical_integration

  USE types

  PRIVATE

  PUBLIC :: is_init_grid_GH, grid_size_GH, grid_points_GH, grid_weights_GH, init_grid_GH, is_init_grid_GLag, grid_size_GLag, alpha_GLag, grid_points_GLag, grid_weights_GLag, init_grid_GLag,is_init_grid_GL, grid_size_GL, grid_points_GL, grid_weights_GL, init_grid_GL

  LOGICAL, protected :: is_init_grid_GH = .false.
  INTEGER, protected :: grid_size_GH
  REAL(KIND=r_kind), ALLOCATABLE, protected :: grid_points_GH(:)
  REAL(KIND=r_kind), ALLOCATABLE, protected :: grid_weights_GH(:)

  INTEGER, parameter :: factorial_table_max = 300
  REAL(KIND=r_kind), protected :: factorial_table(0:factorial_table_max)

  LOGICAL, protected :: is_init_grid_GLag = .false.
  INTEGER, protected :: grid_size_GLag
  REAL(KIND=r_kind), protected :: alpha_GLag
  REAL(KIND=r_kind), ALLOCATABLE, protected :: grid_points_GLag(:)
  REAL(KIND=r_kind), ALLOCATABLE, protected :: grid_weights_GLag(:)
  
  LOGICAL, protected :: is_init_grid_GL = .false.
  INTEGER, protected :: grid_size_GL
  REAL(KIND=r_kind), ALLOCATABLE, protected :: grid_points_GL(:)
  REAL(KIND=r_kind), ALLOCATABLE, protected :: grid_weights_GL(:)
  


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


  SUBROUTINE init_grid_GLag(n,alpha)
    IMPLICIT NONE

    INTEGER, intent(in) :: n
    REAL(kind=r_kind), intent(in) :: alpha

    INTEGER :: II 

    IF(allocated(grid_points_GLag)) DEALLOCATE(grid_points_GLag)
    IF(allocated(grid_weights_GLag)) DEALLOCATE(grid_weights_GLag)


    grid_size_GLag = n
    alpha_GLag = alpha
    ALLOCATE(grid_points_GLag(n),grid_weights_GLag(n))

    write(*,*) 'init_grid_GL: calling gauss_laguerre'
    
    CALL gauss_laguerre(grid_points_GLag,grid_weights_GLag,n,alpha)

    write(*,*) 'init_grid_GL: done'

!!$    !For debug
!!$    WRITE(*,*) 'Generalized Gauss Laguerre'
!!$    WRITE(*,*) 'grid_size_GLag = ',grid_size_GLag
!!$    WRITE(*,*) 'alpha_GLag = ', alpha_GLag
!!$    WRITE(*,'(2A16)') 'grid pts','grid weights'
!!$    DO II = 1,grid_size_GLag
!!$
!!$       WRITE(*,'(2f16.10)') grid_points_GLag(II),grid_weights_GLag(II)
!!$
!!$    END DO
    !end for debug
    
    

    is_init_grid_GLag = .true.


  END SUBROUTINE init_grid_GLag
  



  !  Sets up integration points for
  !  Gaussian quadrature using Laguerre polynomials
  SUBROUTINE gauss_laguerre(x,w,n,alf)
    IMPLICIT NONE
    INTEGER :: MAXIT
    INTEGER, INTENT(IN) :: n
    REAL(KIND=8) :: EPS
    REAL(KIND=8), INTENT(IN) :: alf
    REAL(KIND=8), INTENT(INOUT) :: w(n),x(n)
    INTEGER :: i,its,j
    REAL(KIND=8) :: ai,nn !ai,gammln, nn
    REAL(KIND=8) :: p1,p2,p3,pp,z,z1

    !maxit = 10
    maxit = 50 !increased
    eps = 3.E-14
    nn = n
    DO i=1,n
       IF(i == 1)THEN
          z=(1.+alf)*(3.+.92*alf)/(1.+2.4*n+1.8*alf)
       ELSE IF(i == 2)THEN
          z=z+(15.+6.25*alf)/(1.+9.0*alf+2.5*n)
       ELSE
          ai=i-2
          z=z+((1.+2.55*ai)/(1.9*ai)+1.26*ai*alf/(1.+3.5*ai))* &
               (z-x(i-2))/(1.+.3*alf)
       ENDIF
       DO its=1,MAXIT
          p1=1.d0
          p2=0.d0
          DO  j=1,n
             p3=p2
             p2=p1
             p1=((2*j-1+alf-z)*p2-(j-1+alf)*p3)/j
          ENDDO
          pp=(n*p1-(n+alf)*p2)/z
          z1=z
          z=z1-p1/pp
          IF(ABS(z-z1) <= eps) GOTO 1
       ENDDO
       !PAUSE 'too many iterations in gaulag'
       WRITE(*,*) 'too many iterations in gaulag'
       STOP
1      x(i)=z
       w(i)=-EXP(gammln(alf+nn)-gammln(nn))/(pp*n*p2)
       !write(*,*) x(i), w(i)
    ENDDO

  END SUBROUTINE gauss_laguerre

 
 REAL(KIND =8) FUNCTION gammln(xx)
    REAL(KIND=8), INTENT(IN) :: xx
    INTEGER :: j
    REAL(KIND=8) :: ser,tmp,x,y
    REAL(KIND=8), DIMENSION(6) :: cof(6)
    cof(1)=76.18009172947146; cof(2)=-86.50532032941677
    cof(3)=24.01409824083091; cof(4)=-1.231739572450155
    cof(5)= .1208650973866179E-2; cof(6) = -.5395239384953E-5
    x=xx
    y=x
    tmp=x+5.5
    tmp=(x+0.5)*LOG(tmp)-tmp
    ser=1.000000000190015
    DO j=1,6
       y=y+1.0
       ser=ser+cof(j)/y
    ENDDO
    gammln = tmp+LOG(2.5066282746310005*ser/x)

  END FUNCTION gammln


      SUBROUTINE GAULEG(X1,X2,X,W,N)
      USE types
      IMPLICIT REAL(KIND=r_kind) (A-H,O-Z)
      INTEGER N
      REAL(KIND=r_kind) X1,X2,X(N),W(N)
      REAL(KIND=r_kind), PARAMETER :: EPS=1.E-14_r_kind      
      
      M=(N+1)/2
      XM=0.5_r_kind*(X2+X1)
      XL=0.5_r_kind*(X2-X1)
      DO I=1,M
         Z=COS(3.141592654_r_kind*(I-.25_r_kind)/(N+.5_r_kind))
 1       CONTINUE
         P1=1._r_kind
         P2=0._r_kind
         DO J=1,N
            P3=P2
            P2=P1
            P1=((2._r_kind*J-1._r_kind)*Z*P2-(J-1._r_kind)*P3)/J
         END DO
         PP=N*(Z*P1-P2)/(Z*Z-1._r_kind)
         Z1=Z
         Z=Z1-P1/PP
         IF(ABS(Z-Z1).GT.EPS)GO TO 1
         X(I)=XM-XL*Z
         X(N+1-I)=XM+XL*Z
         W(I)=2._r_kind*XL/((1._r_kind-Z*Z)*PP*PP)
         W(N+1-I)=W(I)
      END DO
      RETURN
    END SUBROUTINE GAULEG

    SUBROUTINE init_grid_GL(n)
    IMPLICIT NONE

    INTEGER, intent(in) :: n
   

    INTEGER :: II 

    IF(allocated(grid_points_GL)) DEALLOCATE(grid_points_GL)
    IF(allocated(grid_weights_GL)) DEALLOCATE(grid_weights_GL)


    grid_size_GL = n
    
    ALLOCATE(grid_points_GL(n),grid_weights_GL(n))

    write(*,*) 'init_grid_GL: calling gauleg'
    
    CALL gauleg(-1.0_r_kind,1.0_r_kind,grid_points_GL,grid_weights_GL,n)

    write(*,*) 'init_grid_GL: done'

!!$    !For debug
!!$    WRITE(*,*) 'Gauss Legendre'
!!$    WRITE(*,*) 'grid_size_GL = ',grid_size_GL
!!$    WRITE(*,'(2A16)') 'grid pts','grid weights'
!!$    DO II = 1,grid_size_GL
!!$
!!$       WRITE(*,'(2f16.10)') grid_points_GL(II),grid_weights_GL(II)
!!$
!!$    END DO
!!$    !end for debug
    
    

    is_init_grid_GL = .true.


  END SUBROUTINE init_grid_GL





END MODULE NUMERICAL_INTEGRATION
