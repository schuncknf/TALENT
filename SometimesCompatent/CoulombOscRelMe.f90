program Hydrogen_Energy
  implicit none 

  real(8) :: coulomb_integral,oscl,E(10),Rnl,overlap_integral
  real(8),allocatable,dimension(:,:) :: HAM ,T,V 
  real(8),allocatable,dimension(:) :: eig,wk
  integer :: n,l,i,j,nmax,info,koscl

  open(unit=31, file='Hyd_energy.dat') 
  l = 0 ! zero partial wave
  
  do koscl = 1,15 ! loop over omega
     oscl = 1./sqrt(koscl*0.2d0)   ! inverse length param
     !print*, koscl
     do nmax =  1,51,5 ! loop over basis size
      
        allocate(HAM(nmax,nmax)) 
        allocate(T(nmax,nmax)) 
        allocate(V(nmax,nmax)) 
        allocate(eig(nmax)) 
        allocate(wk(10*nmax)) 
   
        ! calculate V
        do i = 1, nmax
           do j = i,nmax
              ! see general_pot_integral for other potentials
         
              V(i,j) = coulomb_integral(i-1,j-1,l,oscl) 
              V(j,i) = V(i,j) 
      
           end do
        end do

        ! calculate T
        T = 0.d0
        do i = 1,nmax
           T(i,i) =  0.5*(2*i+l-0.5d0)    

           if (i > 1) then 
              T(i,i-1) =0.5*sqrt((i-1)*(i+l-0.5)) 
           end if
 
           if (i < nmax) then 
              T(i,i+1) =0.5*sqrt((i)*(i+l+0.5))
           end if
        end do
  
        !Scale T by omega
        T = T*koscl*.2d0
        HAM = T+V
  
        call dsyev( 'V','U',nmax,HAM,nmax,eig,wk,10*nmax,info)
   
        ! ground state energy stored
        E(nmax/5+1) = minval(eig) 
   

        deallocate(HAM) 
        deallocate(T) 
        deallocate(V) 
        deallocate(eig) 
        deallocate(wk) 
   
     end do
     !write to file
     write(31,'(11(f14.7))' ) , koscl*0.2d0 ,E 
  end do
  close(31)  
     
end program
!=====================================================
real(8) function overlap_integral(n1,n2,lr,oscl)

  IMPLICIT NONE
  real(8),parameter :: pi = 3.14159
  INTEGER :: n1, n2, lr, i, nrel_max
  REAL(8) :: oscl_r, int_sum, xr, xp, z, factor1, factor2,oscl
  REAL(8), ALLOCATABLE, DIMENSION(:) :: rr, wrr
  REAL(8) :: cx(0:200),fac,dfac,Rnl


  oscl_r=oscl        ! Oscillator parameter for relative

  nrel_max = 400.
  ALLOCATE ( rr(nrel_max ), wrr(nrel_max ))
  CALL gauss_legendre(0.d0, 20.d0, rr, wrr, nrel_max )
  
  int_sum=0.d0
  do i = 1, nrel_max
     int_sum = int_sum +  Rnl(n1,lr,oscl,rr(i)) * &
          Rnl(n2,lr,oscl,rr(i)) * wrr(i) *rr(i) * rr(i) 
  end do 

  overlap_integral = int_sum

END function
!===============================================
real(8) function coulomb_integral(n1,n2,lr,oscl)

  IMPLICIT NONE
  real(8),parameter :: pi = 3.14159
  INTEGER :: n1, n2, lr, i, nrel_max
  REAL(8) :: oscl_r, int_sum, xr, xp, z, factor1, factor2,oscl
  REAL(8), ALLOCATABLE, DIMENSION(:) :: rr, wrr
  REAL(8) :: cx(0:200),fac,dfac,Rnl

  nrel_max = 400
  ALLOCATE ( rr(nrel_max ), wrr(nrel_max ))
  CALL gauss_legendre(0.d0, 20.d0, rr, wrr, nrel_max )

           int_sum = 0.D0
           DO i=1,nrel_max
              
              int_sum=int_sum-wrr(i)*rr(i)*Rnl(n1,lr,oscl,rr(i))*&
                   Rnl(n2,lr,oscl,rr(i)) 
           ENDDO
           !  Coulomb energy in natural units
          coulomb_integral = int_sum
   
  DEALLOCATE ( rr, wrr)

END function
!==============================================
real(8) function general_pot_integral(n1,n2,lr,oscl)

  IMPLICIT NONE
  real(8),parameter :: pi = 3.14159
  INTEGER :: n1, n2, lr, i, nrel_max
  REAL(8) :: oscl_r, int_sum, xr, xp, z, factor1, factor2,oscl
  REAL(8), ALLOCATABLE, DIMENSION(:) :: rr, wrr
  REAL(8) :: cx(0:200),fac,dfac,Rnl,vpot

  nrel_max = 400
  ALLOCATE ( rr(nrel_max ), wrr(nrel_max ))
  CALL gauss_legendre(0.d0, 20.d0, rr, wrr, nrel_max )

           int_sum = 0.D0
           DO i=1,nrel_max
              
              int_sum=int_sum+wrr(i)*rr(i)*rr(i)*Rnl(n1,lr,oscl,rr(i))*&
                   Rnl(n2,lr,oscl,rr(i))*vpot(rr(i))
           ENDDO
           !  Coulomb energy in natural units
          general_pot_integral = int_sum
   
  DEALLOCATE ( rr, wrr)

END function
!==============================================
real(8) function Rnl(n,l,oscl,r) 
  
  implicit none 
  
  real(8),parameter :: pi = 3.14159
  real(8) :: factor2,oscl,r,xi,cx(0:200) ,dfac,fac
  integer :: n,l
  
  ! log of Anl 
  factor2 = 0.5D0*((n+l+1)*LOG(2.D0)+fac(n)-dfac(2*n+2*l+1)-0.5D0*LOG(pi)+&
       log(2.0/oscl**3))
  !Anl 
  factor2 = EXP(factor2)
  
  xi = r/oscl
  
  call laguerre_general(n,l+0.5d0,xi*xi,cx) 
  Rnl = factor2 * exp(-xi*xi*0.5d0) * cx(n) * xi**l 
  
end function
!===========================================  
  SUBROUTINE laguerre_general( n, alpha, x, cx )
    ! calculate laguerre polynomials
    IMPLICIT NONE
    INTEGER, INTENT(IN)  :: n
    REAL (8) ::  alpha
    REAL (8) :: cx(0:n)
    INTEGER :: i
    REAL (8), INTENT(IN) ::  x
    !print*, 'a'
    IF ( alpha <= -1.0D+00 ) THEN
       WRITE ( *, '(a)' ) ' '
       WRITE ( *, '(a)' ) 'LAGUERRE_GENERAL - Fatal error!'
       WRITE ( *, '(a,g14.6)' ) '  The input value of ALPHA is ', alpha
       WRITE ( *, '(a)' ) '  but ALPHA must be greater than -1.'
       STOP
    END IF
    IF ( n < 0 ) THEN
       RETURN
    END IF
    cx(0) = 1.0D+00
    IF ( n == 0 ) THEN
       RETURN
    END IF
    cx(1) = 1.0D+00 + alpha - x
    DO i = 2, n
       cx(i) = ( ( REAL ( 2 * i - 1, kind = 8 ) + alpha - x ) * cx(i-1)   &
            + ( REAL (   - i + 1, kind = 8 ) - alpha     ) * cx(i-2) ) &
            / REAL (     i,     kind = 8 )
    END DO
    !print*, 'a2'
  END SUBROUTINE laguerre_general
!========================================================
subroutine gauss_legendre(x1, x2, x, w, n)
  ! get gauleg quadrature mesh and weights
  implicit none
  integer, intent(in) :: n
  double precision, intent(in) :: x1, x2
  double precision, dimension(n), intent(out) :: x, w
  integer :: i, j, m
  double precision :: p1, p2, p3, pp, xl, xm, z, z1
  double precision, parameter :: eps=3.d-14
  
  m = (n+1)/2
  xm = 0.5d0*(x2+x1)
  xl = 0.5d0*(x2-x1)
  do i=1,m
    z = cos(3.141592654d0*(i-0.25d0)/(n+0.5d0))
    z1 = 0.0
    do while(abs(z-z1) .gt. eps)
      p1 = 1.0d0
      p2 = 0.0d0
      do j=1,n
        p3 = p2
        p2 = p1
        p1 = ((2.0d0*j-1.0d0)*z*p2-(j-1.0d0)*p3)/j
      end do
      pp = n*(z*p1-p2)/(z*z-1.0d0)
      z1 = z
      z = z1 - p1/pp
    end do
    x(i) = xm - xl*z
    x(n+1-i) = xm + xl*z
    w(i) = (2.0d0*xl)/((1.0d0-z*z)*pp*pp)
    w(n+1-i) = w(i)
  end do
  
end subroutine
!==============================
real(8) function dfac(r) 
  ! log of the double factorial
  implicit none 
  
  integer :: r,init,i
  real(8) :: tot 
  !print*, 'b'
  if ( mod(r,2) == 0 ) then 
     init = 2
  else 
     init = 1
  end if 
  
  tot = 0.0
  
   do i = init,r,2
  
     tot = tot+ log(float(i))  
   end do 
  
  dfac = tot
  !print*, 'b2'
end function
!============================
real(8) function fac(r) 
  ! log of the factorial
  implicit none 
  
  integer :: r,i
  real(8) :: tot 
  
  !print*, 'c'
  tot = 0.0
  
  do i = 1,r
     tot = tot + log(float(i))
  end do 
  
  fac = tot
 ! print*, 'c2'
end function
!=============================
real(8) function Vpot(r)  
  ! any general potential may be inserted here
  implicit none 
  
  real(8) :: r
  
  Vpot =  -1.d0/r 

end function 
