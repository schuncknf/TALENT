module ME_minn
  
contains
!========================================================
subroutine calculate_interaction(hbaromega,nmax,fname) 
  ! writes interaction elements to file fname 
  ! hw = oscl 
  ! nmax is max basis state
  implicit none
  
  integer,parameter :: ngauss = 64
  real(8),parameter :: hc = 197.32891, mnuc = 938.9059
  integer :: nmax, qmax,i,j,k,l,n
  integer :: s1,s2,s3,s4
  real(8),dimension(ngauss) ::  wrr,rr
  real(8) :: omega,oscl,hbaromega
  character(*) :: fname 
    
  !omega = hbaromega  ! divide by hbar
  oscl =sqrt(hc*hc/mnuc/hbaromega) ! inverse length parameter 
   
  ! get quadrature mesh and weights 
  call gauss_legendre(0.d0, 30.d0/oscl, rr, wrr, ngauss ) 
 
  open(unit= 37,file = fname//'_twobody.dat') 

  qmax = (nmax + 1)*2 - 1   ! total number of states minus 1 
  
  print*, hbaromega, qmax
  do i = 0,qmax
     do j = i+1,qmax
      
        do k = 0,qmax
           do l = k+1,qmax 
        
              write(37,*) Minn_matrix_element(i,j,k,l,oscl,rr,wrr,ngauss)
   
  end do; end do; end do; end do

  
  close(37)
  
  open(unit = 37, file = fname//'_onebody.dat') 
  do i = 0,qmax
     do j = i , qmax
        
        if (i == j) then 
           n = i/2 
           write(37,*)  (2*n + 1.5d0) * hbaromega
        else 
           write(37,*) 0.d0
        end if
        
     end do 
  end do 
  close(37) 

end subroutine
!===============================================
real(8) function Minn_matrix_element(q1,q2,q3,q4,oscl,rr,wrr,ng) 
  ! q index is just single label of state. 
  ! odd refers to spin down, even to spin up
  implicit none 

  real(8), parameter :: VR = 200.00d0 , VS = 91.85d0
  real(8), parameter :: kapR = 1.487d0,kapS = 0.465d0 
  real(8) :: oscl
  integer :: q1,q2,q3,q4,ng
  integer :: n1,n2,n3,n4,s1,s2,s3,s4
  real(8),dimension(ng) :: wrr,rr
 
  
  ! spins
  s1 = mod(q1,2) 
  s2 = mod(q2,2) 
  s3 = mod(q3,2) 
  s4 = mod(q4,2) 
  
  ! HO number
  n1 = q1/2 
  n2 = q2/2
  n3 = q3/2
  n4 = q4/2
  
  if (kd(s1,s3)*kd(s2,s4) == 1) then
     if (kd(s1,s4)*kd(s2,s3) == 1 )  then 
        Minn_matrix_element = 0.d0 
     else 
  Minn_matrix_element = 0.5d0*(VR*gaussian_integral(n1,n2,n3,n4,oscl,rr,wrr,ng,kapR) - &
       VS * gaussian_integral(n1,n2,n3,n4,oscl,rr,wrr,ng,kapS) + &
       VR * gaussian_integral(n1,n2,n4,n3,oscl,rr,wrr,ng,kapR) - &
       VS * gaussian_integral(n1,n2,n4,n3,oscl,rr,wrr,ng,kapS) )
     end if 
     
  else if (kd(s1,s4)*kd(s2,s3) == 1) then
     Minn_matrix_element = -0.5d0*(VR*gaussian_integral(n1,n2,n3,n4,oscl,rr,wrr,ng,kapR) - &
       VS * gaussian_integral(n1,n2,n3,n4,oscl,rr,wrr,ng,kapS) + &
       VR * gaussian_integral(n1,n2,n4,n3,oscl,rr,wrr,ng,kapR) - &
       VS * gaussian_integral(n1,n2,n4,n3,oscl,rr,wrr,ng,kapS) )
  else 
     Minn_matrix_element = 0.d0 
  end if 
  
end function 
!========================================================================
real(8) function overlap_integral(n1,n2,n3,n4,oscl,rr,wrr,ng)

  ! overlap integral for l = 0 states only
  implicit none
  real(8),parameter :: pi = 3.14159
  integer:: n1, n2, n3, n4, i, j, ng
  real(8) :: oscl_r, int_sum, xr, xp, z, factor1, factor2,oscl
  real(8),dimension(ng) :: rr, wrr

           int_sum = 0.d0
           do i=1,ng 
              do j = 1, ng

         int_sum=int_sum+wrr(i)*wrr(j)*rr(i)*rr(i)*rr(j)*rr(j)* &
              Rnl(n1,0,oscl,rr(i))*Rnl(n2,0,oscl,rr(j)) * &
              Rnl(n3,0,oscl,rr(i))*Rnl(n4,0,oscl,rr(j))
              end do 
           enddo
       
          overlap_integral = int_sum
   
end function
!==============================================
real(8) function gaussian_integral(n1,n2,n3,n4,oscl,rr,wrr,ng,mu)

  ! exp(-r_ij**2/mu) for only l = 0 states 
  implicit none
  integer:: n1, n2, n3, n4, i, j, ng
  real(8) :: oscl_r, int_sum, xr, xp, z, factor1, factor2,oscl
  real(8),dimension(ng) :: rr, wrr
  real(8) :: mu

            int_sum = 0.d0
!$OMP PARALLEL DO PRIVATE(i,j) SHARED(rr,wrr,oscl,n1,n2,n3,n4) REDUCTION(+:int_sum)             
           do i=1,ng 
              do j = 1,ng

         int_sum=int_sum+wrr(i)*wrr(j)*rr(i)*rr(j)* &
              Rnl(n1,0,oscl,rr(i))*Rnl(n2,0,oscl,rr(j)) * &
              Rnl(n3,0,oscl,rr(i))*Rnl(n4,0,oscl,rr(j)) * &
              (exp( -(rr(i)*rr(i) + rr(j)*rr(j) - 2*rr(i)*rr(j))*mu ) - &
              exp( -(rr(i)*rr(i) + rr(j)*rr(j) + 2*rr(i)*rr(j))*mu ) )
         
        
              end do 
           end do
!$OMP END PARALLEL DO          
          gaussian_integral = int_sum * 0.25d0 / mu
   

end function
!==============================================
real(8) function Rnl(n,l,oscl,r) 
  ! verified. This works DFWI 
  implicit none 
  
  real(8),parameter :: pi = 3.14159265359d0
  real(8) :: factor2,oscl,r,xi,cx(0:800)
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
  subroutine laguerre_general( n, alpha, x, cx )
    ! calculate laguerre polynomials
    implicit none
    integer, INTENT(IN)  :: n
    real (8) ::  alpha
    real (8) :: cx(0:n)
    integer:: i
    real (8), INTENT(IN) ::  x
   
    if ( alpha <= -1.0D+00 ) THEN
       write ( *, '(a)' ) ' '
       write ( *, '(a)' ) 'LAGUERRE_GENERAL - Fatal error!'
       write ( *, '(a,g14.6)' ) '  The input value of ALPHA is ', alpha
       write ( *, '(a)' ) '  but ALPHA must be greater than -1.'
       STOP
    end if
    if ( n < 0 ) THEN
       RETURN
    end if
    cx(0) = 1.0D+00
    if ( n == 0 ) THEN
       RETURN
    end if
    cx(1) = 1.0D+00 + alpha - x
    do i = 2, n
       cx(i) = ( ( real ( 2 * i - 1, kind = 8 ) + alpha - x ) * cx(i-1)   &
            + ( real (   - i + 1, kind = 8 ) - alpha     ) * cx(i-2) ) &
            / real (     i,     kind = 8 )
    end do
   
  end subroutine laguerre_general
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
  
end function
!============================
real(8) function fac(r) 
  ! log of the factorial
  implicit none 
  
  integer :: r,i
  real(8) :: tot 
  
  
  tot = 0.0
  
  do i = 1,r
     tot = tot + log(float(i))
  end do 
  
  fac = tot
 
end function
!=============================
real(8) function Vpot(r)  
  ! any general potential may be inserted here
  implicit none 
  
  real(8) :: r
  
  Vpot =  -1.d0/r 

end function 
!========================================================
integer function kd(a,b) 
  implicit none
  
  integer :: a,b
  
  if (a == b ) then 
     kd = 1
  else 
     kd = 0
  end if 
end function 
!========================================================
end module

