program HF

implicit none

real (8) :: omega, sum_v, sum_d, tol, sum_E,sm
integer :: n_max, n1, n2, k, mu, nu, i, j, l, p ,info,Abody
integer , allocatable, dimension(:) :: contract
real (8), allocatable, dimension (:) :: E, E_prev,work
real (8), allocatable, dimension(:,:) :: h, t, delta, gamma
real (8), allocatable, dimension(:,:) :: D, rho, rho_prev 
real (8), allocatable, dimension(:,:,:,:) :: v  


n_max=22
Abody = 2

allocate (D(n_max,n_max),rho(n_max,n_max),rho_prev(n_max,n_max))
allocate (h(n_max,n_max),t(n_max,n_max),v(n_max,n_max,n_max,n_max))
allocate (gamma(n_max,n_max), E(n_max),work(10*n_max),E_prev(n_max))

!Some read comands for the Nathan matrix elements
!==============================

do i=1,n_max
   do j=i,n_max
      
      read(39,*) t(i,j)
      t(j,i) = t(i,j) 
      
   end do
end do      

V = 0.d0
do i=1,n_max
   do j=i+1,n_max
     do k=i,n_max
       do l=k+1,n_max
       
        read(37,*) v(i,j,k,l)
        v(k,l,i,j) = v(i,j,k,l)
        v(l,k,i,j) = -v(i,j,k,l) 
        v(k,l,j,i) = -v(i,j,k,l) 
        v(l,k,j,i) = v(i,j,k,l) 
        v(j,i,k,l) = -v(i,j,k,l) 
        v(i,j,l,k) = -v(i,j,k,l) 
        v(j,i,l,k) = v(i,j,k,l) 
        
       end do
     end do 
   end do
end do  
print*, 'read input',v(1,2,1,2),v(12,13,12,13)
!===========================

!Initial condition

sum_v=0.0d0
tol=1.0d-8

D(:,:)=0.0d0
rho(:,:)=0.0d0
rho_prev(:,:)=0.0d0

do n1=1,n_max

    D(n1,n1)=1.0d0
    
end do
  
do n1=1,Abody

    rho(n1,n1)=1.d0
    rho_prev(n1,n1)=1.d0
   
end do

do n1=1,n_max
   do n2=1,n_max
   
      sum_v = 0.d0
      do k=1,Abody

        sum_v= sum_v + v(n1,k,n2,k)

      end do
 
    gamma(n1,n2)=sum_v   
    !h(n1,n2)=t(n1,n2)+ gamma(n1,n2)
   
   end do
end do 
h = t + gamma
!====================================

sum_E = 10.
do while (sum_E .gt. tol)    

   call dsyev('V','U',n_max,h,n_max,E,work,10*n_max,info)    
    
   D=h
 
   do mu=1,n_max
      do nu=1,n_max   
      
        sum_d=0.0d0
      
        do k=1,Abody
        
!         sum_d = sum_d + D(mu,k)*conjg(D(nu,k))
          sum_d = sum_d +D(mu,k)*D(nu,k)
          
        end do
        
      rho_prev(mu,nu) = sum_d    
            
     end do 
   end do
   
  ! rho=0.5d0*(rho+rho_prev)
   
 rho = rho_prev   
!sm = 0.d0 
!do i = 1, 20
!   sm = sm + rho(i,i)
!end do 
!print*, sm
  do mu=1,n_max    
     do nu=1,n_max
      
      sum_v=0.0d0
      
      do k=1,n_max 
        do p=1,n_max
        
          sum_v= sum_v + v(mu,k,nu,p)*rho(p,k)
        end do
      end do
 
      gamma(mu,nu) = sum_v 
    
     end do
  end do  
  
  h = t + gamma
  
  sum_E=0.0d0
  
  do k=1,n_max
  
     sum_E = sum_E + abs(E(k)-E_prev(k))/n_max
  
  end do   

  E_prev=E
  print*, sum_e,E(1)
end do      

sum_E = 0.d0 

do k = 1, Abody
   sum_E = sum_E + E(k) 
end do 

do k = 1, n_max
   do p = 1,n_max 
      sum_E = sum_E - 0.5*rho(k,p)*gamma(p,k) 
   end do 
end do 

print*, 'E_GS'
print*, sum_E

end program HF
