program HF

implicit none

real (8) :: omega, sum_v, work, info, sum_d, tol, sum_E
integer :: n_max, n1, n2, k, mu, nu, i, j, l, p 
real (8), allocatable, dimension (:) :: E, E_prev
real (8), allocatable, dimension(:,:) :: h, t, delta, gamma
real (8), allocatable, dimension(:,:) :: D, rho, rho_prev 
real (8), allocatable, dimension(:,:,:,:) :: v  


n_max=12

allocate (D(n_max,n_max),rho(n_max,n_max),rho_prev(n_max,n_max))
print*, size(rho_prev(:,1)) , size(rho_prev(1,:))
allocate (h(n_max,n_max),t(n_max,n_max),v(n_max,n_max,n_max,n_max))
allocate (gamma(n_max,n_max), E(n_max))

!Some read comands for the Nathan matrix elements
!==============================

do i=1,n_max
   do j=i,n_max
      
      read(39,*) t(i,j)
      
   end do
end do      


do i=1,n_max
   do j=i+1,n_max
     do k=1,n_max
       do l=k+1,n_max
       
        read(37,*) v(i,j,k,l)
      
       end do
     end do 
   end do
end do  



do i=1,n_max
   do j=1,i
      
      t(i,j)=t(j,i)
      
   end do
end do      



do i=1,n_max
   do j=1,i
     do k=1,n_max
       do l=1,k
       
          v(i,j,k,l)=v(j,i,l,k)
      
       end do
     end do 
   end do
end do


!===========================

!Initial condition


sum_v=0.0d0
tol=1.0d-8

D(:,:)=0.0d0
rho(:,:)=0.0d0
rho_prev(:,:)=0.0d0

do n1=1,n_max

    D(n1,n1)=1.0d0
    rho(n1,n1)=1.0d0
    rho_prev(n1,n1)=1.0d0 
   
end do
    
do n1=1,n_max
   do n2=1,n_max
   
      do k=1,n_max 

        sum_v= sum_v + v(n1,k,n2,k)

      end do
 
    gamma(n1,n2)=sum_v   
    h(n1,n2)=t(n1,n2)+ gamma(n1,n2)
   
   end do
end do 

!====================================

!HF loops

do while (sum_E.lt.tol)
        
   call dsyev('V','U',n_max,h,n_max,E,work,10*n_max,info)    
     
   D=h
 
   do mu=1,n_max
      do nu=1,n_max   
      
        sum_d=0.0d0
      
        do k=1,n_max
        
!         sum_d = sum_d + D(mu,k)*conjg(D(nu,k))
          sum_d = sum_d +D(mu,k)*D(nu,k)
          
        end do
        
      rho_prev(mu,nu) = sum_d    
            
     end do 
   end do
   
   
   rho=0.5d0*(rho+rho_prev)
      
  do mu=1,n_max    
     do nu=1,n_max
      
      sum_v=0.0d0
      
      do k=1,n_max 
        do p=1,n_max
        
          sum_v= sum_v + v(mu,k,nu,p)*rho(p,k)

        end do
      end do
 
      gamma(mu,nu) = sum_v 
      h(mu,nu) = t(mu,nu) + gamma(mu,nu)    
    
     end do
  end do  
  
  sum_E=0.0d0
  
  do k=1,n_max
  
   sum_E = sum_E + (E(k)-E_prev(k))/n_max
  
  end do   

  E_prev=E

end do      

print*, 'E_GS'
print*, E(1)

end program HF
