program HF

implicit none

real (8) :: omega, sum_v, sum_d, tol, sum_E, sum_HF, sum_HF_2
integer :: n_max, k, mu, nu, i, j, l, info, N, A 
real (8), allocatable, dimension (:) :: E, E_prev, work
real (8), allocatable, dimension(:,:) :: h, t, delta, gamma, rho, rho_prev, D, test
real (8), allocatable, dimension(:,:,:,:) :: v  
character (7) :: fname

fname='hw10_n9'
n_max=20
A=2

allocate (D(n_max,n_max),rho(n_max,n_max),rho_prev(n_max,n_max),test(n_max,n_max))
allocate (h(n_max,n_max),t(n_max,n_max),v(n_max,n_max,n_max,n_max))
allocate (gamma(n_max,n_max), E(n_max), work(10*n_max), E_prev(n_max))

!Some read comands for the Nathan matrix elements
!==============================

open(unit = 39, file = fname//'_onebody.dat') 

do i=1,n_max
   do j=i,n_max
      
      read(39,*) t(i,j)
      t(j,i) = t(i,j) 
      
   end do
end do      

open(unit = 37, file = fname//'_twobody.dat') 

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


!======================
!Toy matrix elements

!~ do i=1,n_max

!~    if (i.lt.(n_max/2+1)) then

!~        t(i,i)=1.0d0
!~    else
	
!~ 	   t(i,i)=2.0d0
	   
!~    end if
    
!~ end do

!~ do i=1,n_max
!~   do k=1,n_max
    
!~     if (i.lt.(n_max/2+1) .and. k.lt.(n_max/2+1)) then
    
!~      v(i,k,i,k)=-2.5d0
!~      v(k,i,i,k)=-v(i,k,i,k)
!~      v(i,k,k,i)=-v(i,k,i,k)
!~      
!~     else
    
!~      v(i,k,i,k)=0.5d0
!~      v(k,i,i,k)=-v(i,k,i,k)
!~      v(i,k,k,i)=-v(i,k,i,k)
     
!~     end if
    
!~   end do 
!~ end do      


!===========================
!Initial condition

do i=1,A

	rho(i,i)=1
    rho_prev(i,i)=1
    
end do

do mu=1,n_max
   do nu=1,n_max
   
   sum_v=0.0d0
   
      do k=1,A

        sum_v= sum_v + v(mu,k,nu,k)

      end do
 
    gamma(mu,nu)=sum_v   
       
   end do
end do 

h=t+gamma

!====================================
!HF loops

tol=1.0d-8
sum_E=1.0d0
N=0

print*, 'Energy condition < E-8'

do while (sum_E.gt.tol)
   
   N=N+1
        
   call dsyev('V','U',n_max,h,n_max,E,work,10*n_max,info)    
     
   D=h
 
   do mu=1,n_max
      do nu=1,n_max   
      
        sum_d=0.0d0
      
        do k=1,A
        
!~          sum_d = sum_d + D(mu,k)*conjg(D(nu,k))
          sum_d = sum_d + D(mu,k)*D(nu,k)
          
        end do
        
      rho_prev(mu,nu) = sum_d    
            
     end do 
   end do
   
!~    rho = rho_prev
   rho=0.5d0*(rho+rho_prev)
     
  do mu=1,n_max    
     do nu=1,n_max
      
      sum_v=0.0d0
      
      do k=1,n_max 
        do l=1,n_max
        
          sum_v= sum_v + v(mu,k,nu,l)*rho(l,k)

        end do
      end do
 
      gamma(mu,nu) = sum_v  
    
     end do
  end do  
  
  h=t+gamma
  
!~   test=maxval(rho*rho-rho)
  
  sum_E=0.0d0
  
!~   sum_E=sum(test)
   
  do k=1,n_max

   sum_E = sum_E + abs(E(k)-E_prev(k))/n_max
  
  end do   

  E_prev=E

!~   print*, E,E_prev

 print*, sum_E
 

end do      

sum_HF=0.0d0

do i=1,A

	sum_HF= sum_HF + E(i)

end do

do i=1,A
   do k=1,n_max
     do l=1,n_max
     
         sum_HF= sum_HF + t(k,l)*D(l,i)*D(k,i)

     end do
   end do
end do    

sum_HF=0.5d0*sum_HF

sum_HF_2=0.0d0
 
 do k = 1,A

   sum_HF_2 = sum_HF_2 + E(k) 
   
end do 

do k = 1, n_max
   do l = 1,n_max 
   
      sum_HF_2 = sum_HF_2 - 0.5*rho(k,l)*gamma(l,k) 
   
   end do 
end do 
   
print*, 'steps'
print*, N
    
print*, 'E_HF'
print*, sum_HF, sum_HF_2


!print*, test

end program HF
