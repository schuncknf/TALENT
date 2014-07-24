program HF

implicit none

real (8) :: omega, sum_v, sum_d, tol, sum_E, sum_HF, sum_HF_2
integer :: n_max, l_max k, mu, nu, i, j, l, info, N, A, Rmax, p,q
real (8), allocatable, dimension (:) :: E, E_prev, work !This dimension may change
real (8), allocatable, dimension(:,:,:) :: h, gamm, rho, test !I think gamm wont be necessary
real (8), allocatable, dimension(:,:) :: v, t, D
character (*) :: fname

(n_max,l_max,A,fname)

Rmax = n_max*(n_max-1)/2

allocate (D(n_max,n_max),rho(l_max,n_max,n_max),test(l_max,n_max,n_max))
allocate (h(l_max,n_max,n_max),t(n_max,n_max),v(Rmax,Rmax))
allocate (gamm(l_max,n_max,n_max), E(n_max), work(10*n_max), E_prev(n_max))

!Some read comands for the matrix elements
!==============================

open(unit = 39, file = fname//'_onebody.dat') 

do i=1,n_max
   do j=i,n_max
      
      read(39,*) t(i,j)
      t(j,i) = t(i,j) 
      
   end do
end do      

open(unit = 37, file = fname//'_twobody.dat') 

p = 1
do i = 1, n_max
   do j = i+1, n_max 
      
      q = 1
      do k = 1,n_max
         do l = k+1,n_max
      
            read(37,*) v(p,q) 
            q = q+ 1
            
         end do
      end do
      p = p+ 1
   end do 
end do 

close(37)
close(39) 

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

	rho(:,i,i)=1
    rho_prev(:,i,i)=1
    
end do

do mu=1,n_max
   do nu=1,n_max
   
   sum_v=0.0d0
   
      do k=1,A

        sum_v= sum_v + v_elem(mu,k,nu,k,v,Rmax,n_max)

      end do
 
    gamm(mu,nu)=sum_v   
       
   end do
end do 

h=t+gamm

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
   
   rho = rho_prev
!~    rho=0.5d0*(rho+rho_prev)
     
  do mu=1,n_max    
     do nu=1,n_max
      
      sum_v=0.0d0
      
      do k=1,n_max 
        do p=1,n_max
        
          sum_v= sum_v + v_elem(mu,k,nu,p,v,Rmax,n_max)*rho(p,k)

        end do
      end do
 
      gamm(mu,nu) = sum_v  
    
     end do
  end do  
  
  h=t+gamm
  
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
     do p=1,n_max
     
         sum_HF= sum_HF + t(k,p)*D(p,i)*D(k,i)

     end do
   end do
end do    

sum_HF=0.5d0*sum_HF

sum_HF_2=0.0d0
 
 do k = 1,A

   sum_HF_2 = sum_HF_2 + E(k) 
   
end do 

do k = 1, n_max
   do p = 1,n_max 
   
      sum_HF_2 = sum_HF_2 - 0.5*rho(k,p)*gamm(p,k) 
   
   end do 
end do 
   
print*, 'steps'
print*, N
    
print*, 'E_HF'
print*, sum_HF, sum_HF_2


!print*, test

end subroutine HF

!===============================================
!=========================================
real(8) function v_elem(px,qx,rx,sx,V,x,M) 
  ! grabs the TBME for pqrs
  implicit none 
  
  integer :: p,q,r,s,x,i,j,M,a,b,pre
  integer :: px,qx,rx,sx
  real(8),dimension(x,x) :: V
  real(8) :: out
  
  p = px
  q = qx
  r = rx 
  s = sx
  pre = 1
 
  ! formula only works for p < q 
  ! if this isn't true, switch them, 
  ! and note this with the prefactor
  if (p > q) then 
     a = p 
     p = q 
     q = a
     pre = -1 
  end if 
  
  if (r > s) then 
     a = r 
     r = s 
     s = a
     pre = -1 * pre 
  end if 
   
  if ((q == p ) .or. (s == r )) then 
     v_elem = 0.d0
  else 
     
  ! map from sp indeces to tp indeces
  i = (p-1) * M  - p*(p-1)/2 + q - p
  j = (r-1) * M  - r*(r-1)/2 + s - r
  
  
  v_elem = V(i,j) * pre
  end if
  
end function
!================

end program
