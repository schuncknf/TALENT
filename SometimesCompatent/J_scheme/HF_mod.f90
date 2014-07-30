module HF_mod
  use admin_jscheme
  ! J-scheme Hartree-Fockery. 
  implicit none 
  
  real(8),public,parameter :: al =1.d0 , bet = 0.d0
contains 
!==================================================
!==================================================
subroutine H0_matrix(T,jbas) 
  implicit none 
  
  type(full_sp_block_mat) :: T
  type(spd) :: jbas
  integer :: q,i,r
  
  ! loop over blocks
  do q = 1, T%blocks
     T%blkM(q)%matrix = 0.d0 
     
     r = 1
     ! loop over all states
     do i = 1, jbas%total_orbits 
        
        ! check if in block
        if ( ( jbas%ll(i) .ne. T%blkM(q)%lmda(1) ) .or.&
              (jbas%jj(i) .ne. T%blkM(q)%lmda(2) ) ) cycle
        
        T%blkM(q)%matrix(r,r) = jbas%e(i) 
        if (r == T%map(q)) exit
        r = r + 1
        
     end do
     
  end do 
! i'm just shy. 
end subroutine
!====================================================
!==================================================== 
subroutine gamma_matrix(gam,int,rho,jbas) 
  implicit none
  
  type(full_sp_block_mat) :: gam,rho
  type(twobody_array) :: int
  type(spd) :: jbas
  integer :: q,r,i,jmax,lmax,n1,n2,j,JJ,n3,n4
  integer :: grho,hrho,qrho,jrho,lrho,PAR,jfoc,lfoc
  real(8) :: sm
  
  jmax = jbas%Jtotal_max
  lmax = jbas%lmax
  
  do q = 1, gam%blocks
     gam%blkM(q)%matrix = 0.d0 
     jfoc = rho%blkM(q)%lmda(2)
     lfoc = rho%blkM(Q)%lmda(1)

     ! loop over states in this block
     do i = 1, gam%map(q)
        do j = i,gam%map(q) 
        
           n1 = rho%blkM(q)%states(i) 
           n3 = rho%blkM(q)%states(j) 
           
           ! sum over blocks of the density matrix
           sm = 0.d0 
           do qrho =  1, rho%blocks
           
              jrho = rho%blkM(qrho)%lmda(2) 
              lrho = rho%blkM(qrho)%lmda(1) 
           
              PAR =mod(lfoc + lrho,2) 
                        
              ! sum over elements of the block
              do grho = 1,rho%map(qrho) 
                 do hrho = 1,rho%map(qrho) 
                       
                    n2 = rho%blkM(qrho)%states(grho) 
                    n4 = rho%blkM(qrho)%states(hrho) 
                       
                    ! sum over allowed JJ values
                    
                    do JJ = abs(jrho - jfoc),jrho+jfoc,2
                       
                       sm = sm + rho%blkM(qrho)%matrix(hrho,grho) &
                            * v_elem(JJ,PAR,n1,n2,n3,n4,jbas,int) &
                            * (JJ + 1.d0)/(jrho + 1.d0) * &
                            sqrt( 1.d0 + kd(n1,n2)*(-1)**(JJ/2)) * &
                            sqrt( 1.d0 + kd(n3,n4)*(-1)**(JJ/2))
                      print*, n1,n2,n3,n4,JJ,sm
                    end do 
                    
                 end do 
              end do 
           end do 
   
           gam%blkM(q)%matrix(i,j) = sm/(jfoc + 1.d0)
           gam%blkM(q)%matrix(j,i) = sm/(jfoc + 1.d0)
              
           
        end do 
     end do 

  end do 
end subroutine
!===========================================================
!===========================================================
subroutine density_matrix(rho,D,DX,jbas) 
  implicit none 
  
  type(full_sp_block_mat) :: rho,D,DX 
  type(spd) :: jbas
  integer :: i,j,q,N
  
  
  do q = 1, rho%blocks
     N = rho%map(q)
     if (N == 0) cycle 
    
    do i = 1, N
       DX%blkM(q)%matrix(:,i) = D%blkM(q)%matrix(:,i) *&
            sqrt(float(jbas%con(rho%blkM(q)%states(i))))
    end do 
    
   
    call dgemm('N','T',N,N,N,al,Dx%blkM(q)%matrix,N,&
         Dx%blkM(q)%matrix,N,bet,rho%blkM(q)%matrix,N)
  end do 

end subroutine
!===========================================================
!===========================================================
real(8) function e_HF(rho,gam,F,jbas) 
  implicit none 
  
  real(8) :: sm
  integer :: q,i,j
  type(spd) :: jbas
  type(full_sp_block_mat) :: rho,F,gam 
  
 
  sm = 0.d0
  do q = 1, F%blocks
     do i = 1, F%map(q) 
        sm = sm + F%blkM(q)%eigval(i) &
             * jbas%con(F%blkm(q)%states(i)) 
     end do 
     
     do i = 1, F%map(q)
        do j = 1,F%map(q) 
     
           sm = sm - 0.5d0 * rho%blkM(q)%matrix(i,j) * &
                gam%blkM(q)%matrix(j,i) 
        end do 
     end do 
 end do 
 
 e_HF = sm
     
end function 
!===========================================================
!===========================================================              
end module         
                
              
                 
                 
           
            
