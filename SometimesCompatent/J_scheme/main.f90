program main
  use admin_jscheme
  use HF_mod
  implicit none
  
  type(spd) :: jbas
  type(full_sp_block_mat) :: F,T,Vgam,D,Dx,rho
  type(twobody_array) :: interaction 
  integer :: i,A,q,r
  real(8) :: crit
  
  A = 8
  
  ! set some stuff up brah 
  call allocate_sp_descript(jbas,A) 
  call allocate_sp_mat(jbas,F) 
  call duplicate_sp_mat(F,T)
  call duplicate_sp_mat(F,Vgam) 
  call duplicate_sp_mat(F,D) 
  call duplicate_sp_mat(F,Dx) 
  call duplicate_sp_mat(F,rho)
  call allocate_tp_array(jbas,interaction)   
  call read_tp_matrix_elements(jbas,interaction)

  call H0_matrix(T,jbas) 
  
  print*, v_elem( 0, 0, 1 , 6 , 1, 6 , jbas, interaction) 
  print*, v_elem( 2, 0, 1 , 6 , 1, 6 , jbas, interaction) 

  !initial eigenvectors
  do q = 1, D%blocks
     D%blkM(q)%matrix = 0.d0 
     do i =1,D%map(q)
        D%blkM(q)%matrix(i,i) = 1.d0
     end do 
  end do 
  
  crit = 10.d0 
  r = 0
  do while (crit > 1e-6) 
     
     call density_matrix(rho,D,DX,jbas) 
     call gamma_matrix(Vgam,interaction,rho,jbas)
  
     ! fock matrix
     do q = 1,T%blocks
        F%blkM(q)%matrix = T%blkM(q)%matrix + Vgam%blkM(Q)%matrix
     end do
     
     print*,
     do q = 1,T%blocks
         do i = 1, T%map(q) 
           print*, F%blkM(q)%matrix(i,:) 
        end do 
        print*
     end do 
     
     call diagonalize_blocks(F) 
     
     ! new eigenvectors
     ! calculate conv. criteria
     ! store eigenvalues
 
     crit = 0.d0 
     do q = 1,T%blocks
        D%blkM(q)%matrix = F%blkM(q)%matrix
        crit = crit + sqrt(sum((D%blkM(q)%eigval-F%blkM(q)%eigval)**2))
        D%blkM(q)%eigval = F%blkM(q)%eigval
       
     end do
     print*, crit, F%blkM(1)%eigval(1)
    !r = r + 1 
   
 crit = 1e-12
 end do 
 
 
 print*
 print*, '====================='
 print*, 'Hartree Fock Energy: '
 print*, e_HF(rho,Vgam,F,jbas)
       
end program
