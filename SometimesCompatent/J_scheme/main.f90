program main
  use admin_jscheme
  implicit none
  
  type(spd) :: jbas
  type(full_sp_block_mat) :: F,T,Vgam,D
  type(twobody_array) :: interaction 
  integer :: i,A
  
  A = 2
  call allocate_sp_descript(jbas,A) 
  call allocate_sp_mat(jbas,F) 
  call duplicate_sp_mat(F,T)
  call duplicate_sp_mat(F,Vgam) 
  call duplicate_sp_mat(F,D) 
  
  call allocate_tp_array(jbas,interaction) 
  

  
end program
