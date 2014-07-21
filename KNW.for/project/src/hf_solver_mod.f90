MODULE hf_solver_mod
  USE types
  USE Ho_basis
  USE constants

CONTAINS

  SUBROUTINE hf_solver

  IMPLICIT NONE

  REAL(kind=r_kind) :: b
  INTEGER :: nmax,lmax
  REAL(kind=r_kind) :: delta, tol
  INTEGER :: iteration, iteration_max
  INTEGER :: num_part
  REAL(kind=r_kind) :: e_fermi

  b = 1.0_r_kind !maby have default value based on mass
  iteration_max = 10
  tol = 1e-5_r_kind
  num_part = 2
  lmax = 0
  nmax = 10

  

  CALL init_Ho_basis(b,nmax,lmax) 
  CALL hf_init
  CALL hf_init_dens(num_part) !generates initial density matrix
  
  
  
  WRITE(*,*) 'Starting iteration'
  WRITE(*,'(A5,A14,A14)') 'iter','e_fermi','delta'
  delta = 100.0_r_kind
  iteration = 1
  hf_iteration_loop : DO WHILE (iteration <= iteration_max .and. delta > tol)
     
     CALL hf_update_hamiltonian
     
     CALL hf_diagonalize

     CALL hf_find_fermi_energy(num_part,e_fermi)
     
     CALL hf_update_density_matrix(0.0_r_kind)
     
     !CALL hf_calculate_delta(delta)     

     WRITE(*,'(I5,F14.6,F14.6)') iteration, e_fermi, delta
    
     iteration = iteration + 1
     

  END DO hf_iteration_loop


  CALL print_hf_states

END SUBROUTINE hf_solver



  !CALL hf_calculate_tot_energy

  !CALL hf_print_energies


END MODULE hf_solver_mod
