MODULE hf_solver_mod
  USE types
  USE Ho_basis
  USE constants

CONTAINS

  SUBROUTINE hf_solver(num_part,hbaromega,b,Nmax,lmax,tol,iteration_max,mixing,is_RPA,HF_v_2body_scale)

  IMPLICIT NONE

  REAL(kind=r_kind) :: b, hbaromega
  INTEGER :: nmax,lmax
  REAL(kind=r_kind) :: delta, tol
  INTEGER :: iteration, iteration_max
  INTEGER :: num_part
  REAL(kind=r_kind) :: e_fermi, mixing, E_HF, E_HF_v2
  LOGICAL :: is_RPA

  REAL(kind=r_kind), intent(in) :: HF_v_2body_scale

  
  !b = 1.0_r_kind !maby have default value based on mass
  !hbaromega = 10.0_r_kind
  !iteration_max = 20
  !tol = 1e-15_r_kind
  !num_part = 8
  !lmax = 12
  !Nmax = 12
  !Nmax = 2
  !lmax = 2

  !mixing = 0.0_r_kind


  IF(lmax < 0) THEN
     lmax = Nmax
  END IF

 
  IF(b <=  0.0_r_kind) THEN

     CALL init_Ho_basis(Nmax,lmax,b,hbaromega) !b is not used if hbaromega is supplied

  ELSE
     CALL init_Ho_basis(nmax,lmax,b)

  END IF

  WRITE(*,*) 'hbarc = ',hbarc
  WRITE(*,*) 'neuron mass = ',mnc2
  
  CALL hf_init(is_RPA,HF_v_2body_scale) !generates matrix elements needed etc

  CALL hf_init_dens(num_part) !generates initial density matrix
  
  
  WRITE(*,*) 'Mixing parameter = ',mixing
  WRITE(*,*) 'Starting iteration'
  WRITE(*,'(A5,2A20,A14,A14)') 'iter','E_HF','E_HF_v2','e_fermi','delta'
  delta = 100.0_r_kind
  iteration = 1
  hf_iteration_loop : DO WHILE (iteration <= iteration_max .and. delta > tol)
     
     CALL hf_update_hamiltonian
       
     CALL hf_diagonalize

     CALL hf_find_fermi_energy(num_part,e_fermi)
     
     CALL hf_update_density_matrix(mixing)
     
     CALL hf_total_energy(E_HF)
     CALL hf_total_energy_v2(E_HF_v2)

     CALL hf_calculate_delta(delta)     

     WRITE(*,'(I5,2F20.14,F14.6,E14.6)') iteration, E_HF,E_HF_v2, e_fermi, delta
    
     iteration = iteration + 1
     

  END DO hf_iteration_loop


  CALL print_hf_states

END SUBROUTINE hf_solver



  !CALL hf_calculate_tot_energy

  !CALL hf_print_energies


END MODULE hf_solver_mod
