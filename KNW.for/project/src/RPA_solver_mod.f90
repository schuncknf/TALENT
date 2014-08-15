MODULE RPA_solver_mod
  USE types
  USE Ho_basis

CONTAINS
  
  SUBROUTINE RPA_solver
    IMPLICIT NONE

    REAL(kind=r_kind) :: Ecut, v_scale

    ! cutoff above fermi surface in units of MeV
    Ecut = -1.0_r_kind
    
    CALL RPA_init_ph_space(Ecut)

    v_scale = 1.0_r_kind

    CALL RPA_init_matrix(v_scale)

    CALL RPA_diagonalize

    CALL RPA_print_info
    
    

  END SUBROUTINE RPA_solver


END MODULE RPA_solver_mod

  
