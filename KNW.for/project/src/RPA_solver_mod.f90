MODULE RPA_solver_mod
  USE types
  USE Ho_basis

CONTAINS
  
  SUBROUTINE RPA_solver
    IMPLICIT NONE

    REAL(kind=r_kind) :: Ecut

    ! cutoff above fermi surface in units of MeV
    Ecut = -1.0_r_kind
    
    CALL RPA_init_ph_space(Ecut)

    CALL RPA_init_matrix

    !CALL RPA_diagonalize

    
    

  END SUBROUTINE RPA_solver


END MODULE RPA_solver_mod

  
