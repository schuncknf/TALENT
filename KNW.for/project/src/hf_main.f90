PROGRAM hf_main

  USE types
  USE hf_solver_mod
  USE RPA_solver_mod

  IMPLICIT NONE

  INTEGER :: num_part = 8 !number of particles
  INTEGER :: Nmax = 4 !largest major oscillaotr shell
  INTEGER :: lmax = -1 !largest angular momentum in basis, if <0 lmax = Nmax
  REAL(kind=r_kind) :: tol = 1e-15_r_kind !convergence tolerance for hf iteration
  REAL(kind=r_kind) :: mixing = 0.0_r_kind !mixing parameter for hf iteration
  REAL(kind=r_kind) :: hbaromega = 10.0_r_kind !oscillator width for basis and harmonic trap
  REAL(kind=r_kind) :: b = -1.0_r_kind !oscillator paramter
  INTEGER :: iter_max = 40

  logical :: is_RPA = .false. !to run RPA after hf calc
  REAL(kind=r_kind) :: RPA_v_res_scale = 1.0_r_kind !scaling of residual interaction
  REAL(kind=r_kind) :: HF_v_2body_scale = 1.0_r_kind !scaling of 2body matrix elements used in HF calc

  NAMELIST/hf_data/num_part,hbaromega,Nmax,lmax,tol,iter_max,mixing,is_RPA, RPA_v_res_scale, HF_v_2body_scale
  READ(5,NML=hf_data)

  CALL hf_solver(num_part,hbaromega,b,Nmax,lmax,tol,iter_max,mixing,is_RPA,HF_v_2body_scale)
  
 

  IF(is_RPA) THEN

     CALL RPA_solver(RPA_v_res_scale)

  END IF

  



END PROGRAM hf_main
