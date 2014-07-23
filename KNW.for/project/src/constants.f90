MODULE constants
  
  USE types
  
  IMPLICIT NONE

  !FOR TALENT BENCHMARK
  REAL(KIND=r_kind), PARAMETER :: hbarc = 197.32891 !hbar*c in MeVfm
  REAL(KIND=r_kind), PARAMETER :: mnc2 = 938.9059 !Neutron mass in MeV

  REAL(KIND=r_kind), PARAMETER  :: pi = 3.1415926535897932384626433832795_R_KIND
  !REAL(KIND=r_kind), PARAMETER :: hbarc = 197.3269718 !hbar*c in MeVfm
  !REAL(KIND=r_kind), PARAMETER :: mpc2 =  938.272046 !Proton mass in MeV
  !REAL(KIND=r_kind), PARAMETER :: mnc2 = 939.565379 !Neutron mass in MeV
  !REAL(KIND=r_kind), PARAMETER :: alphaFS = 7.2973525698E-3 !Fine structure constant, elementary charge in Gaussian units  e^2 = hbarc*alphaFS 
  

! for testing 
 !  REAL(KIND=r_kind), PARAMETER :: hbarc = 1.0 !hbar*c in MeVfm
 ! REAL(KIND=r_kind), PARAMETER :: mpc2 =  938.272046 !Proton mass in MeV
 !  REAL(KIND=r_kind), PARAMETER :: mnc2 = 1.0 !Neutron mass in MeV
  
  


END MODULE CONSTANTS
