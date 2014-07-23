program main 
  use ME_minn
  use Hartree_Fock
  implicit none 

  real(8) :: hw
  integer :: A,nmax
  character(5) :: Astr,nmaxstr,hwstr
  logical :: ex1,ex2
  character(50) :: fname
  
  call getarg(1,Astr) !number of neutrons
  call getarg(2,nmaxstr) !maximum n value in basis
  call getarg(3,hwstr) !hw  
  
  read(Astr,'(I5)') A
  read(nmaxstr,'(I5)') nmax
  read(hwstr,'(f12.7)') hw 
  
  hwstr = adjustl(hwstr) 
  nmaxstr = adjustl(nmaxstr) 
  
  fname = adjustl('hw'//trim(hwstr)//'_n'//trim(nmaxstr))
  
  inquire(file=trim(fname)//'_onebody.dat', exist=ex1) 
  inquire(file=trim(fname)//'_twobody.dat', exist=ex2) 
  
  if ( (.not. ex1) .or. (.not. ex2) ) then
      print*, 'calculating interaction files' 
      print*, fname
      call calculate_interaction(hw,nmax,trim(fname)) 
  else 
     print*, 'reading interaction from file' 
     print*, fname
  end if 
  
  call HF(2*(nmax+1),A,trim(fname)) 
  
end program
