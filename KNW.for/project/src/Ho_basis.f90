MODULE Ho_basis
  USE types
  USE numerical_integration
  USE constants

  IMPLICIT NONE

 ! PRIVATE

  !PUBLIC :: RadHO, overlap_ho, Ho_b, init_ho_basis, init_ho_basis_old, h_sp

! not supported in gfortran
!!$  TYPE Ho_block(nmax)
!!$     INTEGER, len :: nmax
!!$     INTEGER :: l,s
!!$     COMPLEX(kind=r_kind) :: block(0:nmax,0:nmax)
!!$
!!$  END type Ho_block

  !Ho_block, allocatable, protected :: hf_transform(:,:)

  REAL(kind=r_kind), protected :: Ho_b, Ho_hbaromega
  INTEGER, protected :: Ho_nmax, Ho_lmax  !TODO change Ho_nmax to array containing different nmax for different l

  REAL(kind=r_kind), allocatable, protected:: h_sp(:,:) !used in testing/coulomb excersize

  ! 2j = 2l + (-1)**jindex

  COMPLEX(kind=r_kind), allocatable, protected :: hf_transform(:,:,:,:) !l,jindex, n,npr
  COMPLEX(kind=r_kind), allocatable, protected :: density_matrix(:,:,:,:) !l,jindex, n,npr
  REAL(kind=r_kind), allocatable, protected :: hf_energies(:,:,:) !l,jindex,n
  REAL(kind=r_kind), allocatable, protected :: hf_energies_sorted(:,:,:)
  INTEGER, allocatable, protected :: fermi_index_hfbasis(:,:) !l,jindex
  COMPLEX(kind=r_kind), allocatable, protected :: hf_hamiltonian(:,:,:,:) !l,jindex,n,npr
  COMPLEX(kind=r_kind), allocatable, protected :: Ho_two_body_matel(:,:,:,:,:,:) !l,jindex,n1,n2,n3,n4


CONTAINS

!!$  !sorts the  columns of the hf_transformation according to hf energies  
!!$  SUBROUTINE sort_hf
!!$    
!!$  END SUBROUTINE sort_hf
  

  INTEGER FUNCTION Ho_degeneracy(l,jindex)
    IMPLICIT NONE
    INTEGER, intent(in) :: l,jindex
    INTEGER :: j2

    IF(jindex < 0 .or. jindex > 1) THEN
       WRITE(*,*) 'In function Ho_degenerace, incorrect jindex'
       STOP
    END IF

    IF(jindex == 0) j2 = 2*l + 1
    IF(jindex == 1) j2 = 2*l - 1

    Ho_degeneracy = j2+1

  END FUNCTION Ho_degeneracy

  SUBROUTINE hf_find_fermi_energy(N)
   
    USE sorting, only : sortrx

     IMPLICIT NONE

    INTEGER :: N !number of particles

   

    

  END SUBROUTINE hf_find_fermi_energy
  

  !rho(nljm,n'l'j'm') = rho(lj;nn')d_ll'd_jj'd_mm'
  !the reduced density rho(lj;nn') is stored, keep track of degeneracy!
  SUBROUTINE hf_update_density_matrix
    IMPLICIT NONE

    INTEGER :: l,jindex,n1,n2,II

    COMPLEX(kind=r_kind), parameter :: zero_c = (0.0_r_kind,0.0_r_kind)
    COMPLEX(kind=r_kind) :: sum_occupied    

    DO l=0,Ho_lmax
       DO jindex = 0,1
          IF(l==0 .and. jindex==1) CYCLE

          
          DO n1 = 0,Ho_nmax
             DO n2 = n1,Ho_nmax

                sum_occupied = zero_c

                DO II = 1,fermi_index_hfbasis(l,jindex)  !hf_transform should be sorted according to hf sp-energies in each block
                   
                   sum_occupied = sum_occupied + hf_transform(l,jindex,n1,II)*CONJG(hf_transform(l,jindex,n2,II))

                END DO
                
                density_matrix(l,jindex,n1,n2) = density_matrix(l,jindex,n1,n2) + sum_occupied !Ho_degeneracy(l,jindex)*sum_occupied
                density_matrix(l,jindex,n2,n1) = CONJG(density_matrix(l,jindex,n2,n1))
                
             END DO
          END DO
       END DO
    END DO

  END SUBROUTINE hf_update_density_matrix


  SUBROUTINE hf_update_hamiltonian
    IMPLICIT NONE

    INTEGER :: l,jindex,n1,n2,n3,n4
    COMPLEX(kind=r_kind), parameter :: zero_c = (0.0_r_kind,0.0_r_kind)
    COMPLEX(kind=r_kind) :: Gamma

    DO l=0,Ho_lmax
       DO jindex = 0,1
          IF(l==0 .and. jindex==1) CYCLE
          DO n1 = 0,Ho_nmax
             DO n2 = n1,Ho_nmax

                Gamma = zero_c
                DO n3 = 0,Ho_nmax
                   DO n4 = 0,Ho_nmax !check symmetries of Gamma

                      Gamma = Gamma + Ho_degeneracy(l,jindex)*Ho_two_body_matel(l,jindex,n1,n4,n2,n3)*density_matrix(l,jindex,n3,n4) 

                   END DO
                END DO

                !hf_hamiltonian(l,jind,n1,n2) = Ho_one_body_matel(l,jind,n1,n2) + Gamma
                hf_hamiltonian(l,jindex,n1,n2) = Ho_hbaromega  + Gamma

             END DO
          END DO
       END DO
    END DO

  END SUBROUTINE hf_update_hamiltonian


  SUBROUTINE hf_diagonalize
    IMPLICIT NONE

    INTEGER :: l,jindex
    COMPLEX(kind=r_kind), ALLOCATABLE :: work(:)
    REAL(kind=r_kind), ALLOCATABLE :: rwork(:)
    INTEGER :: lwork, info

    INTERFACE 
       SUBROUTINE zheev(jobz, uplo, n, a, lda, w, work, lwork, rw, info)
         CHARACTER(len=1) ::  jobz, uplo
         INTEGER          ::  info, lda, lwork, n
         REAL(KIND=8)     ::  rw(*), w(*)
         COMPLEX(8)   ::  a(lda, *), work(*)
       END SUBROUTINE zheev
    END INTERFACE

    

    DO l=0,Ho_lmax
       DO jindex = 0,1
          IF(l==0 .and. jindex==1) CYCLE

                ALLOCATE(work(1))
                ALLOCATE(rwork( max(1,3*(Ho_nmax+1)-2) ) )
                lwork = -1

                CALL zheev('V','U',Ho_nmax+1, hf_transform(l,jindex,:,:),Ho_nmax+1,hf_energies(l,jindex,:),work,lwork,rwork,info)

                lwork = work(1)
                DEALLOCATE(work)
                ALLOCATE(work(lwork))

                CALL zheev('V','U',Ho_nmax+1, hf_transform(l,jindex,:,:),Ho_nmax+1,hf_energies(l,jindex,:),work,lwork,rwork,info)

                DEALLOCATE(work)
                DEALLOCATE(rwork)

         
       END DO
    END DO

  END SUBROUTINE hf_diagonalize


  SUBROUTINE init_Ho_basis(b,nmax,lmax)
    IMPLICIT NONE

    REAL(kind=r_kind), intent(in) :: b
    INTEGER, intent(in) :: nmax, lmax

    INTEGER :: n1,n2,l1,l2

    WRITE(*,*) 'Initializing ho basis'

    !currently only allowing for lmax = 0, need to consider better structure for block diagonality for lmax>0
    IF(lmax > 0) THEN
       WRITE(*,*) 'Only lmax = 0 supportet ATM, quitting'
       STOP
    END IF

    CALL set_Ho_b(b)


  END SUBROUTINE init_Ho_basis



  SUBROUTINE init_ho_basis_old(b,nmax,lmax)
    IMPLICIT NONE
    
    REAL(kind=r_kind), intent(in) :: b
    INTEGER, intent(in) :: nmax, lmax

    INTEGER :: n1,n2,l1,l2

    WRITE(*,*) 'Initializing ho basis'

    !currently only allowing for lmax = 0, need to consider better structure for block diagonality for lmax>0
    IF(lmax > 0) THEN
       WRITE(*,*) 'Only lmax = 0 supportet ATM, quitting'
       STOP
    END IF

    CALL set_Ho_b(b)
    
    !TODO fix this
    Ho_hbaromega = 1.0_r_kind
    !


    Ho_nmax = nmax
    Ho_lmax = lmax

    IF(ALLOCATED(h_sp)) DEALLOCATE(h_sp)

    ALLOCATE(h_sp(0:nmax,0:nmax))    

    !
    l1 = 0
    l2 = 0
    !

   
    WRITE(*,*) 'Calculating matrix elements'
    WRITE(*,*) 'Ho_nmax =', Ho_nmax
    
    DO n1=0,Ho_nmax
       DO n2=n1,Ho_nmax

          h_sp(n1,n2) = one_body_radial_matel_GH(n1,l1,n2,l2,coulomb_pot)+kinetic_matel_analytic(n1,l1,n2,l2)
          h_sp(n2,n1) = h_sp(n1,n2)
          !WRITE(*,*) 'n1,n2 = ',n1,n2
       END DO
    END DO



  END SUBROUTINE init_ho_basis_old


  REAL(kind=r_kind) FUNCTION coulomb_pot(r)
    IMPLICIT NONE
    REAL(kind=r_kind), intent(in) :: r
    !fix units    
    coulomb_pot = -1/ABS(r)
  END FUNCTION coulomb_pot



  SUBROUTINE set_Ho_b(b)
    IMPLICIT NONE
    REAL(kind=r_kind) :: b

    IF(b<=0) THEN
       WRITE(*,*) 'b<=0 not allowed'
       STOP
    END IF
    Ho_b = b

    Ho_hbaromega = hbarc**2/(Ho_b**2*mnc2)

  END SUBROUTINE set_Ho_b

  !The prefactor for radial ho wave funciton with osc. lenght b = 1
  REAL(kind=r_kind) FUNCTION radial_wf_prefactor(n,l)

    IMPLICIT NONE

    INTEGER, intent(in) :: n,l
    REAL(kind=r_kind) :: lnfac

    IF (n == 0) THEN
      lnfac = 0.0_r_kind
    ELSE
      lnfac = gammln(DBLE(n)+1.0_r_kind)
    END IF

    !WRITE(*,*) 'lfac in function = ' , lnfac
    
    radial_wf_prefactor = sqrt(2.0_r_kind)&
         * EXP( 0.5_r_kind*(lnfac - gammln(DBLE(n)+DBLE(l)+ 1.5_r_kind)))
    
     !radial_wf_prefactor = 2.0_r_kind * EXP( 0.5_r_kind*(lnfac + lnfac - gammln(DBLE(n)+DBLE(l)+ 1.5_r_kind) - gammln(DBLE(n)+DBLE(l)+ 1.5_r_kind)) )
     !radial_wf_prefactor = sqrt(radial_wf_prefactor)

    

  END FUNCTION radial_wf_prefactor
  
  REAL(kind=r_kind) FUNCTION overlap_ho(n1,l1,n2,l2)
    

    IMPLICIT NONE
    INTEGER, intent(in) :: n1,l1,n2,l2
    INTEGER :: II

    REAL(kind=r_kind), ALLOCATABLE :: f1(:),f2(:),fw(:)
    REAL(kind=r_kind) :: lnfac1, lnfac2


    overlap_ho = 0.0

    IF(.not. is_init_grid_GH) THEN
       WRITE(*,*) 'GH grid is not initialized'
       STOP
    END IF

    ALLOCATE(f1(grid_size_GH),f2(grid_size_GH),fw(grid_size_GH))
    
    CALL LaguerreL2(n1, l1, grid_points_GH**2, f1, fw, grid_size_GH)
    CALL LaguerreL2(n2, l2, grid_points_GH**2, f2, fw, grid_size_GH)

    DO II = 1,grid_size_GH
       overlap_ho = overlap_ho + grid_weights_GH(II)*grid_points_GH(II)**(l1+l2+2)*f1(II)*f2(II)  
       !overlap_ho = overlap_ho + grid_points_GH(II)**2*grid_weights_GH(II)
    END DO

!!$    IF (n1 == 0) THEN
!!$      lnfac1 = 0.0_r_kind
!!$    ELSE
!!$      lnfac1 = gammln(DBLE(n1)+1.0_r_kind)
!!$    END IF
!!$    IF (n2 == 0) THEN
!!$      lnfac2 = 0.0_r_kind
!!$    ELSE
!!$      lnfac2 = gammln(DBLE(n2)+1.0_r_kind)
!!$    END IF

    
    !overlap_ho = 2.0_r_kind * EXP( 0.5_r_kind*(lnfac1 + lnfac2 - gammln(DBLE(n1)+DBLE(l1)+ 1.5_r_kind) - gammln(DBLE(n2)+DBLE(l2)+ 1.5_r_kind)) ) * overlap_ho

    overlap_ho = radial_wf_prefactor(n1,l1)*radial_wf_prefactor(n2,l2)* overlap_ho
    
    !WRITE(*,*) 'lnfac1 = ', lnfac1, 'lnfac2= ',lnfac2

    !WRITE(*,*) 2.0_r_kind * EXP( 0.5_r_kind*(lnfac1 + lnfac2 - gammln(DBLE(n1)+DBLE(l1)+ 1.5_r_kind) - gammln(DBLE(n2)+DBLE(l2)+ 1.5_r_kind)) )
    !WRITE(*,*) radial_wf_prefactor(n1,l1)*radial_wf_prefactor(n2,l2)


    DEALLOCATE(f1,f2,fw)

    RETURN 

  END FUNCTION overlap_ho


  REAL(kind=r_kind) FUNCTION one_body_radial_matel_GH(n1,l1,n2,l2,potential_function)
    
    IMPLICIT NONE
    
    INTERFACE       
       !REAL(kind=r_kind) FUNCTION potential_function(r) !cannot acces r_kind in this scope
       REAL(8) FUNCTION potential_function(r)
         !REAL(kind=r_kind), intent(in) :: r
         REAL(8), intent(in) :: r
       END FUNCTION potential_function
    END INTERFACE

    INTEGER, intent(in) :: n1, l1, n2, l2  
    INTEGER :: II
    REAL(kind=r_kind) :: matel

    REAL(kind=r_kind), ALLOCATABLE :: f1(:),f2(:),fw(:)
    REAL(kind=r_kind) :: lnfac1, lnfac2

    IF(.not. is_init_grid_GH) THEN
       WRITE(*,*) 'GH grid is not initialized'
       STOP
    END IF

    ALLOCATE(f1(grid_size_GH),f2(grid_size_GH),fw(grid_size_GH))

    CALL LaguerreL2(n1, l1, grid_points_GH**2, f1, fw, grid_size_GH)
    CALL LaguerreL2(n2, l2, grid_points_GH**2, f2, fw, grid_size_GH)

   
    matel = 0.0_r_kind

    DO II = 1,grid_size_GH
       matel = matel + grid_weights_GH(II)*potential_function(Ho_b*ABS(grid_points_GH(II)))*ABS(grid_points_GH(II))**(2+l1+l2)*f1(II)*f2(II)
       
       !matel = matel + ABS(grid_points_GH(II))**(2+l1+l2)*f1(II)*f2(II)
       !for debug
       !WRITE(*,*) matel
       !end for debug

    END DO

    matel = radial_wf_prefactor(n1,l1)*radial_wf_prefactor(n2,l2)*matel 

    DEALLOCATE(f1,f2,fw)

    one_body_radial_matel_GH = matel

    RETURN 

    
  END FUNCTION one_body_radial_matel_GH

!!$  SUBROUTINE test
!!$
!!$    IMPLICIT NONE
!!$    REAL(kind=r_kind) :: x
!!$
!!$    WRITE(*,*) is_init_grid_GH
!!$
!!$  END SUBROUTINE test


  REAL(kind=r_kind) FUNCTION kinetic_matel_analytic(n1,l1,n2,l2)
    IMPLICIT NONE
    INTEGER, intent(in) :: n1,l1,n2,l2
    REAL(kind=r_kind) :: matel

    IF(l1 /= l2) THEN
      kinetic_matel_analytic = 0.0_r_kind
       RETURN
    END IF
    
    SELECT CASE(n1-n2)
    CASE(0)
       matel = 2*n1 + l1 + 1.5_r_kind
    CASE(-1)
       matel = sqrt(n2*(n2+l1+0.5_r_kind))
    CASE(1)
       matel = sqrt(n1*(n1+l1+0.5_r_kind))
    CASE DEFAULT
       matel = 0.0_r_kind
    END SELECT

   kinetic_matel_analytic = 0.5_r_kind * Ho_hbaromega * matel 


  END FUNCTION kinetic_matel_analytic



  ! log(Gamma(xx)) From numerical recipies
  FUNCTION gammln(xx)
    IMPLICIT NONE
    DOUBLE PRECISION gammln,xx
    ! Returns the value ln[GAM(xx)] for xx > 0.
    INTEGER j
    DOUBLE PRECISION ser,stp,tmp,x,y,cof(6)
    SAVE cof,stp
    DATA cof,stp/76.18009172947146d0,-86.50532032941677d0,&
          24.01409824083091d0,-1.231739572450155d0,.1208650973866179d-2,&
          -.5395239384953d-5,2.5066282746310005d0/
    x=xx
    y=x
    tmp=x+5.5d0
    tmp=(x+0.5d0)*log(tmp)-tmp
    ser=1.000000000190015d0
    do j=1,6
       y=y+1.d0
       ser=ser+cof(j)/y
    enddo
    gammln=tmp+log(stp*ser/x)
    return
  END FUNCTION gammln

  
!  L = L_(n)^(l+1/2) ! Using recursion relations

! Originally returned sign(L) and log(abs(L)) /Daniel W
  Subroutine LaguerreL2(n, l, RVEC, FVEC, FVEC_S, dimr) 
    IMPLICIT NONE
    INTEGER :: n,l,Ik,II, dimr
    DOUBLE PRECISION :: RVEC(dimr), FVEC(dimr),FVEC_S(dimr),bin,alpha
    DOUBLE PRECISION :: L0(dimr), L1(dimr)

    IF(n.lt.0) THEN
       FVEC   = 0.D0
       FVEC_S = 0.D0
       RETURN
    END IF

    alpha = l + 0.5d0
    
    L1 = 0.d0; FVEC = 1.d0;
    DO ii = 1 , n 
       L0 = L1
       L1 = FVEC
       FVEC = ((2.d0 * ii- 1.d0 + alpha- RVEC)* L1- (ii- 1.d0 + alpha)* L0)/ dble(ii)
    END DO
  
    FVEC_S = sign(1.D0,FVEC)
    !Commented away statement below I want the value not the log /Daniel W
    !FVEC   = log(abs(FVEC)) 
    !End Comment /Daniel W
    RETURN 
  END Subroutine LaguerreL2

! Calculates a vector of values of g_nl(r), where g_nl(r) is the radial part of a HO wave function with size parameter b as defined
! on page 49 of Jouni Sohonen, From Nucleons to Nucleus, Springer
  Subroutine RadHO(n, l, b, RVEC, FVEC, dimr)
    IMPLICIT NONE
    INTEGER :: n,l,dimr
    DOUBLE PRECISION :: nR,lR,b, RVEC(dimr), FVEC(dimr), FVEC_S(dimr), FVECtmp(dimr), lnfac
    nR = DBLE(n) 
    lR = DBLE(l)
    CALL LaguerreL2(n, l, (RVEC/b)**2, FVEC, FVEC_S, dimr)
     
    ! gamma(n+1) = n!
    IF (n == 0) THEN
      lnfac = 0d0
    ELSE
      lnfac = gammln(nR+1d0)
    END IF

   FVEC = SQRT(2d0/b**3)* EXP(0.5d0* ( lnfac - gammln(nR+lR+1.5d0) ) )* (RVEC/b)**l* EXP(-RVEC**2/2/b**2)* FVEC 
   
  END Subroutine RadHO

! Computes overlap certain of radial HO functions < R_n0^(bN) | R_00^(ba) > with a closed expression
! n is nubmer of nodes of radial ho-wave function with l=0 and oscillator length bN, ba is oscillator length of n=0,l=0 wave function
  FUNCTION olrho(n,bN,ba)
    ! Returns the value of < R_n0^(bN) | R_00^(ba) >
    IMPLICIT NONE
    INTEGER :: n
    DOUBLE PRECISION :: nDbl, olrho, bN, ba  
    IF (n<0 .OR. bN<=0d0 .OR. ba<=0d0) THEN
      olrho = 0d0
      RETURN
    END IF

    IF (n==0) THEN
      olrho = ( 2d0*bN*ba/(bN**2+ba**2) )**(3d0/2d0)
      RETURN
    END IF
    
    nDbl = DBLE(n)
    olrho = EXP( 0.5d0*gammln(2*nDbl + 2) - gammln(nDbl + 1) )*2d0**(3d0/2d0 - nDbl)*(bN*ba)**(3d0/2d0)*(bN**2 - ba**2)**nDbl/(bN**2 + ba**2)**(3d0/2d0+nDbl)
    

  END FUNCTION olrho

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  FUNCTION ylm(l,m,theta, phi )
    IMPLICIT NONE
  INTEGER :: l, m
  DOUBLE PRECISION :: theta, phi, lD, mD 
  COMPLEX :: ylm
  DOUBLE PRECISION , PARAMETER ::  pi = 2d0*ACOS(0.0)
  
  IF(ABS(m)>l) THEN
    ylm = 0d0
    RETURN
  END IF
  
  lD = DBLE(l)
  IF(m >= 0) THEN
    mD = DBLE(m)
    ylm = sqrt( (2*l+1d0)/4d0/pi ) * exp( .5d0*( gammln(lD-mD+1d0) - gammln(lD+mD+1d0) ) )*plgndr(l,m,COS(theta))*exp( (0d0,1d0)*m*phi )
    RETURN
  ELSE !Y_l(-m) = (-1)^m Y^*_lm
    mD = DBLE(-1d0*m)
    m = -1*m
    ylm = (-1d0)**m*sqrt( (2*l+1d0)/4d0/pi ) * exp( .5d0*( gammln(lD-mD+1d0) - gammln(lD+mD+1d0) ) )*plgndr(l,m,COS(theta))*exp( (0d0,-1d0)*m*phi )
  END IF  

  END FUNCTION ylm



!FROM NUMERICAL RECIPES IN FORTRAN77 and FORTRAN 90 PressW.H. et al. Cambridge University Publishing 1997
!page 247
! CHANGED REAL TO DOUBLE PRECISION and removed do labels

  FUNCTION plgndr(l,m,x)
    IMPLICIT NONE
  INTEGER l,m
  DOUBLE PRECISION :: plgndr,x
  !Computes the associated Legendre polynomial Plm (x). Here m and l are integers satisfying
  !0 <= m <= l, while x lies in the range -1 <= x <= 1.
  INTEGER i,ll
  DOUBLE PRECISION :: fact,pll,pmm,pmmp1,somx2
  if(m.lt.0.or.m.gt.l.or.abs(x).gt.1.) THEN !pause ?bad arguments in plgndr?
    plgndr = 0d0
    RETURN
  end if
  
  pmm=1.0d0  !Compute P_m^m
  if(m.gt.0) then
      somx2=sqrt((1.0d0-x)*(1.0d0+x))
      fact=1.0d0
      do i=1,m
	pmm=-pmm*fact*somx2
	fact=fact+2.d0
      enddo 
  endif
  if(l.eq.m) then
      plgndr=pmm
  else
    pmmp1=x*(2d0*m+1d0)*pmm !Compute P_m+1^m 
    if(l.eq.m+1) then
      plgndr=pmmp1
    else !Compute P_l^m , l > m + 1.
      do ll=m+2,l
	pll=(x*(2d0*ll-1d0)*pmmp1-(ll+m-1d0)*pmm)/(ll-m)
	pmm=pmmp1
	pmmp1=pll
      enddo 
      plgndr=pll
    endif
  endif
  return
  END FUNCTION plgndr
    

END MODULE
