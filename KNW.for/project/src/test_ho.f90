PROGRAM test_ho

USE types
USE Ho_basis
USE numerical_integration
USE symmat_mod

IMPLICIT NONE

INTEGER :: gpts
!ho stuff
REAL(kind=r_kind) :: bosc, ol
INTEGER :: n, npr, l, lpr
REAL(kind=r_kind), ALLOCATABLE :: wf1(:), wf2(:)

INTEGER :: nmax, II
REAL(kind=r_kind), allocatable :: eigs(:)


gpts = 70 !larger than 70 causes numerical errors or floating poit exeptions

CALL init_grid_GH(gpts)

ALLOCATE(wf1(grid_size_GH),wf2(grid_size_GH))

bosc = 1.0_r_kind

n = 50
npr = 50
l = 1
lpr = 1

CALL RadHO(n,l,bosc,grid_points_GH,wf1,grid_size_GH)

CALL RadHO(npr,lpr,bosc,grid_points_GH,wf2,grid_size_GH)

ol  = overlap_ho(n,l,npr,lpr)

WRITE(*,*) 'ol = ', ol




WRITE(*,*) 'Input nmax'
READ(*,*) nmax
WRITE(*,*) 'nmax = ',nmax


CALL init_ho_basis(bosc,nmax,0)

!CALL print_matrix(nmax+1,h_sp)

ALLOCATE(eigs(0:nmax))

CALL calculate_eigs_symmetric_matrix(nmax+1,h_sp,eigs)


WRITE(*,*) 'First 3 eigenvalues'
DO II=0,min(ABS(nmax),3)
   WRITE(*,'(f16.10)') eigs(II)
END DO


END PROGRAM test_ho




