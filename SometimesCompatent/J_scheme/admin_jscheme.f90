module admin_jscheme
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~  
TYPE :: spd    ! single particle discriptor
     INTEGER :: total_orbits,Jtotal_max,lmax
     INTEGER, ALLOCATABLE,DIMENSION(:) :: nn, ll, jj, itzp, nshell, mvalue,con
     ! for clarity:  nn, ll, nshell are all the true value
     ! jj is j+1/2 (so it's an integer) 
     ! likewise itzp is 2*tz 
     CHARACTER (LEN=10), ALLOCATABLE,DIMENSION(:) :: orbit_status, model_space,basis
     REAL(8), ALLOCATABLE,DIMENSION(:) :: e, evalence, e_original
END TYPE spd 
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
type sp_block_mat
   real(8),allocatable,dimension(:,:) :: matrix!,eigvec
   real(8),allocatable,dimension(:) :: Eigval,extra
   integer,allocatable,dimension(:,:) :: labels
   integer,dimension(2) :: lmda
end type sp_block_mat
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
type full_sp_block_mat
   type(sp_block_mat),allocatable,dimension(:) :: blkM
   integer,allocatable,dimension(:) :: map
   integer :: blocks
end type full_sp_block_mat
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
type array_elem
   integer,dimension(2) :: lam !(J,parity)  
   real(8),allocatable,dimension(:,:) :: V
   integer,allocatable,dimension(:,:) :: qn
end type array_elem
   
type twobody_array
   type(array_elem),allocatable,dimension(:) :: mat
   integer :: Abody,nblocks 
end type twobody_array 

contains
!=========================================================
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!=========================================================
subroutine allocate_tp_array(jbas,int) 
  implicit none 
  
  type(spd) :: jbas
  type(twobody_array) :: int
  integer :: JJ,par,q,parnow,Jnow,r
  integer :: i,j
  ! parity = 0 is even
  ! parity = 1 is odd 
  
  int%nblocks = 2 * jbas%Jtotal_max + 2
  
  q = 1
  do par = 0,1
     do JJ = 0,jbas%Jtotal_max
        
        int%mat(q)%lam(1) =  JJ 
        int%mat(q)%lam(2) = par
        
        r = 0
        do i = 1,jbas%total_orbits
           do j = i,jbas%total_orbits
              
              parnow = mod(jbas%ll(i) + jbas%ll(j) ,2) 

              if ( (parnow == par) .and.&
                   triangle(jbas%jj(i),jbas%jj(j),JJ) ) then 
                 r = r + 1
              end if 
           end do 
        end do 
        
        allocate(int%mat(q)%V(r,r)) 
        allocate(int%mat(q)%qn(r,2)) 
        
        r = 1
        do i = 1,jbas%total_orbits
           do j = i,jbas%total_orbits
              
              parnow = mod(jbas%ll(i) + jbas%ll(j) ,2) 

              if ( (parnow == par) .and.&
                   triangle(jbas%jj(i),jbas%jj(j),JJ) ) then 
                 
                 int%mat(q)%qn(r,1) =  i
                 int%mat(q)%qn(r,2) =  j 
                 r = r + 1
              end if 
           end do 
        end do 
                 
        q = q + 1
     end do 
  end do 

end subroutine   
!=========================================================
!========================================================= 
subroutine allocate_sp_descript(jbas,A) 
  implicit none 
  
  type(spd) :: jbas 
  character(8) :: num,part,inside
  character(45) :: dumb
  integer :: ist,i,label,ni,li,ji,tzi,ix,holes,e1,A
  real(8) :: e,eval
  
  open(unit=31,file = 'spJ.dat') 
  
  ! skip the bull
  read(31,*); read(31,*); read(31,*); read(31,*); read(31,*)
  read(31,*); read(31,*); read(31,*)
  
  read(31,'(A,I10)') dumb, ix
 
  ix = ix/2 ! only neutrons
  read(31,*) 
 
  ! build the jscheme descriptor
  jbas%total_orbits=ix
  allocate(jbas%nn(ix)) ! n 
  allocate(jbas%ll(ix)) ! l
  allocate(jbas%jj(ix)) ! j*2
  allocate(jbas%itzp(ix)) !isospin*2
  allocate(jbas%nshell(ix)) !shell number
  allocate(jbas%e(ix)) ! sp energies
  allocate(jbas%con(ix)) ! hole or particle (1 or 0) 
 
  ! read in info
  do i = 1,2*ix
    
     read(31,*) num,label,ni,li,ji,tzi,e1,e,eval,part,inside

     if (tzi == -1) cycle ! only neutrons
     jbas%nn(i/2) = ni 
     jbas%ll(i/2) = li
     jbas%jj(i/2) = ji
     jbas%itzp(i/2) = tzi 
     jbas%nshell(i/2) = e1
     jbas%e(i/2) =  e 

  end do 
  close(31)
  
  jbas%Jtotal_max = maxval(jbas%jj) 
  jbas%lmax = maxval(jbas%ll) 
  call find_holes(jbas,A) 
  
end subroutine 
!==============================================
!==============================================
subroutine find_holes(jbas,holes) 
  implicit none 
  
  type(spd) :: jbas
  integer :: holes,i,minpos(1)
  real(8),dimension(jbas%total_orbits) :: temp
  
  temp = jbas%e
  jbas%con = 0 ! zero if particle, one if hole
  
  do i = 1,holes
     minpos = minloc(temp) !find lowest energy 
     jbas%con(minpos(1)) = 1 ! thats a hole
     temp(minpos(1)) = 9.e9  ! replace it with big number
  end do 
     
end subroutine 
!=======================================================
!=======================================================
subroutine allocate_sp_mat(jbas,H) 
  ! allocate block sizes for sp matrix
  implicit none 
  
  type(spd) :: jbas
  type(full_sp_block_mat) :: H 
  integer :: jmax,lmax,m,blcks,d
  integer :: q,r,i,j ,l
  
  ! how many blocks are there 
  jmax = jbas%Jtotal_max
  lmax = jbas%lmax
  
  blcks = (jmax+1)*(lmax+1)/2 
  H%blocks = blcks
  
  allocate(H%blkM(blcks))
  allocate(H%map(blcks)) 
  
  ! find how many states in each block, allocate
  q = 1
  do l = 0,lmax
     do j = 1,jmax,2
        H%blkM(q)%lmda(1) = l
        H%blkM(q)%lmda(2) = j 
        
        r = 0
        do i = 1, jbas%total_orbits
           if (l == jbas%ll(i)) then
              if (j == jbas%jj(i)) then 
                 
                 r = r + 1
              end if
           end if 
        end do 
        
        allocate(H%blkM(q)%matrix(r,r)) 
        allocate(H%blkM(q)%eigval(r))
        allocate(H%blkM(q)%extra(10*r)) 
        
        H%map(q) = r 
        q = q + 1
     end do 
  end do 
  
end subroutine  
!==========================================
!==========================================  
subroutine duplicate_sp_mat(H,T) 
  ! use all of the information about H to make 
  ! an empty block matrix of the same size
  implicit none 
  
  type(full_sp_block_mat) :: H,T 
  integer :: q,r
  
  T%blocks = H%blocks
  
  allocate(T%map(H%blocks)) 
  allocate(T%blkM(H%blocks)) 
  
  do q = 1, H%blocks
     T%blkM(q)%lmda = H%blkM(q)%lmda
     T%map(q) = H%map(q) 
     r = H%map(q) 
     allocate(T%blkM(q)%matrix(r,r)) 
     allocate(T%blkM(q)%eigval(r)) 
     allocate(T%blkM(q)%extra(10*r)) 
  end do 

end subroutine
!=====================================================
!===================================================== 
!fuck
!==================================================================  
!==================================================================
logical function triangle(j1,j2,J) 
  ! all entries are multiplied by 2 (to eliminate half integers) 
  implicit none 
  
  integer :: j1,j2,J
  

  if ( ( J .le. j1+j2 ) .and. (J .ge. abs(j1-j2) ) )then
     triangle = .true. 
  else 
     triangle = .false.
  end if 
end function 
end module
