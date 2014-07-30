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
   real(8),allocatable,dimension(:,:) :: matrix
   real(8),allocatable,dimension(:) :: Eigval,extra
   integer,allocatable,dimension(:,:) :: labels
   integer,allocatable,dimension(:) :: states
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
   integer,allocatable,dimension(:) :: qn
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
  integer :: i,j,N

  ! parity = 0 is even
  ! parity = 1 is odd 

! total number of blocks in tp basis  
  N = jbas%total_orbits
  int%nblocks = 2 * jbas%Jtotal_max + 2
  
  q = 1
! allocate block storage
  allocate(int%mat(int%nblocks)) 

  ! go through blocks and find the dimension of each
  ! allocate space
  do par = 0,1
     do JJ = 0,jbas%Jtotal_max

        
        ! store the conserved quantities
        int%mat(q)%lam(1) =  JJ 
        int%mat(q)%lam(2) = par
        
        r = 0
        do i = 1,jbas%total_orbits
           do j = i,jbas%total_orbits
              
              parnow = mod(jbas%ll(i) + jbas%ll(j) ,2) 

              if ( (parnow == par) .and.&
                   triangle(jbas%jj(i),jbas%jj(j),2*JJ) ) then 
                 r = r + 1 ! dimension of this block
              end if 
           end do 
        end do 
        
        allocate(int%mat(q)%V(r,r)) 
        allocate(int%mat(q)%qn(r)) 
        
        r = 1
        do i = 1,jbas%total_orbits
           do j = i,jbas%total_orbits
              
              parnow = mod(jbas%ll(i) + jbas%ll(j) ,2) 

              if ( (parnow == par) .and.&
                   triangle(jbas%jj(i),jbas%jj(j),2*JJ) ) then 
                 
                 int%mat(q)%qn(r) =  (i-1)*N + j   ! sp quantum numbers stored
          
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
  integer :: holes,i,minpos(1),q
  real(8),dimension(jbas%total_orbits) :: temp
  
  temp = jbas%e
  jbas%con = 0 ! zero if particle, one if hole
  
  i = 1
  q = 0
  do while (q < holes)
     minpos = minloc(temp) !find lowest energy 
     jbas%con(minpos(1)) = jbas%jj(i) + 1 ! thats a hole
     temp(minpos(1)) = 9.e9  ! replace it with big number
     q = q + jbas%jj(i) + 1
     i = i + 1
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
        
        H%map(q) = r
        allocate(H%blkM(q)%matrix(r,r)) 
        allocate(H%blkM(q)%eigval(r))
        allocate(H%blkM(q)%extra(10*r)) 
        allocate(H%blkM(q)%states(r)) 
        
        r = 1
        do i = 1, jbas%total_orbits
           if (l == jbas%ll(i)) then
              if (j == jbas%jj(i)) then 
                 H%blkM(q)%states(r) = i 
                 r = r + 1
              end if
           end if 
        end do 
        

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
  
  T%map = H%map 
  do q = 1, H%blocks
     T%blkM(q)%lmda = H%blkM(q)%lmda
     T%map(q) = H%map(q) 
     r = H%map(q) 
     allocate(T%blkM(q)%matrix(r,r)) 
     allocate(T%blkM(q)%eigval(r)) 
     allocate(T%blkM(q)%extra(10*r)) 
     allocate(T%blkM(q)%states(r))
     T%blkM(q)%states = H%blkM(q)%states
  end do 

end subroutine
!=====================================================
!===================================================== 
subroutine read_tp_matrix_elements(jbas,H) 
  ! reads TBME from morten's file
  implicit none
  
  type(spd) :: jbas
  type(twobody_array) :: H
  integer :: TZ,PAR,JJ,j1,j2,j3,j4,ist
  integer :: p1,p2,Q,pre,I1,I2,NR,N
  integer :: n1,n2,n3,n4
  real(8) :: v1,v2,v3,v4
  
  N = jbas%total_orbits
  
  open(unit=32,file='VJ-scheme.dat') 
  
  do 
     
     read(32,*,iostat=ist) TZ,PAR,JJ,n1,n2,n3,n4, &
          v1,v2,v3,v4
      
    
     if (ist < 0 ) exit
     
     if (TZ .ne. 1) cycle
    
     JJ = JJ/2
     n1 = n1/2
     n2 = n2/2
     n3 = n3/2
     n4 = n4/2  !neutrons only   
     
     ! find sp mapping to tp for the bra and ket
     pre = 1 ! prefactor if the states aren't ordered right
     if (n1 > n2) then 
        I1 = (n2-1)*N +n1 ! unique integer label for 2 particle bra 
        pre = (-1)**( 1+ (jbas%jj(n1) + jbas%jj(n2))/2 - JJ ) 
     else 
        I1 = (n1-1)*N +n2
     end if 
     
     if (n3 > n4) then 
        I2 = (n4-1)*N +n3 
        pre = pre * (-1)**(1+ (jbas%jj(n3) + jbas%jj(n4))/2 -JJ ) 
     else 
        I2 = (n3-1)*N +n4
     end if 
  
     ! get the block index
     Q = (Jbas%Jtotal_max+1)*par + JJ + 1
     
     ! find the tp indeces for the bra and the ket
     p1 = 1 

     do   
        ! search for the integer associated with the 2 particle bra
      
        if (H%mat(Q)%qn(p1) == I1) exit
        p1 = p1 + 1
     end do 

     p2 = 1
     do 
        if (H%mat(Q)%qn(p2) == I2) exit
        p2 = p2 + 1
        
     end do 
     
     ! assign matrix element
     H%mat(Q)%V(p1,p2) = v1*pre
     H%mat(Q)%V(p2,p1) = v1*pre
  end do
  
  close(32) 
  
end subroutine
!==================================================================
!==================================================================
real(8) function v_elem(JJtimestwo,PAR,n1,n2,n3,n4,jbas,H) 
  ! reads TBME from morten's file
  implicit none
  
  type(spd) :: jbas
  type(twobody_array) :: H
  integer :: TZ,PAR,JJtimestwo,j1,j2,j3,j4,ist
  integer :: p1,p2,Q,pre,I1,I2,NR,N
  integer :: n1,n2,n3,n4,PB,PK
  real(8) :: v1,v2,v3,v4
  
  N = jbas%total_orbits
  
  PB = mod(jbas%ll(n1) + jbas%ll(n2),2) 
  PK = mod(jbas%ll(n3) + jbas%ll(n4),2) 
  
  j1 = jbas%jj(n1)
  j2 = jbas%jj(n2)
  j3 = jbas%jj(n3)
  j4 = jbas%jj(n4)

  if (PB .ne. PAR) then 
     v_elem = 0.d0 
  else if (PK .ne. PAR) then
     v_elem = 0.d0 
     
  else if (( j1 == j2 ) .and. (j3 == j4) .and. &
       (mod(JJtimestwo/2,2) == 1)) then
     v_elem = 0.d0  ! this is a minnesota thing
     ! it shouldn't be necessary but I think there is a bug somewhere

  else if ((triangle(j1,j2,JJtimestwo)) .and.& 
       (triangle(j3,j4,JJtimestwo))) then 
  
     ! find sp mapping to tp for the bra and ket
     pre = 1 ! prefactor if the states aren't ordered right
     if (n1 > n2) then 
        I1 = (n2-1)*N +n1 ! unique integer label for 2 particle bra 
        pre = (-1)**( 1 + (jbas%jj(n1) + jbas%jj(n2) -JJtimestwo)/2 ) 
     else 
        I1 = (n1-1)*N +n2
     end if 
     

     if (n3 > n4) then 
        I2 = (n4-1)*N +n3 
        pre = pre *(-1)**( 1 + (jbas%jj(n3) + jbas%jj(n4) -JJtimestwo)/2 ) 
     else 
        I2 = (n3-1)*N +n4
     end if 
  
     ! get the block index
     Q = (Jbas%Jtotal_max+1)*par + JJtimestwo/2 + 1
     
     ! find the tp indeces for the bra and the ket
     p1 = 1
     
     do       
        ! search for the integer associated with the 2 particle bra
        if (H%mat(Q)%qn(p1) == I1) exit
        p1 = p1 + 1
     end do 

     p2 = 1
     do 
        if (H%mat(Q)%qn(p2) == I2) exit
        p2 = p2 + 1
     end do 
     
     ! assign matrix element
     v_elem = H%mat(Q)%V(p1,p2)*pre

     else 
        v_elem = 0.d0
     end if
  
end function
!==================================================================  
!==================================================================
logical function triangle(j1,j2,J) 
  ! all entries are multiplied by 2 (to eliminate half integers) 
  implicit none 
  
  integer :: j1,j2,J
  
  ! is a J = j1 + j2 possible? 
  if ( ( J .le. j1+j2 ) .and. (J .ge. abs(j1-j2) ) )then
     triangle = .true. 
  else 
     triangle = .false.
  end if 
end function
!================================================
subroutine diagonalize_blocks(R)
  implicit none 
  
  type(full_sp_block_mat) :: R
  integer :: a,q,info
  
  do q=1,R%blocks
     
     a=R%map(q)
     if (a == 0) cycle
  
     call dsyev('V','U',a,R%blkM(q)%matrix,a, &
          R%blkM(q)%eigval,R%blkM(q)%extra,10*a,info)
   
  end do 
end subroutine
!=====================================================
integer function kd(a,b) 
  implicit none 
  
  integer :: a,b
  
  if (a == b) then 
     kd = 1
  else 
     kd = 0
  end if 
end function 
!====================================================
end module
