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

  TYPE hf_state
     REAL(kind=r_kind) :: e
     INTEGER :: l,jindex
     INTEGER :: n
     REAL(kind=r_kind) :: occ
  END type hf_state

  REAL(kind=r_kind), protected :: Ho_b, Ho_hbaromega
  INTEGER, protected :: Ho_Nmax, Ho_lmax, Ho_size_all  
  ! Ho_Nmax means the major osc. qn from this version on!



  REAL(kind=r_kind), allocatable, protected:: h_sp(:,:) !used in testing/coulomb excersize

  ! 2j = 2l + (-1)**jindex

  COMPLEX(kind=r_kind), allocatable, protected :: hf_transform(:,:,:,:) !l,jindex, n,npr
  COMPLEX(kind=r_kind), allocatable, protected :: density_matrix(:,:,:,:) !l,jindex, n,npr
  COMPLEX(kind=r_kind), allocatable, protected :: density_matrix_old(:,:,:,:) !l,jindex, n,npr
  REAL(kind=r_kind), allocatable, protected :: hf_energies(:,:,:) !l,jindex,n 
  REAL(kind=r_kind), allocatable, protected :: hf_energies_old(:,:,:) !l,jindex,n 
  !INTEGER, allocatable, protected :: hf_energies_sort(:,:,:)

  type(hf_state), allocatable, protected :: hf_states_all(:) !1..total number of states
  REAL(kind=r_kind), allocatable, protected :: hf_energies_all(:)  
  REAL(kind=r_kind), allocatable, protected :: hf_energies_all_old(:) 
  INTEGER, allocatable, protected :: hf_energies_all_sort(:)

  INTEGER, allocatable, protected :: fermi_index_hfbasis(:,:) !l,jindex
  COMPLEX(kind=r_kind), allocatable, protected :: hf_hamiltonian(:,:,:,:) !l,jindex,n,npr
  

  COMPLEX(kind=r_kind), allocatable, protected :: Ho_two_body_matel(:,:,:,:,:,:) !l,jindex,n1,n2,n3,n4
  !Check the symmtetries of the two-body matrix elements, this form might be wrong. Should work in the case of l=0 j=0.5.

  !The symmetries of the rot invariant interaction gives that we need the following structure
  COMPLEX(kind=r_kind), allocatable, protected :: Ho_two_body_matels(:)
  INTEGER, protected :: Ho_size_TBME
  !Ho_two_body_matels(:,:,:,:,:,:,:,:) !l,jindex, lpr, jindexpr, n1,n2,n3,n4


  COMPLEX(kind=r_kind), allocatable, protected :: hf_gamma(:,:,:,:)
  INTEGER, protected :: hf_num_part_request

CONTAINS

!!$  !sorts the  columns of the hf_transformation according to hf energies  
!!$  SUBROUTINE sort_hf
!!$    
!!$  END SUBROUTINE sort_hf


  !The following assumed that mortens matels where antisymmetric, big mistake!
  ! I found that triplet-triplet scattering of two particles of same species by minnesota was non-zero!!!!!
  SUBROUTINE read_mortens_matels(Nmax,lmax,filename_sp_in,filename_matels_in)
    IMPLICIT NONE

    INTEGER, intent(in) :: Nmax, lmax
    CHARACTER(100), optional :: filename_sp_in, filename_matels_in
    CHARACTER(100) :: filename_sp, filename_matels

    CHARACTER(100) :: slask
    INTEGER :: II,morten_index, n, l, j2, tz, Nmaj, jindex
    INTEGER :: ios, par, a,b,c,d, line, a_temp, c_temp
    REAL(kind=r_kind) :: e_ho, matel_val
    INTEGER :: n1,n2,n3,n4,l1,l2,l3,l4,jindex1,jindex2,jindex3,jindex4,j2max
    INTEGER :: lpr,jindexpr,deg,deg_pr, j12, j12_2, Ho_size_all_morten, Ho_size_TBME_morten
    INTEGER :: j2_1,j2_2


    INTEGER, allocatable :: lookup_mi(:,:,:)

    REAL(kind=r_kind), allocatable :: matels_jcoup(:,:)
    REAL(kind=r_kind) :: norm_L, norm_R
    INTEGER :: sign_ket_perm, sign_bra_perm, ket_perm, bra_perm
    
    TYPE nlj
       INTEGER :: n,l,jindex, j2
    END type nlj
    
    TYPE(nlj), allocatable :: lookup_nlj(:)

    IF(.not. present(filename_sp_in)) THEN
       filename_sp= 'spJ.dat'
    ELSE
       filename_sp = filename_sp_in
    END IF

    IF(.not. present(filename_matels_in)) THEN
       filename_matels= 'VJ-scheme.dat'
    ELSE
       filename_matels = filename_matels_in
    END IF


    !construct lookup table

    OPEN(unit=1,file=TRIM(filename_sp), status="old")

    !first ten lines are legend etc
    DO II = 1,10
       READ(1,*)
    END DO

    

    Ho_size_all_morten = 0

    DO l=0,lmax
        DO jindex = 0,1
           IF(l==0 .and. jindex==1) CYCLE

           Ho_size_all_morten = Ho_size_all_morten + (Nmax-l)/2 + 1

        END DO
     END DO


    ALLOCATE(lookup_mi(0:Nmax/2,0:lmax,0:1))
    ALLOCATE(lookup_nlj(2*Ho_size_all_morten))

    !WRITE(*,*) 'Ho_size_all = ',Ho_size_all

    DO II=1,2*Ho_size_all_morten
       READ(1,*) slask, morten_index, n, l, j2, tz, Nmaj, e_ho
       !WRITE(*,*) morten_index,n,l,j2,tz,Nmaj,e_ho
       
       IF(ABS(e_ho - Ho_hbaromega*(Nmaj + 1.5_r_kind))>1e-10) THEN
          WRITE(*,*) 'Oscillator paramters in sp file:',TRIM(filename_sp),' are inconsistent.'
          STOP
       END IF

       IF(tz == -1) THEN !assuming this means neutron
          IF(j2-2*l == 1) jindex = 0
          IF(j2-2*l == -1) jindex = 1
          lookup_mi(n,l,jindex) = morten_index
          lookup_nlj(morten_index)%n = n
          lookup_nlj(morten_index)%l = l
          lookup_nlj(morten_index)%jindex = jindex
          lookup_nlj(morten_index)%j2 = j2

          WRITE(*,*) 'n,l,jindex = ',n,l,jindex
          WRITE(*,*) 'lookup_mi(n,l,jindex) = ',lookup_mi(n,l,jindex)
       END IF      

    END DO  
       
    CLOSE(1)

    !read matels
    WRITE(*,*) 'Reading J scheme matrix elements from file: ',TRIM(filename_matels)

    OPEN(unit=1,file=TRIM(filename_matels), status="old")

    j2max = 2*lmax + 1
    Ho_size_TBME_morten = (lmax+1)*2*(lmax+1)*2*(Nmax/2+1)**4


    WRITE(*,*) 'Basis for matrix elements read'
    WRITE(*,*) 'Nmax = ',Nmax,'lmax = ',lmax
    WRITE(*,*) 'size_TBME = ',HO_size_TBME_morten
    WRITE(*,*) 'j2max = ',j2max
    WRITE(*,*) 'Basis used in hf'
    WRITE(*,*) 'Ho_Nmax = ',Ho_nmax,'Ho_lmax = ',Ho_lmax
    WRITE(*,*) 'Ho_size_TBME = ',Ho_size_TBME

    
    ALLOCATE(matels_jcoup(0:j2max,Ho_size_TBME_morten))
    matels_jcoup = 0.0_r_kind

    OPEN(unit=3,file="matels_jcoup.dat")


    line = 1

    DO 

       READ(1,*,iostat=ios) tz, par, j12_2, a, b, c, d, matel_val

       !WRITE(*,*) tz, par, j2, a, b, c, d, matel_val

       IF(tz == -1) THEN

          !TODO
          !derive correct normalization factors for j-scheme states
          norm_L = 1.0_r_kind
          norm_R = 1.0_r_kind
          IF(a == b) norm_L = sqrt(2.0_r_kind)
          IF(c == d) norm_R = sqrt(2.0_r_kind)
          !end TODO

          !need to do permutations in a,b,c,d indices, 4! = 4*3*2*1 = 24 

          ! ab|cd  ab|dc  ba|cd  ba|dc  4 perms
          ! cd|ab  dc|ab  cd|ba  dc|ba  4 perms


          DO ket_perm = 1,2
             sign_ket_perm = 1
             IF(ket_perm == 2) THEN
                j2_1 = lookup_nlj(a)%j2
                j2_2 = lookup_nlj(b)%j2
                sign_ket_perm = (-1)**(1+(j12_2-j2_1-j2_2)/2)
                a_temp = a
                a = b
                b = a_temp
             END IF
             DO bra_perm = 1,2
                sign_bra_perm = 1
                IF(bra_perm == 2) THEN
                   j2_1 = lookup_nlj(c)%j2
                   j2_2 = lookup_nlj(d)%j2
                   sign_bra_perm = (-1)**(1+(j12_2-j2_1-j2_2)/2)
                   c_temp = c
                   c = d
                   d = c_temp
                END IF




                l1 = lookup_nlj(a)%l
                jindex1 = lookup_nlj(a)%jindex
                n1 = lookup_nlj(a)%n
                l2 = lookup_nlj(b)%l
                jindex2 = lookup_nlj(b)%jindex
                n2 = lookup_nlj(b)%n
                l3 = lookup_nlj(c)%l
                jindex3 = lookup_nlj(c)%jindex
                n3 = lookup_nlj(c)%n
                l4 = lookup_nlj(d)%l
                jindex4 = lookup_nlj(d)%jindex
                n4 = lookup_nlj(d)%n

                IF(l1 /= l3 .or. jindex1 /= jindex3 .or. l2 /= l4 .or. jindex2 /= jindex4) THEN
!!$             WRITE(*,*) 'Unused matel'
!!$             WRITE(*,*) 'a,b,c,d = ',a,b,c,d
!!$             WRITE(*,*) 'l1, l3, jindex1, jindex3 = ',l1,l3,jindex1,jindex3
!!$             WRITE(*,*) 'l2, l4, jindex2, jindex4 = ',l2,l4,jindex2,jindex4
                   CYCLE
                END IF

                IF(MOD(j12_2,2) /= 0) THEN
                   WRITE(*,*) 'Incorrect jj coupling'
                   STOP
                END IF

                j12 = j12_2 / 2

!!$          !for debug
!!$          WRITE(*,*) 'line = ',line
!!$          WRITE(*,*) 'a,b,c,d = ',a,b,c,d
!!$          WRITE(*,*) 'j12,l1,jindex1,l2,jindex2 = ',j12,l1,jindex1,l2,jindex2
!!$          WRITE(*,*) 'n1,n2,n3,n4 =',n1,n2,n3,n4
!!$
!!$          !end for debug


               
                matels_jcoup(j12,TBME_index_v2(Nmax,lmax,l1,jindex1,l2,jindex2,n1,n2,n3,n4)) = sign_ket_perm*sign_bra_perm*norm_L*norm_R*matel_val
                matels_jcoup(j12,TBME_index_v2(Nmax,lmax,l1,jindex1,l2,jindex2,n3,n4,n1,n2)) = sign_ket_perm*sign_bra_perm*norm_L*norm_R*matel_val

                WRITE(3,'(5I3,F14.8)') j12,a,b,c,d,sign_ket_perm*sign_bra_perm*matel_val
                



             END DO
          END DO

          !TODO derive correct symmetry relations
!!$          matels_jcoup(j12,TBME_index_v2(Nmax,lmax,l1,jindex1,l2,jindex2,n3,n4,n1,n2)) = matels_jcoup(j12,TBME_index_v2(Nmax,lmax,l1,jindex1,l2,jindex2,n1,n2,n3,n4))
!!$
!!$          matels_jcoup(j12,TBME_index_v2(Nmax,lmax,l1,jindex1,l2,jindex2,n2,n1,n4,n3)) = matels_jcoup(j12,TBME_index_v2(Nmax,lmax,l1,jindex1,l2,jindex2,n1,n2,n3,n4))
!!$
!!$          j2_1 = 2*l1 + (-1)**jindex1
!!$          j2_2 = 2*l2 + (-1)**jindex2
!!$          
!!$          matels_jcoup(j12,TBME_index_v2(Nmax,lmax,l1,jindex1,l2,jindex2,n2,n1,n3,n4)) = &
!!$               (-1)**(1+j12-(j2_1+j2_2)/2)*matels_jcoup(j12,TBME_index_v2(Nmax,lmax,l1,jindex1,l2,jindex2,n1,n2,n3,n4))
!!$
!!$          matels_jcoup(j12,TBME_index_v2(Nmax,lmax,l1,jindex1,l2,jindex2,n1,n2,n4,n3)) = &
!!$               (-1)**(1+j12-(j2_1+j2_2)/2)*matels_jcoup(j12,TBME_index_v2(Nmax,lmax,l1,jindex1,l2,jindex2,n1,n2,n3,n4))
!!$
!!$          matels_jcoup(j12,TBME_index_v2(Nmax,lmax,l1,jindex1,l2,jindex2,n4,n3,n1,n2)) = &
!!$               (-1)**(1+j12-(j2_1+j2_2)/2)*matels_jcoup(j12,TBME_index_v2(Nmax,lmax,l1,jindex1,l2,jindex2,n3,n4,n1,n2))
!!$
!!$          matels_jcoup(j12,TBME_index_v2(Nmax,lmax,l1,jindex1,l2,jindex2,n3,n4,n2,n1)) = &
!!$               (-1)**(1+j12-(j2_1+j2_2)/2)*matels_jcoup(j12,TBME_index_v2(Nmax,lmax,l1,jindex1,l2,jindex2,n3,n4,n1,n2))

          !matels_jcoup(j12,TBME_index_v2(Nmax,lmax,l1,jindex1,l2,jindex2,n2,n1,n3,n4)) = matels_jcoup(j12,TBME_index_v2(Nmax,lmax,l1,jindex1,l2,jindex2,n1,n2,n3,n4))          




          !end TODO


       END IF


       IF(ios<0) EXIT
       IF(ios>0) THEN
          WRITE(*,*) 'Errror reading file:',TRIM(filename_matels)
          STOP
       END IF

       line = line +1 

    END DO


    CLOSE(1)

    CLOSE(3)
    WRITE(*,*) 'Done'

    WRITE(*,*) 'Converting matrix elements'
    OPEN(unit=2, file="test_read_morten.txt")


    
   

    DO l=0,MIN(Ho_lmax,lmax)
       DO jindex = 0,1
          IF(l==0 .and. jindex==1) CYCLE
          deg = Ho_degeneracy(l,jindex)
          DO lpr = 0,MIN(Ho_lmax,lmax)
             DO jindexpr = 0,1
                IF(lpr==0 .and. jindexpr==1) CYCLE
                deg_pr = Ho_degeneracy(lpr,jindexpr)

                DO n1 = 0,(MIN(Ho_Nmax,Nmax)-l)/2
                   DO n2 = 0,(MIN(Ho_Nmax,Nmax)-lpr)/2
                      DO n3 = 0,(MIN(Ho_Nmax,Nmax)-l)/2
                         DO n4 = 0,(MIN(Ho_Nmax,Nmax)-lpr)/2 

                            matel_val = 0.0_r_kind
                            
                            j2max = (2*l+(-1)**jindex + 2*lpr+(-1)**jindex)/2
                            !for test
                            !j2max = 0
                            !end for test

                            DO j12 = 0,j2max
                               matel_val = matel_val + REAL((2*j12+1),kind=r_kind)/REAL(deg*deg_pr,kind=r_kind)*matels_jcoup(j12,TBME_index_v2(Nmax,lmax,l,jindex,lpr,jindexpr,n1,n2,n3,n4))

                               !for debug
                               IF(ABS(matels_jcoup(j12,TBME_index_v2(Nmax,lmax,l,jindex,lpr,jindexpr,n1,n2,n3,n4))) > 1.0e-20) THEN

                                  WRITE(2,'(9I2,F12.8)') j12,l,2*l+(-1)**jindex,lpr,2*lpr+(-1)**jindex,n1,n2,n3,n4,REAL((2*j12+1),kind=r_kind)/REAL(deg*deg_pr,kind=r_kind)*matels_jcoup(j12,TBME_index_v2(Nmax,lmax,l,jindex,lpr,jindexpr,n1,n2,n3,n4))
                                  !WRITE(2,*) 2*j12, lookup_mi(n1,l,jindex), lookup_mi(n2,lpr,jindexpr), lookup_mi(n3,l,jindex), lookup_mi(n4,lpr,jindexpr)
                                  !WRITE(2,*) (2*j12+1), deg, deg_pr
                                  !WRITE(2,*) matel_val !, matels_jcoup(j12,TBME_index_v2(Nmax,lmax,l,jindex,lpr,jindexpr,n1,n2,n3,n4))
                                  
                               END IF
                               !end for debug

                            END DO
                            Ho_two_body_matels(TBME_index(l,jindex,lpr,jindexpr,n1,n2,n3,n4)) =COMPLEX(matel_val,0.0_r_kind) 

!!$                            !for debug
!!$                            IF(ABS(matel_val) > 1.0e-20) THEN
!!$                               WRITE(*,*) l,jindex,lpr,jindexpr,n1,n2,n3,n4,Ho_two_body_matels(TBME_index(l,jindex,lpr,jindexpr,n1,n2,n3,n4))
!!$                            END IF
!!$                            !end for debug
                            
                         END DO
                      END DO
                   END DO
                END DO
             END DO
          END DO
       END DO
    END DO

   CLOSE(2)


    
    DEALLOCATE(lookup_mi,lookup_nlj)

  END SUBROUTINE read_mortens_matels


  SUBROUTINE read_mortens_matels_v2(Nmax,lmax,filename_sp_in,filename_matels_in)
    IMPLICIT NONE

    INTEGER, intent(in) :: Nmax, lmax
    CHARACTER(100), optional :: filename_sp_in, filename_matels_in
    CHARACTER(100) :: filename_sp, filename_matels

    CHARACTER(100) :: slask
    INTEGER :: II,morten_index, n, l, j2, tz, Nmaj, jindex
    INTEGER :: ios, par, a,b,c,d, line, a_temp, c_temp
    REAL(kind=r_kind) :: e_ho, matel_val
    INTEGER :: n1,n2,n3,n4,l1,l2,l3,l4,jindex1,jindex2,jindex3,jindex4,j2max
    INTEGER :: lpr,jindexpr,deg,deg_pr, j12, j12_2, Ho_size_all_morten, Ho_size_TBME_morten
    INTEGER :: j2_1,j2_2, length


    INTEGER, allocatable :: lookup_mi(:,:,:)

    REAL(kind=r_kind), allocatable :: matels_jcoup(:,:,:,:,:)
    REAL(kind=r_kind) :: norm_L, norm_R
    INTEGER :: sign_ket_perm, sign_bra_perm, ket_perm, bra_perm
    
    TYPE nlj
       INTEGER :: n,l,jindex, j2
    END type nlj
    
    TYPE(nlj), allocatable :: lookup_nlj(:)

    IF(.not. present(filename_sp_in)) THEN
       filename_sp= 'spJ.dat'
    ELSE
       filename_sp = filename_sp_in
    END IF

    IF(.not. present(filename_matels_in)) THEN
       filename_matels= 'VJ-scheme.dat'
    ELSE
       filename_matels = filename_matels_in
    END IF


    !construct lookup table

    OPEN(unit=1,file=TRIM(filename_sp), status="old")

    !first ten lines are legend etc
    DO II = 1,10
       READ(1,*)
    END DO

    

    Ho_size_all_morten = 0

    DO l=0,lmax
        DO jindex = 0,1
           IF(l==0 .and. jindex==1) CYCLE

           Ho_size_all_morten = Ho_size_all_morten + (Nmax-l)/2 + 1

        END DO
     END DO

     WRITE(*,*) 'Ho_size_all_morten = ',Ho_size_all_morten


    ALLOCATE(lookup_mi(0:Nmax/2,0:lmax,0:1))
    
    !TODO fix this, morten uses som funny triangular truncation
    !ALLOCATE(lookup_nlj(2*Ho_size_all_morten))

    ALLOCATE(lookup_nlj(4*Ho_size_all_morten))

    


    WRITE(*,*) 'Ho_size_all = ',Ho_size_all

    DO II=1,2*Ho_size_all_morten
       READ(1,*) slask, morten_index, n, l, j2, tz, Nmaj, e_ho
       !WRITE(*,*) morten_index,n,l,j2,tz,Nmaj,e_ho
       
       IF(ABS(e_ho - Ho_hbaromega*(Nmaj + 1.5_r_kind))>1e-10) THEN
          WRITE(*,*) 'Oscillator paramters in sp file:',TRIM(filename_sp),' are inconsistent.'
          STOP
       END IF       

       IF(tz == -1) THEN !assuming this means neutron
          IF(j2-2*l == 1) jindex = 0
          IF(j2-2*l == -1) jindex = 1
          lookup_mi(n,l,jindex) = morten_index
          lookup_nlj(morten_index)%n = n
          lookup_nlj(morten_index)%l = l
          lookup_nlj(morten_index)%jindex = jindex
          lookup_nlj(morten_index)%j2 = j2

          !WRITE(*,*) 'n,l,jindex = ',n,l,jindex
          !WRITE(*,*) 'lookup_mi(n,l,jindex) = ',lookup_mi(n,l,jindex)
       END IF      

    END DO  
       
    CLOSE(1)

    !read matels
    WRITE(*,*) 'Reading J scheme matrix elements from file: ',TRIM(filename_matels)

    OPEN(unit=1,file=TRIM(filename_matels), status="old")

    j2max = 2*lmax + 1
    Ho_size_TBME_morten = (lmax+1)*2*(lmax+1)*2*(Nmax/2+1)**4


    WRITE(*,*) 'Basis for matrix elements read'
    WRITE(*,*) 'Nmax = ',Nmax,'lmax = ',lmax
    WRITE(*,*) 'size_TBME = ',HO_size_TBME_morten
    WRITE(*,*) 'j2max = ',j2max
    WRITE(*,*) 'Basis used in hf'
    WRITE(*,*) 'Ho_Nmax = ',Ho_nmax,'Ho_lmax = ',Ho_lmax
    WRITE(*,*) 'Ho_size_TBME = ',Ho_size_TBME

    
    !TODO fix this
    length = MIN(2*SIZE(lookup_nlj),4*Ho_size_all_morten)
    !TODO

    WRITE(*,*) 'length = ',length

    ALLOCATE(matels_jcoup(0:j2max,length,length,length,length))
    matels_jcoup = 0.0_r_kind

    OPEN(unit=3,file="matels_jcoup.dat")

    
    line = 1

    DO 

       READ(1,*,iostat=ios) tz, par, j12_2, a, b, c, d, matel_val

       !WRITE(*,*) tz, par, j2, a, b, c, d, matel_val

       IF(tz == -1) THEN

          Norm_R = 1.0_r_kind
          Norm_L = 1.0_r_kind

          IF(a/=b) Norm_L = 1.0_r_kind/sqrt(2.0_r_kind)
          IF(c/=d) Norm_R = 1.0_r_kind/sqrt(2.0_r_kind)

          DO ket_perm = 1,2

             
             IF(ket_perm == 2) THEN
                j2_1 = lookup_nlj(a)%j2
                j2_2 = lookup_nlj(b)%j2
                sign_ket_perm = (-1)**((j12_2-j2_1-j2_2)/2)
                a_temp = a
                a = b
                b = a_temp
             ELSE
                sign_ket_perm = 1
             END IF
             DO bra_perm = 1,2

                IF(bra_perm == 2) THEN
                   j2_1 = lookup_nlj(c)%j2
                   j2_2 = lookup_nlj(d)%j2
                   sign_bra_perm = (-1)**((j12_2-j2_1-j2_2)/2)
                   c_temp = c
                   c = d
                   d = c_temp
                ELSE
                   sign_bra_perm = 1
                END IF



                matels_jcoup(j12_2/2,a,b,c,d) = Norm_L*Norm_R*matel_val
                matels_jcoup(j12_2/2,c,d,a,b) = Norm_L*Norm_R*matel_val

                WRITE(3,'(2I3)') ket_perm, bra_perm
                WRITE(3,'(5I3,F14.8)') j12_2/2,a,b,c,d,sign_ket_perm*sign_bra_perm*matel_val 

             END DO
          END DO

       END IF


       IF(ios<0) EXIT
       IF(ios>0) THEN
          WRITE(*,*) 'Errror reading file:',TRIM(filename_matels)
          STOP
       END IF

       line = line +1 

    END DO


    CLOSE(1)

    CLOSE(3)
    WRITE(*,*) 'Done'

    WRITE(*,*) 'Converting matrix elements'
    OPEN(unit=2, file="test_read_morten.txt")


    
   

    DO l=0,MIN(Ho_lmax,lmax)
       DO jindex = 0,1
          IF(l==0 .and. jindex==1) CYCLE
          deg = Ho_degeneracy(l,jindex)
          DO lpr = 0,MIN(Ho_lmax,lmax)
             DO jindexpr = 0,1
                IF(lpr==0 .and. jindexpr==1) CYCLE
                deg_pr = Ho_degeneracy(lpr,jindexpr)

                DO n1 = 0,(MIN(Ho_Nmax,Nmax)-l)/2
                   DO n2 = 0,(MIN(Ho_Nmax,Nmax)-lpr)/2
                      DO n3 = 0,(MIN(Ho_Nmax,Nmax)-l)/2
                         DO n4 = 0,(MIN(Ho_Nmax,Nmax)-lpr)/2 

                            matel_val = 0.0_r_kind
                            
                            j2max = (2*l+(-1)**jindex + 2*lpr+(-1)**jindexpr)/2
                            !for test
                            !j2max = 0
                            !end for test

                            DO j12 = 0,j2max
                               a = lookup_mi(n1,l,jindex)
                               b = lookup_mi(n2,lpr,jindexpr)
                               c = lookup_mi(n3,l,jindex)
                               d = lookup_mi(n4,lpr,jindexpr)

                               j2_1 = 2*l+(-1)**jindex
                               j2_2 = 2*lpr+(-1)**jindexpr


                               matel_val = matel_val + REAL((2*j12+1),kind=r_kind)/REAL(deg*deg_pr,kind=r_kind) * (matels_jcoup(j12,a,b,c,d) -  (-1)**(j12-(j2_1+j2_2)/2)*matels_jcoup(j12,a,b,d,c))

                               WRITE(2,'(5I3,F14.10)') j12,a,b,c,d,matels_jcoup(j12,a,b,c,d)
                               WRITE(2,'(5I3,F14.10)') j12,a,b,c,d,matels_jcoup(j12,a,b,d,c)
                               WRITE(2,*) 
                               !matel_val = REAL((2*j12+1),kind=r_kind)/REAL(deg*deg_pr,kind=r_kind)*(matels_jcoup(j12,a,b,c,d)+matels_jcoup(j12,a,b,d,c))

                               !for debug                              
                               !end for debug

                            END DO
                            Ho_two_body_matels(TBME_index(l,jindex,lpr,jindexpr,n1,n2,n3,n4)) =COMPLEX(matel_val,0.0_r_kind) 

!!$                            !for debug
!!$                            IF(ABS(matel_val) > 1.0e-20) THEN
!!$                               WRITE(*,*) l,jindex,lpr,jindexpr,n1,n2,n3,n4,Ho_two_body_matels(TBME_index(l,jindex,lpr,jindexpr,n1,n2,n3,n4))
!!$                            END IF
!!$                            !end for debug
                            
                         END DO
                      END DO
                   END DO
                END DO
             END DO
          END DO
       END DO
    END DO

   CLOSE(2)


    
    DEALLOCATE(lookup_mi,lookup_nlj)
    WRITE(*,*) 'Done'


  END SUBROUTINE read_mortens_matels_v2


  SUBROUTINE calculate_matels_full

    USE talmi, only: TMBR, INIT_TALMI
    USE geometric, only : sixj

    IMPLICIT NONE

    INTEGER :: l12,l1,l2,n1,n2,E1,E2,Etot,N_CM,L_CM,n_rel,l_rel,E_CM,E_rel  
    INTEGER :: l3,l4,n3,n4,E3,E4,Etot34,n_rel34,E_rel34
    INTEGER :: jindex1,j1_2, jindex2, j2_2, jindex3, j3_2, jindex4, j4_2
    INTEGER :: n_rel12,E_rel12,Etot12

    INTEGER :: deg1,l12max, II, term, JJ

    REAL(kind=r_kind), ALLOCATABLE :: T(:), matels_j(:), R(:,:,:)
    REAL(kind=r_kind), ALLOCATABLE :: P1(:), P2(:), scaled_gp(:)

    REAL(kind=r_kind) :: T12,T34,RadInt,geom12,geom34,matel_val,r_val

    REAL(kind=r_kind) :: kappa, V0, r_test
    REAL(kind=r_kind), parameter :: Minnesota_kappas(2) =  (/ 1.487_r_kind, 0.465_r_kind /) !(/ 1.0e-15, 1.0e-15 /) !
    REAL(kind=r_kind), parameter :: Minnesota_Vs(2) = (/ 200.0_r_kind, -91.85_r_kind /) ! (/0.5 , 0.5 /)!


    WRITE(*,*) 'Calcualting matrix elements'


    WRITE(*,*) '   Initializing Talmi module'
    CALL INIT_TALMI
    WRITE(*,*) '   Done'

    WRITE(*,*) 'Precalucalting Talmi-Moshinsky brackets'

    WRITE(*,*) 'talmi_index_max = ',talmi_index(2*Ho_Nmax+2,Ho_Nmax/2,Ho_Nmax,Ho_Nmax/2,Ho_Nmax,Ho_Nmax,2*Ho_Nmax,Ho_Nmax,2*Ho_nmax)

    ALLOCATE( T(talmi_index(2*Ho_Nmax+1,Ho_Nmax/2,Ho_Nmax,Ho_Nmax/2,Ho_Nmax,Ho_Nmax,2*Ho_Nmax,Ho_Nmax,2*Ho_nmax)) )
    T = 0.0_r_kind

    DO l12 = 0, 2*Ho_Nmax
       DO l1 = 0, Ho_Nmax
          DO l2 = 0, Ho_Nmax
             DO n1 = 0, (Ho_Nmax - l1)/2
                DO n2 = 0, (Ho_Nmax - l2)/2
                   E1 = 2*n1 + l1
                   E2 = 2*n2 + l2
                   Etot = E1 + E2


                   DO l_rel = 0, Etot, 2
                      DO L_CM = 0, Etot - l_rel
                         DO N_CM = 0, (Etot - L_CM - l_rel)/2
                            DO n_rel = 0, (Etot - 2*N_CM - L_CM - l_rel)/2

                               E_CM = 2*N_CM + L_CM
                               E_rel = 2*n_rel + l_rel

                               T(talmi_index(l12,n1,l1,n2,l2,N_CM,L_CM,n_rel,l_rel)) = TMBR(E_CM,L_CM,E_rel,l_rel,E1,l1,E2,l2,l12)
                            END DO
                         END DO
                      END DO
                   END DO
                END DO
             END DO
          END DO
       END DO
    END DO

    WRITE(*,*) '   Done'


    WRITE(*,*) 'Precalcalcualting radial overlaps'


    

    ALLOCATE(R(0:2*Ho_Nmax,0:Ho_Nmax,0:Ho_Nmax))

IF(.true.) THEN

    WRITE(*,*) 'Using Gauss-Hermite quadrature'
    CALL init_grid_GH(50) !larger than 70 causes errors
    
    ALLOCATE(P1(grid_size_GH),P2(grid_size_GH),scaled_gp(grid_size_GH))

    R = 0.0_r_kind

    DO term = 1,2

       kappa = Minnesota_kappas(term)
       V0 = Minnesota_Vs(term)

       DO l_rel = 0, 2*Ho_Nmax
          DO n_rel12 = 0, (2*Ho_Nmax - l_rel)/2
             DO n_rel34 = 0, (2*Ho_Nmax - l_rel)/2

                scaled_gp = grid_points_GH/sqrt(1.0_r_kind+Ho_b**2*2.0_r_kind*kappa)
                CALL RadHO_poly(n_rel12,l_rel,1.0_r_kind,scaled_gp,P1,grid_size_GH)
                CALL RadHO_poly(n_rel34,l_rel,1.0_r_kind,scaled_gp,P2,grid_size_GH)

                r_val = 0.0_r_kind
                DO II=1,grid_size_GH
                   r_val = r_val + grid_weights_GH(II)*grid_points_GH(II)**2*P1(II)*P2(II)
                END DO
                R(l_rel,n_rel12,n_rel34) = R(l_rel,n_rel12,n_rel34) + V0*r_val/(1.0_r_kind + Ho_b**2*2.0_r_kind*kappa)**1.5_r_kind
             END DO
          END DO
       END DO
    END DO


    DEALLOCATE(P1,P2,scaled_gp)

    WRITE(*,*) '    Done'

ELSE


    WRITE(*,*) 'Using Gauss-Legendre quadrature'

    CALL init_grid_GL(400)

    ALLOCATE(P1(grid_size_GL),P2(grid_size_GL),scaled_gp(grid_size_GL))

    scaled_gp = tan(pi/4.0_r_kind*(grid_points_GL+1.0_r_kind))

    R = 0.0_r_kind

    DO term = 1,2

       kappa = Minnesota_kappas(term)
       V0 = Minnesota_Vs(term)

       DO l_rel = 0, 2*Ho_Nmax
          DO n_rel12 = 0, (2*Ho_Nmax - l_rel)/2
             DO n_rel34 = 0, (2*Ho_Nmax - l_rel)/2

                CALL RadHO(n_rel12,l_rel,Ho_b,scaled_gp,P1,grid_size_GL)
                CALL RadHO(n_rel34,l_rel,Ho_b,scaled_gp,P2,grid_size_GL)

                r_val = 0.0_r_kind
                DO II=1,grid_size_GL
                   r_val = r_val + grid_weights_GL(II)*(scaled_gp(II)**2+1.0_r_kind)*scaled_gp(II)**2*P1(II)*P2(II)*V0*EXP(-2.0_r_kind*kappa*scaled_gp(II)**2)
                END DO
                R(l_rel,n_rel12,n_rel34) = R(l_rel,n_rel12,n_rel34) + pi/4.0_r_kind*r_val
             END DO
          END DO
       END DO
    END DO

    WRITE(*,*) 'R(0,0,0)= ', R(0,0,0) 

    !test
    
    CALL RadHO(0,0,Ho_b,scaled_gp,P1,grid_size_GL)
    CALL RadHO(0,0,Ho_b,scaled_gp,P2,grid_size_GL)
    r_test = 0.0
    DO term = 1,2

       kappa = Minnesota_kappas(term)
       V0 = Minnesota_Vs(term)


       r_val = 0.0_r_kind
       DO II = 1,grid_size_GL
          DO JJ = 1,grid_size_GL

             r_val = r_val + grid_weights_GL(II)*(scaled_gp(II)**2+1.0_r_kind)*grid_weights_GL(JJ)*(scaled_gp(JJ)**2+1.0_r_kind)&
                  *scaled_gp(II)*P1(II)**2*scaled_gp(JJ)*P2(JJ)**2*(EXP(-kappa*(scaled_gp(II)**2 + scaled_gp(JJ)**2) + 2.0*kappa*scaled_gp(II)*scaled_gp(JJ)) - EXP(-kappa*(scaled_gp(II)**2 + scaled_gp(JJ)**2) -2.0*kappa*scaled_gp(II)*scaled_gp(JJ)))/(4.0*kappa)

          END DO
       END DO
       r_test = r_test + V0*pi**2/16.0*r_val

    END DO

    WRITE(*,*) 'r_test = ', r_test

    !end test


    DEALLOCATE(P1,P2,scaled_gp)

    WRITE(*,*) '    Done'

END IF


    WRITE(*,*) 'Calculating jj-coupled matrix elments'


    WRITE(*,*) 'TBME_j_index_max = ', TBME_j_index(2*Ho_Nmax + 1, Ho_Nmax/2 + 1, Ho_Nmax + 1, 2, Ho_Nmax/2 + 1, Ho_Nmax + 1, 2, Ho_Nmax/2 + 1, Ho_Nmax + 1, 2, Ho_Nmax/2 + 1, Ho_Nmax + 1, 2) 

    ALLOCATE(matels_j( TBME_j_index(2*Ho_Nmax + 1, Ho_Nmax/2 + 1, Ho_Nmax + 1, 2, Ho_Nmax/2 + 1, Ho_Nmax + 1, 2, Ho_Nmax/2 + 1, Ho_Nmax + 1, 2, Ho_Nmax/2 + 1, Ho_Nmax + 1, 2)))

    matels_j = 0.0_r_kind

    !The following calcualates the two body matrix elements of the Minnesota in the most general case
    !For the HF only a small subset of the matrix elements are used

    DO l12 = 0, 2*Ho_Nmax 
       DO l1 = 0, Ho_Nmax
          DO l2 = 0, Ho_Nmax
             DO n1 = 0, (Ho_Nmax - l1)/2
                DO n2 = 0, (Ho_Nmax - l2)/2
                   E1 = 2*n1 + l1
                   E2 = 2*n2 + l2
                   Etot12 = E1 + E2

                   DO l_rel = 0, Etot12, 2
                      DO L_CM = 0, Etot12 - l_rel
                         DO N_CM = 0, (Etot12 - L_CM - l_rel)/2
                            DO n_rel12 = 0, (Etot12 - 2*N_CM - L_CM - l_rel)/2

                               E_CM = 2*N_CM + L_CM
                               E_rel12 = 2*n_rel12 + l_rel

                               T12 = TMBR(E_CM,L_CM,E_rel12,l_rel,E1,l1,E2,l2,l12)
                               !T(talmi_index(l12,n1,l1,n2,l2,N_CM,L_CM,n_rel12,l_rel))

                               DO l3 = 0, Ho_Nmax
                                  DO l4 = 0, Ho_Nmax
                                     DO n3 = 0, (Ho_Nmax - l3)/2
                                        DO n4 = 0, (Ho_Nmax - l4)/2
                                           E3 = 2*n3 + l3
                                           E4 = 2*n4 + l4
                                           Etot34 = E3 + E4

                                           DO n_rel34 = 0, (Etot34 - 2*N_CM - L_CM - l_rel)/2
                                              E_rel34 = 2*n_rel34 + l_rel
                                              T34 = TMBR(E_CM,L_CM,E_rel34,l_rel,E3,l3,E4,l4,l12) 
                                              !T(talmi_index(l12,n3,l3,n4,l4,N_CM,L_CM,n_rel34,l_rel))

                                              
                                              RadInt = R(l_rel,n_rel12,n_rel34)
                                              

                                              DO jindex1 = 0,1
                                                 IF(l1==0 .and. jindex1==1) CYCLE
                                                 j1_2 = twoj(l1,jindex1)
                                                 DO jindex2 = 0,1
                                                    IF(l2==0 .and. jindex2==1) CYCLE
                                                    j2_2 = twoj(l2,jindex2)

                                                    geom12 = sqrt(REAL(j1_2+1,kind=r_kind)*REAL(j2_2+1,kind=r_kind))&
                                                         *sixj(2*l2,2*l1,2*l12,j1_2,j2_2,1)

                                                    DO jindex3 = 0,1
                                                       IF(l3==0 .and. jindex3==1) CYCLE
                                                       j3_2 = twoj(l3,jindex3)
                                                       DO jindex4 = 0,1
                                                          IF(l4==0 .and. jindex4==1) CYCLE
                                                          j4_2 = twoj(l4,jindex4)

                                                          geom34 = sqrt(REAL(j3_2+1,kind=r_kind)*REAL(j4_2+1,kind=r_kind))&
                                                               *sixj(2*l4,2*l3,2*l12,j3_2,j4_2,1)

                                                          matels_j(TBME_j_index(l12,n1,l1,jindex1,n2,l2,jindex2,n3,l3,jindex3,n4,l4,jindex4))= &
                                                                matels_j(TBME_j_index(l12,n1,l1,jindex1,n2,l2,jindex2,n3,l3,jindex3,n4,l4,jindex4))&
                                                                +REAL((-1)**(l1+l3+(j2_2+j4_2)/2+1),kind=r_kind)*geom12*geom34*T12*T34*RadInt
                                                                

                                                          !for debug
!!$                                                          IF( ABS(matels_j(TBME_j_index(l12,n1,l1,jindex1,n2,l2,jindex2,n3,l3,jindex3,n4,l4,jindex4))) > 1.0e-10 ) THEN 
!!$                                                          WRITE(*,'(3I3,A1,3I3,A1,3I3,A1,3I3,A1,I3)') n1,l1,jindex1,'|',n2,l2,jindex2,'|',n3,l3,jindex3,'|',n4,l4,jindex4,'|',l12
!!$                                                          WRITE(*,'(2F16.10)') matels_j(TBME_j_index(l12,n1,l1,jindex1,n2,l2,jindex2,n3,l3,jindex3,n4,l4,jindex4)), RadInt
!!$                                                          WRITE(*,'(4F12.6)') geom12, geom34, T12, T34
!!$                                                          
!!$                                                          END IF

                                                          IF(n1 == 0 .and. n2 == 0 .and. n3 == 0 .and. n4 == 0 .and. l1 == 0 .and. l2 == 0 .and. l3 == 0 .and. l4 == 0) THEN
                                                             WRITE(*,'(3I3,A1,3I3,A1,3I3,A1,3I3,A1,I3)') n1,l1,jindex1,'|',n2,l2,jindex2,'|',n3,l3,jindex3,'|',n4,l4,jindex4,'|',l12
                                                          WRITE(*,'(2F16.10)') matels_j(TBME_j_index(l12,n1,l1,jindex1,n2,l2,jindex2,n3,l3,jindex3,n4,l4,jindex4)), RadInt
                                                          WRITE(*,'(4F12.6)') geom12, geom34, T12, T34

                                                          END IF


                                                       END DO
                                                    END DO
                                                 END DO
                                              END DO
                                           END DO
                                        END DO
                                     END DO
                                  END DO
                               END DO
                            END DO
                         END DO
                      END DO
                   END DO
                END DO
             END DO
          END DO
       END DO
    END DO

    WRITE(*,*) '   Done'

    WRITE(*,*) 'Calcualting matrix elements needed for HF'


    !calculates sum_m2 <n1 l1 j1 m1, n2 l2 j2 m2 | v | n3 l1 j1 m1, n4 l2 j2 m2> 

    DO l1=0,Ho_lmax
       DO jindex1 = 0,1
          IF(l1==0 .and. jindex1==1) CYCLE
          deg1 = Ho_degeneracy(l1,jindex1)
          DO l2 = 0,Ho_lmax
             DO jindex2 = 0,1
                IF(l2==0 .and. jindex2==1) CYCLE               

                DO n1 = 0,(Ho_Nmax-l1)/2
                   DO n2 = 0,(Ho_Nmax-l2)/2
                      DO n3 = 0,(Ho_Nmax-l1)/2
                         DO n4 = 0,(Ho_Nmax-l2)/2 

                            matel_val = 0.0_r_kind
                            
                            l12max = l1 + l2
                            

                            DO l12 = 0,l12max
                               !j2_1 = twoj(l1,jindex1)!2*l1+(-1)**jindex1
                               !j2_2 = 2*lpr+(-1)**jindexpr


                               matel_val = matel_val + REAL((2*l12+1),kind=r_kind)/REAL(deg1,kind=r_kind) *matels_j(TBME_j_index(l12,n1,l1,jindex1,n2,l2,jindex2,n3,l1,jindex1,n4,l2,jindex2))                            
                               

                            END DO
                            !to conform to the definition in the talent notes
                            matel_val = matel_val / REAL(twoj(l2,jindex2)+1,kind=r_kind)
                            !

                            Ho_two_body_matels(TBME_index(l1,jindex1,l2,jindex2,n1,n2,n3,n4)) =COMPLEX(matel_val,0.0_r_kind) 

!!$                            !for debug
!!$                            IF(ABS(matel_val) > 1.0e-20) THEN
!!$                               WRITE(*,*) l,jindex,lpr,jindexpr,n1,n2,n3,n4,Ho_two_body_matels(TBME_index(l,jindex,lpr,jindexpr,n1,n2,n3,n4))
!!$                            END IF
!!$                            !end for debug
                            
                         END DO
                      END DO
                   END DO
                END DO
             END DO
          END DO
       END DO
    END DO

    WRITE(*,*) '   Done'




    DEALLOCATE(T)
    DEALLOCATE(matels_j)


    WRITE(*,*) 'Finished calculating matrix elements'

  END SUBROUTINE calculate_matels_full
  
  INTEGER FUNCTION talmi_index(l12,n1,l1,n2,l2,N_CM,L_CM,n_rel,l_rel)
    IMPLICIT NONE
    INTEGER, intent(in) :: n1,l1,n2,l2,N_CM,L_CM,n_rel,l_rel,l12

    INTEGER :: n_lab_dim, l_lab_dim, n_rel_dim, l_rel_dim !, l12_dim

    n_lab_dim = Ho_Nmax/2 + 1
    l_lab_dim = Ho_Nmax + 1
    n_rel_dim = Ho_Nmax + 1
    l_rel_dim = 2*Ho_Nmax + 1
    !l12_dim = 2*Ho_Nmax + 2
   
    talmi_index = 1 &
         + l12*l_lab_dim**2*n_lab_dim**2*l_rel_dim**2*n_rel_dim**2 &
         + (l_lab_dim*l1 + l2)*n_lab_dim**2*l_rel_dim**2*n_rel_dim**2 & 
         + (n_lab_dim*n1 + n2)*l_rel_dim**2*n_rel_dim**2 &
         + (l_rel_dim*L_CM + l_rel)*n_rel_dim**2 &
         + n_rel_dim*N_CM +  n_rel

  END FUNCTION talmi_index

  INTEGER FUNCTION TBME_j_index(j12,n1,l1,jindex1,n2,l2,jindex2,n3,l3,jindex3,n4,l4,jindex4)
    IMPLICIT NONE
    INTEGER, intent(in) :: j12,n1,l1,jindex1,n2,l2,jindex2,n3,l3,jindex3,n4,l4,jindex4

    INTEGER :: dim_j12, dim_n, dim_l, dim_ji

    dim_j12 = 2*Ho_Nmax + 3
    dim_n = Ho_Nmax/2 +1
    dim_l = Ho_Nmax + 1
    dim_ji = 2

    TBME_j_index = 1 + j12*dim_ji**4*dim_l**4*dim_n**4 &
    + (dim_ji*jindex3 + jindex4)*dim_ji**2*dim_l**4*dim_n**4 &
    + (dim_ji*jindex1 + jindex2)*dim_l**4*dim_n**4 &
    + (dim_l*l3 + l4)*dim_l**2*dim_n**4 &
    + (dim_l*l1 + l2)*dim_n**4 &
    + (dim_n*n1 + n2)*dim_n**2 &
    + dim_n*n3  + n4    
    
  END FUNCTION TBME_j_index




  !as fortran does not support matrices with more than 7 dimensions this function
  !returns the index in the collapsed 8-dim array to where the matrix element is stored
  INTEGER FUNCTION TBME_index(l,jindex,lpr,jindexpr,n1,n2,n3,n4)
    IMPLICIT NONE
    INTEGER, intent(in) :: l,jindex,lpr,jindexpr,n1,n2,n3,n4

    INTEGER :: size_block, size_block_2, size_block_4

    TBME_index = 1

    !n1,n3 belong to l,j
    !n2,n4 belong to lpr,jpr

    size_block = Ho_Nmax/2+1
    size_block_2 = size_block**2
    size_block_4 = size_block_2**2

    TBME_index = TBME_index + l*2*(Ho_lmax+1)*2*size_block_4 &
         + jindex*(Ho_lmax+1)*2*size_block_4 &
         + lpr*2*size_block_4 + jindexpr*size_block_4 &
         + (n1 + size_block*n3)*size_block_2 &
         + n2 + size_block*n4
         
    


  END FUNCTION TBME_index


  INTEGER FUNCTION TBME_index_v2(Nmax,lmax,l,jindex,lpr,jindexpr,n1,n2,n3,n4)
    IMPLICIT NONE
    INTEGER, intent(in) :: Nmax,lmax,l,jindex,lpr,jindexpr,n1,n2,n3,n4

    INTEGER :: size_block, size_block_2, size_block_4

    TBME_index_v2 = 1

    !n1,n3 belong to l,j
    !n2,n4 belong to lpr,jpr

    size_block = Nmax/2+1
    size_block_2 = size_block**2
    size_block_4 = size_block_2**2

    TBME_index_v2 = TBME_index_v2 + l*2*(lmax+1)*2*size_block_4 &
         + jindex*(lmax+1)*2*size_block_4 &
         + lpr*2*size_block_4 + jindexpr*size_block_4 &
         + (n1 + size_block*n3)*size_block_2 &
         + n2 + size_block*n4
         
    


  END FUNCTION TBME_index_v2



  SUBROUTINE hf_calculate_delta(delta)
    IMPLICIT NONE

    REAL(kind=r_kind), intent(out) :: delta
    
    INTEGER :: II
    REAL(kind=r_kind) :: e, e_old

    delta = 0.0_r_kind

     DO II = 0,Ho_size_all-1
       e = hf_energies_all(II)
       e_old = hf_energies_all_old(II)

       delta = delta + ABS(e-e_old)/hf_num_part_request

    END DO


  END SUBROUTINE hf_calculate_delta
  

  INTEGER FUNCTION twoj(l,jindex)
    IMPLICIT NONE
    INTEGER, intent(in) :: l,jindex
    
    IF(jindex == 0) twoj = 2*l + 1
    IF(jindex == 1) twoj = 2*l - 1
    
  END FUNCTION twoj

  

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

  SUBROUTINE hf_find_fermi_energy(N_request,e_fermi)

    USE sorting, only : sortrx

    IMPLICIT NONE

    INTEGER, intent(in) :: N_request !number of particles
    REAL(kind=r_kind), intent(out) :: e_fermi
    INTEGER :: particle_number, II, index !, fermi_index_tot
    INTEGER :: n,l,jindex
    REAL(kind=r_kind) :: e


    hf_energies_all_old = hf_energies_all


    !stores information on the eigenvectors
    II = 0

    DO l=0,Ho_lmax
       DO jindex = 0,1
          IF(l==0 .and. jindex==1) CYCLE

          DO n=0,(Ho_Nmax-l)/2
             hf_states_all(II)%l = l
             hf_states_all(II)%jindex = jindex
             hf_states_all(II)%n = n
             hf_states_all(II)%e = hf_energies(l,jindex,n)

             hf_energies_all(II) = hf_energies(l,jindex,n)
             II = II + 1

          END DO
       END DO

    END DO

    !The following assumes that the diagonalization sorts according to increasing energies in each lj block

    !sorts the energies
    CALL sortrx(Ho_size_all,hf_energies_all,hf_energies_all_sort)
    !indexing starts from 0
    hf_energies_all_sort = hf_energies_all_sort -1


    !finds the fermi indicies

    hf_states_all(:)%occ = 0.0_r_kind
    
!    fermi_index_hfbasis(:,:) = 0
        fermi_index_hfbasis(:,:) = -1

   particle_number = 0
    !fermi_index_tot = 0
    II = 0
    DO WHILE(particle_number < N_request)

       index = hf_energies_all_sort(II)
       
       n = hf_states_all(index)%n
       l = hf_states_all(index)%l
       jindex = hf_states_all(index)%jindex
       e = hf_states_all(index)%e
       
       particle_number = particle_number + Ho_degeneracy(l,jindex)
       hf_states_all(index)%occ = REAL(Ho_degeneracy(l,jindex),kind=r_kind)

       fermi_index_hfbasis(l,jindex) = n

       e_fermi = e
       !fermi_index_tot = fermi_index_tot + 1
       II = II + 1

       !WRITE(*,*) 'Find fermi: l = ',l,'n= ',n,'jindex = ',jindex,'fermi_index = ',n 


    END DO

    IF(N_request-particle_number /= 0) THEN
       WRITE(*,*) 'ERROR: Cannot get correct particle number by filling degenerate shells'
       WRITE(*,*) 'Requested particle number',N,'Actual particle number',particle_number 
       STOP
    END IF

      


  END SUBROUTINE hf_find_fermi_energy

  SUBROUTINE print_hf_states
    IMPLICIT NONE

    INTEGER :: II, n, l , jindex, j2, index
    REAL(kind=r_kind) :: e, occ

    
    WRITE(*,*) 'Ho_b = ',Ho_b, 'Ho_hbaromega = ',Ho_hbaromega
    WRITE(*,*) 'Ho_lmax = ',Ho_lmax,'Ho_nmax = ',Ho_nmax
    WRITE(*,*) 'Ho_size_all = ',Ho_size_all


    WRITE(*,'(3A3,2A14)') 'n','l','2j','e','occ'
    DO II = 0,Ho_size_all-1
       index = hf_energies_all_sort(II)

       n = hf_states_all(index)%n
       l = hf_states_all(index)%l
       jindex = hf_states_all(index)%jindex
       j2 = l + (-1)**jindex
       e = hf_states_all(index)%e
       occ = hf_states_all(index)%occ

       WRITE(*,'(3I3,2F14.6)') n,l,j2,e,occ


    END DO
    

  END SUBROUTINE print_hf_states


  SUBROUTINE hf_init_dens(N_request)
    IMPLICIT NONE

    INTEGER, intent(in) :: N_request
    INTEGER :: n, l,jindex, particle_number
    INTEGER :: N_major, N_major_max
    

    INTEGER :: n1,n2

    hf_num_part_request = N_request


    particle_number = 0
    
    !N_major_max = 2*Ho_nmax + Ho_lmax !changed def of Ho_Nmax

    density_matrix = (0.0_r_kind,0.0_r_kind)

     
    WRITE(*,*) 'hf_init_dens: filling states'
    WRITE(*,'(A3,A3,A3)') 'n','l','2j'

    Ho_fill_shells :  DO l = 0,Ho_lmax
       DO jindex = 0,1
          IF(l==0 .and. jindex==1) CYCLE

          DO n = 0,(Ho_Nmax - l)/2

             !IF(n > Ho_nmax .or. l > Ho_lmax) CYCLE !unnencecary with changed def of Nmax

             density_matrix(l,jindex,n,n) = (1.0_r_kind,0.0_r_kind)
             particle_number = particle_number + Ho_degeneracy(l,jindex)

             WRITE(*,'(I3,I3,I3)') n,l,l+(-1)**jindex


             IF(particle_number >= N_request) THEN
                EXIT Ho_fill_shells
             END IF
          END DO

       END DO

    END DO Ho_fill_shells

    WRITE(*,*) 'hf_init_dens:'
    WRITE(*,*) 'Requested particle number',N_request
    WRITE(*,*) 'Particle number',particle_number

!!$
!!$    DO l=0,Ho_lmax
!!$       DO jindex = 0,1
!!$          IF(l==0 .and. jindex==1) CYCLE
!!$
!!$
!!$          DO n1 = 0,Ho_nmax
!!$             DO n2 = n1,Ho_nmax
!!$
!!$                WRITE(*,*) 'n1 = ',n1,'n2=',n2,'density_matrix(l,jindex,n1,n2)=',density_matrix(l,jindex,n1,n2)
!!$
!!$
!!$             END DO
!!$          END DO
!!$       END DO
!!$    END DO
    
    

  END SUBROUTINE hf_init_dens



  !rho(nljm,n'l'j'm') = rho(lj;nn')d_ll'd_jj'd_mm'
  !the reduced density rho(lj;nn') is stored, keep track of degeneracy!
  SUBROUTINE hf_update_density_matrix(mixing)
    IMPLICIT NONE

    REAL(kind=r_kind), intent(in) :: mixing
    
    INTEGER :: l,jindex,n1,n2,II

    COMPLEX(kind=r_kind), parameter :: zero_c = (0.0_r_kind,0.0_r_kind)
    COMPLEX(kind=r_kind) :: sum_occupied    
    

    density_matrix_old = density_matrix
    

    DO l=0,Ho_lmax
       DO jindex = 0,1
          IF(l==0 .and. jindex==1) CYCLE

          
          DO n1 = 0,(Ho_Nmax-l)/2
             DO n2 = n1,(Ho_Nmax-l)/2
             !DO n2 = 0,Ho_nmax
             
                sum_occupied = zero_c
!!$                
!!$                WRITE(*,*) 'Update density: l=',l,'n1=',n1,'n2=',n2
!!$                WRITE(*,*) 'Update density: fermi_index =',fermi_index_hfbasis(l,jindex)  
!!$                WRITE(*,*) 'Update density: hf_transform(l,jindex,n1,n2) =',hf_transform(l,jindex,n1,n2) 

                DO II = 0,fermi_index_hfbasis(l,jindex)  !hf_transform should be sorted according to hf sp-energies in each block

                   sum_occupied = sum_occupied + hf_transform(l,jindex,n1,II)*CONJG(hf_transform(l,jindex,n2,II))

                END DO

!!$                WRITE(*,*) 'Update density: sum_occupied = ',sum_occupied
                
                density_matrix(l,jindex,n1,n2) = mixing*density_matrix_old(l,jindex,n1,n2) + (1.0_r_kind-mixing)*sum_occupied !Ho_degeneracy(l,jindex)*sum_occupied
                density_matrix(l,jindex,n2,n1) = CONJG(density_matrix(l,jindex,n1,n2))
                
             END DO
          END DO
       END DO
    END DO

    

  END SUBROUTINE hf_update_density_matrix


  SUBROUTINE hf_update_hamiltonian
    IMPLICIT NONE

    INTEGER :: l,jindex,n1,n2,n3,n4,lpr,jindexpr,deg
    COMPLEX(kind=r_kind), parameter :: zero_c = (0.0_r_kind,0.0_r_kind)
    COMPLEX(kind=r_kind) :: Gamma


    hf_hamiltonian = zero_c
    hf_gamma = zero_c

    DO l=0,Ho_lmax
       DO jindex = 0,1
          IF(l==0 .and. jindex==1) CYCLE

          DO n1 = 0,(Ho_Nmax-l)/2
             DO n2 = n1,(Ho_Nmax-l)/2
                !DO n2 = 0,Ho_nmax

                Gamma = zero_c
                DO lpr = 0,Ho_lmax
                   DO jindexpr = 0,1
                      IF(lpr==0 .and. jindexpr==1) CYCLE
                      deg = Ho_degeneracy(lpr,jindexpr)

                      DO n3 = 0,(Ho_Nmax-lpr)/2
                         DO n4 = 0,(Ho_Nmax-lpr)/2 

                            !for the l=0 subspace
                            !Gamma = Gamma + Ho_two_body_matel(l,jindex,n1,n4,n2,n3)*density_matrix(l,jindex,n3,n4) 
                            !General
                            Gamma = Gamma + deg*Ho_two_body_matels(TBME_index(l,jindex,lpr,jindexpr,n1,n4,n2,n3))*density_matrix(lpr,jindexpr,n3,n4) 

!!$                            !for debug
!!$                            WRITE(*,*) 'n1,n2,n3,n4 = ',n1,n2,n3,n4
!!$                            WRITE(*,*) TBME_index(l,jindex,lpr,jindexpr,n1,n2,n3,n4)
!!$                            WRITE(*,*) Ho_two_body_matel(l,jindex,n1,n2,n3,n4)
!!$                            WRITE(*,*) Ho_two_body_matels(TBME_index(l,jindex,lpr,jindexpr,n1,n2,n3,n4))
!!$                            !end for debug

                            !





                         END DO
                      END DO
                   END DO
                END DO

                !WRITE(*,*) 'Gamma = ',Gamma

                !
                hf_gamma(l,jindex,n1,n2) = Gamma
                hf_gamma(l,jindex,n2,n1) = CONJG(Gamma)
                !


                !hf_hamiltonian(l,jind,n1,n2) = Ho_one_body_matel(l,jind,n1,n2) + Gamma                
                hf_hamiltonian(l,jindex,n1,n2) = Gamma 
                hf_hamiltonian(l,jindex,n2,n1) = CONJG(Gamma)!CONJG(hf_hamiltonian(l,jindex,n1,n2))


             END DO

          END DO
       END DO
    END DO

    DO l=0,Ho_lmax
       DO jindex = 0,1
          IF(l==0 .and. jindex==1) CYCLE
          DO n1 = 0,(Ho_Nmax-l)/2
             hf_hamiltonian(l,jindex,n1,n1) = hf_hamiltonian(l,jindex,n1,n1) + Ho_hbaromega*(2*n1+l+1.5_r_kind)
          END DO
       END DO
    END DO

  END SUBROUTINE hf_update_hamiltonian


  SUBROUTINE hf_total_energy(Etot)
    IMPLICIT NONE

    REAL(kind=r_kind) :: Etot
    INTEGER :: l,jindex,n1,n2,j2    

    Etot = 0.0_r_kind

    !WRITE(*,*) 'hf_total_energy: Ho_lmax = ',Ho_lmax

    DO l=0,Ho_lmax
       DO jindex = 0,1
          IF(l==0 .and. jindex==1) CYCLE
          j2 = 2*l + (-1)**jindex
          DO n1 = 0,(Ho_Nmax-l)/2
             Etot = Etot + Ho_hbaromega*(2.0_r_kind*n1+l+1.5_r_kind)*(j2+1)*density_matrix(l,jindex,n1,n1)

             !WRITE(*,*) 'l',l,'n1=',n1,'density_matrix(l,jindex,n1,n1) = ',density_matrix(l,jindex,n1,n1)

             DO n2 = 0,(Ho_Nmax-l)/2

                Etot = Etot + 0.5_r_kind*(j2+1.0_r_kind)*density_matrix(l,jindex,n2,n1)*hf_gamma(l,jindex,n1,n2)


             END DO
          END DO
       END DO
    END DO


  END SUBROUTINE hf_total_energy

   SUBROUTINE hf_total_energy_v2(Etot)
    IMPLICIT NONE

    REAL(kind=r_kind) :: Etot
    INTEGER :: l,jindex,n1,n2,j2,n3,n4,lpr,jindexpr,j2pr    

    Etot = 0.0_r_kind
    
    DO l=0,Ho_lmax
       DO jindex = 0,1
          IF(l==0 .and. jindex==1) CYCLE
          j2 = 2*l + (-1)**jindex
          DO n1 = 0,(Ho_Nmax-l)/2
             Etot = Etot + Ho_hbaromega*(2.0_r_kind*n1+l+1.5_r_kind)*(j2+1)*density_matrix(l,jindex,n1,n1)
             !!$             WRITE(*,*) 'n1=',n1,'density_matrix(l,jindex,n1,n1) = ',density_matrix(l,jindex,n1,n1)
             DO n3 = 0,(Ho_Nmax-l)/2
                DO lpr=0,Ho_lmax
                   DO jindexpr = 0,1
                      IF(lpr==0 .and. jindexpr==1) CYCLE
                      j2pr = 2*lpr + (-1)**jindexpr


                      DO n2 = 0,(Ho_Nmax-lpr)/2
                         DO n4 = 0,(Ho_Nmax-lpr)/2
                            ! l = 0 subspace
                            ! Etot = Etot + 0.5_r_kind*(j2+1.0_r_kind)*density_matrix(l,jindex,n3,n1)*Ho_two_body_matel(l,jindex,n1,n2,n3,n4)*density_matrix(l,jindex,n4,n2)

                            !general
                            Etot = Etot + 0.5_r_kind*(j2+1)*(j2pr+1)&
                                 *density_matrix(l,jindex,n3,n1)&
                                 *Ho_two_body_matels(TBME_index(l,jindex,lpr,jindexpr,n1,n2,n3,n4))&
                                 *density_matrix(lpr,jindexpr,n4,n2)


                         END DO
                      END DO
                   END DO
                END DO
             END DO
          END DO
       END DO
    END DO



  END SUBROUTINE hf_total_energy_v2


  SUBROUTINE hf_diagonalize
    IMPLICIT NONE

    INTEGER :: l,jindex
    COMPLEX(kind=r_kind), ALLOCATABLE :: work(:)
    REAL(kind=r_kind), ALLOCATABLE :: rwork(:)
    COMPLEX(kind=r_kind), ALLOCATABLE :: h_block(:,:)
    REAL(kind=r_kind), ALLOCATABLE :: e_block(:)
    INTEGER :: lwork, info, size_block
    !INTEGER :: count

    INTERFACE 
       SUBROUTINE zheev(jobz, uplo, n, a, lda, w, work, lwork, rw, info)
         CHARACTER(len=1) ::  jobz, uplo
         INTEGER          ::  info, lda, lwork, n
         REAL(KIND=8)     ::  rw(*), w(*)
         COMPLEX(8)   ::  a(lda, *), work(*)
       END SUBROUTINE zheev
    END INTERFACE

    !count = 1

    DO l=0,Ho_lmax
       DO jindex = 0,1
          IF(l==0 .and. jindex==1) CYCLE

          size_block = (Ho_Nmax-l)/2+1

          !ALLOCATE(h_block(Ho_nmax+1,Ho_nmax+1))
          !ALLOCATE(e_block(Ho_nmax+1))

          ALLOCATE(h_block(size_block,size_block),e_block(size_block))

          h_block(:,:) = hf_hamiltonian(l,jindex,0:size_block-1,0:size_block-1)
          

                ALLOCATE(work(1))
                ALLOCATE(rwork( max(1,3*size_block-2) ) )
                lwork = -1

                CALL zheev('V','U',size_block,h_block ,size_block,e_block,work,lwork,rwork,info)

                lwork = work(1)
                DEALLOCATE(work)
                ALLOCATE(work(lwork))

                CALL zheev('V','U',size_block, h_block,size_block,e_block,work,lwork,rwork,info)

                DEALLOCATE(work)
                DEALLOCATE(rwork)

                hf_transform(l,jindex,0:size_block-1,0:size_block-1) = h_block(:,:)
                hf_energies(l,jindex,0:size_block-1) = e_block(:)
                DEALLOCATE(h_block,e_block)

!!$                !debug
!!$                WRITE(*,*) 'In hf_diagonalize, the hf transformation matrix:'
!!$                CALL print_matrix(Ho_nmax+1,REAL(hf_transform(l,jindex,0:Ho_nmax,0:Ho_nmax),kind=r_kind))
!!$                !
                        
       END DO

       !count = count + 1
    END DO

  END SUBROUTINE hf_diagonalize


  SUBROUTINE init_Ho_basis(Nmax,lmax,b,hbaromega)
    IMPLICIT NONE

    REAL(kind=r_kind), intent(in) :: b
    REAL(kind=r_kind), optional :: hbaromega
    INTEGER, intent(in) :: nmax, lmax

    WRITE(*,*) 'Initializing ho basis'
   
    IF(lmax > Nmax) THEN
       WRITE(*,*) "In init_ho_basis: lmax > Nmax"
       STOP
    END IF
    
    Ho_lmax = lmax
    Ho_Nmax = Nmax
 
    IF(.not. present(hbaromega)) THEN

       IF(b<=0) THEN
          WRITE(*,*) 'b<=0 not allowed'
          STOP
       END IF
       Ho_b = b
       Ho_hbaromega = hbarc**2/(Ho_b**2*mnc2)

    END IF

    IF(present(hbaromega)) THEN
        IF(hbaromega<=0) THEN
          WRITE(*,*) 'hbaromega<=0 not allowed'
          STOP
       END IF

       Ho_hbaromega  = hbaromega
       Ho_b = hbarc/sqrt(hbaromega*mnc2)
    END IF
       
    WRITE(*,*) 'IN init_Ho_basis'
    WRITE(*,*) 'mnc2 = ',mnc2
    WRITE(*,*) 'hbarc = ',hbarc
    WRITE(*,*) 'Ho_b = ',Ho_b
    WRITE(*,*) 'Ho_hbaromega = ',Ho_hbaromega
    WRITE(*,*) 'init_Ho_basis, done'



  END SUBROUTINE init_Ho_basis

  SUBROUTINE hf_init
    IMPLICIT NONE

    INTEGER :: l,jindex,n1,n2,n3,n4,max_n

    COMPLEX(kind=r_kind) :: zero_c = (0.0_r_kind,0.0_r_kind)
    REAL(kind=r_kind) :: zero_r = 0.0_r_kind
    
    REAL(kind=r_kind), parameter :: kappa_R = 1.487_r_kind
    REAL(kind=r_kind), parameter :: V_0R = 200.00_r_kind
    REAL(kind=r_kind), parameter :: kappa_S = 0.465_r_kind
    REAL(kind=r_kind), parameter :: V_0S = -91.85_r_kind

    CHARACTER(100) :: file1,file2


    WRITE(*,*) 'Allocating matricies for HF'
    WRITE(*,*) 'Ho_lmax = ',Ho_lmax
    WRITE(*,*) 'Ho_nmax = ',Ho_nmax

    max_n = Ho_nmax/2

    ALLOCATE(hf_transform(0:Ho_lmax,0:1,0:max_n,0:max_n))
    hf_transform = zero_c
    ALLOCATE(density_matrix(0:Ho_lmax,0:1,0:max_n,0:max_n))
    density_matrix = zero_c
    ALLOCATE(density_matrix_old(0:Ho_lmax,0:1,0:max_n,0:max_n))
    density_matrix_old = zero_c
    ALLOCATE(hf_energies(0:Ho_lmax,0:1,0:max_n))
    hf_energies = zero_r
    ALLOCATE(hf_energies_old(0:Ho_lmax,0:1,0:max_n))
    hf_energies_old = zero_r

    Ho_size_all = 0

    DO l=0,Ho_lmax
        DO jindex = 0,1
           IF(l==0 .and. jindex==1) CYCLE

           Ho_size_all = Ho_size_all + (Ho_Nmax-l)/2 + 1

        END DO
     END DO

     ALLOCATE(hf_states_all(0:Ho_size_all-1)) 
     ALLOCATE(hf_energies_all(0:Ho_size_all-1))
     hf_energies_all = zero_r
     ALLOCATE(hf_energies_all_old(0:Ho_size_all-1))
     hf_energies_all_old = zero_r

     ALLOCATE(hf_energies_all_sort(0:Ho_size_all-1))

     ALLOCATE(fermi_index_hfbasis(0:Ho_lmax,0:1))
     ALLOCATE(hf_hamiltonian(0:Ho_lmax,0:1,0:max_n,0:max_n))
     hf_hamiltonian = zero_c
     ALLOCATE(hf_gamma(0:Ho_lmax,0:1,0:max_n,0:max_n))
     hf_gamma = zero_c
     ALLOCATE(Ho_two_body_matel(0:Ho_lmax,0:1,0:max_n,0:max_n,0:max_n,0:max_n))
     Ho_two_body_matel = zero_c



     Ho_size_TBME = (Ho_lmax+1)*2*(Ho_lmax+1)*2*(max_n+1)**4

     ALLOCATE(Ho_two_body_matels(Ho_size_TBME))
     Ho_two_body_matels = zero_c

    

     IF(.false.) THEN     
        

        !CALL calculate_two_body_matels_GL(kappa_R,0.5_r_kind*V_0R)
        !CALL calculate_two_body_matels_GL(kappa_S,0.5_r_kind*V_0S)

        !should give 0,1,and 2 based on orthogonality of radial wfs
        CALL calculate_two_body_matels_GL(1e-15_r_kind,1.0_r_kind)


!!$     Ho_two_body_matel = (0.0_r_kind,0.0_r_kind)
!!$      
!!$
!!$     CALL calculate_two_body_matels(kappa_R,0.5_r_kind*V_0R)
!!$     CALL calculate_two_body_matels(kappa_S,0.5_r_kind*V_0S)
!!$
!!$     !should give 0,1,and 2 based on orthogonality of radial wfs
!!$     !CALL calculate_two_body_matels(1e-5_r_kind,1.0_r_kind)


        !DO l=0,Ho_lmax
        !DO jindex = 0,1
        !IF(l==0 .and. jindex==1) CYCLE
        l = 0
        jindex = 0
        DO n1 = 0,Ho_nmax/2
           DO n2 = 0,Ho_nmax/2
              DO n3 = 0,Ho_nmax/2
                 DO n4 = 0,Ho_nmax/2


                    Ho_two_body_matels(TBME_index(l,jindex,l,jindex,n1,n2,n3,n4)) = 0.5_r_kind*Ho_two_body_matel(l,jindex,n1,n2,n3,n4)


                 END DO
              END DO
           END DO
        END DO

        !END DO
        !END DO
     END IF

     IF(.false.) THEN
!!$        file1 = 'spJ.swave.dat'
!!$        file2 = 'VJ-scheme.swave.dat'
!!$        CALL read_mortens_matels_v2(8, 0,file1,file2 )


        file1 = 'spJ.Nmax3.dat'
        file2 = 'VJ-scheme.Nmax3.dat'
        !Nmax lmax
        CALL read_mortens_matels_v2(3, 3,file1,file2 )
        !CALL read_mortens_matels(2,2)
     END IF


     CALL  calculate_matels_full


     CALL print_TBMEs

     !for test
     !Ho_two_body_matels = zero_c
     !end for test

     



  END SUBROUTINE hf_init


  SUBROUTINE print_TBMEs
    IMPLICIT NONE

    INTEGER :: l,jindex,lpr,jindexpr
    INTEGER :: n1,n2,n3,n4

     OPEN(unit=1,file='matels_hf_v2.dat')

     WRITE(1,*) 'Two-body matrix elements'

     DO l=0,Ho_lmax
        DO jindex = 0,1
           IF(l==0 .and. jindex==1) CYCLE

           DO lpr = 0,Ho_lmax
              DO jindexpr = 0,1
                 IF(lpr==0 .and. jindexpr==1) CYCLE

                 DO n1 = 0,(Ho_Nmax-l)/2
                    DO n2 = 0,(Ho_Nmax-lpr)/2
                       DO n3 = 0,(Ho_Nmax-l)/2
                          DO n4 = 0,(Ho_Nmax-lpr)/2                             

                             WRITE(1,'(8I2,2F20.12)') l,jindex,lpr,jindexpr,n1,n2,n3,n4,REAL(Ho_two_body_matels(TBME_index(l,jindex,lpr,jindexpr,n1,n2,n3,n4)),kind=r_kind),AIMAG(Ho_two_body_matels(TBME_index(l,jindex,lpr,jindexpr,n1,n2,n3,n4)))

                          END DO
                       END DO
                    END DO
                 END DO
              END DO
           END DO
        END DO
     END DO

     CLOSE(1)


  END SUBROUTINE print_TBMEs


   SUBROUTINE calculate_two_body_matels(mu,V_0)
    IMPLICIT NONE

    REAL(kind=r_kind), intent(in) :: mu, V_0

    INTEGER :: n
    REAL(kind=r_kind) :: alpha

    REAL(kind=r_kind), ALLOCATABLE :: wf_1(:), wf_3(:)
    REAL(kind=r_kind), ALLOCATABLE :: wf_2(:), wf_4(:)
    REAL(kind=r_kind), ALLOCATABLE :: scaled_gp(:)
    INTEGER :: II, JJ, dim_dg
    INTEGER :: l, jindex

    INTEGER :: n1,n2,n3,n4
    COMPLEX(kind=r_kind) :: matel

    

    INTEGER :: term
    INTEGER, parameter :: no_terms = 1
    REAL(kind=r_kind) :: prefactor(no_terms)
    
    
    REAL(kind=r_kind) :: xII,xJJ,pf,norm

    

    n = 100
    alpha = 0.0_r_kind

    WRITE(*,*) 'Calculating matix elements'

    WRITE(*,*) 'Initializing grid'
    CALL init_grid_GLag(n,alpha)
    WRITE(*,*) '  Done'

    ALLOCATE(wf_1(n),wf_2(n),wf_3(n),wf_4(n),scaled_gp(n))

    !only s-wave implemented
    l = 0
    jindex = 0

    !scaling   
    prefactor = (/ 1.0_r_kind/(Ho_b**2*16.0_r_kind*mu*(Ho_b**2*mu+1.0_r_kind)**2) /)
    prefactor = V_0*prefactor
    
    
    WRITE(*,*) 'Perfoming sums'

    WRITE(*,*) 'mu = ',mu,'V_0 = ',V_0,'Ho_b = ', Ho_b

    scaled_gp = (Ho_b**2*mu+1.0_r_kind)**(-0.5_r_kind)*sqrt(grid_points_GLag)
  
    DO term = 1,no_terms
       pf = prefactor(term)
       DO n1 = 0,Ho_Nmax/2
          DO n2 = 0,Ho_Nmax/2
             DO n3 = 0,Ho_Nmax/2
                DO n4 = 0,Ho_Nmax/2

!!$                   !for debug
!!$                   WRITE(*,*) 'n1=',n1,'n2=',n2,'n3=',n3,'n4=',n4
!!$                   
!!$
!!$                    norm = 0.0_r_kind
!!$                    CALL RadHO_poly(n1,0,1.0_r_kind,sqrt(grid_points_GLag),wf_1,n)
!!$                    CALL RadHO_poly(n2,0,1.0_r_kind,sqrt(grid_points_GLag),wf_2,n)
!!$                     DO II = 1,n
!!$                        xII = sqrt(grid_points_GLag(II))
!!$                        norm = norm + grid_weights_GLag(II)*xII*wf_1(II)*wf_2(II)                       
!!$
!!$                     END DO
!!$                     norm = norm/2.0
!!$                     !norm = Ho_b**3/2.0*norm
!!$                     WRITE(*,*) 'norm n1,n2 =',norm
!!$
!!$                     !end for debug
                   
                   CALL RadHO_poly(n1,0,1.0,scaled_gp,wf_1,n)
                   CALL RadHO_poly(n3,0,1.0,scaled_gp,wf_3,n)
                   CALL RadHO_poly(n2,0,1.0,scaled_gp,wf_2,n)
                   CALL RadHO_poly(n4,0,1.0,scaled_gp,wf_4,n)

!!$                   wf_1 = 0.0
!!$                   wf_2 = 0.0
!!$                   wf_3 = 0.0
!!$                   wf_4 = 0.0
                   
                   matel = (0.0_r_kind,0.0_r_kind) 


                  

                   DO II = 1,n
                      DO JJ = 1,n
                         xII = scaled_gp(II)
                         xJJ = scaled_gp(JJ)
                         
                         !WRITE(*,*) 'II = ',II,'JJ = ',JJ,'xII = ',xII,'xJJ = ',xJJ 
                         
!!$                         matel = matel + grid_weights_GLag(II)*grid_weights_GLag(JJ)&
!!$                              *(wf_1(II)*wf_3(II)*wf_2(JJ)*wf_4(JJ) + wf_1(II)*wf_4(II)*wf_2(JJ)*wf_3(JJ))&
!!$                              *(EXP(2.0_r_kind*Ho_b**2*mu*xII*xJJ)-EXP(-2.0_r_kind*Ho_b**2*mu*xII*xJJ))   

                         matel = matel + EXP(LOG(grid_weights_GLag(II)*grid_weights_GLag(JJ)) + 2.0_r_kind*Ho_b**2*mu*xII*xJJ)*(wf_1(II)*wf_3(II)*wf_2(JJ)*wf_4(JJ) + wf_1(II)*wf_4(II)*wf_2(JJ)*wf_3(JJ))  - grid_weights_GLag(II)*grid_weights_GLag(JJ)&
                              *(wf_1(II)*wf_3(II)*wf_2(JJ)*wf_4(JJ) + wf_1(II)*wf_4(II)*wf_2(JJ)*wf_3(JJ))&
                              *EXP(-2.0_r_kind*Ho_b**2*mu*xII*xJJ)


                       

                         
!!$                         !
!!$                         WRITE(*,*) xII,grid_weights_GLag(II),xJJ,grid_weights_GLag(JJ),2.0_r_kind*Ho_b**2*mu*xII*xJJ
!!$                         WRITE(*,*) matel
!!$                         WRITE(*,*) pf*matel
!!$                         !
                         

                      END DO
                   END DO
                   
!!$                   !for debug
!!$                   WRITE(*,*) 'matel =',matel,'Ho_two_body_matel(l,jindex,n1,n2,n3,n4)=',Ho_two_body_matel(l,jindex,n1,n2,n3,n4)
!!$                   STOP
!!$                   !end for debug
                   
                   Ho_two_body_matel(l,jindex,n1,n2,n3,n4) = Ho_two_body_matel(l,jindex,n1,n2,n3,n4) +  pf*matel
                   !Ho_two_body_matel(l,jindex,n1,n2,n3,n4) =  pf*matel
                   
                END DO
             END DO
          END DO
       END DO
    END DO

    WRITE(*,*) 'Done'

    WRITE(*,*) 'Two-body matrix elements'
    
    DO n1 = 0,Ho_Nmax/2
       DO n2 = 0,Ho_Nmax/2
          DO n3 = 0,Ho_Nmax/2
             DO n4 = 0,Ho_Nmax/2
                
                WRITE(*,'(4I2,2F20.14)') n1,n2,n3,n4,REAL(Ho_two_body_matel(l,jindex,n1,n2,n3,n4),kind=r_kind),AIMAG(Ho_two_body_matel(l,jindex,n1,n2,n3,n4))

             END DO
          END DO
       END DO
    END DO



  END SUBROUTINE calculate_two_body_matels



  




   SUBROUTINE calculate_two_body_matels_GL(mu,V_0)
    IMPLICIT NONE

    REAL(kind=r_kind), intent(in) :: mu, V_0

    INTEGER :: n
    

    REAL(kind=r_kind), ALLOCATABLE :: wf_1(:), wf_3(:)
    REAL(kind=r_kind), ALLOCATABLE :: wf_2(:), wf_4(:)
    REAL(kind=r_kind), ALLOCATABLE :: scaled_gp(:)
    INTEGER :: II, JJ
    INTEGER :: l, jindex

    INTEGER :: n1,n2,n3,n4
    COMPLEX(kind=r_kind) :: matel

    

    INTEGER :: term
    INTEGER, parameter :: no_terms = 1
    REAL(kind=r_kind) :: prefactor(no_terms)
    
    
    REAL(kind=r_kind) :: xII,xJJ,pf,norm

    

    n = 100
    
    
    WRITE(*,*) 'Calculating matix elements'

    IF(.not. is_init_grid_GL .or. grid_size_GL /= n) THEN
       WRITE(*,*) 'Initializing grid'
       CALL init_grid_GL(n)
       WRITE(*,*) '  Done'
    END IF

    ALLOCATE(wf_1(n),wf_2(n),wf_3(n),wf_4(n),scaled_gp(n))

    !only s-wave implemented
    l = 0
    jindex = 0

    !scaling   
    prefactor = (/ pi**2/16.0_r_kind/(4.0_r_kind*mu) /)
    prefactor = V_0*prefactor
    
    
    WRITE(*,*) 'Perfoming sums'

    WRITE(*,*) 'mu = ',mu,'V_0 = ',V_0,'Ho_b = ', Ho_b

    scaled_gp = tan(pi/4.0_r_kind*(grid_points_GL+1.0_r_kind))
  
    DO term = 1,no_terms
       pf = prefactor(term)
       DO n1 = 0,Ho_Nmax/2
          DO n2 = 0,Ho_Nmax/2
             DO n3 = 0,Ho_Nmax/2
                DO n4 = 0,Ho_Nmax/2
!!$
!!$                   !for debug
!!$                   WRITE(*,*) 'n1=',n1,'n2=',n2,'n3=',n3,'n4=',n4
!!$                   
!!$
!!$                    norm = 0.0_r_kind
!!$                    CALL RadHO(n1,0,Ho_b,scaled_gp,wf_1,n)
!!$                    CALL RadHO(n2,0,Ho_b,scaled_gp,wf_2,n)
!!$                     DO II = 1,n
!!$                        xII = scaled_gp(II) ! tan(pi*grid_points_GL(II)/2.0_r_kind)
!!$                        norm = norm + grid_weights_GL(II)*(xII**2+1.0_r_kind)&
!!$                             *xII**2*wf_1(II)*wf_2(II)       
!!$
!!$                        !norm = norm + grid_weights_GL(II)*(xII**2+1.0_r_kind)*EXP(-xII)
!!$
!!$                        !xJJ = pi/4.0_r_kind*(grid_points_GL(II)+1.0_r_kind)
!!$                        !norm = norm + grid_weights_GL(II)/(cos(xJJ)**2)*EXP(-xII)
!!$
!!$                     END DO
!!$                     norm = pi*norm/4.0_r_kind
!!$                     !norm = Ho_b**3/2.0*norm
!!$                     WRITE(*,*) 'norm n1,n2 =',norm
!!$
!!$                     !end for debug
                   
                   CALL RadHO(n1,0,Ho_b,scaled_gp,wf_1,n)
                   CALL RadHO(n3,0,Ho_b,scaled_gp,wf_3,n)
                   CALL RadHO(n2,0,Ho_b,scaled_gp,wf_2,n)
                   CALL RadHO(n4,0,Ho_b,scaled_gp,wf_4,n)

!!$                   wf_1 = 0.0
!!$                   wf_2 = 0.0
!!$                   wf_3 = 0.0
!!$                   wf_4 = 0.0
                   
                   matel = (0.0_r_kind,0.0_r_kind) 


                  

                   DO II = 1,n
                      DO JJ = 1,n
                         xII = scaled_gp(II)
                         xJJ = scaled_gp(JJ)



                         matel = matel + grid_weights_GL(II)*grid_weights_GL(JJ)*(xII**2+1.0_r_kind)*(xJJ**2+1.0_r_kind)&
                              *(wf_1(II)*wf_3(II)*wf_2(JJ)*wf_4(JJ) + wf_1(II)*wf_4(II)*wf_2(JJ)*wf_3(JJ))&
                              *xII*xJJ*(EXP(-mu*(xII**2+xJJ**2)+2.0_r_kind*mu*xII*xJJ)-EXP(-mu*(xII**2+xJJ**2)-2.0_r_kind*mu*xII*xJJ))




                         !EXP(-mu*(xII**2+xJJ**2))*(EXP(2*mu*xII*xJJ)-EXP(-2*mu*xII*xJJ))



                         
!!$                         !
!!$                         WRITE(*,*) xII,grid_weights_GLag(II),xJJ,grid_weights_GLag(JJ),2.0_r_kind*Ho_b**2*mu*xII*xJJ
!!$                         WRITE(*,*) matel
!!$                         WRITE(*,*) pf*matel
!!$                         !
                         

                      END DO
                   END DO
                   
!!$                   !for debug
!!$                   WRITE(*,*) 'matel =',matel,'Ho_two_body_matel(l,jindex,n1,n2,n3,n4)=',Ho_two_body_matel(l,jindex,n1,n2,n3,n4)
!!$                   STOP
!!$                   !end for debug
                   
                   Ho_two_body_matel(l,jindex,n1,n2,n3,n4) = Ho_two_body_matel(l,jindex,n1,n2,n3,n4) +  pf*matel
                   !Ho_two_body_matel(l,jindex,n1,n2,n3,n4) =  pf*matel
                   
                END DO
             END DO
          END DO
       END DO
    END DO

    WRITE(*,*) 'Done'

!!$    OPEN(unit=1,file='matels_hf.dat')
!!$
!!$    WRITE(1,*) 'Two-body matrix elements'
!!$    
!!$    DO n1 = 0,Ho_Nmax/2!MIN(1,Ho_nmax)
!!$       DO n2 = 0,Ho_Nmax/2!MIN(1,Ho_nmax)
!!$          DO n3 = 0,Ho_Nmax/2!MIN(1,Ho_nmax)
!!$             DO n4 = 0,Ho_Nmax/2!MIN(1,Ho_nmax)
!!$                
!!$                WRITE(1,'(4I2,2F20.12)') n1,n2,n3,n4,REAL(Ho_two_body_matel(l,jindex,n1,n2,n3,n4),kind=r_kind),AIMAG(Ho_two_body_matel(l,jindex,n1,n2,n3,n4))
!!$
!!$             END DO
!!$          END DO
!!$       END DO
!!$    END DO
!!$
!!$    CLOSE(1)


  END SUBROUTINE calculate_two_body_matels_GL











  !the following was a failed attempt at calculating the two-body matrix elements
  !where the integration variables where changed. Might work with some more thought
  SUBROUTINE calculate_two_body_matels_incorrect(mu,V_0)
    IMPLICIT NONE

    REAL(kind=r_kind), intent(in) :: mu, V_0

    INTEGER :: n
    REAL(kind=r_kind) :: alpha

    REAL(kind=r_kind), ALLOCATABLE :: double_grid_1(:), wf_1(:), wf_3(:)
    REAL(kind=r_kind), ALLOCATABLE :: double_grid_2(:), wf_2(:), wf_4(:)
    INTEGER :: II, JJ, dim_dg
    INTEGER :: l, jindex

    INTEGER :: n1,n2,n3,n4
    COMPLEX(kind=r_kind) :: matel

    

    INTEGER :: term
    INTEGER, parameter :: no_terms = 2
    REAL(kind=r_kind) :: scaling_1(no_terms), scaling_2(no_terms), prefactor(no_terms)
    
    
    REAL(kind=r_kind) :: xII,xJJ,pf,s1,s2

    

    n = 100
    dim_dg = n**2
    alpha = -0.5_r_kind

    WRITE(*,*) 'Calculating matix elements'

    WRITE(*,*) 'Initializing grid'
    CALL init_grid_GLag(n,alpha)
    WRITE(*,*) '  Done'

    ALLOCATE(double_grid_1(n**2),double_grid_2(n**2),wf_1(n**2),wf_2(n**2),wf_3(n**2),wf_4(n**2))

    DO II = 1,n
       DO JJ = 1,n
          
          double_grid_1((II-1)*n + JJ) = grid_points_GLag(II)+0.5_r_kind*grid_points_GLag(JJ)
          double_grid_2((II-1)*n + JJ) = grid_points_GLag(II)-0.5_r_kind*grid_points_GLag(JJ)

       END DO
    END DO

  

    !only s-wave implemented
    l = 0
    jindex = 0

    !scaling
    scaling_1 = (/ 2.0_r_kind**(-0.5_r_kind),(4.0_r_kind*mu*Ho_b**2+2.0_r_kind)**(-0.5_r_kind)/)
    scaling_2 = (/ (mu*Ho_b**2+0.5_r_kind)**(-0.5_r_kind),2.0_r_kind**(-0.5_r_kind) /)
    prefactor = (/ Ho_b**4/(16.0_r_kind*mu*sqrt(2.0_r_kind*mu*Ho_b**2+1.0_r_kind)) &
         ,  -1.0_r_kind*Ho_b**4/(16.0_r_kind*mu*sqrt(8.0_r_kind*mu*Ho_b**2+4.0_r_kind)) /)
    prefactor = V_0*prefactor
    
    
    WRITE(*,*) 'Perfoming sums'

  
    DO term = 1,no_terms
       s1 = scaling_1(term)
       s2 = scaling_2(term)
       pf = prefactor(term)
       DO n1 = 0,Ho_Nmax/2
          DO n2 = 0,Ho_Nmax/2
             DO n3 = 0,Ho_Nmax/2
                DO n4 = 0,Ho_Nmax/2

                   WRITE(*,*) 'n1=',n1,'n2=',n2,'n3=',n3,'n4=',n4

                   CALL RadHO_poly(n1,0,1.0,s1*double_grid_1**(0.5_r_kind),wf_1,dim_dg)
                   CALL RadHO_poly(n3,0,1.0,s1*double_grid_1**(0.5_r_kind),wf_3,dim_dg)
                   CALL RadHO_poly(n2,0,1.0,s2*double_grid_2**(0.5_r_kind),wf_2,dim_dg)
                   CALL RadHO_poly(n4,0,1.0,s2*double_grid_2**(0.5_r_kind),wf_4,dim_dg)

!!$                   wf_1 = 0.0
!!$                   wf_2 = 0.0
!!$                   wf_3 = 0.0
!!$                   wf_4 = 0.0
                   
                   matel = (0.0_r_kind,0.0_r_kind) 

                   DO II = 1,n
                      DO JJ = 1,n
                         xII = s1*grid_weights_GLag(II)**(0.5_r_kind)
                         xJJ = s2*grid_weights_GLag(JJ)**(0.5_r_kind)
                         matel = matel + grid_weights_GLag(II)*grid_weights_GLag(JJ)&
                              *(xII**2 - xJJ**2/4.0_r_kind)&
                              *wf_1((II-1)*n + JJ)*wf_3((II-1)*n + JJ)&
                              *wf_2((II-1)*n + JJ)*wf_4((II-1)*n + JJ)

                      END DO
                   END DO
                   Ho_two_body_matel(l,jindex,n1,n2,n3,n4) = Ho_two_body_matel(l,jindex,n1,n2,n3,n4) +  pf*matel

                END DO
             END DO
          END DO
       END DO
    END DO

    WRITE(*,*) 'Done'




  END SUBROUTINE calculate_two_body_matels_incorrect



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

    WRITE(*,*) 'Ho_b = ',Ho_b
    WRITE(*,*) 'Ho_hbaromega = ',Ho_hbaromega

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
    REAL(kind=r_kind) :: RVEC(dimr), FVEC(dimr),FVEC_S(dimr),bin,alpha
    REAL(kind=r_kind) :: L0(dimr), L1(dimr)

    IF(n.lt.0) THEN
       FVEC   = 0.0_r_kind
       FVEC_S = 0.0_r_kind
       RETURN
    END IF

    alpha = l + 0.5_r_kind
    
    L1 = 0.0_r_kind; FVEC = 1.0_r_kind;
    DO ii = 1 , n 
       L0 = L1
       L1 = FVEC
       FVEC = ((2.0_r_kind * ii- 1.0_r_kind + alpha- RVEC)* L1- (ii- 1.0_r_kind + alpha)* L0)/ REAL(ii,kind=r_kind)
    END DO
  
    FVEC_S = sign(1.0_r_kind,FVEC)
    !Commented away statement below I want the value not the log /Daniel W
    !FVEC   = log(abs(FVEC)) 
    !End Comment /Daniel W
    RETURN 
  END Subroutine LaguerreL2

  !Calculates the polynomial part of the radial HO-wf
  Subroutine RadHO_poly(n, l, b, RVEC, FVEC, dimr)
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

   FVEC = SQRT(2d0/b**3)* EXP(0.5d0* ( lnfac - gammln(nR+lR+1.5d0) ) )* (RVEC/b)**l* FVEC 
   
 END Subroutine RadHO_poly


! Calculates a vector of values of g_nl(r), where g_nl(r) is the radial part of a HO wave function with size parameter b as defined
! on page 49 of Jouni Sohonen, From Nucleons to Nucleus, Springer
  Subroutine RadHO(n, l, b, RVEC, FVEC, dimr)
    IMPLICIT NONE
    INTEGER :: n,l,dimr
    REAL(kind=r_kind) :: nR,lR,b, RVEC(dimr), FVEC(dimr), FVEC_S(dimr), FVECtmp(dimr), lnfac
    nR = REAL(n,kind=r_kind) 
    lR = REAL(l,kind=r_kind)
    CALL LaguerreL2(n, l, (RVEC/b)**2, FVEC, FVEC_S, dimr)
     
    ! gamma(n+1) = n!
    IF (n == 0) THEN
      lnfac = 0.0_r_kind
    ELSE
      lnfac = gammln(nR+1.0_r_kind)
    END IF

   FVEC = SQRT(2.0_r_kind/b**3)* EXP(0.5_r_kind* ( lnfac - gammln(nR+lR+1.5_r_kind) ) )* (RVEC/b)**l* EXP(-RVEC**2/2.0_r_kind/b**2)* FVEC 
   
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


  SUBROUTINE print_matrix(size,matrix)

    USE types

    IMPLICIT NONE

    INTEGER, intent(in) :: size
    REAL(kind=r_kind) :: matrix(size,size)

    INTEGER :: II, JJ

    DO II = 1,size
       DO JJ = 1,size
          WRITE(*,'(E18.8)',advance = "no") matrix(II,JJ)
       END DO
       WRITE(*,*)
    END DO



  END SUBROUTINE PRINT_MATRIX

    

END MODULE
