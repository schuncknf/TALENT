program rpa

!     NRMAX  - maximum number of points in spatial grid
!     NHMAX  - maximum number of hole wave functions
!     MATMAX - maximum dimensionality of matrix
  integer, parameter :: nrmax = 50
  integer, parameter :: nhmax = 50
  integer, parameter :: matmax = 40

! mathematical constants
  real(kind = 8), parameter :: half = 0.5
  real(kind = 8), parameter :: third =1.0/3.0
  real(kind = 8), parameter :: pi = 3.1415926535
  real(kind = 8), parameter :: pi4 = 4.0*pi

! global variables

! DEL - mesh size for Sch. equation and Green's function
! NGRID - total number of mesh points
  integer :: ngrid
  real(kind = 8) :: del

! mesh for the response function
  integer :: n
  real(kind = 8) :: del2

! interaction parameters
  real(kind = 8) :: t0
  real(kind = 8) :: t3
  real(kind = 8) :: x
  real(kind = 8) :: vscal

! quantum numbers, energies and eigenvectors
  integer, dimension(nhmax) :: lh
  integer, dimension(nhmax) :: jh
  integer, dimension(nhmax) :: nq
  real(kind = 8) :: nmass,npro ! nuclear mass and proton number
  real(kind = 8), dimension(nhmax) :: eh
  real(kind = 8), dimension(nhmax,nrmax) :: uh

! potentials, central, spin-orbit and Coulomb
  real(kind = 8), dimension(nrmax) :: v
  real(kind = 8), dimension(nrmax) :: vs
  real(kind = 8), dimension(nrmax) :: vc  

! residual interaction
  real(kind = 8), dimension(matmax) :: vnn
  real(kind = 8), dimension(matmax) :: vnp
  real(kind = 8), dimension(matmax) :: vext

! densities
  real(kind = 8), dimension(nrmax) :: rhop
  real(kind = 8), dimension(nrmax) :: rhon  

! nuclear response
  integer :: l ! multipolarity
  real(kind = 8) :: ex ! starting energy  
  real(kind = 8) :: exm ! maximum energy
  real(kind = 8) :: dex ! delta Ex
  real(kind = 8) :: gam ! gamma, width of the transitions
  real(kind = 8) :: iso ! isospin of the external field

! RPA variables
  complex(kind = 8), dimension(matmax,matmax) :: a
  complex(kind = 8), dimension(matmax,matmax) :: b
  complex(kind = 8), dimension(matmax,matmax) :: g
  complex(kind = 8), dimension(nrmax) :: up,up1
  complex(kind = 8), dimension(nrmax) :: wp,wp1
  complex(kind = 8), dimension(matmax) :: drhof
  complex(kind = 8), dimension(matmax) :: drhrpa
  complex(kind = 8) :: sfree,srpa

  real(kind = 8), dimension(130) :: faclog

  integer :: nh,n2,nx
  real(kind = 8) :: h,fl

  integer :: ll
  integer :: i2,m1,m2,jp
  integer :: jmin,jmax,lp
  integer :: ma,mb,mx
  real(kind = 8) :: temp,sumf,sumrpa,value

  call factor

  call static
  print *,' # mass = ',nmass,'z = ',npro
  read(*,*) n,del2
  write(*,'(''# ''i5,'' points in response mesh, with mesh size='',f4.2)') n,del2
  nx = del2/del

  call vresid(vnn,vnp,nmass)

!  Multipolarity,starting energy, ending energy, and energy step for the response function
  read(*,*) l,ex,exm,dex,gam
  ll = 2*l
  print*,'# l= ',l,' gam= ', gam

!  make multipole external field
  call extf(vext,l)

! loop over excitation energy
  do
    if(ex.gt.exm) exit
    n2 = 2*n
    g(:,:) = 0.0

    do i2 = 1,nh
      jmin = abs(ll - jh(i2))
      jmax = ll + jh(i2)
      do jp = jmin,jmax,2
        if(mod((jp + 1)/2 + lh(i2) + l,2).ne.0) then
          lp = (jp - 1)/2
        else
          lp = (jp + 1)/2
        end if
        call clebsh(jh(i2),jp,ll,1,-1,0,temp)

        temp = (jp + 1)*(jh(i2) + 1)*(temp**2)/(pi4*(ll + 1.0))

        call green(up,wp,ex+eh(i2),gam,nq(i2),lp,jp)
        call green(up1,wp1,-ex+eh(i2),-gam,nq(i2),lp,jp)

        do m1 = 1,n ! independent particle response
          ma = m1*nx
          do m2 = m1,n
            mb = m2*nx
            if(nq(i2).eq.1) then ! p and n in separate blocks
              mx = 0
            else
              mx = n
            end if
            g(m1+mx,m2+mx) = g(m1+mx,m2+mx) + temp*uh(i2,ma)*uh(i2,mb)*(up(ma)*wp(mb)+up1(ma)*wp1(mb))
          end do ! m2
        end do ! m1
      end do ! jp
    end do ! i2

! multiply G with the residual interaction V
    do m1 = 1,n2
      do m2 = 1,n2
        g(m2,m1) = g(m1,m2)
        if(m1.le.n) then
         if(m2.le.n) then
            a(m1,m2) = g(m1,m2)*vnn(m2)
         else
            a(m1,m2) = g(m1,m2-n)*vnp(m2-n)
         end if
        else
          if(m2.le.n) then
            a(m1,m2) = g(m1,m2+n)*vnp(m2)
          else
            a(m1,m2) = g(m1,m2)*vnn(m2-n)
          end if
        end if
        if(m1.eq.m2) a(m1,m2) = 1.0 + a(m1,m2)
      end do ! m2
    end do ! m1    

! invert 1 + VG
    call matr(n2,a,b)

! evaluate response to field VEXT
    sfree = 0.0
    do m1 = 1,n2
      drhof(m1) = 0.0
      do m2 = 1,n2
        drhof(m1) = drhof(m1) + g(m1,m2)*vext(m2)*del2
      end do ! m2
      sfree = sfree + vext(m1)*drhof(m1)*del2
    end do ! m1

    srpa = 0.0
    do m1=1,n2
      drhrpa(m1) = 0.0
      do m2 = 1,n2
        drhrpa(m1) = drhrpa(m1) + b(m1,m2)*drhof(m2)
      end do ! m2
      srpa = srpa + vext(m1)*drhrpa(m1)*del2
    end do ! m1

    write(*,'(f11.5,4e12.4)') ex,sfree,srpa
    sumf = sumf + aimag(sfree)*ex*dex/pi
    sumrpa = sumrpa + aimag(srpa)*ex*dex/pi
    ex = ex + dex
  end do

  call sumrle(vext,value,l)
  print*,'# total strength in free response, rpa response, and sum rule'
  write(*,'(''#''21x,e12.4,e14.4,e15.4)') sumf,sumrpa,value


contains

  subroutine factor

!    real(kind = 8), dimension(130), intent(out) :: faclog

    integer :: m
    real(kind = 8) :: fn

    faclog(1) = 0.0
    faclog(2) = 0.0
    fn = 1.0
    do m = 3,130
      fn = fn + 1.0
      faclog(m) = faclog(m - 1) + log(fn)
    end do ! m

  end subroutine factor

  subroutine static

! Woods-Saxon parameters
    real(kind = 8), parameter :: vr = -53.0
    real(kind = 8), parameter :: vt = 20.0
    real(kind = 8), parameter :: vso = -15.5
    real(kind = 8), parameter :: a0 = 0.65
    real(kind = 8), parameter :: r0 = 1.25

    integer :: i,j

    real(kind = 8) :: r,temp,ex
    real(kind = 8) :: rr,dd,all
    real(kind = 8) :: emin,emax,etrial
    real(kind = 8) :: norm

    real(kind = 8), dimension(nrmax) :: u

    read(*,*) del,ngrid
    print *,'# DEL = ',del,'  NGRID = ',ngrid
    read(*,*) nmass,npro
    print *,'# A = ',nmass,' Z = ',npro

    rr = r0*((nmass - 1)**third)
    h = 20.75
    dd = (del**2)/h
    if(npro.gt.0) npro = npro - 1

! set up the Woods-Saxon potential    
    do i = 1,ngrid
      r=del*i
      temp=(r - rr)/a0
      if(temp.gt.30) temp = 30
      ex=exp(temp)
      v(i)=(vr + vt*(nmass - 2*npro - 1)/nmass)/(1 + ex)
      vs(i)=vso*ex/(((1. + ex)**2)*r*a0)
      vc(i)=1.44*npro/r
      if(r.lt.rr) vc(i)=1.44*npro/rr*(1.5 - 0.5*(r/rr)**2)
      vc(i)=vc(i) - 2*vt*(nmass - 2*npro - 1)/nmass/(1 + ex)
    end do

!  determine occupied wave functions
    nmass = 0.0
    npro = 0.0
    nh = 1
    do i=1,ngrid ! initialize densities
      rhop(i) = 0.0
      rhon(i) = 0.0
    end do

    do 
      read(*,*) lh(nh),jh(nh),nq(nh),node
      if(lh(nh).lt.0) then
        nh = nh - 1
        exit
      end if

      all = lh(nh)*(lh(nh) + 1)*h
      if(2*lh(nh).lt.jh(nh)) then
        fl = lh(nh)
      else
        fl = -(lh(nh) + 1)
      end if

      emin = -50.0
      emax = 0.0
      do i = 1,25
!  integrate schroedinger equation
        etrial = (emin + emax)/2.0
        u(1) = (0.1)**(lh(nh) + 1)
        u(2) = (0.2)**(lh(nh) + 1)
        nd = 0
        do j = 2,ngrid-1
          r = del*j
          vv = v(j) + nq(nh)*vc(j) + fl*vs(j) + all/r**2
          sx = dd*(vv - etrial)
          u(j+1) = (2 + sx)*u(j) - u(j-1)
          if(u(j+1)*u(j).lt.0) nd = nd + 1
        end do ! j
        if(nd.gt.node) then
          emax=etrial
        else
          emin=etrial
        end if
      end do ! i 

      norm = sqrt(del*sum(u(:)*u(:)))
      do i = 1,ngrid
        r = del*i
        uh(nh,i) = u(i)/norm
        if(nq(nh).eq.1) then
         rhop(i) = rhop(i) + (jh(nh) + 1)*((uh(nh,i)/r)**2)/pi4
        else
         rhon(i) = rhon(i) + (jh(nh) + 1)*((uh(nh,i)/r)**2)/pi4
        end if
      end do    

      eh(nh) = etrial
      nmass = nmass + jh(nh) + 1
      if(nq(nh).eq.1) npro = npro + jh(nh) + 1
      write(*,'(''# l,2j,q,node,e = '',4i3,f10.5)') lh(nh),jh(nh),nq(nh),node,etrial
      nh = nh + 1
    end do ! reading occupied states and solving Sch. eq.

  end subroutine static

! calculate residual interaction

  subroutine vresid(vnn,vnp,nmass)

    real(kind = 8), intent(in) :: nmass

    real(kind = 8), dimension(matmax), intent(out) :: vnn
    real(kind = 8), dimension(matmax), intent(out) :: vnp

    real(kind = 8), parameter :: rho = 0.16
    real(kind = 8), parameter :: r0 = 1.20
    real(kind = 8), parameter :: a0 = 0.50

    real(kind = 8) :: r,temp

    read(*,*) t0,t3,x,vscal
    write(*,'(''# t0='',f10.2,''   t3='',f8.0,''   x='',f4.2,''   vscal='',f4.2)') t0,t3,x,vscal
    do i=1,n
      r=i*del2
      temp=0.5*t3*rho/(1.0 + exp((r - r0*(nmass**third))/a0))
      vnp(i)=del2*vscal*(t0*(1 + x*half)+temp)/((i*del2)**2)  !  pn interaction
      vnn(i)=del2*vscal*(t0*(1 - x) + temp)/(((i*del2)**2)*2.0) !  pp interaction
    end do

  end subroutine vresid

  subroutine extf(f,l)

    integer, intent(in) :: l

    real(kind = 8), dimension(:), intent(out) :: f

    integer :: i,k

    read(*,*) iso
    print *,'# type of field: ',iso

    if(l.eq.0) then
      k = 2
    else
      k = l
    end if

    do i = 1,n
      r = del2*i
      f(i) = r**k
      if(iso.eq.0) then
        f(i+n) = f(i)
      else if (iso.eq.1) then
        f(i+n) = -f(i)
      else
        f(i+n) = 0.0
      endif
    end do

  end subroutine extf

  subroutine clebsh(j1,j2,j3,m1,m2,m3,ans)

    integer, intent(in) :: j1,j2,j3
    integer, intent(in) :: m1,m2,m3
    real(kind = 8), intent(out) :: ans

    integer :: numin1,numax1,kb,kc,num1,nu
    integer, dimension(10) :: ffi
    real(kind = 8) :: first,part1,ff,second,con
!    real(kind = 8), dimension(10) :: ffi

    if(m1 + m2.ne.m3) then
      ans = 0.0
    else if(j3.le.j1+j2.and.j3.ge.abs(j1-j2)) then
      ffi(1) = half*(j1 + j2 - j3 + 2)
      ffi(2) = half*(j3 + j1 - j2 + 2)
      ffi(3) = half*(j3 + j2 - j1 + 2)
      ffi(4) = half*(j1 + j2 + j3 + 4)
      ffi(5) = half*(j1 + m1 + 2)
      ffi(6) = half*(j1 - m1 + 2)
      ffi(7) = half*(j2 + m2 + 2)
      ffi(8) = half*(j2 - m2 + 2)
      ffi(9) = half*(j3 + m3 + 2)
      ffi(10) = half*(j3 - m3 + 2)

      first = 0.0
      do num1 = 1,10
        first = first + faclog(ffi(num1))
      end do
      first = (first - 2.0*faclog(ffi(4)) + log(j3 + 1.0))*half
      part1 = exp(first)

      numin1 = abs(min(((j3 - j2 + m1)*half),((j3 - j1 - m2)*half),0.0)) + 1.0
      numax1 = min(((j1 + j2 - j3)*half),((j1 - m1)*half),((j2 + m2)*half)) + 1
      kb = (j3 - j2 + m1)*half + 1.0
      kc = (j3 - j1 - m2)*half + 1.0
      ff = 0.0
      do num1 = numin1,numax1
        nu = num1 - 1
        second = -(faclog(num1) + faclog(ffi(1)-nu) + faclog(ffi(6)-nu) + faclog(ffi(7)-nu) + faclog(kb+nu) + faclog(kc+nu))
        con = ((-1)**nu)*exp(second)
        ff = ff + con
      end do
      ans = part1*ff
    else
      ans = 0.0
    end if

  end subroutine clebsh

  subroutine green(up,wp,e,gam,nt,l,j)

    integer, intent(in) :: nt,l,j
    real(kind = 8), intent(in) :: e,gam
    complex(kind = 8), dimension(:), intent(out) :: up,wp

    integer :: nmax,nmax1,nl,k,i,nq

    real(kind = 8) :: dd,all,z

    complex(kind = 8) :: one,ei,ww,ecomp
    complex(kind = 8), dimension(nrmax) :: s

    dd = (del**2)/h
    one = (1.0,0.0)
    ei = (0.0,1.0)
    ecomp = e*one + gam*ei*half
    all = l*(l + 1)*h
    fl = -(l + 1)
    if(2*l.lt.j) fl = l
    nmax = ngrid
    nmax1 = nmax + 1
    nl = nmax - 1

! initialize regular solution
    up(1) = (0.1)**(l+1)
    up(2) = (0.2)**(l+1)
    do i = 2,nmax1
      r = del*(i-1)
      s(i-1) = dd*(v(i-1) + nt*vc(i-1) + fl*vs(i-1) + all/r**2 - ecomp)
      if(i.lt.3.or.i.gt.nmax) cycle
      up(i) = (2.0 + s(i-1))*up(i-1) - up(i-2)
    end do
    z = -s(nl)
    if(z.gt.0.0) then !  initialize irregular solution to outgoing wave
      pk = sqrt(z)
      wp(nl) = (1.0,0.0)
      wp(nmax) = one*(1. - half*(pk**2.)) + ei*pk
    else
      wp(nl) = (0.0,0.0)
      wp(nmax) = 0.001
    endif

    do i = nmax1 - 3,nmax1 - nmax,-1
      wp(i) = (2. + s(i+1))*wp(i+1) - wp(i+2)
    end do

!  Wronskian
    nq = nmax/2
    ww = -(up(nq)*wp(nq+1) - up(nq+1)*wp(nq))/del
    wp(:) = wp(:)/(h*ww)

  end subroutine green

  subroutine matr(nmax,c,d)

    integer, intent(in) :: nmax
    complex(kind = 8), dimension(:,:), intent(inout) :: c,d

    integer :: j,k,l,m,n

    complex(kind = 8) :: t,a,b
    complex(kind = 8), dimension(matmax) :: u,v

    d = 0.0
    do m = 1,nmax
      d(m,m) = 1.0
    end do

    do n = 1,nmax
      t = c(n,n)
      if(abs(t).le.1.e-10) then 
        j=n
        do
          if(j.gt.nmax) then
            print *,'Matrix not invertible'
            exit
          end if
          j=j+1
          t=c(n,j)
          if(abs(t).gt.1.e-10) exit
        end do
        u(:)=c(n,:) ! swap n and j
        v(:)=d(n,:)
        c(n,:)=c(j,:)
        d(n,:)=d(j,:)
        c(j,:)=u(:)
        d(j,:)=v(:)
      end if
      do k = 1,nmax
        if(k.eq.n) cycle
        a = c(k,n)/c(n,n)
        c(k,:) = c(k,:) - a*c(n,:)
        d(k,:) = d(k,:) - a*d(n,:)
      end do
      b = c(n,n)
      c(n,:) = c(n,:)/b
      d(n,:) = d(n,:)/b
    end do

  end subroutine matr

  subroutine sumrle(f,value,j)

    integer, intent(in) :: j
    real(kind = 8), intent(out) :: value
    real(kind = 8), dimension(:), intent(in) :: f

    integer :: i
    real(kind = 8) :: fg,ss
    real(kind = 8), dimension(matmax) :: fp

    do  i = 1,n2
      if(i.eq.1.or.i.eq.n+1) then
        fl = 0.0
      else
        fl = f(i-1)
      endif
      if(i.eq.n.or.i.eq.2*n) then
        fg = (2*f(i) - f(i-1))
      else
        fg = f(i+1)
      endif
      fp(i) = (fg-fl)/(2.0*del2)
    end do

    do i = 1,n
      r = i*del2
      ss = ss + ((fp(i)**2 + (f(i)/r)**2*j*(j+1))*rhop(i*nx) + &
                (fp(i+n)**2 + (f(i+n)/r)**2*j*(j+1))*rhon(i*nx))*(r**2)*del2*pi4
    end do
    value = h*ss/pi4

  end subroutine sumrle

end program rpa