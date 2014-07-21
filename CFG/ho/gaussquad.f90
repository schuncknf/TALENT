! gaussquad.f90: Module for gaussian quadrature rules
! http://infty.us/gaussquad/gaussquad.html
! v1.0
!
! Copyright (c) 2012 Christopher N. Gilbreth
!
! Permission is hereby granted, free of charge, to any person obtaining a copy
! of this software and associated documentation files (the "Software"), to deal
! in the Software without restriction, including without limitation the rights
! to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
! copies of the Software, and to permit persons to whom the Software is
! furnished to do so, subject to the following conditions:
!
! The above copyright notice and this permission notice shall be included in all
! copies or substantial portions of the Software.
!
! THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
! IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
! FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
! AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
! LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
! OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
! SOFTWARE.

module gaussquad
  implicit none

  ! References:
  !  1. G. H. Golub and J. H. Welsch, Math. Comp. 23, 221 (1969)
  !  2. Digital Library of Mathematical Functions, Ch. 3.5, 18.9

  integer, parameter :: rk = selected_real_kind(p=15)

contains


  subroutine cpquad(n,alpha,name,w,x)
    ! Compute the weights w_i and abscissas x_i for an n-point Gaussian
    ! quadrature rule for classical polynomials of type "name".
    ! Input:
    !   n:      Degree of quadrature (= number of points)
    !   alpha:  Parameter associated with the polynomials
    !   name:   Name of polynomials: "Legendre", "Laguerre", "Hermite"
    ! Output:
    !   w:      Weights of the quadrature rule
    !   x:      Abscissas of the quadrature rule.
    !
    ! Notes:
    !   The weights w(i) and abscissas x(i) are computed so that
    !     ∫W(x) f(x) dx = ∑_i w(i) f(x(i))
    !   exactly for polynomials f of degree ≤ 2n-1.
    !
    !   Here W(x) is the wieght function associated with the orthogonal
    !   polynomials p_i(x), i.e.
    !     ∫W(x) p_i(x) p_j(x) dx = C_i δ_{i,j}.
    !   Each type of polynomial corresponds to a particlar weight function W.
    !
    !      name         weight function W(x)
    !      ----         --------------------
    !      Legendre       1
    !      Laguerre       x^{alpha} exp(-x)
    !      Hermite        exp(-x^2)
    !
    !   The parameter alpha only affects the Laguerre calculation.
    implicit none
    integer,          intent(in) :: n
    real(rk),         intent(in) :: alpha
    character(len=*), intent(in) :: name
    real(rk),         intent(out) :: w(n), x(n)

    real(rk) :: diag(n), subdiag(n-1), work(max(1,2*n-2)), mu0
    integer  :: info, i
    real(rk), allocatable :: P(:,:)

    allocate(P(n,n))

    ! Compute the symmetric tridiagonal matrix "J" for the given type of polynomial
    call cpj(n,alpha,name,diag,subdiag)
    ! Diagonalize (using LAPACK)
    call DSTEV('V', n, diag, subdiag, P, n, work, info)
    if (info .ne. 0) stop "Error diagonalizing in cpquad"
    ! Eigenvalues are the abscissas, i.e. the zeroes of p_n(x)
    x = diag
    ! Calculate the zeroth moment of the weight function W(x)
    call cpzm(alpha,name,mu0)
    ! Weights (Eq. 2.6 in Golub & Welsch)
    do i=1,n
       w(i) = P(1,i)**2 * mu0
    end do
  end subroutine cpquad


  subroutine cpzm(alpha,name,mu0)
    ! Calculate
    !   μ_0 = ∫W(x)dx,
    ! the zeroth moment of the wieght function for a given type of classical
    ! orthgonal polynomial.
    ! Input:
    !   alpha:  Parameter associated with the polynomials
    !   name:   Name of the polynomials
    ! Output:
    !   mu0:    The zeroth moment, as above.
    implicit none
    real(rk),         intent(in) :: alpha
    character(len=*), intent(in) :: name
    real(rk),         intent(out) :: mu0

    if (name == "Legendre") then
       mu0 = 2
    else if (name == "Laguerre") then
       mu0 = gamma(alpha + 1)
    else if (name == "Hermite") then
       mu0 = sqrt(2*asin(1._rk))
    else
       stop "Unknown polynomial in cpzm"
    end if
  end subroutine cpzm


  subroutine cpj(n,alpha,name,diag,subdiag)
    ! Calculate the symmetric tridiagonal matrix J associated with the classical
    ! orthogonal polynimals given by 'name'.
    ! Input:
    !   n:       Degree (= number of roots)
    !   alpha:   Parameter associated with the polynomials
    !   name:    Name of polynimials: "legendre", "laguerre", "hermite"
    ! Output:
    !   diag:    Diagonal matrix elements of J
    !   subdiag: Subdiagonal matrix elements of J
    ! Notes:
    !   Reference: G. H. Golub and J. H. Welsch, Math. Comp. 23, 221 (1969)
    implicit none
    integer,          intent(in) :: n
    real(rk),         intent(in) :: alpha
    character(len=*), intent(in) :: name
    real(rk),         intent(out) :: diag(n), subdiag(n-1)

    integer :: i

    ! In general, if
    !   P_n(x) = (a_n x + b_n) P_{n-1}(x) - c_n P_{n-2}(x)
    ! Then the n diagonal elements of the matrix are
    !   α_i = -b_i/a_i,   i=1,...,n
    ! And the n-1 subdiagonal elements are
    !   β_i = Sqrt(c_{i+1}/(a_i a_{i+1})),  i=1,...,n-1.

    if (name == "Legendre") then
       ! Legendre polynomials:
       !   P_n(x) = (a_n x + b_n) P_{n-1}(x) - c_n P_{n-2}(x)
       ! where
       !   a_n = (2n-1)/n
       !   b_n = 0
       !   c_n = (n-1)/n.
       ! So
       !   α_i = -b_i/a_i = 0
       !   β_i = Sqrt(i^2/[(2i-1)(2i+1)])
       diag = 0._rk
       do i=1,n-1
          subdiag(i) = sqrt(i**2 / ((2*i-1._rk)*(2*i+1._rk)))
       end do
    else if (name == "Laguerre") then
       ! Associated Laguerre polynomials:
       !   L^{α}_n(x) = (a_n x + b_n) L^{α}_{n-1}(x) - c_n L^{α}_{n-2}(x)
       ! where
       !   a_n = -1/n
       !   b_n = (α + 2n-1)/n
       !   c_n = (α + n-1)/n
       ! So
       !   α_i = (α + 2i-1)
       !   β_i = sqrt((α+i)i)
       do i=1,n
          diag(i) = alpha+2*i-1
       end do
       do i=1,n-1
          subdiag(i) = sqrt((alpha+i)*i)
       end do
    else if (name == "Hermite") then
       ! Hermite polynomials:
       !   H_n(x) = (a_n x + b_n) H_{n-1}(x) - c_n H_{n-2}(x)
       ! where
       !   a_n = 2
       !   b_n = 0
       !   c_n = 2*(n-1)
       ! So
       !   α_i = 0
       !   β_i = sqrt(i/2)
       diag = 0._rk
       do i=1,n-1
          subdiag(i) = sqrt(i/2._rk)
       end do
    else
       stop "Unknown polynomial in cpj"
    end if
  end subroutine cpj

end module gaussquad

