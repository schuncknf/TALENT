!=========================================================
  SUBROUTINE laguerre_general( n, alpha, x, cx )
    USE PARS
    IMPLICIT NONE
    INTEGER, INTENT(IN)  :: n
    REAL (dp ) ::  alpha
    REAL ( dp ) :: cx(0:n)
    INTEGER :: i
    REAL ( dp ), INTENT(IN) ::  x

    IF ( alpha <= -1.0D+00 ) THEN
       WRITE ( *, '(a)' ) ' '
       WRITE ( *, '(a)' ) 'LAGUERRE_GENERAL - Fatal error!'
       WRITE ( *, '(a,g14.6)' ) '  The input value of ALPHA is ', alpha
       WRITE ( *, '(a)' ) '  but ALPHA must be greater than -1.'
       STOP
    END IF
    IF ( n < 0 ) THEN
       RETURN
    END IF
    cx(0) = 1.0D+00
    IF ( n == 0 ) THEN
       RETURN
    END IF
    cx(1) = 1.0D+00 + alpha - x
    DO i = 2, n
       cx(i) = ( ( REAL ( 2 * i - 1, kind = 8 ) + alpha - x ) * cx(i-1)   &
            + ( REAL (   - i + 1, kind = 8 ) - alpha     ) * cx(i-2) ) &
            / REAL (     i,     kind = 8 )
    END DO

  END SUBROUTINE laguerre_general
