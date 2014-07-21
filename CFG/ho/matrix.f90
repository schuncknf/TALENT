!=========================================================
SUBROUTINE ANLS(LOGANL)!calculate log(Anl) at each n,l point,B^1.5 factor is not included
USE PARS
IMPLICIT NONE
    REAL(dp)::alpha
    REAL(dp)::LOGANL(0:NMAX,0:LMAX)
    INTEGER::n,L    
    DO n=0,NMAX
       DO L=0,LMAX
          LOGANL(n,L)=0.5D0*LOG(2.0D0)&
                   +0.5D0*(LOG_GAMMA(REAL(n+1,kind=DP))&
                   -LOG_GAMMA(REAL(n+l+1.5,kind=DP)))
       ENDDO
    ENDDO
    RETURN
END SUBROUTINE

SUBROUTINE LNLS(X,LPOLY,NAM)! CONSTRUCT L_{nl}^{alpha} at each l, n, x point
USE PARS
IMPLICIT NONE
    INTEGER::I,L,N,NR
    REAL(dp)::alpha,XI2
    REAL(dp)::X(1:NGAUSS),LPOLY(0:NMAX,0:LMAX,1:NGAUSS)
    CHARACTER(LEN=*)::NAM
    LPOLY(:,:,:)=0.0D0
!    
    DO N=0,NMAX
       DO L=0,LMAX
          ALPHA=L+0.5D0
          DO I=1,NGAUSS
             IF(NAM=="Hermite")THEN
                XI2=X(I)**2
                CALL laguerre_general(N,ALPHA,XI2,LPOLY(0:N,L,I))
             ELSEIF(NAM=="Legendre")THEN  
                XI2=TAN(PI4*(X(I)+1.0d0))**2
                CALL laguerre_general(N,ALPHA,XI2,LPOLY(0:N,L,I))
             ELSEIF(NAM=="Laguerre")THEN
                XI2=X(I)
                CALL laguerre_general(N,ALPHA,XI2,LPOLY(0:N,L,I))
             ELSE
                STOP"WRONG QUAD"
             ENDIF
          ENDDO
       ENDDO
    ENDDO      
    RETURN
END SUBROUTINE
!=========================================================
SUBROUTINE TMATRX(NR,L,TMAT)!calculate T matrix element with fixed L 
USE PARS
IMPLICIT NONE
INTEGER::L,NR,N
REAL(DP)::HBOM2
REAL(DP)::TMAT(0:NMAX,0:NMAX)
INTEGER::I,J

   TMAT(:,:)=0
   HBOM2=0.5D0*HBOMG

   TMAT(0,0)=HBOM2*(L+1.5D0)
   DO I=1,NR
      N=2*I+L
      TMAT(I,I)=HBOM2*(N+1.5D0)
      J=I-1
      TMAT(J,I)=HBOM2*SQRT(I*(I+L+0.5D0))
      TMAT(I,J)=TMAT(J,I)!the t matrix is symmetric
   ENDDO
   RETURN
END SUBROUTINE
!=========================================================
SUBROUTINE INTACT(X,V,NAMI,NAMQ)!calculate a given v(x) at each integral point
USE PARS
IMPLICIT NONE
 INTEGER::I
 REAL(DP)::X(1:NGAUSS),V(1:NGAUSS)
 CHARACTER(LEN=*)::NAMI,NAMQ
 REAL(DP)::XREAL
    DO I=1,NGAUSS
       IF(NAMQ=="Hermite")THEN
          XREAL=X(I)
          IF(NAMI=="COL")THEN
             V(I)=1/XREAL!X ARE VALUES OF br RATHER TAHN r
          ELSEIF(NAMI=="OSC")THEN
             V(I)=0.5d0*HBAR**2/MASS*XREAL**2!1/2*m*omega^2*x^2
          ELSE
             STOP "WRONG INTERACTION"
          ENDIF
       ELSEIF(NAMQ=="Legendre")THEN
          XREAL=TAN(PI4*(X(I)+1.0d0))
          IF(NAMI=="COL")THEN
             V(I)=1/XREAL
          ELSEIF(NAMI=="OSC")THEN
             V(I)=0.5d0*HBAR**2/MASS*XREAL**2
          ELSE
             STOP "WRONG INTERACTION"
          ENDIF
       ELSEIF(NAMQ=="Laguerre")THEN
          XREAL=X(I)
          IF(NAMI=="COL")THEN
             V(I)=1/SQRT(XREAL)
          ELSEIF(NAMI=="OSC")THEN
             V(I)=0.5d0*HBAR**2/MASS*XREAL
          ELSE
             STOP "WRONG INTERACTION"
          ENDIF
       ELSE
          STOP"WRONG QUAD"
       ENDIF
    ENDDO
    RETURN
END SUBROUTINE
!=========================================================
SUBROUTINE VMATRX(NR,L,X,W,V,LPOLY,LOGANL,NAM,VMAT)!CALCULATE <nl|V|n'l> WITH L FIXED, n=0~NR
USE PARS
IMPLICIT NONE
 INTEGER::NR,L
 REAL(DP)::X(1:NGAUSS),W(1:NGAUSS),V(1:NGAUSS)
 REAL(DP)::LPOLY(0:NMAX,1:NGAUSS),LOGANL(0:NMAX),VMAT(0:NMAX,0:NMAX)
 INTEGER::NI,NJ,K,KSTART
 REAL(DP)::SUMS,XREAL,XREAL2
 CHARACTER(LEN=*)::NAM
 KSTART=1
 IF(NAM=="Hermite") KSTART=NGAUSS/2+1
DO NI=0,NR
   DO NJ=NI,NR
      SUMS=0.0D0
      DO K=KSTART,NGAUSS
         IF(NAM=="Hermite")THEN
            XREAL=X(K)
            XREAL2=XREAL**2
            SUMS=SUMS+W(K)*V(K)*LPOLY(NI,K)*LPOLY(NJ,K)&
             *XREAL**(2*L+2)&
             *EXP(LOGANL(NI)+LOGANL(NJ))
         ELSEIF(NAM=="Legendre")THEN
            XREAL=TAN(PI4*(X(K)+1.0d0))
            XREAL2=XREAL**2
            SUMS=SUMS+W(K)*V(K)*LPOLY(NI,K)*LPOLY(NJ,K)&
             *XREAL**(2*L+2)&
             *EXP(-XREAL2+LOGANL(NI)+LOGANL(NJ))&
             *PI4/COS(PI4*(X(K)+1.0d0))**2
         ELSEIF(NAM=="Laguerre")THEN
            XREAL2=X(K)
            SUMS=SUMS+W(K)*V(K)*LPOLY(NI,K)*LPOLY(NJ,K)&
             *EXP(LOGANL(NI)+LOGANL(NJ))&
             *0.5D0
         ELSE
            STOP"WRONG QUAD"
         ENDIF
      ENDDO
      VMAT(NI,NJ)=SUMS!in V there could be B factor, but B from r^2dr and Rnl*Rn'l would cancel eachother 
      VMAT(NJ,NI)=VMAT(NI,NJ)
   ENDDO
ENDDO
RETURN
END SUBROUTINE 
!=========================================================
SUBROUTINE CMATRX(NR,L,X,W,LPOLY,LOGANL,INTAC,QUAD,VMAT)!CALCULATE MATRIX ELEMENT FOR COLUMB WITH L FIXED
USE PARS
IMPLICIT NONE
 INTEGER::NR,L
 REAL(DP)::X(1:NGAUSS),W(1:NGAUSS)
 REAL(DP)::LPOLY(0:NMAX,1:NGAUSS),LOGANL(0:NMAX),VMAT(0:NMAX,0:NMAX)
 INTEGER::NI,NJ,K,KSTART
 REAL(DP)::SUMS,XREAL,XREAL2
 CHARACTER(LEN=*)::INTAC,QUAD
FACTOR=1.0D0
KSTART=1
IF(QUAD=="Hermite") KSTART=NGAUSS/2+1
DO NI=0,NR
   DO NJ=NI,NR
      SUMS=0.0D0
      DO K=KSTART,NGAUSS
         IF(QUAD=="Hermite")THEN
            IF(INTAC=="COL")THEN
               XREAL=X(K)
               XREAL2=XREAL**2
               SUMS=SUMS+W(K)*LPOLY(NI,K)*LPOLY(NJ,K)&
                *XREAL**(2*L+1)&
                *EXP(LOGANL(NI)+LOGANL(NJ))
            ELSEIF(INTAC=="OSC")THEN
               XREAL=X(K)
               XREAL2=XREAL**2
               SUMS=SUMS+W(K)*LPOLY(NI,K)*LPOLY(NJ,K)&
                *XREAL**(2*L+4)&
                *EXP(LOGANL(NI)+LOGANL(NJ))&
                *0.5D0
            ELSE
               STOP "WRONG INTERACTION"
            ENDIF
         ELSEIF(QUAD=="Legendre")THEN
            IF(INTAC=="COL")THEN
               XREAL=TAN(PI4*(X(K)+1.0d0))
               XREAL2=XREAL**2
               SUMS=SUMS+W(K)*LPOLY(NI,K)*LPOLY(NJ,K)&
                *XREAL**(2*L+1)&
                *EXP(-XREAL2+LOGANL(NI)+LOGANL(NJ))&
                *PI4/COS(PI4*(X(K)+1.0d0))**2
            ELSEIF(INTAC=="OSC")THEN
               XREAL=TAN(PI4*(X(K)+1.0d0))
               XREAL2=XREAL**2
               SUMS=SUMS+W(K)*LPOLY(NI,K)*LPOLY(NJ,K)&
                *XREAL**(2*L+4)&
                *EXP(-XREAL2+LOGANL(NI)+LOGANL(NJ))&
                *PI4/COS(PI4*(X(K)+1.0d0))**2*0.5d0
            ELSE
               STOP "WRONG INTERACTION"
            ENDIF
         ELSEIF(QUAD=="Laguerre")THEN
            IF(INTAC=="COL")THEN!ALPHA=L
               XREAL2=X(K)
               SUMS=SUMS+W(K)*LPOLY(NI,K)*LPOLY(NJ,K)&
                *EXP(LOGANL(NI)+LOGANL(NJ))&
                *0.5D0
            ELSEIF(INTAC=="OSC")THEN!ALPHA=L+1.5D0
               XREAL2=X(K)
               SUMS=SUMS+W(K)*LPOLY(NI,K)*LPOLY(NJ,K)&
                *EXP(LOGANL(NI)+LOGANL(NJ))&
                *0.25D0
            ELSE
               STOP "WRONG INTERACTION"
            ENDIF
         ELSE
            STOP"WRONG QUAD "
         ENDIF
      ENDDO
      VMAT(NI,NJ)=SUMS*BOSC!from coloumb interaction there is still one b
      VMAT(NJ,NI)=VMAT(NI,NJ)
   ENDDO
ENDDO
RETURN
END SUBROUTINE 
!=========================================================
