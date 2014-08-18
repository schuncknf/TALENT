! Modified from f77 to f90 by B.G. Carlsson
MODULE talmi

  PRIVATE

  PUBLIC :: TMBR, INIT_TALMI, TMSH, TRI, TMB
! Addded TRI and TMB to the list of public stuff /Daniel W
  INTERFACE TMBR
     MODULE PROCEDURE TMBR
  END INTERFACE

  INTEGER, parameter :: max_bin = 200 !/added
  INTEGER, parameter :: max_tin = 200 !/added

!  REAL(8), SAVE :: BIN(0:99,0:99),TIN(0:99,0:99,0:99) /changed
  REAL(8), SAVE :: BIN(0:max_bin,0:max_bin),TIN(0:max_tin,0:max_tin,0:max_tin)


CONTAINS

  ! TMB(EE,LL,ER,LR,E1,L1,E2,L2,LM,D)
  ! 
  ! Input: EE-centre of mass energy, LL-centre of mass angular moment
  !        ER-relative energy, LR-relative angular moment
  !        E1-energy of the first particle, L1-angular moment of the first particle
  !        E2-energy of the second particle, L2-angular moment of the second particle
  !        LM-total angular moment   
  !        D = 1d0 gives transformation R = 1/sqrt(2) (r1 + r2) and r = 1/sqrt(2) (r1 - r2)

! B.G. Carlsson
         REAL*8 FUNCTION TMBR(EE,LL,ER,LR,E1,L1,E2,L2,LM)
!   	   TALMI-MOSHINSKY BRACKET
!	    (EE,LL;ER,LR:LM/E1,L1;E2,L2:LM)D
	IMPLICIT NONE
	REAL*8 D
	INTEGER EE,LL,ER,LR,E1,L1,E2,L2,LM 

  !     Compensates for the different phases in the 
  !     definitions of the basis functions in the 
  !     hosphe code compared to those used to derive
  !     the TMB routine. The difference is 
  !     G_nl = (-1)^n g_nl where g_nl is hosphe's
  !     radial part of the wavefunction. 


        D = 1d0

!        TMBR = 1.d0
        
        TMBR = TMB(EE,LL,ER,LR,E1,L1,E2,L2,LM,D) &              
             * (-1)**( (L1 + L2 - LL - LR)/2 )
 
 
        END FUNCTION TMBR


!  B.G. Carlsson
  SUBROUTINE INIT_TALMI
    IMPLICIT NONE
    CALL BINOM
    CALL TRINOM
    RETURN
  END SUBROUTINE INIT_TALMI

!** HARMONIC OSCILLATOR BRACKETS CALCULATION PROGRAM 
!** Version 1.0: May 2001 
!** E-mail: Gintautas_Kamuntavicius@fc.vdu.lt
!** Reference: nucl-th/0105009
!** WWW: http://www.nuclear.physics.vdu.lt  
!** COMMON-blocks:
!** BIN-array of binomial coefficients
!** TIN-array of trinomial coefficients
!** TMB-real*8 function of HOBs:
!** Input: EE-centre of mass energy, LL-centre of mass angular moment
!**        ER-relative energy, LR-relative angular moment
!**        E1-energy of the first particle, L1-angular moment of the first particle
!**        E2-energy of the second particle, L2-angular moment of the second particle
!**        LM-total angular moment   
!** Output: TMB-HARMONIC OSCILLATOR BRACKET
      SUBROUTINE TMSH
      IMPLICIT NONE
      REAL tim1,tim2,timd
      REAL*8 T1,skirtm !,TMB
      INTEGER LAST
!	COMMON /BINO/ BIN
!	COMMON /TINO/ TIN
!c	LAST=0,1,2,3,4...
      LAST=3
	CALL CPU_TIME(tim1)
      CALL BINOM
	CALL TRINOM
      T1=TMB(2,2,2,0,1,1,3,3,2,1.D0)
      write(6,2) T1
      T1=TMB(1,1,5,5,1,1,5,5,6,1.D0)
      write(6,2) T1
      T1=TMB(1,1,3,3,2,2,2,2,4,1.D0)
      write(6,2) T1
      T1=TMB(5,3,1,1,2,2,4,4,3,1.D0)
      write(6,2) T1
      T1=TMB(5,5,2,2,2,2,5,5,4,1.D0)
      write(6,2) T1
      T1=TMB(3,3,8,6,6,2,5,3,4,1.D0)
      write(6,2) T1
      T1=TMB(2,0,9,5,6,2,5,3,5,1.D0)
      write(6,2) T1
      T1=TMB(2,2,10,2,6,2,6,4,2,1.D0)
      write(6,2) T1
      T1=TMB(8,2,4,4,6,2,6,4,4,1.D0)
      write(6,2) T1
2     FORMAT(' TMB =',D22.15)
      CALL ORTNTMB(LAST,skirtm)
	CALL CPU_TIME(tim2)
      write(6,1) skirtm
1     FORMAT(' tikslumas =',D22.15)
	timd=tim2-tim1
	write(6,3) timd
3	FORMAT(' laikas =' ,E9.3) 
      END SUBROUTINE TMSH
!C
      SUBROUTINE BINOM
!C 	THE ARRAY OF BINOMIAL COEFFICIENTS
!C     BIN(I,J)= = I!/J!/(I-J)! 
      IMPLICIT NONE
!      REAL*8 BIN(0:99,0:99)
      INTEGER I,K
      !COMMON /BINO/ BIN
!	DO I=0,99 /changed
        DO I=0,max_bin
	BIN(I,0)=1.D0
	BIN(I,I)=1.D0
	DO K=1,I/2
      BIN(I,K)=DNINT(BIN(I,K-1)/DFLOAT(K)*DFLOAT(I-K+1))
	BIN(I,I-K)=BIN(I,K)
	END DO
      END DO
	RETURN
       END SUBROUTINE BINOM
!C
      INTEGER FUNCTION TRI(I,J,K)
!C     TRIADIC CONDITION FOR MOMENTS I/2,J/2,K/2:
!C     I+J>=K, I+K>=J, J+K>=I,
!C     I/2+J/2+K/2 = INTEGER.
!C     TRI=1, WHEN TRIADIC CONDITION IS FULFILLED, TRI=0 OTHERWISE	 
      IMPLICIT NONE
      INTEGER I,J,K,L
      TRI=0
	L=I+J+K
      IF(L/2*2.NE.L) RETURN
	L=L/2
      IF((L-I)*(L-J)*(L-K).LT.0) RETURN
      TRI=1
      RETURN
    END FUNCTION TRI
!C
      REAL*8 FUNCTION C6J(I,J,K,L,M,N)
!C     6J - COEFFICIENT
!C     ( I/2  J/2  K/2 )
!C     ( L/2  M/2  N/2 ) 
!C     [JB 65] (22.1.4)
	IMPLICIT NONE
      INTEGER I,J,K,L,M,N,I1,I2,I3,I4,I5,IZ,JZ,KZ !,TRI
 	REAL*8 T,DZ !,BIN(0:99,0:99)
!      COMMON /BINO/ BIN
	C6J=0.D0
	IF(TRI(I,J,K)*TRI(I,M,N)*TRI(J,L,N)*TRI(K,L,M).EQ.0) RETURN
	I1=(I+J+K)/2
      I2=(I+M+N)/2
      I3=(J+L+N)/2
      I4=(K+L+M)/2
      I5=(I+K+L+N)/2
	T=DSQRT(DFLOAT((I1+1)*(I4+1))/DFLOAT((I2+1)*(I3+1))*&
      BIN(I2,I)*BIN(I,I2-N)*BIN(I4,L)/&
      BIN(I3,N)*BIN(L,I4-K)*BIN(I1,K)*BIN(K,I1-J)/&
      BIN(N,I3-L))/DFLOAT(I2-N+1)
	JZ=MAX0(0,(I+L-J-M)/2)
	DZ=1.D0
	IF((JZ+I5)/2*2.NE.(JZ+I5)) DZ=-1.D0
      KZ=MIN0(I2-M,I3-J,I5-K)
	DO IZ=JZ,KZ
	C6J=C6J+DZ*T*BIN(I2-M,IZ)/BIN(I2,I5-K-IZ)*&
      BIN(I2-I,I3-J-IZ)/BIN(I4-I3+J+IZ+1,I2-N+1)
      DZ=-DZ
      END DO
	RETURN
END FUNCTION C6J
!C
      REAL*8 FUNCTION C9J(J1,J2,J3,L1,L2,L3,K1,K2,K3)
!C     9J COEFICIENT
!C     (J1/2 J2/2 J3/2)
!C     (L1/2 L2/2 L3/2)
!C	(K1/2 K2/2 K3/2)  	 
!C     [JB 65] (24.33)
      IMPLICIT NONE
! 	REAL*8 C6J
      INTEGER J1,J2,J3,L1,L2,L3,K1,K2,K3,I,J,K,L !,TRI
      C9J=0.D0
      L=TRI(J1,J2,J3)*TRI(L1,L2,L3)*TRI(K1,K2,K3)*&
           TRI(J1,L1,K1)*TRI(J2,L2,K2)*TRI(J3,L3,K3)
      IF(L.EQ.0) RETURN
      J=MAX0(IABS(J1-K3),IABS(J2-L3),IABS(L1-K2))
      K=MIN0(J1+K3,J2+L3,L1+K2)
      DO I=J,K,2
	C9J=C9J+DFLOAT(I+1)*C6J(J1,J2,J3,L3,K3,I)*&
     C6J(L1,L2,L3,J2,I,K2)*C6J(K1,K2,K3,I,J1,L1)
      END DO
	IF(J/2*2.NE.J) C9J=-C9J
 	RETURN
END FUNCTION C9J

!C
	REAL*8 FUNCTION KL0(I,J,K)
!C	KLEBS-GORDAN COEFFICIENT WITH ZERO PROJECTIONS OF MOMENTA
!C	(I, J, K)
!C	(0, 0, 0)  
!C	I,J,K - MOMENTA = INTEGER NUMBERS
!C	[JB,65] (15.10)
      IMPLICIT NONE
	REAL*8 T !,BIN(0:99,0:99)
      INTEGER I,J,K,L,M !,TRI
!      COMMON /BINO/ BIN
      KL0=0.D0
      IF(TRI(I,J,K).EQ.0) RETURN
 	L=(I+J+K)/2
	M=L-K
	T=1.D0
	IF(M/2*2.NE.M) T=-1.D0
      KL0=T*BIN(K,L-J)*BIN(L,K)/&
           DSQRT(BIN(2*K,2*(L-J))*BIN(2*L+1,2*K+1))
      RETURN
    END FUNCTION KL0
!C
 	SUBROUTINE TRINOM
!C	THE ARRAY OF TRINOMIAL COEFFICIENTS
!C	TIN(I,J,K)=I!!/J!!/K!!
	IMPLICIT NONE
!	REAL*8 !TIN(0:99,0:99,0:99)
	INTEGER I,J,K,M,N
!	COMMON /TINO/ TIN
	TIN(0,0,0)=1.D0
	TIN(1,1,1)=1.D0
!	DO I=2,99 /changed
        DO I=2,max_tin
	M=I-I/2*2
	TIN(I,I,M)=1.D0
	TIN(I,M,I)=1.D0
	N=M+2
	DO J=I,N,-2
	DO K=N,J,2
	TIN(I,J,K)=TIN(I,J,K-2)/DFLOAT(K)
	TIN(I,K,J)=TIN(I,J,K)
	END DO
	TIN(I,J-2,M)=TIN(I,J,M)*DFLOAT(J)
	TIN(I,M,J-2)=TIN(I,J-2,M)
	END DO
	END DO
	RETURN
END SUBROUTINE TRINOM
!C	
	REAL*8 FUNCTION G(E1,L1,EA,LA,EB,LB)
	IMPLICIT NONE
!	REAL*8 KL0 !, TIN(0:99,0:99,0:99)
	INTEGER E1,L1,EA,LA,EB,LB
!	COMMON /TINO/ TIN
	G=KL0(LA,LB,L1)*DSQRT(DFLOAT((2*LA+1)*(2*LB+1))*&
        TIN(E1-L1,EA-LA,EB-LB)*TIN(E1+L1+1,EA+LA+1,EB+LB+1))
	RETURN
END FUNCTION G
!C
	REAL*8 FUNCTION TMB(EE,LL,ER,LR,E1,L1,E2,L2,LM,D)
!C     	   TALMI-MOSHINSKY BRACKET
!C	    (EE,LL;ER,LR:LM/E1,L1;E2,L2:LM)D
	IMPLICIT NONE
	REAL*8 S,D,T !G,C9J
	INTEGER EE,LL,ER,LR,E1,L1,E2,L2,LM !,TRI
	INTEGER L,M,ED,LD,EB,LB,EC,LC,EA,LA
 	TMB=0.D0
	IF(EE+ER.NE.E1+E2) RETURN
	IF(TRI(2*LL,2*LR,2*LM)*TRI(2*L1,2*L2,2*LM).EQ.0) RETURN
	T=DSQRT((D**(E1-ER))/((1.D0+D)**(E1+E2)))
	M=MIN0(ER,E2)
	S=1.D0
	DO 1 ED=0,M
	EB=ER-ED
	EC=E2-ED
	EA=E1-ER+ED
	DO 2 LD=ED,0,-2 
	DO 3 LB=EB,0,-2
	IF(TRI(LD,LB,LR).EQ.0) GO TO 3
	DO 4 LC=EC,0,-2
	IF(TRI(LD,LC,L2).EQ.0) GO TO 4
	DO 5 LA=EA,0,-2
	IF((TRI(LA,LB,L1).EQ.0).OR.(TRI(LA,LL,LC).EQ.0)) GO TO 5
	TMB=TMB+S*T*&
          C9J(2*LA,2*LB,2*L1,2*LC,2*LD,2*L2,2*LL,2*LR,2*LM)*&
          G(E1,L1,EA,LA,EB,LB)*G(E2,L2,EC,LC,ED,LD)*&
          G(EE,LL,EA,LA,EC,LC)*G(ER,LR,EB,LB,ED,LD)
5 	CONTINUE
4 	CONTINUE
3 	CONTINUE
2 	CONTINUE
	S=S*(-D)
1 	CONTINUE
	RETURN
        END FUNCTION TMB

	SUBROUTINE ORTNTMB(LAST,skirtm)
      IMPLICIT NONE 
      REAL*8 D,skirt,skirtm,atsak,check,tarp1,tarp2
      INTEGER LAST,LMB,EE,EE_,LL,LL_,e,e_,l,l_,l1,l2,e1,e2
      INTEGER*4 skaitl
      skaitl=0
      D=1.D0
      DO LMB=0,LAST
      DO EE=0,LAST
       LL=EE
       DO WHILE(LL.GE.0)
        DO EE_=0,LAST
         LL_=EE_
         DO WHILE(LL_.GE.0)
          DO e=0,LAST
           l=e
           DO WHILE(l.GE.0)
            IF(IABS(LL-l).GT.LMB) GOTO 42
            IF(LL+l.LT.LMB) GOTO 42
            DO e_=0,LAST
             l_=e_
             DO WHILE(l_.GE.0)
              IF(IABS(LL_-l_).GT.LMB) GOTO 43
              IF((LL_+l_).LT.LMB) GOTO 43
              atsak=0.D0
              DO e1=0,EE+e
               l1=e1
               DO WHILE(l1.GE.0)
                DO e2=0,EE+e
                 IF((EE+e).NE.(e1+e2)) CYCLE
                 IF((EE_+e_).NE.(e1+e2)) CYCLE
                 l2=e2
                 DO WHILE(l2.GE.0)
                  IF(IABS(l1-l2).GT.LMB) GOTO 41
                  IF(l1+l2.LT.LMB) GOTO 41
                  tarp1=TMB(EE,LL,e,l,e1,l1,e2,l2,LMB,D)
                  tarp2=TMB(e1,l1,e2,l2,EE_,LL_,e_,l_,LMB,D)
                  atsak=atsak+tarp1*tarp2
                  skaitl=skaitl+2
41                l2=l2-2
                 END DO
                END DO
                l1=l1-2
               END DO
              END DO
              IF((EE.EQ.EE_).AND.(LL.EQ.LL_).AND.&
                 (e.EQ.e_).AND.(l.EQ.l_)) THEN
               check=1.0D0
              ELSE
               check=0.0D0
              END IF
              skirt=DABS(atsak-check)
              IF(skirtm<skirt) skirtm=skirt
43            l_=l_-2
             END DO
            END DO
42          l=l-2
           END DO
          END DO
          LL_=LL_-2
         END DO
        END DO
        LL=LL-2
       END DO
      END DO
      END DO
      write(6,4) skaitl
4     FORMAT(' TMB viso =',I10)
      RETURN
    END SUBROUTINE ORTNTMB

  END MODULE talmi
