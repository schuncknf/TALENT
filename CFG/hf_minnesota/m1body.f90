SUBROUTINE TPLUSU(L)
USE PARS
USE TVUMAT
IMPLICIT NONE
 INTEGER::L
 INTEGER::I
 REAL(DP)::FACTOR
 FACTOR=HBROMG
 TUMAT(:,:)=0.0D0
 DO I=0,NUSE
    TUMAT(I,I)=FACTOR*(2*I+1.5D0)
 ENDDO
END SUBROUTINE TPLUSU
