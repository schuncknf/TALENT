!spin 对角化且两块相等时，轨道计半，然后乘2，能量同理
PROGRAM MAIN
USE PARS
USE GAUSS
USE TVUMAT
USE EIGENS
USE HRHOMAT

IMPLICIT NONE
 INTEGER::NREAL,NFIX,L,ISPINC
 INTEGER::ERROR,ITNOW,ITMAX,FLAGMAX,IFLAG
 REAL(DP)::RMIX,EOLD,ENEW,VFACTOR,RSURF,EPSI
 REAL(DP)::ALPHA
 L=0;NUSE=2;GUSEX=80;QUAD="Laguerre" !BASIS SETTINGS
 NREAL=2                              !PARTICLE NUMBER
 ERROR=1.0D0;EPSI=1.0D-14;RMIX=0.5D0  !ITERATION OPTIONS
 ITNOW=0;ITMAX=50;FLAGMAX=5;IFLAG=0
 ISPINC=1                             !TREAT ONLY SPIN UPPER PART
! 
 IF(ISPINC.NE.1) STOP "NOT IMPLEMENTED YET"
 IF((ISPINC.EQ.1).AND.MOD(NREAL,2)==1) STOP 'NFIX NOW CAN ONLY BE EVEN'
 NFIX=NREAL/(ISPINC+1)!WHEN IPSINC=1, ONE ONLY TREAT THE UPPER BLOCK OF H MATRIX 
!
 HBAR2M=20.73D0;HBROMG=10.0D0;BOSC=SQRT(0.5D0*HBROMG/HBAR2M)!CONSTANTS  
 V0R=2.0D2;KAPPAR=1.487D0
 V0S=-9.185D1;KAPPAS=0.465D0
!
 CALL INITIAL(NFIX,L,EOLD,ENEW)!INITIALIZE T,V,RHO,E! CALL V2TEST(L,RSURF,VFACTOR)
 CALL ITERS(NFIX,ISPINC,EPSI,ITMAX,ITNOW,IFLAG,FLAGMAX,RMIX,EOLD,ENEW)
        
END PROGRAM
