      SUBROUTINE PRICC (ND,NF,WCHARM,GAMMAC,DELTAC,MODHEAD,JOBNUM)
C***********************************************************************
C***  PRINTOUT OF SCHARMER CONTINUUM CORES (WEIGHT FUNCTION)
C***    ( FIRST DECIMAL DIGIT, ROUNDED TO THE NEAREST INTEGER )
C***********************************************************************
  
      IMPLICIT NONE
  
      INTEGER, INTENT(IN) :: ND, NF, JOBNUM
      REAL, INTENT(IN) :: GAMMAC, DELTAC

      REAL, DIMENSION(ND,NF) :: WCHARM
      CHARACTER(100) :: MODHEAD

      INTEGER :: K, K1, K2, L
 
      PRINT 1,MODHEAD,JOBNUM
    1 FORMAT (1X,  A  ,20X,'JOB NO.',I7,//)

      PRINT 2,GAMMAC,DELTAC
    2 FORMAT (10X,'CONTINUUM SCHARMER WEIGHTS (FIRST DECIMAL DIGIT) -',
     $      '   GAMMAC=',F5.1,'   DELTAC=',F5.1,//)

C***  BLOCKS OF 100 FREQUENCY POINTS
      K1 = 1
   10 K2 = MIN0(NF, K1+99)

      PRINT 21,(K/100          , K=K1, K2)
      PRINT 21,(K/10 - K/100*10, K=K1, K2)
      PRINT 21,(K    - K/10 *10, K=K1, K2)
   21 FORMAT (10X,100I1)
      PRINT 23
   23 FORMAT (1X)

      DO L=1, ND
        PRINT 22,L,(INT(WCHARM(L,K)*10.+0.5), K=K1, K2)
   22   FORMAT (I5,5X,100I1)
      ENDDO

      IF (K2 .LT. NF) THEN
         PRINT 11
   11    FORMAT (//, ' ... CONTINUED ...', /)
         K1 = K1+100
         GOTO 10
      ENDIF

      RETURN
      END
