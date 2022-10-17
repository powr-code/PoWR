      SUBROUTINE PRIGAHIST(GAHIST, MAXGAHIST, LEVEL, N, AG, BAG, 
     >                     BPGAHISTE)
C*****************************************************************
C***  Output of GAMMA HISTORY
C*****************************************************************

      DIMENSION GAHIST(26,MAXGAHIST), AG(6)
      CHARACTER*10 LEVEL(N), NAME*2, ACTLEVEL(4), NAME2*1
      LOGICAL BAG, BPGAHISTE

      WRITE (*,*) 
      WRITE (*,*) 
      IF (BPGAHISTE) THEN
        WRITE (*,*) '          G A M M A   H I S T O R Y   Extended'
        WRITE (*,*) '          ===================================='
      ELSE
        WRITE (*,*) '          G A M M A   H I S T O R Y'
        WRITE (*,*) '          ========================='
      ENDIF
      WRITE (*,*) 
      IF (BAG) THEN
        WRITE (*,*) 'Automatic GAMMA Adjustment Active'
        WRITE (*,'(A,2(1X,F6.1),A,2(1X,F6.1))') 
     >    'Tresh:', AG(3), AG(4), '   Gamma:', AG(1), AG(2)
        WRITE (*,*) 
      ENDIF
      WRITE (*,'(A4, 3X, A, A,1X, A, A, 2X,A)')
     >  'JOB', '-------GAMMA(CLRD)-----', 
     >  ' --MAX--', 'T: L, MAX--', 
     >' ------ CORRECTIONS: L, Level, Corr ---------------------', 
     >  'N: A, M,  C, D'
      DO I=MAXGAHIST, 1, -1
        IF (GAHIST(1,I) .LT. 0. .OR. 
     >      (.NOT. BPGAHISTE .AND. I .GT. 5)) CYCLE
        IF (GAHIST(2,I) .EQ. 0.) THEN
          NAME = 'ST'
        ELSE
          NAME = 'EX'
        ENDIF
        IF (GAHIST(26,I) .EQ. 0.) THEN
          NAME2 = 'F'
        ELSE
          NAME2 = 'T'
        ENDIF
        IF (GAHIST(19,I) .GT. 0.) THEN
          CORRMAX = ALOG10(GAHIST(19,I))
        ELSE 
          CORRMAX = -999.
        ENDIF
        WRITE (*,'(I3,1X,A2,1X,4(F5.0,1X),F7.2,1X,A1,1X,I2,1X,F6.3,1X,
     >         4(I2,1X,I3,1X,F6.2,1X),
     >         F6.2,1X,I2,1X,I3,1X,I2)') 
     >    INT(GAHIST(1,I)), NAME, 
     >    GAHIST(3,I), GAHIST(4,I), GAHIST(5,I), GAHIST(6,I), 
     >    CORRMAX, NAME2, 
     >    INT(GAHIST(20,I)), GAHIST(21,I), 
     >    INT(GAHIST(7,I)),  INT(GAHIST(8,I)), GAHIST(9,I), 
     >    INT(GAHIST(10,I)), INT(GAHIST(11,I)), GAHIST(12,I), 
     >    INT(GAHIST(13,I)), INT(GAHIST(14,I)), GAHIST(15,I), 
     >    INT(GAHIST(16,I)), INT(GAHIST(17,I)), GAHIST(18,I), 
     >    GAHIST(22,I), INT(GAHIST(23,I)), 
     >    INT(GAHIST(24,I)), INT(GAHIST(25,I))
      ENDDO

      RETURN
      END
