
      SUBROUTINE FLGRID (NFLDIM,NFL,PHI,PWEIGHT,DELTAX,XMAX)
C***********************************************************************
C***  FREQUENCY SCALE ETC. FOR PROGRAM STEAL (CALLED FROM LINPOP)
C***  NFL          NUMBER OF FREQUENCY POINTS >>>>>> DEFINED IN "STEAL" <<<<<<
C***  XMAX         MAXIMUM FREQUENCY (DEFINED IN STEAL)
C***  XK           FREQUENCY POINTS (FALLING, DOPPLER UNITS) -  NOT STORED
C***  PHI(K)       NORMALIZED LINE PROFILE
C***  PWEIGHT(K)   RENORMALIZED INTEGRATION WEIGHTS
C***********************************************************************

      DIMENSION PHI(NFLDIM),PWEIGHT(NFLDIM)
C***  WPIINV = PI ** (-1/2)
      DATA WPIINV / 0.5641896 /
 
      IF (NFL .GT. NFLDIM) THEN
         CALL REMARK ('FLGRID: DIMENSION NFLDIM TOO SMALL')
         STOP 'ERROR'
         ENDIF

      DELTAX = 2. * XMAX / FLOAT(NFL-1)
      WS = .0
      DO 1 K=1, NFL
      XK = XMAX - (K-1) * DELTAX
      EX = EXP(-XK*XK)
      PWEIGHT(K) = EX
      WS = WS + EX
      PHI(K) = EX * WPIINV
    1 CONTINUE
 
C***  RENORMALIZATION
      DO 2 K=1, NFL
    2 PWEIGHT(K) = PWEIGHT(K) / WS
      RETURN
      END
