      SUBROUTINE BFCROSS (SIGMAKI,NF,N,ELEVEL,EION,EINST,NDIM,
     $                    XLAMBDA,ALPHA,SEXPO,
     $                    ADDCON1, ADDCON2, ADDCON3, 
     $                    IGAUNT,
     $                    KONTNUP,KONTLOW,LASTKON)
C***********************************************************************
C***  THIS ROUTINE PREPARES AN ARRAY SIGMAKI WITH THE BOUND-FREE CROSS SECTIONS
C***  ( IN CM**2) TO AVOID UNNECCESSARY MULTIPLE CALCULATIONS
C***********************************************************************
 
      DIMENSION EION(NDIM),ELEVEL(NDIM),EINST(NDIM,NDIM)
      DIMENSION KONTNUP(LASTKON),KONTLOW(LASTKON)
      DIMENSION XLAMBDA(NF)
      DIMENSION SIGMAKI(NF,LASTKON)
 
C***  LOOP OVER ALL CONTINUUM TRANSITIONS
      DO 8 KON=1,LASTKON
      NUP=KONTNUP(KON)
      LOW=KONTLOW(KON)
C***  EINST = THRESHOLD CROSS SECTION IN 10**-18 CM**2
      SIGMATH=EINST(LOW,NUP)*1.E-18
C***  EDGE = THRESHOLD ENERGY IN KAYSER *****
      EDGE=ELEVEL(NUP)+EION(LOW)-ELEVEL(LOW)

C***  LOOP OVER ALL CONTINUUM FREQUENCY POINTS
      DO 9 K=1,NF
      WAVENUM=1.E8/XLAMBDA(K)
      IF (WAVENUM .LT. EDGE) THEN
         SIGMAKI(K,KON)=.0
      ELSE
         CALL PHOTOCS (SIGMAKI(K,KON),SIGMATH,EDGE,WAVENUM,ALPHA,
     $                 SEXPO,
     $                 ADDCON1, ADDCON2, ADDCON3, 
     $                 IGAUNT,KON)
      ENDIF
    9 CONTINUE

    8 CONTINUE
 
      RETURN
      END
