      SUBROUTINE DRDAT( NDIM, MAXNDR, MAXAUTO, LOWAUTO, EION, LEVEL,
     $                  NCHARG, NOM, IONGRND, NDR, ELEVDR, NCHARGDR, 
     $                  INDNUPDR, LEVELDR, NOMDR, IONDR, WEIGHTDR)
C******************************************************************************
C***  COLLECTS DATA FOR UPPER DR-LEVELS (READS FROM TAPE4=DATOM)
C***  ONLY FOR OUTPUT  (CALLED FROM SUBROUTINE PRIDAT)
C******************************************************************************
 
      DIMENSION NCHARGDR(MAXNDR),ELEVDR(MAXNDR),NOMDR(MAXNDR)
      DIMENSION WEIGHTDR(MAXNDR)
      DIMENSION NOM(NDIM),NCHARG(NDIM),EION(NDIM),IONGRND(NDIM)
      DIMENSION INDNUPDR(MAXAUTO),LOWAUTO(MAXAUTO)
      CHARACTER KARTE*80
      CHARACTER*10 LEVEL(NDIM)
      CHARACTER*10 IONDR(MAXNDR),LEVELDR(MAXNDR),LEVUP,LEVLOW,LEVION

      NDR = 1
      LEVONE = 1
      INDONE = 1
      INDDR = 0

      OPEN (4, FILE='DATOM', STATUS='UNKNOWN')
 
C  READING INPUT DATA FROM DATOM-FILE  ------------------------
***************************************************************

    1 READ(4,2,END=3) KARTE
    2 FORMAT(A)

      IF (KARTE(:10) .NE. 'DRTRANSIT ' ) GOTO 1

      INDDR = INDDR + 1 

      IF (NDR .GT. MAXNDR) THEN
         CALL REMARK ('TOO MANY UPPER DR-LEVELS (MAXNDR TOO SMALL)')
         STOP 'DRDAT'
         ENDIF   

      IF (INDDR .GT. MAXAUTO) THEN
         CALL REMARK ('TOO MANY DR-TRANSITIONS (MAXAUTO TOO SMALL)')
         STOP 'DRDAT'
         ENDIF   

      READ (KARTE,10) LEVLOW,LEVUP,NW,EUP,ADUMMY,LEVION
   10 FORMAT(10X,A10,2X,A10,1X,I4,1X,F10.0,1X,G10.0,2X,A10)


C ******   SET INDEX, IF NEW ION OR NEW ELEMENT IS DETECTED
***********************************************************

      IF (INDDR .GT. 1) THEN
         IF ((NOM(LOWAUTO(INDDR)) .NE. NOM(LOWAUTO(INDDR-1))) .OR.
     >      (NCHARG(LOWAUTO(INDDR)) .NE. NCHARG(LOWAUTO(INDDR-1)))) 
     >      THEN
              LEVONE = NDR
              INDONE = INDDR
            ENDIF
         ENDIF


      EUPGES = EUP + EION(LOWAUTO(INDDR))

C*******  NEW UPPER DR-LEVEL ?
C*****************************

      DO 20 I = LEVONE,NDR-1
         IF (LEVUP .EQ. LEVELDR(I)) THEN
            INDNUPDR(INDDR) = I
            GOTO 1
            ENDIF
   20 CONTINUE

C*******    IF YES, SEARCHING FOR THE RIGHT POSITION (INCREASING ENERGIES)
C*******               AND ENLARGING THE ARRAYS INDICATED BY LEVEL-NUMBERS
C*************************************************************************

      DO 21 K = LEVONE,NDR-1
         IF (EUPGES .LT. ELEVDR(K)) THEN
            JK = K
            CALL SHIFT(ELEVDR,K,NDR-1)
            CALL SHIFT(NCHARGDR,K,NDR-1)
            CALL SHIFT(NOMDR,K,NDR-1)
            CALL SHIFT(WEIGHTDR,K,NDR-1)
            DO 27 LK = NDR,JK+1,-1
            LEVELDR(LK) = LEVELDR(LK-1)
   27       IONDR(LK) = IONDR(LK-1)
            DO 28 LIND =INDONE,INDDR-1
   28       IF (INDNUPDR(LIND) .GE. JK) 
     >          INDNUPDR(LIND) = INDNUPDR(LIND)+1
            GOTO 22
            ENDIF

   21 CONTINUE

      JK = NDR

   22 CONTINUE


C ***   NEW VALUES ARE INSERTED
*******************************

      INDNUPDR(INDDR) = JK
      LEVELDR(JK) = LEVUP
      ELEVDR(JK) = EUPGES
      NCHARGDR(JK) = NCHARG(LOWAUTO(INDDR))
      NOMDR(JK) = NOM(LOWAUTO(INDDR))
      WEIGHTDR(JK) = NW
      IF (LEVION .EQ. '          ') THEN
         IONDR(JK) = LEVEL(IONGRND(LOWAUTO(INDDR)))
      ELSE
         IONDR(JK) = LEVION
      ENDIF

      NDR = NDR + 1

C ***    READING NEXT DATA
**************************

      GOTO 1


C***  END OF DR-DATA REACHED  ---------------------------------------
*********************************************************************

    3 CLOSE (4)

      NDR = NDR - 1

      RETURN
      END 
