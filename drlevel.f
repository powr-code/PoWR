      SUBROUTINE DRLEVEL (N, NDIM, MAXIND, MAXAUTO, NAUTO, KRUDAUT, 
     $                    LOWAUTO, IONAUTO, AAUTO, EAUTO, ELEVEL, 
     $                    LEVEL, EINST, EION, WEIGHT, NCHARG, INDLOW, 
     $                    INDNUP, LASTIND, NHIGH, LAINDHI,
     $                    ND, T, ENTOT, RNE, POPNUM, NOM, DENSCON)
C******************************************************************************
C***  COLLECTS DATA FOR UPPER DR-LEVELS (READS FROM TAPE4=DATOM),
C***  ENLARGES THE ARRAYS INDICATED BY LEVEL- OR LINENUMBERS AND
C***  CALCULATES THE SAHA-POPNUMBERS OF THE UPPER DR-LEVELS
C***  CALLED ONLY FROM CMF
C******************************************************************************

cccccccccccccccccccccccccccccccccccccccc
ccc   Diese Subroutine fuehrt zum Absturz wegen WEIGHT(NION) .eq. 0
ccc   Ueberhaupt scheint mir das Programm ziemlicher Quatsch zu sein.
ccc   Sehr unschoen, dass File DATOM erneut gelesen wird!! 
ccc   Vermutlich ist die Bestimmung des naechsten Grundniveaus hier
ccc   falsch. 
ccc   Die Absicht dieser Routine ist offenbar, die DRTRANSITS konsistent im
ccc   Strahlungstransport zu beruecksichtigen.
ccc   Als Notmassnahme habe ich versucht, die ganze Routine einfach 
ccc   abzuschalten, was aber zu Folgefehlern fuehrte. 
ccc   wrh 25-Sep-2010 14:30:15
ccccccccccccccccccccccccccccccccccc


      DIMENSION NCHARG (NDIM),WEIGHT(NDIM),ELEVEL(NDIM)
      DIMENSION EION(NDIM),EINST(NDIM,NDIM),NOM(NDIM)
      DIMENSION INDNUP(MAXIND),INDLOW(MAXIND)
      DIMENSION LOWAUTO(MAXAUTO),IONAUTO(MAXAUTO)
      DIMENSION AAUTO(MAXAUTO),EAUTO(MAXAUTO),KRUDAUT(MAXAUTO)
      DIMENSION T(ND),ENTOT(ND),RNE(ND),POPNUM(ND,NDIM)
      CHARACTER KARTE*80
      CHARACTER*10 LEVEL(NDIM),LEVUP,LEVLOW, IDUMMY
      CHARACTER*3 RDUMMY
      DIMENSION DENSCON(ND)

C***  CI : FACTOR IN SAHA EQUATION (MIHALAS, P. 113)
      DATA CI / 2.07E-16 /
C***  C1 = H * C / K    ( CM*KELVIN )
      DATA C1 / 1.4388 /

      NHIGH = N + 1
      LEVONE = N + 1
      INDONE = LASTIND + 1
      INDDR = 0

      OPEN (4, FILE='DATOM', STATUS='UNKNOWN')
 
C  READING INPUT DATA FROM DATOM-FILE  ------------------------
***************************************************************

    1 READ(4,2,END=3) KARTE
    2 FORMAT(A)

      IF (KARTE(:10) .NE. 'DRTRANSIT ' ) GOTO 1

      INDDR = INDDR + 1 

      IF (NHIGH .GT. NDIM) THEN
         CALL REMARK ('TOO MANY UPPER DR-LEVELS (NDIM TOO SMALL)')
         STOP 'DRLEVEL'
         ENDIF   

      IF (INDDR .GT. MAXAUTO) THEN
         CALL REMARK ('TOO MANY DR-TRANSITIONS (MAXAUTO TOO SMALL)')
         STOP 'DRLEVEL'
         ENDIF   

      READ (KARTE,10) LEVLOW,LEVUP,NW,EUP,ADUMMY,IDUMMY,RDUMMY
   10 FORMAT(10X,A10,2X,A10,1X,I4,1X,F10.0,1X,G10.0,2X,A10,1X,A3)

C *****     TESTING, IF LOWER LEVEL-NAME AND UPPER ENERGY IS UNCHANGED
**********************************************************************

      IF ((LEVEL(LOWAUTO(INDDR)) .NE. LEVLOW) .OR. 
     $             (EAUTO(INDDR) .NE. EUP)) THEN
         CALL REMARK ('CHANGED DR-TRANSITIONS IN DATOM-FILE')
         STOP 'DRLEVEL'
         ENDIF   

C ******   SET INDEX, IF NEW ION OR NEW ELEMENT IS DETECTED
***********************************************************

      IF (INDDR .GT. 1) THEN
         IF ((NOM(LOWAUTO(INDDR)) .NE. NOM(LOWAUTO(INDDR-1))) .OR.
     >      (NCHARG(LOWAUTO(INDDR)) .NE. NCHARG(LOWAUTO(INDDR-1)))) 
     >      THEN
              LEVONE = NHIGH
              INDONE = LASTIND + INDDR
            ENDIF
         ENDIF

C*****   ENLARGING THE ARRAYS INDICATED BY LINE NUMBERS (INDDR)
C**************************************************************

      INDLOW(LASTIND+INDDR) = LOWAUTO(INDDR)
      EUPGES = EUP + EION(LOWAUTO(INDDR))

C*******  NEW UPPER DR-LEVEL ?
C*****************************

      DO 20 I = LEVONE,NHIGH-1
         IF (LEVUP .EQ. LEVEL(I)) THEN
            INDNUP(LASTIND+INDDR) = I
            EINST(I,LOWAUTO(INDDR)) = AAUTO(INDDR)
            IF (KRUDAUT(INDDR) .EQ. 1) EINST(LOWAUTO(INDDR),I) = -2.
            GOTO 1
            ENDIF
   20 CONTINUE

C*******    IF YES, SEARCHING FOR THE RIGHT POSITION (INCREASING ENERGIES)
C*******               AND ENLARGING THE ARRAYS INDICATED BY LEVEL-NUMBERS
C*************************************************************************

      DO 21 K = LEVONE,NHIGH-1
         IF (EUPGES .LT. ELEVEL(K)) THEN
            JK = K
            CALL SHIFT(ELEVEL,K,NHIGH-1)
            CALL SHIFT(WEIGHT,K,NHIGH-1)
            CALL SHIFT(NCHARG,K,NHIGH-1)
            CALL SHIFT(EION,K,NHIGH-1)
            DO 23 LN = 1,ND
             DO 23 LK = NHIGH,JK+1,-1
   23       POPNUM(LN,LK) = POPNUM(LN,LK-1)
            DO 25 LK1 = NHIGH,JK+1,-1
             DO 25 LK2 = 1,NHIGH
             EINST(LK1,LK2) = EINST(LK1-1,LK2)
             EINST(LK2,LK1) = EINST(LK2,LK1-1)
   25       CONTINUE
            DO 26 LK = NHIGH,JK+1,-1
   26       LEVEL(LK) = LEVEL(LK-1)
            DO 27 LIND =INDONE,LASTIND+INDDR-1
   27       IF (INDNUP(LIND) .GE. JK) INDNUP(LIND) = INDNUP(LIND)+1
            GOTO 22
            ENDIF

   21 CONTINUE

      JK = NHIGH

   22 CONTINUE

C  ******** PRESET FOR NEW EINST-VALUES  = -11.
***********************************************
      DO 28 J=1,NHIGH
      EINST(JK,J) = -11.
      EINST(J,JK) = -11.
   28 CONTINUE

C ***   NEW VALUES ARE INSERTED
*******************************

      INDNUP(LASTIND+INDDR) = JK
      LEVEL(JK) = LEVUP
      EINST(JK,LOWAUTO(INDDR)) = AAUTO(INDDR)
      IF (KRUDAUT(INDDR) .EQ. 1) EINST(LOWAUTO(INDDR),JK) = -2.
      ELEVEL(JK) = EUPGES
      WEIGHT(JK) = NW
      EION(JK) = EION(LOWAUTO(INDDR))
      NCHARG(JK) = NCHARG(LOWAUTO(INDDR))

C*****    CALCULATING THE LTE POPNUMBER OF THE NEW UPPER DR-LEVEL JK
********************************************************************

      DO 30 L=1,ND
        TL = T(L)
        SQRTL = SQRT(TL)
C***  Density increased by DENSCON
        ENTOTL = ENTOT(L)*DENSCON(L)
        RNEL = RNE(L)
        NION = IONAUTO(INDDR)
        POPNUP = POPNUM(L,NION)
ccc     test stop !!!
        if (WEIGHT(NION) .eq. 0) then
           write (0,*) 'WEIGHT(NION) .eq. 0'
           write (0,*) 'INDDR =', inddr
           write (0,*) 'NION  =', nion
           stop '*** teststop in SUBROUTINE DRLEVEL'
        endif


        SAHAPHI = WEIGHT(JK) / WEIGHT(NION) * CI / TL / SQRTL
     *           * EXP(-C1 * (EUP-ELEVEL(NION)) / TL)
        POPNUM(L,JK) = POPNUP * ENTOTL * RNEL * SAHAPHI

   30 CONTINUE

      NHIGH = NHIGH + 1

C ***    READING NEXT DATA
**************************

      GOTO 1


C***  END OF DR-DATA REACHED  ---------------------------------------
*********************************************************************

    3 CLOSE (4)

      IF (NAUTO .NE. INDDR) THEN
         CALL REMARK ('CHANGED NUMBER OF DR-TRANSITIONS IN DATOM-FILE')
         STOP 'DRLEVEL'
         ENDIF   

      NHIGH = NHIGH - 1

      RETURN
      END 
