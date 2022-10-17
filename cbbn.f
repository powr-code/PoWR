      SUBROUTINE CBBN (OMEGA,IND,NUP,LOW,TL,TROOT,T32,NDIM,N,NCHARG,
     $                 EINST,COCO,KEYCBB,WEIGHT,WAVENUM,WN2,WN3)
C***********************************************************************
C***  COLLISIONAL BOUND-BOUND TRANSITIONS OF NITROGEN
C***  CALCULATION OF OMEGA(UP-LOW)
C***  OMEGA(UP-LOW) IS CALCULATED DEPENDING ON KEYWORD KEYCBB  *********
C***********************************************************************

      DIMENSION NCHARG(N),WEIGHT(N)
      DIMENSION EINST(NDIM,NDIM),COCO(4,IND)
      CHARACTER*4 KEYCBB(IND)
      
      REAL, EXTERNAL :: EXPINT1EXP

C***  C1 = H * C / K    ( CM * KELVIN )
      DATA C1 / 1.4388 /
 
C-----------------------------------------------------------------------
C***  HANDLING OF FORMULA "KB23" (PROGRAM 'DETAIL', K. BUTLER)
      PARAMETER ( KB23MAX = 36 ) 
      DIMENSION AI(8,KB23MAX)
C***  TRANSITION 1:   N 32P2D2.3  -->  N III2P2.1
      DATA (AI(I,1),I=1,8) / 5.39184E-6, 3.58593E-7, 3.08826E-7, 0.,
     >      0., 0., 0., 0. /
C***  TRANSITION 2:   N 32P2S2.4  -->  N III2P2.1
      DATA (AI(I,2),I=1,8) / 1.02795E-5, 8.75551E-7, 8.80781E-7, 0.,
     >      0., 0., 0., 0. /
C***  TRANSITION 3:   N 32P2S2.4  -->  N 32P2D2.3
      DATA (AI(I,3),I=1,8) / 3.91835E-6, 6.99741E-10, 1.76365E-8, 0.,
     >      0., 0., 0., 0. /
C***  TRANSITION 4:   N 32P2P2.5  -->  N III2P2.1
      DATA (AI(I,4),I=1,8) / 9.87075E-6, 1.76246E-7, 1.08741E-6,
     >      9.87657E-7, 0., 0., 0., 0. /
C***  TRANSITION 5:   N 32P2P2.5  -->  N 32P2D2.3
      DATA (AI(I,5),I=1,8) / 3.48691E-6, -5.47237E-8, -1.68796E-8, 0.,
     >      0., 0., 0., 0. /
C***  TRANSITION 6:   N 32P2P2.5  -->  N 32P2S2.4
      DATA (AI(I,6),I=1,8) / 5.98409E-7, 1.0402E-8, -6.26268E-9, 0.,
     >      0., 0., 0., 0. /
C***  TRANSITION 7:   N 32P3D2.7  -->  N III2P2.1
      DATA (AI(I,7),I=1,8) / 4.07432E-7, 2.14008E-8, 1.14144E-8, 0.,
     >      0., 0., 0., 0. /
C***  TRANSITION 8:   N 32P3D2.7  -->  N 32P2D2.3
      DATA (AI(I,8),I=1,8) / 1.06453E-5, 1.00047E-6, 9.4256E-7, 0.,
     >      0., 0., 0., 0. /
C***  TRANSITION 9:   N 32P3D2.7  -->  N 32P2S2.4
      DATA (AI(I,9),I=1,8) / 4.08307E-8, -3.62028E-9, -2.98355E-9, 0.,
     >      0., 0., 0., 0. /
C***  TRANSITION 10:   N 32P3D2.7  -->  N 32P2P2.5
      DATA (AI(I,10),I=1,8) / 1.5764E-5, -3.20815E-6, -1.97302E-6,
     >      2.53603E-6, 2.50697E-6, 0., 0., 0. /
C***  TRANSITION 11:   N III3S2.8  -->  N III2P2.1
      DATA (AI(I,11),I=1,8) / 2.29679E-6, 5.65358E-7, 7.02663E-7,
     >      -5.05955E-7, -7.24263E-7, 0., 0., 0. /
C***  TRANSITION 12:   N III3S2.8  -->  N 32P2D2.3
      DATA (AI(I,12),I=1,8) / 3.90011E-7, 1.21137E-7, 1.50587E-7,
     >      -8.25036E-8, -1.30093E-7, 0., 0., 0. /
C***  TRANSITION 13:   N III3S2.8  -->  N 32P2S2.4
      DATA (AI(I,13),I=1,8) / 1.22032E-7, 5.79084E-8, 1.3954E-7,
     >      7.4935E-8, -6.72064E-8, -5.46334E-8, 0., 0. /
C***  TRANSITION 14:   N III3S2.8  -->  N 32P2P2.5
      DATA (AI(I,14),I=1,8) / 5.88388E-7, 3.84043E-7, 8.4958E-7,
     >      -1.94291E-7, -1.36535E-6, -7.83788E-8, 5.93208E-7, 0. /
C***  TRANSITION 15:   N III3S2.8  -->  N 32P3D2.7
      DATA (AI(I,15),I=1,8) / 2.62954E-8, 1.308E-8, 2.71565E-8,
     >      -7.41727E-9, -4.2968E-8, -1.20766E-9, 1.89689E-8, 0. /
C***  TRANSITION 16:   N 32P3P2.9  -->  N III2P2.1
      DATA (AI(I,16),I=1,8) / 1.12487E-6, 5.68246E-7, 6.72242E-7,
     >      -5.79867E-7, -1.25604E-6, 1.37004E-7, 5.71814E-7, 0. /
C***  TRANSITION 17:   N 32P3P2.9  -->  N 32P2D2.3
      DATA (AI(I,17),I=1,8) / 3.62802E-6, 1.85155E-6, 1.80305E-6,
     >      6.51599E-7, 0., 0., 0., 0. /
C***  TRANSITION 18:   N 32P3P2.9  -->  N 32P2S2.4
      DATA (AI(I,18),I=1,8) / 3.56698E-6, 1.42097E-7, 1.18743E-7, 0.,
     >      0., 0., 0., 0. /
C***  TRANSITION 19:   N 32P3P2.9  -->  N 32P2P2.5
      DATA (AI(I,19),I=1,8) / 1.3623E-5, 2.89568E-6, 2.99736E-6,
     >      1.20512E-6, 0., 0., 0., 0. /
C***  TRANSITION 20:   N 32P3P2.9  -->  N 32P3D2.7
      DATA (AI(I,20),I=1,8) / 3.13572E-6, 3.99969E-7, 2.99522E-7, 0.,
     >      0., 0., 0., 0. /
C***  TRANSITION 21:   N 32P3P2.9  -->  N III3S2.8
      DATA (AI(I,21),I=1,8) / 4.32887E-7, 6.86322E-7, 5.19394E-7,
     >      -5.40623E-7, -8.46581E-7, 1.94588E-7, 4.02574E-7, 0. /
C***  TRANSITION 22:   N III3P210  -->  N III2P2.1
      DATA (AI(I,22),I=1,8) / 2.20374E-6, 1.61833E-6, -7.87773E-7,
     >      -2.36562E-6, 1.53515E-7, 1.15393E-6, 0., 0. /
C***  TRANSITION 23:   N III3P210  -->  N 32P2D2.3
      DATA (AI(I,23),I=1,8) / 1.05007E-6, 4.08775E-7, -6.27183E-7,
     >      -8.59982E-7, 2.03092E-7, 4.1169E-7, 0., 0. /
C***  TRANSITION 24:   N III3P210  -->  N 32P2S2.4
      DATA (AI(I,24),I=1,8) / 1.52085E-7, 1.42043E-7, -1.54179E-7,
     >      -2.61798E-7, 1.2092E-7, 2.49035E-7, -4.16792E-8,
     >      -9.01979E-8 /
C***  TRANSITION 25:   N III3P210  -->  N 32P2P2.5
      DATA (AI(I,25),I=1,8) / 1.78494E-6, 7.45508E-7, -1.32394E-6,
     >      -1.93529E-6, 1.17438E-6, 1.01048E-6, -4.77597E-7, 0. /
C***  TRANSITION 26:   N III3P210  -->  N 32P3D2.7
      DATA (AI(I,26),I=1,8) / 1.32525E-7, 1.04535E-7, -1.1821E-7,
     >      -1.98249E-7, 9.03186E-8, 1.87273E-7, -3.26434E-8,
     >      -6.85308E-8 /
C***  TRANSITION 27:   N III3P210  -->  N III3S2.8
      DATA (AI(I,27),I=1,8) / 6.13693E-5, 8.17511E-6, 1.10851E-5,
     >      6.12946E-6, 0., 0., 0., 0. /
C***  TRANSITION 28:   N III3P210  -->  N 32P3P2.9
      DATA (AI(I,28),I=1,8) / 6.81862E-6, -2.21207E-6, -2.65653E-6,
     >      -1.72966E-7, 8.22934E-7, 0., 0., 0. /
C***  TRANSITION 29:   N III3D211  -->  N III2P2.1
      DATA (AI(I,29),I=1,8) / 1.31576E-6, 1.30853E-7, 1.26168E-7, 0.,
     >      0., 0., 0., 0. /
C***  TRANSITION 30:   N III3D211  -->  N 32P2D2.3
      DATA (AI(I,30),I=1,8) / 1.66095E-7, 1.66731E-8, 1.60763E-8, 0.,
     >      0., 0., 0., 0. /
C***  TRANSITION 31:   N III3D211  -->  N 32P2S2.4
      DATA (AI(I,31),I=1,8) / 2.66142E-8, 2.25754E-9, 3.53866E-9,
     >      2.16874E-9, 0., 0., 0., 0. /
C***  TRANSITION 32:   N III3D211  -->  N 32P2P2.5
      DATA (AI(I,32),I=1,8) / 6.71233E-8, -4.95973E-9, -4.78218E-9, 0.,
     >      0., 0., 0., 0. /
C***  TRANSITION 33:   N III3D211  -->  N 32P3D2.7
      DATA (AI(I,33),I=1,8) / 7.14224E-8, 7.80894E-9, 7.52939E-9, 0.,
     >      0., 0., 0., 0. /
C***  TRANSITION 34:   N III3D211  -->  N III3S2.8
      DATA (AI(I,34),I=1,8) / 7.95835E-6, 7.59789E-8, 7.3259E-8, 0.,
     >      0., 0., 0., 0. /
C***  TRANSITION 35:   N III3D211  -->  N 32P3P2.9
      DATA (AI(I,35),I=1,8) / 1.7438E-7, -8.23105E-9, -7.93639E-9, 0.,
     >      0., 0., 0., 0. /
C***  TRANSITION 36:   N III3D211  -->  N III3P210
      DATA (AI(I,36),I=1,8) / 7.74475E-5, 1.20402E-5, 1.88729E-5,
     >      1.15666E-5, 0., 0., 0., 0. /
C-----------------------------------------------------------------------

C***  GLOBAL KEYWORDS  *************************************************
C***  'JEFF': OPTICALLY PERMITTED  TRANSITIONS. JEFFERIES P. 119 (EQ. 6.25)
      IF (KEYCBB(IND) .EQ. 'JEFF') THEN
C                           ====
              OMEGA=4.06*EINST(NUP,LOW)/WN3/TROOT
              GOTO 99

C***  'KB..': DATA FROM PROGRAM "DETAIL" OF K.BUTLER (MUNICH)
C***          '..22': FORMULA 22, I.E. THE VAN REGEMORTER (1962, APJ 136, 906)
C***                              FORMULA WITH  G = COCO(1,IND),
C***                                        GAMMA = MAX(G,.276*EXP(U0)*E1(U0))
C***                                           (E1: FIRST EXPONENTIAL INTEGRAL)
      ELSE IF (KEYCBB(IND) .EQ. 'KB22') THEN
C                                ====
              G=COCO(1,IND)
              U0=C1*WAVENUM/TL
              GAMMA=AMAX1(G,0.276*EXPINT1EXP(U0))
              OMEGA=20.56*GAMMA*EINST(NUP,LOW)/WN3/TROOT
              GOTO 99

C***          '..24': FORMULA 24, I.E. A SEMI-EMPIRICAL CROSS SECTION
C***                              (CONSTANT UPSILON) FROM ALLEN (1973),
C***                              ASTROPHYSICAL QUANTITIES
      ELSE IF (KEYCBB(IND) .EQ. 'KB24') THEN
C                                ====
              OMEGA=8.6287E-6*COCO(1,IND)/WEIGHT(NUP)/TROOT
              GOTO 99
      ENDIF
      
C***  N II  ************************************************************
      IF (NCHARG(LOW) .EQ. 1) THEN
          CALL REMARK (' N II - LINE TRANSITION ')
          STOP 'ERROR'

C***  N III  ***********************************************************
      ELSE IF (NCHARG(LOW) .EQ. 2) THEN
C***  'KB23': FORMULA 23 FROM PROGRAM "DETAIL" (K. BUTLER, PRIVATE
C***          COMMUNICATION), I.E. A FIT FORMULA TO SOME UNPUBLISHED 
C***          RESULTS BY HUMMER (PRIVATE COMMUNICATION)
          IF (KEYCBB(IND) .EQ. 'KB23') THEN
C                               ====
              T4=ALOG10(TL/1.E4)
              KBL=IFIX(COCO(1,IND))
              GAMMA=AI(8,KBL)
              DO 9 I=7,1,-1
              GAMMA=AI(I,KBL)+GAMMA*T4
    9         CONTINUE
              OMEGA=GAMMA/TROOT

          ELSE
C***  NO/UNKNOWN KEYWORD 'KEYCBB' DECODED
              CALL REMARK (' N III: WRONG KEYWORD KEYCBB ')
              STOP 'ERROR'
          ENDIF
 
C***  N IV  ************************************************************
      ELSE IF (NCHARG(LOW).EQ.3) THEN
C***  'UPS.': OPTICALLY FORBIDDEN TRANSITIONS WITH CONSTANT EFFECTIVE
C***          COLLISIONAL STRENGTH (UPSILON)
          IF (KEYCBB(IND) .EQ. 'UPS0') THEN
C                               ====
              OMEGA=0.0
 
          ELSE IF (KEYCBB(IND) .EQ. 'UPS1') THEN
C                                    ====
              UPSILON=1.
              OMEGA=8.6287E-6*UPSILON/WEIGHT(NUP)/TROOT
 
          ELSE
C***  NO/UNKNOWN KEYWORD 'KEYCBB' DECODED
              CALL REMARK (' N IV: WRONG KEYWORD KEYCBB ')
              STOP 'ERROR'
          ENDIF
 
C***  N V  *************************************************************
      ELSE IF (NCHARG(LOW).EQ.4) THEN
C***  '.CMW': COLLISIONAL EXCITATION RATE COEFFICIENTS: 
C***          COCHRANE + MCWHIRTER (1983), PHYSICA SCRIPTA 28, 25-44
C***          'A...': ALLOWED TRANSITIONS
          IF (KEYCBB(IND) .EQ. 'ACMW') THEN
C                               ====
              GFIT=COCO(1,IND)+COCO(2,IND)*ALOG(TL/WAVENUM/C1+
     +             COCO(3,IND))
              OMEGA=20.56*GFIT*EINST(NUP,LOW)/WN3/TROOT
 
C***          'F...': NON-ALLOWED TRANSITIONS
          ELSE IF (KEYCBB(IND) .EQ. 'FCMW') THEN
C                                    ====
              GFFIT=COCO(1,IND)+COCO(2,IND)*ALOG((TL/WAVENUM/C1+
     +              COCO(3,IND))/(TL/WAVENUM/C1+COCO(4,IND)))
              OMEGA=13.71/WAVENUM/TROOT*GFFIT*WEIGHT(LOW)/WEIGHT(NUP)
 
          ELSE
C***  NO/UNKNOWN KEYWORD 'KEYCBB' DECODED
              CALL REMARK (' N V: WRONG KEYWORD KEYCBB ')
              STOP 'ERROR'
          ENDIF
 
C***  N VI  ************************************************************
      ELSE
          CALL REMARK (' N VI - LINE TRANSITION ')
          STOP 'ERROR'
      ENDIF
 
   99 CONTINUE

      RETURN
      END
