      SUBROUTINE CBBMORE(OMEGA,IND,NUP,LOW,TL,TROOT,T32,NDIM,N,NCHARG,
     $                 EINST,COCO,KEYCBB,WEIGHT,WAVENUM,WN2,WN3)
C***********************************************************************
C***  COLLISIONAL BOUND-BOUND TRANSITIONS OF CARBON, OXYGEN, NEON
C***  MAGNESIUM, ALUMINIUM AND SILICON
C***  CALCULATION OF OMEGA(UP-LOW)
C***  OMEGA(UP-LOW) IS CALCULATED DEPENDING ON KEYWORD KEYCBB  *********
C***********************************************************************

      DIMENSION NCHARG(N),WEIGHT(N)
      DIMENSION EINST(NDIM,NDIM),COCO(4,IND)
      CHARACTER*4 KEYCBB(IND)
      
      REAL, EXTERNAL :: EXPINT1EXP

C***  C1 = H * C / K    ( CM * KELVIN )
      DATA C1 / 1.4388 /

C***  ONLY GLOBAL KEYWORDS  *************************************************
C***  'JEFF': OPTICALLY PERMITTED  TRANSITIONS. JEFFERIES P. 119 (EQ. 6.25)
      IF (KEYCBB(IND) .EQ. 'JEFF') THEN
C                           ====
              OMEGA=4.06*EINST(NUP,LOW)/WN3/TROOT


C***  'KB..': DATA FROM PROGRAM "DETAIL" OF K.BUTLER (MUNICH)
C***          '..22': FORMULA 22, I.E. THE VAN REGEMORTER (1962, APJ 136, 906)
C***                              FORMULA WITH  G = COCO(1,IND),
C***                                        GAMMA = MAX(G,.276*EXP(U0)*E1(U0))
C***                                           (E1: FIRST EXPONENTIAL INTEGRAL)
      ELSE IF (KEYCBB(IND) .EQ. 'KB22') THEN
C                                ====
              G=COCO(1,IND)
              U0=C1*WAVENUM/TL
C***          prevent overflow for very high energies: restrict U0
C              EXPU0 = EXP(AMIN1(U0,500.)) !no longer needed due to new function EXPINT1EXP
              GAMMA=AMAX1(G, 0.276 * EXPINT1EXP(U0))
              OMEGA=20.56*GAMMA*EINST(NUP,LOW)/WN3/TROOT


C***          '..24': FORMULA 24, I.E. A SEMI-EMPIRICAL CROSS SECTION
C***                              (CONSTANT UPSILON) FROM ALLEN (1973),
C***                              ASTROPHYSICAL QUANTITIES
      ELSE IF (KEYCBB(IND) .EQ. 'KB24') THEN
C                                ====
              OMEGA=8.6287E-6*COCO(1,IND)/WEIGHT(NUP)/TROOT


C***  WRONG KEYWORD KEYCBB

      ELSE
              CALL REMARK ('CBBMORE: WRONG KEYWORD KEYCBB ')
              STOP 'ERROR'

      ENDIF


      RETURN
      END
