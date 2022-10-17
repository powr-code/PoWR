      SUBROUTINE CBBHE (OMEGA,IND,NUP,LOW,TL,TROOT,T32,NDIM,N,NCHARG,
     $                  EINST,COCO,KEYCBB,WEIGHT,WAVENUM,WN2,WN3)
C***********************************************************************
C***  COLLISIONAL BOUND-BOUND TRANSITIONS OF HELIUM
C***  CALCULATION OF OMEGA(UP-LOW)
C***  OMEGA(UP-LOW) IS CALCULATED DEPENDING ON KEYWORD KEYCBB  *********
C***********************************************************************
 
      DIMENSION NCHARG(N),WEIGHT(N)
      DIMENSION EINST(NDIM,NDIM),COCO(4,IND)
      CHARACTER*4 KEYCBB(IND)
 
C***  C1 = H * C / K    ( CM * ANGSTROEM )
      DATA C1 / 1.4388 /
 
C***  HE I  *******************************************************************
      IF (NCHARG(LOW).EQ.0) THEN
C***  'JEFF': OPTICALLY PERMITTED  TRANSITIONS. JEFFERIES P. 118 (EQ. 6.24)
         IF (KEYCBB(IND) .EQ. 'JEFF') THEN
C                              ====
              OMEGA=3.24*EINST(NUP,LOW)/WN2/T32/(C1*WAVENUM/TL)**1.68
 
C***  'BFK ': OPTICALLY FORBIDDEN TRANSITIONS BETWEEN N=2, N=1 OR WITHIN N=2: 
C***          BERRINGTON, FON + KINGSTON ('BFK') 1982, MNRAS 200, 347
         ELSE IF (KEYCBB(IND) .EQ. 'BFK1') THEN
C                                   ====
              PBFK=COCO(1,IND)/TROOT
     +                          +COCO(2,IND)+COCO(3,IND)*TROOT
     +                          +COCO(4,IND)*T32
              OMEGA=PBFK*WEIGHT(LOW)/WEIGHT(NUP)
 
         ELSE IF (KEYCBB(IND) .EQ. 'BFK2') THEN
C                                   ====
              PBFK=COCO(1,IND)/TL
     +                          +COCO(2,IND)/TROOT+COCO(3,IND)
     +                          +COCO(4,IND)*TL
              OMEGA=PBFK*WEIGHT(LOW)/WEIGHT(NUP)
 
C***  'BKMS' OR 'BKGR': 
C***            FORBIDDEN TRANSITIONS BETWEEN N.GT.2, N=1 OR N.GT.3, N.NE.1: 
C***        BENSON + KULANDER ('BK..') 1972, SOLAR PHYSICS 27, 305 (FORMULA 3)
C***        '..MS': MIHALAS + STONE (REF. 11)
C***        '..GR': GREEN (REF. 7)
C***        ATTENTION: COCO(3,IND) := 1.-ALPHA
          ELSE IF (KEYCBB(IND) .EQ. 'BKMS'
C                                    ====
     $        .OR. KEYCBB(IND) .EQ. 'BKGR') THEN
C                                    ====
              OMEGA=COCO(1,IND)*TL**COCO(2,IND)*
     *               EXP(COCO(3,IND)*C1*WAVENUM/TL)*WEIGHT(LOW)/
     /               WEIGHT(NUP)
 
C***  'UPS.': OPTICALLY FORBIDDEN TRANSITIONS BETWEEN N.GT.2, N.NE.1: 
C***        EFFECTIVE COLLISIONAL STRENGTH (UPSILON) FROM WERNER SCHMUTZ
          ELSE IF (KEYCBB(IND) .EQ. 'UPS0') THEN
C                                    ====
            OMEGA=.0
 
          ELSE IF (KEYCBB(IND) .EQ. 'UPS1') THEN
C                                    ====
              UPSILON=0.05
              OMEGA=8.6287E-6*UPSILON/WEIGHT(NUP)/TROOT
 
          ELSE IF (KEYCBB(IND) .EQ. 'UPS2') THEN
C                                    ====
              UPSILON=1.
              OMEGA=8.6287E-6*UPSILON/WEIGHT(NUP)/TROOT
 
          ELSE
C***  NO/UNKNOWN KEYWORD 'KEYCBB' DECODED
              CALL REMARK(' HE I: WRONG KEYWORD KEYCBB ')
              STOP 'ERROR'
          ENDIF
 
C***  HE II  ******************************************************************
      ELSE IF (NCHARG(LOW).EQ.1) THEN
C***  OMEGA (UP-LOW) ACCORDING TO JEFFERIES P. 119 (EQ. 6.25)
         IF (KEYCBB(IND) .EQ. '    ') THEN
C                              ====
          OMEGA=4.06*EINST(NUP,LOW)/WN3/TROOT
         ELSE
C***  NO/UNKNOWN KEYWORD 'KEYCBB' DECODED
          CALL REMARK (' HE II: WRONG KEYWORD KEYCBB ')
          STOP 'ERROR'
         ENDIF
 
C***  HE III  *****************************************************************
C      ELSE
C          CALL REMARK('  HE III - LINE TRANSITION  ')
C          STOP 'ERROR'
      ENDIF
C!!!      OMEGA=0.
 
      RETURN
      END
