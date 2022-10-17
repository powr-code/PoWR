      SUBROUTINE CBBH (OMEGA,IND,NUP,LOW,TL,TROOT,T32,NDIM,N,NCHARG,
     $                 EINST,COCO,KEYCBB,WEIGHT,WAVENUM,WN2,WN3)
C***********************************************************************
C***  COLLISIONAL BOUND-BOUND TRANSITIONS OF HYDROGEN
C***  CALCULATION OF OMEGA(UP-LOW)
C***  OMEGA(UP-LOW) IS CALCULATED DEPENDING ON KEYWORD KEYCBB  *********
C***********************************************************************
 
      DIMENSION NCHARG(N),WEIGHT(N)
      DIMENSION EINST(NDIM,NDIM),COCO(4,IND)
      CHARACTER*4 KEYCBB(IND)
 
C***  C1 = H * C / K    ( CM * ANGSTROEM )
      DATA C1 / 1.4388 /
 
C***  H I  *************************************************************
      IF (NCHARG(LOW) .EQ. 0) THEN
C***  OMEGA (UP-LOW) ACCORDING TO JEFFERIES P. 118 (EQ. 6.24)
         IF (KEYCBB(IND) .EQ. '    ') THEN
C                              ====
          OMEGA=3.24*EINST(NUP,LOW)/WN2/T32/(C1*WAVENUM/TL)**1.68
         ELSE
C***  NO/UNKNOWN KEYWORD 'KEYCBB' DECODED
          CALL REMARK (' H I: WRONG KEYWORD KEYCBB ')
          STOP 'ERROR'
         ENDIF
 
C***  H II  ************************************************************
      ELSE
          CALL REMARK ('H II - LINE TRANSITION')
          STOP 'ERROR'
      ENDIF
 
      RETURN
      END
