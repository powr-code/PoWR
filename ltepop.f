      SUBROUTINE LTEPOP (N,ENLTE,TL,ENE,WEIGHT,NCHARG,EION,ELEVEL,NOM,
     $                  ABXYZ,NFIRST,NLAST,NATOM)
C***********************************************************************
C***  POPULATION NUMBERS IN THERMODYNAMIC EQUILIBRIUM FOR ALL ELEMENTS
C***********************************************************************
 
      DIMENSION ENLTE(N), WEIGHT(N),NCHARG(N),EION(N),ELEVEL(N),NOM(N)
      DIMENSION ABXYZ(NATOM),NFIRST(NATOM),NLAST(NATOM)

C***  C1 = H * C / K    ( CM * KELVIN )
      DATA C1/1.4388/
C***  C2 = FACTOR IN SAHA EQ.   ( C.F. MIHALAS P.113 )
      DATA C2/2.07E-16/
 
      T32=TL*SQRT(TL)
C***  LOOP FOR EACH ELEMENT  -------------------------------------------
      DO 9 NA=1,NATOM
      NFIRNA=NFIRST(NA)
      NLANA=NLAST(NA)
      ENLTE(NFIRNA)=1.
      DO 1 J=NFIRNA+1,NLANA
      IF(NOM(J) .NE. NOM(NFIRNA)) STOP'LTEPOP:WRONG ELEMENT MEMBERSHIP'
      IF (NCHARG(J) .EQ. NCHARG(J-1) ) THEN
C***     BOLTZMANN FACTOR
         ENLTE(J) = EXP(C1*(ELEVEL(J-1)-ELEVEL(J))/TL) 
     >             * WEIGHT(J)/WEIGHT(J-1) * ENLTE(J-1)

      ELSE IF (NCHARG(J) .EQ. NCHARG(J-1)+1 ) THEN
C***     SAHA FACTOR
         ENLTE(J) = EXP(C1*(ELEVEL(J-1)-ELEVEL(J)-EION(J-1))/TL) 
     >       * T32/ENE/C2 * WEIGHT(J)/WEIGHT(J-1) * ENLTE(J-1)
      ELSE
         STOP 'LTEPOP: INVALID CHARGE DIFFERENCE'
      ENDIF
    1 CONTINUE
 
C***  NORMALIZATION
      SUM=.0
      DO 2 J=NFIRNA,NLANA
        ENLTE(J) = MAX(1.E-200,ENLTE(J))
    2   SUM = SUM+ENLTE(J)
      SUM = SUM/ABXYZ(NA)
      DO 3 J=NFIRNA,NLANA
    3   ENLTE(J) = ENLTE(J)/SUM
 
    9 CONTINUE
C***  ENDLOOP  ---------------------------------------------------------
 
      RETURN
      END
