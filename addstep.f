      SUBROUTINE ADDSTEP (XNEW,I,LTOT,ZRAY,XCMF,OPARAY,ETARAY,ETACRAY,
     $                   OPALRAY,ETALRAY,S,SC,RRAY,PJPJ,NDADDIM,NBLINE,
     >                   POROLENGTHRAY)
C***********************************************************************
C***  CALLED FROM OBSFRAM
C***  INSERTION OF AN ADDITIONAL DEPTH POINT INTO A RAY
C***********************************************************************
 
      DIMENSION ZRAY(LTOT),XCMF(LTOT),RRAY(LTOT)
      DIMENSION OPARAY(LTOT),ETARAY(LTOT),ETACRAY(LTOT)
      DIMENSION S(LTOT),SC(LTOT)
      DIMENSION OPALRAY(NDADDIM,NBLINE),ETALRAY(NDADDIM,NBLINE)
      DIMENSION POROLENGTHRAY(LTOT)
      DATA EPS / 1.E-10 /
      SAVE EPS
 
C***  FIND INDEX I SO THAT XCMF(I-1) .LT. XNEW .LT. XCMF(I)
      I = ISRCHFGT (LTOT-1,XCMF,1,XNEW)

C***  CUBIC INTERPOLATION OF ZNEW(XNEW)
      CALL CUBIC (I,LTOT,ZRAY,XCMF,P1,P2,P3,P4)

C***  EVALUATION OF THE CUBIC INTERPOLATION FOR Z
      DXA = XNEW-XCMF(I-1)
      DXB = XCMF(I)-XNEW
      DXA3 = DXA*DXA*DXA
      DXB3 = DXB*DXB*DXB
      ZNEW = P1*DXA3 + P2*DXA + P3*DXB3 + P4*DXB

C***  Suppress insertion if new point is non-monotonic or almost identical 
C***     with one interval endpoint. This is neccessary, as the cubic 
C***     interpolation does not guarrantee a monotonic function!
C***     W.-R. H., 17-Aug-1998
      IF (ZNEW .LE. ZRAY(I)+EPS) RETURN
      IF (ZNEW .GE. ZRAY(I-1)-EPS) THEN
         I=I-1
         RETURN
         ENDIF

      CALL SHIFT (ZRAY   ,I,LTOT)
      CALL SHIFT (RRAY   ,I,LTOT)
      CALL SHIFT (XCMF   ,I,LTOT)
      CALL SHIFT (OPARAY ,I,LTOT)
      CALL SHIFT (ETARAY ,I,LTOT)
      CALL SHIFT (ETACRAY,I,LTOT)
      CALL SHIFT (S      ,I,LTOT)
      CALL SHIFT (SC     ,I,LTOT)
      CALL SHIFT (POROLENGTHRAY,I,LTOT)
      DO 10 NBL=1,NBLINE
        CALL SHIFT (OPALRAY(1,NBL),I,LTOT)
        CALL SHIFT (ETALRAY(1,NBL),I,LTOT)
   10 CONTINUE
 
      ZRAY (I) = ZNEW
      XCMF(I)=XNEW
C***  INTERPOLATION OF OPACITIES ETC.: LINEAR IN RADIUS R
      RNEW = SQRT(PJPJ + ZRAY(I)*ZRAY(I))
      RRAY(I) = RNEW
      P = (RNEW-RRAY(I-1)) / (RRAY(I+1)-RRAY(I-1))
      Q = 1. - P
C***  INTERPOLATION WEIGHTS MODIFIED: LINEAR INTERPOLATION IN R**2
      Q = Q * RRAY(I-1) * RRAY(I-1) / (RNEW*RNEW)
      P = P * RRAY(I+1) * RRAY(I+1) / (RNEW*RNEW)
      OPARAY (I) = Q*OPARAY (I-1) + P*OPARAY (I+1)
      ETARAY (I) = Q*ETARAY (I-1) + P*ETARAY (I+1)
      ETACRAY(I) = Q*ETACRAY(I-1) + P*ETACRAY(I+1)
      S      (I) = Q*S      (I-1) + P*S      (I+1)
      SC     (I) = Q*SC     (I-1) + P*SC     (I+1)
      POROLENGTHRAY(I) = Q*POROLENGTHRAY(I-1) + P*POROLENGTHRAY(I+1)
      DO 20 NBL=1,NBLINE
      OPALRAY(I,NBL) = Q*OPALRAY(I-1,NBL) + P*OPALRAY(I+1,NBL)
      ETALRAY(I,NBL) = Q*ETALRAY(I-1,NBL) + P*ETALRAY(I+1,NBL)
   20 CONTINUE
 
      LTOT = LTOT + 1
      IF (LTOT .GT. NDADDIM) THEN
        WRITE (0,*) 'WARNING : LTOT .GT. NDADDIM'
        WRITE (*,*) 'WARNING : LTOT .GT. NDADDIM'
        WRITE (0,*) 'LTOT, NDADDIM=', LTOT, NDADDIM
        WRITE (*,*) 'LTOT, NDADDIM=', LTOT, NDADDIM
        STOP 'ERROR IN SUBR. ADDSTEP'
C!!!            8.1.96, LARS
      ENDIF

      RETURN
      END
