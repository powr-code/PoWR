      SUBROUTINE BRTPUP (AINV, S, SP, Y, SCRATCH, SCVEC, N, NDIM, 
     >                   TWOPNT, BRDIVFAILED)
C********************************************************************
C***  UPDATING OF THE INVERSE DERIVATIVE MATRIX FOR THE NEWTON-LIKE
C***    ITERATION, USING THE BROYDEN ALGORITHM.
C***  DIE AUFDATIERUNG ERFOLGT MIT DER INVERSEN AUDATIERUNGSFORMEL
C***    VON BROYDEN (KOSMOL SEITE 147) TRANSPONED FOR ROW VEKTORS
C***  AINV = INVERSE MATRIX 
C***  S    = X(K+1) - X(K)  (LAST-STEP CORRECTIONS)
C***  Y    = F[X(K+1)] - F[X(K)]
C********************************************************************

      LOGICAL TWOPNT, BRDIVFAILED

      DIMENSION SCVEC(1), S(1), SP(1)

      IF (TWOPNT) THEN
C***  CALCULATE V(K)
        FAC = SDOT (N,SP,1,S,1) / SDOT (N,SP,1,SP,1)
        DO 1, I=1, N
          SCVEC(I) = S(I) - FAC * SP(I)
1       CONTINUE
C***  CALCULATE THETA
        THETA = ABS (SDOT (N,S,1,SCVEC,1)) / 
     >          ( SQRT(SDOT(N,S,1,S,1)) * SQRT(SDOT(N,SCVEC,1,SCVEC,1)))
        IF (THETA .LT. 0.0175) THEN
C!!!        IF (THETA .LT. 0.0) THEN
          TWOPNT = .FALSE.
          DO 2, I=1, N
            SP(I) = S(I)
2         CONTINUE
        ELSE
          DO 3, I=1, N
            SP(I) = SCVEC(I)
3         CONTINUE
        ENDIF
      ENDIF

C***  Compute the term y(k) * A(k) ( row-Vektor times Matrix ) -> Y
C***  TM is used as workspace
      CALL BRVM (Y, AINV, SCRATCH, N, NDIM)

C***  Compute the denominator of broyden's update formula
C***        RNENN <-- y(k) * A(k) * s(k)
      RNENN=SDOT(N,SP,1,Y,1)

C***  Calculate the negativ of the term
C***  (s(k)-y(k)*A(k)) of Broyden's formula -> Y
      CALL VSUB (Y, S, N)

C***  Divide Vector Y by RNENN
      BRDIVFAILED = .FALSE.
      CALL BRVDIVS (Y, RNENN, N, BRDIVFAILED)

C***  Multiply A(k) by s(k)  (Matrix by column vektor)
      CALL BRMV (SP, AINV, SCRATCH, N, NDIM)

C***  Calculate the dyadic product of S and Y (column by row Vector),
C***   and substract the result from matrix AINV
      CALL BRVVDY (AINV, SP, Y, N, NDIM)

      RETURN
      END
