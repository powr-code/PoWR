      SUBROUTINE INVTRI (A,B,C,Q,N)
C***********************************************************************
C***  TRIDIAGONALE MATRIX   -A(L)   B(L)  -C(L)
C***  RECHTE SEITE Q(L)
C***  LOESUNG AUF VEKTOR Q(L)
C***  ACHTUNG -- AUCH C(L) WIRD VERAENDERT --
C***********************************************************************

      DIMENSION A(N),B(N),C(N),Q(N)
      DIMENSION AHELP(200),BHELP(200)

      IF (N .GT. 200) STOP '*** ERROR: dimension overflow in INVTRI'
      CI=C(1)/B(1)
      C(1)=CI
      QI=Q(1)/B(1)
      Q(1)=QI
      NM=N-1
      IF (N.EQ.1) RETURN
      IF (N.EQ.2) GOTO 3
      DO 1 I=2,NM
      AI=A(I)
      H=B(I)-AI*CI
      CI=C(I)/H
      C(I)=CI
      QI=(Q(I)+QI*AI)/H
    1 Q(I)=QI
    3 QI=(Q(N)+QI*A(N))/(B(N)-CI*A(N))
      Q(N)=QI
C**  THE BACKWARD ELIMINATION MAY BE SPEEDED UP BY CRAY VECTOR ROUTINE FOLR
      DO 4 L=2,N
      AHELP(L)=-C(N+1-L)
    4 BHELP(L)= Q(N+1-L)
      BHELP(1)=Q(N)
      CALL FOLR(N,AHELP,1,BHELP,1)
      DO 5 L=1,N
    5 Q(N+1-L)=BHELP(L)
      RETURN
C***  DEAD BRANCH: RECURSION BY HAND (CAN NOT BE AUTO-VECTORIZED)
      DO 2 I=1,NM
      L=N-I
      QI=Q(L)+C(L)*QI
    2 Q(L)=QI
      RETURN
      END
