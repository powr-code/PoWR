      SUBROUTINE LIPO (F,X,FI,XI,N)
C***********************************************************************
C***  LINEAR INTERPOLATION: FOR GIVEN X, FIND F(X) FROM A TABLE FI(XI)
C ------ INDEXSUCHE DURCH BISECTION ------
C***********************************************************************

      DIMENSION FI(N),XI(N)

      NA=1
      A=XI(1)
      NB=N
      B=XI(N)
      IF((X-A)*(X-B).GT..0) STOP 'ERROR IN SUBR. LIPO'
    1 IF((NB-NA).EQ.1) GOTO 2
      NH=(NA+NB)/2
      H=XI(NH)
      IF((X-A)*(X-H).GT..0) GOTO 3
      NB=NH
      B=H
      GOTO 1
    3 NA=NH
      A=H
      GOTO 1
 
    2 P=(X-A)/(B-A)
      F=P*FI(NB)+(1.-P)*FI(NA)
      RETURN
      END
