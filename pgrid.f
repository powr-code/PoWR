      SUBROUTINE PGRID (NPDIM,NP,ND,R,P, NC)
C***********************************************************************
C***  GRID OF IMPACT-PARAMETER POINTS
C***********************************************************************

      DIMENSION  P(NPDIM),R(ND)

C***  NC = NUMBER OF CORE-INTERSECTING RAYS
      NP=ND+NC
      IF (NP.GT.NPDIM) STOP 'PGRID: TOO MANY IMPACT-PARAMETER POINTS'

C***  CORE RAYS EQUALLY SPACED
      D=1./FLOAT(NC)
      DO 1 J=1,NC
    1 P(J)=(J-1)*D
      DO 2 L=1,ND
      J=NP+1-L
    2 P(J)=R(L)

      RETURN
      END
