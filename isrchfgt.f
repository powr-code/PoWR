      FUNCTION ISRCHFGT(N,X,INCX,TARGET)

C***  NAME
C***       ISRCHFGT - Searches a vector for the first element greater 
C***       to a target

C***  SYNOPSIS
C***  SEE ISRCHEQ

      REAL X(*), TARGET

      J=1
      ISRCHFGT=0
      IF(N.LE.0) RETURN
      IF(INCX.LT.0) J=1-(N-1)*INCX
      DO 100 I=1,N
        IF(X(J).GT.TARGET) GOTO 200
          J=J+INCX
  100 CONTINUE
  200 ISRCHFGT=I

      RETURN
      END

