      FUNCTION ISRCHFGE(N,X,INCX,TARGET)

C***  NAME
C***       ISRCHFGE - Searches a vector for the first element greater 
C***       than or equal to a target

C***  SYNOPSIS
C***  SEE ISRCHEQ

      REAL X(*), TARGET

      J=1
      ISRCHFGE=0
      IF(N.LE.0) RETURN
      IF(INCX.LT.0) J=1-(N-1)*INCX
      DO 100 I=1,N
        IF(X(J).GE.TARGET) GOTO 200
          J=J+INCX
  100 CONTINUE
  200 ISRCHFGE=I

      RETURN
      END

