      FUNCTION ISRCHFLT(N,X,INCX,TARGET)

C***  NAME
C***     ISRCHFLT searches a real vector for the first element that is less
C***     than a real target.

C***  SYNOPSIS
C***  SEE ISRCHLT

      REAL X(*), TARGET

      J=1
      ISRCHFLT=0
      IF(N.LE.0) RETURN
      IF(INCX.LT.0) J=1-(N-1)*INCX
      DO 100 I=1,N
        IF(X(J).LT.TARGET) GOTO 200
          J=J+INCX
  100 CONTINUE
  200 ISRCHFLT=I

      RETURN
      END

