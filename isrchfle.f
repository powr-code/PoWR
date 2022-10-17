      FUNCTION ISRCHFLE(N,X,INCX,TARGET)

C***  NAME
C***     ISRCHFLE searches a real vector for the first element that is less
C***     than or equal to a real target.

C***  SYNOPSIS
C***  SEE ISRCHLE

      REAL X(*), TARGET

      J=1
      ISRCHFLE=0
      IF(N.LE.0) RETURN
      IF(INCX.LT.0) J=1-(N-1)*INCX
      DO 100 I=1,N
        IF(X(J).LE.TARGET) GOTO 200
          J=J+INCX
  100 CONTINUE
  200 ISRCHFLE=I

      RETURN
      END

