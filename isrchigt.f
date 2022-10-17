      FUNCTION ISRCHIGT(N,IAR,INCX,ITARGET)

C***  NAME
C***     ISRCHIGT searches an integer vector for the first element that is
C***     greater than an integer target.

C***  SYNOPSIS
C***  SEE ISRCHEQ

      DIMENSION IAR(N)

      J=1
      ISRCHIGT=0
      IF(N.LE.0) RETURN
      IF(INCX.LT.0) J=1-(N-1)*INCX
      DO 100 I=1,N
        IF(IAR(J).GT.ITARGET) GOTO 200
          J=J+INCX
  100 CONTINUE
  200 ISRCHIGT=I

      RETURN
      END

