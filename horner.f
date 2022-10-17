      SUBROUTINE HORNER (X, Y, KPLUS1, Z)
C***  Polynomial Y(X) 
C***  The first coefficient belongs to the highest exponent!
      DIMENSION Z(KPLUS1)

      Y = 0.

      DO I=1, KPLUS1-1
        Y = (Y + Z(I)) * X
      ENDDO

      Y = Y + Z(KPLUS1)

      RETURN
      END
