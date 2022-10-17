      SUBROUTINE MOMENT1 (R,NP,P,U,H)
C***********************************************************************
C***  CALCULATES AN ANGLE INTEGRAL H OF THE RADIATION FIELD U
C***  BESIDES OF THE OUTER BOUNDARY, THIS IS NOT THE 1. MOMENT H,
C***  BUT RATHER AN INTENSITY-LIKE QUANTITY
C***  INTEGRATION WITH TRAPEZOIDAL RULE, WEIGHTS P * DP
C***    -  VECTORIZING VERSION, REPLACED 26-MARCH-1991
C***********************************************************************
 
      DIMENSION U(NP),P(NP)

C***  FIRST POINT
      A = P(1)
      B = P(2)
      W = (B - A) * (B + 2. * A)
      H = W * U(1)

C***  INNER POINT
      DO 2 J=2,NP-1
      A = P(J-1)
      B = P(J)
      C = P(J+1)
      W = (A + B + C) * (C - A)
      H = H + W * U(J)
    2 CONTINUE

C***  LAST POINT
      A = P(NP-1)
      B = P(NP)
      W = (B - A) * (2. * B + A)
      H = H + W * U(NP)
      H = H / (R * R * 6.)

      RETURN
      END
