      SUBROUTINE GENWP1 (NP, P, WP1, WP1LAST)
C***********************************************************************
C***  GENERATE INTEGRATION WEIGHTS FUR FLUX-INTEGRATION, P*DP
C***  Quadrature sum: H = H + V(J) * WP1(J)   J=1...JMAX-1
C***   The additional weight vector WP1LAST allows the integration 
C***   along interstice radius shells: 
C***   The last integration step for those rays with p >=1. is 
C***   between P(JMAX) and the intermesh point 0.5*(P(JMAX)+P(JMAX+1); 
C***   WP1LAST replaces the weight at JMAX acccordingly:
C***     H = H + V(JMAX) * WP1LAST(LZ) 
C***     V does not exist at the additional point JMAX+0.5, but is 
C***     implicitely assumed to be ZERO (symmetry!)
C***********************************************************************
 
      DIMENSION P(NP), WP1(NP), WP1LAST(NP)

C***  Normal integration weights J=1 ... NP 

C***  FIRST POINT
      A = P(1)
      B = P(2)
      WP1(1) = (B - A) * (B + 2. * A) / 6.

C***  INNER POINT
      DO J=2,NP-1
        A = P(J-1)
        B = P(J)
        C = P(J+1)
        WP1(J) = (A + B + C) * (C - A) / 6.
      ENDDO

C***  LAST POINT
      A = P(NP-1)
      B = P(NP)
      WP1(NP) = (B - A) * (2. * B + A) / 6.

C***  Special integration weight for J=JMAX if the integration interval ends
C***   at the interstice (J, J+1)

      DO J=1, NP-1
         IF (P(J) .LT. 1.) THEN
            WP1LAST(J) = WP1(J)
         ELSE
           A = P(J-1)
           B = P(J)
           C = 0.5 * (P(J) + P(J+1))
           WP1LAST(J) = (A + B + C) * (C - A) / 6.
         ENDIF
      ENDDO

      RETURN
      END
