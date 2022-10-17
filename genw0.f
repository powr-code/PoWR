      SUBROUTINE GENW0 (JP,LMAX,ND,NP,Z,RADIUS,P,W0)
C***********************************************************************
C***  AT GIVEN IMPACT-PARAMETER INDEX JP, THIS SUBROUTINE CALCULATES THE
C***  WEIGHTS FOR THE ANGLE INTEGRATION OF THE 0. MOMENT, W0(L),
C***  AT ALL DEPTH POINTS L.
C***  THE ANGLE INTEGRATION IS PERFORMED IN THE VARIABLE "Z".
C***   - CALLED FROM: SUBROUTINE CMFRAY
C***********************************************************************
 
      DIMENSION W0(ND), Z(ND,NP),P(NP),RADIUS(ND)
 
C***  LOOP OVER ALL DEPTH POINTS
      DO 10 L=1,LMAX

      RL=RADIUS(L)
      RL2=RL+RL
      IF (JP.GT.1) GOTO 8

C***  FIRST STEP IF JP=1
      B=Z(L,1)
      A=Z(L,2)
      W0(L)=(B-A)/RL2
      GOTO 10

    8 IF (L.EQ.LMAX .AND. JP.GT.(NP-ND) ) GOTO 9

C***  INTERMEDIATE STEP
      A=Z(L,JP+1)
      B=Z(L,JP)
      C=Z(L,JP-1)
      W0(L)=(C-A)/RL2
      GOTO 10

C***  LAST STEP, IMPLYING Z(L,JMAX)=0
    9 B=Z(L,JP-1)
      W0(L)=B/RL2

   10 CONTINUE

      RETURN
      END
