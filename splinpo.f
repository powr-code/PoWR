      SUBROUTINE SPLINPO (F, X, FI, XI, N)
C***********************************************************************
C***  CUBIC SPLINE INTERPOLATION, READY-FOR-USE
C***  XI(I), FI(I)  TABLE WHICH DEFINES THE FUNCTION TO BE INTERPOLATED
C***  X             ARGUMENT FOR WHICH THE FUNCTION VALUE IS REQUIRED
C***  FX            RESULTING FUNCTION VALUE
C***  THE RESULTING FUNCTION IS A PIECEWISE CUBIC INTERPOLATION 
C***     POLYNOMIAL WITH CONTINUOUS DERIVATIVE. EXTREMA CAN ONLY OCCUR 
C***     AT GIVEN MESHPOINTS (MONOTONIC VERSION AFTER M. STEFFEN)
C***********************************************************************

      DIMENSION XI(N), FI(N)

C***  CHECK FOR STRICTLY MONOTONIC ORDER
      DN = XI(N) - XI(1)
      DO 1 L=2, N
      DX = XI(L) - XI(L-1)
      IF (DX*DN .LE. .0) THEN
         WRITE (0, 3) XI(1), L-1, XI(L-1), L, XI(L), N, XI(N)
    3    FORMAT (' *** BAD USE OF SUBROUTINE SPLINPO:', 
     >           ' X-VALUES NOT IN STRICLY MONOTONIC ORDER!', /, 
     >           ' X(1)=',        G12.5, 5X,
     >           ' X(', I3, ')=', G12.5, 5X, 
     >            'X(', I3, ')=', G12.5, 5X,
     >            'X(N=', I3, ')=', G12.5)
         CALL TRBK
         STOP 'ERROR'
         ENDIF
    1 CONTINUE

C***  FIND THE INTERVAL XI(L-1), XI(L) WHICH INCLUDES X
      DO 4 I=2, N
      L = I
      IF ( (X-XI(L-1)) * (X-XI(L)) .LE. .0) GOTO 2
    4 CONTINUE

      CALL REMARK('BAD USE OF SUBR. SPLINPO - X OUTSIDE TABLE')
      WRITE (0,'(A,G12.5,5X,A,G12.5,5X,A,I3,A,G12.5)') 
     >         ' X=', X,' X(1)=', XI(1),' X(', N, ')=', XI(N)
      CALL TRBK
      STOP 'ERROR'

    2 CONTINUE

C***  DETERMINATION OF THE COEFFICIENTS P1, P2, P3, P4 (CF. SUBR. CUBIC)
 
C***  SET UP THE COEFFICIENT MATRIX
      D1=1./(XI(L)-XI(L-1))
      D2=D1*D1
      D3=D1*D2
      D23=D2/3.
      H11=D3
      H12=-D3
      H13=D23
      H14=2.*D23
      H21=-D1
      H22=2.*D1
      H23=-0.333333333333333
      H24=-0.666666666666666
      H31=-D3
      H32=D3
      H33=-2.*D23
      H34=-D23
      H41=2.*D1
      H42=-D1
      H43=0.666666666666666
      H44=0.333333333333333
C***  FOR THE BOUNDARY INTERVALS THE DERIVATIVE CANNOT EXTEND OVER THE BOUNDARY
      LA=MAX0(L-2,1)
      LB=MIN0(L+1,N)
C***  FUNCTION TO BE INTERPOLATED: FI
      F1 = FI(L-1)
      F2 = FI(L)

      IF (.TRUE.) THEN
C***     Standard version: zentrierte Ableitungen an den Stuetzstellen. Das 
C***     ist bei nicht aequidistanten Stuetzstellen fragwuerdig, verringert 
C***     aber andererseits das Ueberschwingen insbesondere wenn man nicht 
C***     MONO verwendet 
         F3 = (FI(L) - FI(LA)) / (XI(L) - XI(LA))
         F4 = (FI(LB) - FI(L-1)) / (XI(LB) - XI(L-1))

      ELSE
C***  Alternative Version wrh 16-May-2002 13:39:29
C***  Statt der zentrierten Ableitung werden gewichtete Mittel der
C***  Ableitungen der angrenzenden Intervalle genommen. Das entspricht 
C***  der Ableitung eines Parabel-Fits durch drei Punkte.  
         FS0 = (FI(L) - FI(L-1)) / (XI(L) - XI(L-1))
         IF (L .EQ. 2) THEN 
            F3 = FS0
         ELSE
            P = (XI(L) - XI(L-1)) / (XI(L) - XI(L-2))
            FSM = (FI(L-1) - FI(L-2)) / (XI(L-1) - XI(L-2))
            F3 = P * FSM + (1.-P) * FS0
         ENDIF
         IF (L .EQ. N) THEN 
            F4 = FS0
         ELSE
            P = (XI(L) - XI(L-1)) / (XI(L+1) - XI(L-1))
            FSP = (FI(L+1) - FI(L)) / (XI(L+1) - XI(L))
            F4 = P * FSP + (1.-P) * FS0
         ENDIF
      ENDIF

C***  SET TRUE FOR MONO OPTION
      IF (.TRUE.) THEN

ccc   Diese bis heute (3-Sep-2002) verwendete Version erscheint mir 
ccc   merckwuerdig und an mehreren Stellen fehlerhaft! Ich lasse sie 
ccc   aus Dokumentationsgruenden hier stehen. Nachfolgend dann eine 
ccc   Version nach heutiger Erkenntnis. wrh  
c       S4 = ( FI(L) - FI(L-1) ) / ( XI(L) - XI(L-1) )
c       IF (LA .NE. L-2 .AND. LB .NE. L) THEN
c         S3 = S4
c         S5 = S4
c       ELSE
c         IF (LA .EQ. L-2) THEN
c           S3 = ( FI(L-1) - FI(L-2) ) / ( XI(L-1) - XI(L-2) )
c         ELSE
c           S3 = 1.4 * S4 - 0.5 * F3
c         ENDIF
c         IF (LB .EQ. L+1) THEN
c           S5 = ( FI(L+1) - FI(L) ) / ( XI(L+1) - XI(L) )
c         ELSE
c           S5 = 1.5 * S4 - 0.5 * F4
c         ENDIF
c       ENDIF

          S4 = ( FI(L) - FI(L-1) ) / ( XI(L) - XI(L-1) )
C***   We are not in the first interval:
          IF (LA .NE. L-2) THEN
            S3 = S4
          ELSE
            S3 = ( FI(L-1) - FI(L-2) ) / ( XI(L-1) - XI(L-2) )
          ENDIF
C***   We are not in the last interval:
          IF (LB .NE. L+1) THEN
             S5 = S4
          ELSE
             S5 = ( FI(L+1) - FI(L) ) / ( XI(L+1) - XI(L) )
          ENDIF

       F3 = (SIGN(1.0,S3)+SIGN(1.0,S4))*MIN(ABS(S3),ABS(S4),0.5*ABS(F3))
       F4 = (SIGN(1.0,S4)+SIGN(1.0,S5))*MIN(ABS(S4),ABS(S5),0.5*ABS(F4))

      ENDIF
 
C***  CALCULATE POLYNOMIAL COEFFICIENTS: P(VECTOR) = H(MATRIX) * F(VECTOR)
      P1=H11*F1+H12*F2+H13*F3+H14*F4
      P2=H21*F1+H22*F2+H23*F3+H24*F4
      P3=H31*F1+H32*F2+H33*F3+H34*F4
      P4=H41*F1+H42*F2+H43*F3+H44*F4
 

C***  EVALUATION OF THE INTERPOLATION POLYNOMIAL
      DXM = X - XI(L-1)
      DX  = XI(L) - X
      F = (P1 * DXM * DXM + P2 ) * DXM
     >  + (P3 * DX  * DX  + P4 ) * DX

      RETURN
      END
