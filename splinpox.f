      SUBROUTINE SPLINPOX(F, X, FI, XI, N, SAFE, LPARAM, DFDX, D2FD2X)
C***********************************************************************
C***  CUBIC SPLINE INTERPOLATION, READY-FOR-USE
C***  XI(I), FI(I)  TABLE WHICH DEFINES THE FUNCTION TO BE INTERPOLATED
C***  X             ARGUMENT FOR WHICH THE FUNCTION VALUE IS REQUIRED
C***  FX            RESULTING FUNCTION VALUE
C***  THE RESULTING FUNCTION IS A PIECEWISE CUBIC INTERPOLATION 
C***     POLYNOMIAL WITH CONTINUOUS DERIVATIVE. EXTREMA CAN ONLY OCCUR 
C***     AT GIVEN MESHPOINTS (MONOTONIC VERSION AFTER M. STEFFEN)
C
C     Unified version implementing all features from classic routines
C       SPLINPO, SPLINPO_FAST and SPLINP from Goetz branch (A. Sander, Jan 2012)
C
C     Due to the optional arguments (e.g. for the derivatives), all main routines
C     need to have the following interface block:
C
C      INTERFACE SPLINPO
C        SUBROUTINE SPLINPO(F, X, FI, XI, N, DFDX, D2FD2X)
C          INTEGER, INTENT(IN) :: N          
C          REAL, DIMENSION(N), INTENT(IN) :: XI, FI
C          REAL, INTENT(OUT) :: F
C          REAL, INTENT(IN) :: X
C          LOGICAL, INTENT(IN), OPTIONAL :: SAFE
C          INTEGER, INTENT(IN), OPTIONAL :: LPARAM
C          REAL, INTENT(OUT), OPTIONAL :: DFDX, D2FD2X
C        END SUBROUTINE
C      END INTERFACE SPLINPO
C
C***********************************************************************

      IMPLICIT NONE

      INTEGER, INTENT(IN) :: N
      
      REAL, DIMENSION(N), INTENT(IN) :: XI, FI
      REAL, INTENT(OUT) :: F
      REAL, INTENT(IN) :: X
      LOGICAL, INTENT(IN), OPTIONAL :: SAFE
      INTEGER, INTENT(IN), OPTIONAL :: LPARAM
      REAL, INTENT(OUT), OPTIONAL :: DFDX, D2FD2X

      REAL :: DN, DX, DXM, FS0,
     >        D1, D2, D3, D23, H11, H12, H13, H14,
     >        H21, H22, H23, H24, H31, H32, H33, H34,
     >        H41, H42, H43, H44,
     >        F1, F2, F3, F4, FSM, FSP, S3, S4, S5,
     >        P, P1, P2, P3, P4

      INTEGER :: L, I, LA, LB

      LOGICAL :: bFoundInterval, bSafeMode

      !File and channel handles
      INTEGER, PARAMETER :: hCPR = 0        !write to wruniqXX.cpr

      SAVE      !this is required due to the _SAME_X entry point

      IF (present(SAFE)) THEN
        bSafeMode = SAFE
      ELSE
        bSafeMode = .TRUE.              !safe mode is default
      ENDIF

C***  CHECK FOR STRICTLY MONOTONIC ORDER (safe mode only)
      IF (bSafeMode) THEN
        DN = XI(N) - XI(1)
        DO L=2, N
          DX = XI(L) - XI(L-1)
          IF (DX*DN <= .0) THEN
            WRITE (hCPR, 3) XI(1), L-1, XI(L-1), L, XI(L), N, XI(N)
    3       FORMAT (' *** BAD USE OF SUBROUTINE SPLINPOX:', 
     >           ' X-VALUES NOT IN STRICLY MONOTONIC ORDER!', /, 
     >           ' X(1)=',        G12.5, 5X,
     >           ' X(', I3, ')=', G12.5, 5X, 
     >            'X(', I3, ')=', G12.5, 5X,
     >            'X(N=', I3, ')=', G12.5)
            CALL TRBK
            STOP 'ERROR'
          ENDIF
        ENDDO

C***    FIND THE INTERVAL XI(L-1), XI(L) WHICH INCLUDES X
        bFoundInterval = .FALSE.
        DO I=2, N
          L = I
          IF ( (X-XI(L-1)) * (X-XI(L)) <= .0) THEN
            bFoundInterval = .TRUE.
            EXIT
          ENDIF
        ENDDO

        IF (.NOT. bFoundInterval) THEN
          CALL REMARK('BAD USE OF SUBR. SPLINPOX - X OUTSIDE TABLE')
          WRITE (hCPR,'(A,G12.5,5X,A,G12.5,5X,A,I3,A,G12.5)') 
     >         ' X=', X,' X(1)=', XI(1),' X(', N, ')=', XI(N)
          CALL TRBK
          STOP 'ERROR'
        ENDIF
        IF (present(LPARAM)) THEN
          IF ((LPARAM > 0) .AND. (LPARAM /= L)) THEN
            !If LPARAM has been specified (and is > 0), it must have the correct value
            WRITE (hCPR,'(A, I5)') '**** LPARAM =', LPARAM
            WRITE (hCPR,'(A, I5)') '**** CORRECT L =', L
            STOP ' *** ERROR IN SUBR. SPLINPO - WRONG INDEX LPARAM'
          ENDIF
        ENDIF

      ENDIF

      !Note: Unlike in other languages, the following two IF statements
      !      cannot be combined, because the second check would be made
      !      even if present(LPARAM) already returns FALSE. Therere it
      !      would cause a crash if LPARAM is not set 
      IF (present(LPARAM)) THEN
        IF (LPARAM >= 0) THEN
          !Preset L if LPARAM has been specified and is inside 2..N
          L = LPARAM
        ENDIF
      ENDIF

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

C      WRITE (hCPR,*) 'Debug: L-2=', L-2, '  LA=', LA

C***  Entry point for a subsequent call with same x point, 
C***     but different function 
C      ENTRY SPLINPOX_SAME_X (F, X, FI, XI, N, SAFE)

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
         IF (L == 2) THEN 
            F3 = FS0
         ELSE
            P = (XI(L) - XI(L-1)) / (XI(L) - XI(L-2))
            FSM = (FI(L-1) - FI(L-2)) / (XI(L-1) - XI(L-2))
            F3 = P * FSM + (1.-P) * FS0
         ENDIF
         IF (L == N) THEN 
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
          IF (LA /= L-2) THEN
            S3 = S4
          ELSE
            S3 = ( FI(L-1) - FI(L-2) ) / ( XI(L-1) - XI(L-2) )
          ENDIF
C***   We are not in the last interval:
          IF (LB /= L+1) THEN
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

C***  Calculation of derivatives (optional)
      !added on 06.10.2011 to provide the same functionality as goetz
      IF (present(DFDX)) THEN
        DFDX = 3. * P1 *  DXM * DXM + P2 - 3. * P3 * DX * DX - P4
      ENDIF
      IF (present(D2FD2X)) THEN
        D2FD2X = 6. * (P1 *  DXM + P3 * DX)
      ENDIF

      RETURN
      END
