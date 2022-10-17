      SUBROUTINE CLSAVEWC (IFF_WCHARM, IFF_MAX, IFF_N, ND, DJDSMO, 
     >                     K, XLAMK)

C*********************************************************************
C***  Stores DJDSMO in special (short) array IFF_WCHARM
C***  To avoid logarithms, array T is used
C*********************************************************************
      INTEGER, PARAMETER :: TINYINT = SELECTED_INT_KIND(2)

      INTEGER (KIND=TINYINT), DIMENSION(IFF_MAX, ND) :: IFF_WCHARM
      DIMENSION DJDSMO(ND), T(0:15)

      DATA T / 1.0,
     >         0.1,
     >         0.01,
     >         0.001,
     >         0.0001,
     >         0.00001,
     >         0.000001,
     >         0.0000001,
     >         0.00000001,
     >         0.000000001,
     >         0.0000000001,
     >         0.00000000001,
     >         0.000000000001,
     >         0.0000000000001,
     >         0.00000000000001,
     >         0.000000000000001 /

      DO L=1, ND
        WC = AMAX1(0.,DJDSMO(L))
        WCI = 1. - WC
        WCI = AMAX1(0.,WCI)
        IF (WCI .EQ. 0.) THEN
          I = 16
          J = 1
        ELSE
          DO I=0, 15
            IF (WCI .GE. T(I)) EXIT
          ENDDO
          J = INT(WCI/T(I))
        ENDIF
        INDEX = 10*I + J - 101
        IF (INDEX .GT. 100) THEN
          WRITE (0,*) 'Trouble in storing DJDSMO, INDEX too large'
          STOP 'ERROR in Subr. CLSAVEWC'
        ENDIF
        IFF_WCHARM(IFF_N,L) = INDEX

C***  TEST-OUTPUT
        IF (.false.) THEN
          IF (L .EQ. 10) THEN
            OPEN (UNIT=140, FILE='wcharm_coli.dat', STATUS='UNKNOWN')
            WRITE (140, '(I5,1X,E20.10,1X,3(I4,1X),1X,E20.10)')
     >            K, WCI, I, J, INDEX+101, XLAMK
          ENDIF
        ENDIF
      ENDDO

      RETURN
      END
