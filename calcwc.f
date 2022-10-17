      SUBROUTINE CALCWC(IFF_WCHARM, IFF_N, L, WCHARM_FINE, 
     >                  INDSTART, GAMMAC)
C*******************************************************************
C***  Calculates WCHARM_FINE from the compreessed Integer*1-Array 
c***    IFF_WCHARM
C***  Called by SETXJFINE
C*******************************************************************

      INTEGER, PARAMETER :: TINYINT = SELECTED_INT_KIND(2)

      INTEGER (KIND=TINYINT), DIMENSION(IFF_N) :: IFF_WCHARM
      REAL, DIMENSION(IFF_N) :: WCHARM_FINE
      DIMENSION T(0:15)

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

      LTESTOUT = 0
C!!!      LTESTOUT = 10

c      if (l .ne. ltestout) return

      IF (LTESTOUT .EQ. L) THEN
        OPEN (UNIT=141, FILE='wcharm_setxjfine.dat', STATUS='UNKNOWN')
        write (141,*) '* iff_n=', iff_n
      ENDIF

      DO K=1, IFF_N
        INDP = IFF_WCHARM(K) + 101
        I = INT(INDP/10)
        J = INDP - 10*I

        if (i.gt.16 .or. j .eq. 0) then
          write (0,*) 'Problems in CALCWC'
          write (0,*) 'indp, i, j=', indp, i, j
        endif

        IF (I .EQ. 16) THEN
          IF (J .EQ. 1) THEN
            WCI = 0.
          ELSE
            STOP 'STOP: Problems in CALCWC'
          ENDIF
        ELSE
          WCI = AMIN1(1., FLOAT(J+1) * T(I))
        ENDIF

C***  Apply Gamma damping factor GAMMAC
        IF (WCI .GT. 0.) THEN
          WCI = EXP(ALOG(WCI) / GAMMAC)
        ENDIF

        WCHARM_FINE(K) = 1. - WCI

C***  TEST-OUTPUT
        IF (LTESTOUT .EQ. L) THEN
          WRITE (141, '(I5,1X,2(E20.10,1X),3(I4,1X))')
     >          K+INDSTART-1, WCI, WCHARM_FINE(K), I, J, INDP
        ENDIF

      ENDDO

      RETURN
      END
