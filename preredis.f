      SUBROUTINE PREREDIS (WREDI, NREDI, DKR, NREDMAX, ND, T, 
     >                     BWES, DXMAX, VDOP, BPLOT, WREDISUM)
C**********************************************************************
C***  Calculates the Weigths for the Redistribution (WREDI)
C**********************************************************************

      DIMENSION WREDI(0:NREDMAX, ND), T(ND), NREDI(ND), WREDISUM(ND)
      LOGICAL BPLOT

      ATMASS = 4.

cC***  Initialize WREDI for Plot purposes
      IF (BPLOT) THEN
        DO L=1, ND
          WREDISUM(L) = 0.
        ENDDO
      ENDIF

C***  Initialize WREDI
      DO L=1, ND
        DO I=0, NREDMAX
          WREDI(I,L) = 0.
        ENDDO
      ENDDO

C***  Maximum width of the Redistribution
      DKR = 0.

      DO L=1, ND
        VTHHE2 = 0.01651 * T(L) / ATMASS 
        VMICRO2 = VDOP * VDOP - VTHHE2
        VTHEL2 = 30.3165 * T(L)
        VDUEL = SQRT (VMICRO2 + VTHEL2) / VDOP
C***  Widths of the Redistribution
        NREDI(L) = NINT (BWES * 2. * VDUEL / DXMAX)
        IF (NREDI(L) .GT. NREDMAX) THEN
          WRITE (0,*) 'Dimension NREDMAX of WREDI insufficient'
          WRITE (0,*) 'L, NREDI(L), NREDMAX=', NREDI(L), NREDMAX
          STOP 'ERROR in Subr. PREREDIS'
        ENDIF
        FNREDI = FLOAT(NREDI(L))
        IF (FNREDI .GT. DKR) THEN
          DKR = FNREDI
        ENDIF
        SUM = 0.
C***  Calculate the Weigths
        DO I=0, NREDI(L)
          ARG = I * DXMAX / (2. * VDUEL)
          WREDI(I,L) = FIERFC(ARG)
          IF (I .EQ. 0) THEN
            SUM = SUM + WREDI(I,L)
          ELSE
            SUM = SUM + 2. * WREDI(I,L)
          ENDIF
        ENDDO
c        write (0,'(a,i3,1x,3(f10.2,1x),i6)') 'PREREDIS : L, SUM=', 
c     >  L, T(L), VDUEL, SUM, NREDI(L)
C***  Normalization
        SUMI = 1. / SUM
        DO I=0, NREDI(L)
          WREDI(I,L) = WREDI(I,L) * SUMI
        ENDDO
      ENDDO

C***  Testplot facility
      IF (BPLOT) THEN
        DO L=1, ND
          DO I=0, NREDI(L)
            IF (I .EQ. 0) THEN
              WREDISUM(L) = WREDISUM(L) + WREDI(I,L)
            ELSE
              WREDISUM(L) = WREDISUM(L) + 2. * WREDI(I,L)
            ENDIF
          ENDDO
        ENDDO
        OPEN (UNIT=50, FILE='coli_redis.dat')
        WRITE (50, '(A3,1X,6(A9,1X))')
     >    '*L', 'DEPTH1', 'D2', 'D3', '...', '...', '...'
        DO I=0, NREDMAX
          WRITE (50,'(I3,3X,150(F9.7,1X))') I, (WREDI(I,LL), LL=1, ND)
        ENDDO
        WRITE (50,'(A,1X,150(F9.7,1X))') 
     >        '* SUM', (WREDISUM(LL), LL=1, ND)
        CLOSE(50)
      ENDIF

      RETURN
      END
