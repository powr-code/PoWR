      SUBROUTINE POLYFIT (XFIT, YFIT, WFIT, M, KPLUS1, 
     >                    X, A, B, SCRATCH, ATEST, BTEST, DTEST)
C******************************************************************************
C***  THIS SUBROUTINE CALCULATES THE POLYNOM KOEFFICIENTS FOR AN
C***  POLYNOMINAL INTERPOLATION OF RANK (KPLUS1-1) FOR M POINTS
C***  GIVEN BY XFIT AND YFIT AND WEIGHTED WITH WFIT.
C***  SCRATCH, ATEST, BTEST AND DTEST ARE USED FOR THE ROUTINE LINSOL
C******************************************************************************

      DIMENSION XFIT(M),YFIT(M),WFIT(M)
      DIMENSION X(KPLUS1), A(KPLUS1,KPLUS1), B(KPLUS1)
      DIMENSION SCRATCH(2*KPLUS1), ATEST(KPLUS1,KPLUS1)
      DIMENSION BTEST(KPLUS1), DTEST(KPLUS1)
      CHARACTER*4 CKEY

      CKEY = 'OWN'

      IF (M .LT. KPLUS1) THEN
         CALL REMARK ('KPLUS1 EXCEEDS M')
         STOP 'ERROR IN POLYFIT'
         ENDIF

C***  NEW BRANCH TO AVOID NAG-CALL E02ADF
C***  Probably (!) this is only made for cubic degree
      IF (KPLUS1 .NE. 4) THEN
         WRITE (0,*) 'INVALID DEGREE OF POLYNOMIAL FIT: .NE. 4'  
         STOP '****** ERROR IN POLYFIT **********'
      ENDIF

      DO I=0, 6
        SCRATCH(I+1) = 0.
      ENDDO
      DO I=0, 3
        DTEST(I+1) = 0.
      ENDDO
      DO J=1, M
        XF = 1.
        DO I=0, 6
          SCRATCH(I+1) = SCRATCH(I+1) + WFIT(J) * XF
          IF (I .LT. KPLUS1) DTEST(I+1) = DTEST(I+1) + WFIT(J) * 
     >                       XF * YFIT(J)
          XF = XF * XFIT(J)
        ENDDO
      ENDDO

      DO I=1, KPLUS1
        B(I) = DTEST(KPLUS1+1-I)
      ENDDO

      A(1,1) = SCRATCH(7)
      A(2,1) = SCRATCH(6)
      A(1,2) = SCRATCH(6)
      A(3,1) = SCRATCH(5)
      A(2,2) = SCRATCH(5)
      A(1,3) = SCRATCH(5)
      A(4,1) = SCRATCH(4)
      A(3,2) = SCRATCH(4)
      A(2,3) = SCRATCH(4)
      A(1,4) = SCRATCH(4)
      A(4,2) = SCRATCH(3)
      A(3,3) = SCRATCH(3)
      A(2,4) = SCRATCH(3)
      A(4,3) = SCRATCH(2)
      A(3,4) = SCRATCH(2)
      A(4,4) = SCRATCH(1)

      CALL INV (KPLUS1, KPLUS1, A, CKEY)
      CALL VMF (X, B, A, KPLUS1, KPLUS1)

C!!!      CALL LINSOL (X, A, B, KPLUS1, KPLUS1, 
C!!!     >             SCRATCH, ATEST, BTEST, DTEST)

      RETURN
      END
