      SUBROUTINE ADJUST_LX (XDATA, MAXXDAT, EMCOLI, EMCOLI_LAMBDA, 
     >                      XLAMBDA, NF, XLOGLXLBOL, NO_UPDATE)
C*********************************************************************
C***  If logLxLbol is requested, adjust XFILL accordingly 
C***  If NO_UPDATE: XLOGLXLBOL is calculated, but XFILL is not updated
C*********************************************************************

      IMPLICIT NONE

      REAL EMCOLI(NF), XLAMBDA(NF), EMCOLI_LAMBDA(NF) 
      REAL XDATA(MAXXDAT)
      INTEGER MAXXDAT, NF, K
      LOGICAL NO_UPDATE

      REAL XFILL, XFILL_OLD, XMIN, XMAX, XLX, XLBOL
      REAL XLOGLXLBOL, XLOGLXLBOL_AIM, RATIO

      XFILL_OLD = XDATA(1)
      IF (XFILL_OLD .LE. .0) THEN
         WRITE (0,*) '*** A non-zero starting value for XFILL is'
         WRITE (0,*) '*** required for an automatic adjustment of L_X'
         STOP 'FATAL ERROR detected by subroutine ADJUST_LX'
      ENDIF

C***  Calculate current LX/LBOL ratio

C***  Convert f_nu to f_lambda (clight factor omitted)
      DO K=1, NF
         EMCOLI_LAMBDA(K) = EMCOLI(K) / XLAMBDA(K)**2      
      ENDDO

C***  Note: both integrals are missing factors which cancel out
      XMIN = XLAMBDA(1)
      XMAX = XLAMBDA(NF)
      CALL INTEGRATE (XLBOL, XLAMBDA, EMCOLI_LAMBDA, NF, XMIN, XMAX)
      XMIN = MAX (XDATA(9), XLAMBDA(1))
      XMAX = MAX (XDATA(10), XLAMBDA(1))
      CALL INTEGRATE (XLX, XLAMBDA, EMCOLI_LAMBDA, NF, XMIN, XMAX)

      XLOGLXLBOL = ALOG10(XLX/XLBOL)

C***  Exit if only XLOGLXLBOL is wanted for output
      IF (NO_UPDATE) RETURN

C**** Branch in order to update XFILL
      XLOGLXLBOL_AIM = XDATA(8)

C***  Difference in log -> dex
      RATIO = 10.**(XLOGLXLBOL_AIM - XLOGLXLBOL)      

C***  restrict correction factor to 2dex
      RATIO = MIN (RATIO, 100.)
      RATIO = MAX (RATIO, 0.01)

C***  Damp corrections
      RATIO = 1 + (RATIO-1.0)*0.8
      
C***  Adjust XFILL (and XFILL2)
      XFILL = XFILL_OLD * RATIO
      XDATA(1) = XFILL
      XDATA(5) = XDATA(5) * RATIO

ccc   test output
      write (0,'(A,F10.3)') 'Current log(Lx/Lbol) =', XLOGLXLBOL
      write (0,'(A,F10.3)') 'Aim     log(Lx/Lbol) =', XLOGLXLBOL_AIM
      write (0,'(A,1PG14.3)') 'Current XFILL =', XFILL_OLD
      write (0,'(A,1PG14.3)') 'New     XFILL =', XFILL

      RETURN
      END
