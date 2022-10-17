      SUBROUTINE ARIPAR (ACTPAR, VALUE, IERR)
C**********************************************************************
C***  Decodes a floting-point type argument ACTPAR
C***  that is enclosed in parantheses and thus shall be avaluated 
C***  as an arithmetic expression
C**********************************************************************
      CHARACTER ACTPAR*(*), CALCLINE*1000, AREX*20
      PARAMETER (MAXVAR=2)
      CHARACTER*132 VAR(MAXVAR,2)
      LOGICAL CALC_DEBUG

      IERR = 0
      LAST = IDX(ACTPAR)
      CALCLINE = 'KASDEF CALC AREX = ' // ACTPAR(2:LAST-1)
      NVAR = 0
      DUMMY=0.
      CALC_DEBUG = .FALSE.
      CALL KD_CALC (CALCLINE, DUMMY, DUMMY, 
     >                 VAR, MAXVAR, NVAR, CALC_DEBUG)   
      READ (VAR(2,2), '(F20.0)', ERR=99) VALUE
 
      RETURN

   99 PRINT *, '*** ERROR WHEN EVALUATING ARITHMETIX EXPRESSION: ***'
      PRINT *, '>>> ', ACTPAR(:IDX(ACTPAR)), 
     >         ' <<<< CANNOT BE DECODED AS FLOATING POINT NUMBER!'

      IERR = 1
      RETURN

      END

