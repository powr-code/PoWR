      SUBROUTINE REGULA(F,X,Y,X1,X2,EPS)
C***********************************************************************
C***  THIS ROUTINE CALCULATES THE SOLUTION X OF F(X)=Y IN THE INTERVAL
C***  (X1,X2) ! METHOD: 
C***  -------- REGULA FALSI --------
C***  USING BISECTION STEPS TO GUARANTEE CONVERGENCE, PRECISION IN X: EPS
C***********************************************************************
      LOGICAL BI
 
      BI=.TRUE.
      A=X1
      B=X2
      FA=F(A)-Y
      IF(FA.EQ..0) GOTO 3
      FB=F(B)-Y
      IF(FB.EQ..0) GOTO 4
      IF (FA*FB.GT..0) THEN
         WRITE (0,*) '*** INVALID ARGUMENTS WHEN CALLING REGULA ***'
         WRITE (0,10) A, B
   10    FORMAT ('INTERVAL:     A  =', E15.5, 5X, '  B  =', E15.5)
         WRITE (0,11) FA, FB
   11    FORMAT ('FUNCTION:   F(A) =', E15.5, 5X, 'F(B) =', E15.5)
         CALL TRBK
         STOP 'ERROR'
         ENDIF

    1 D=A-B
      IF(ABS(D).LT.EPS) GOTO 5
      BI=.NOT.BI
      IF(BI) X=A-FA*D/(FA-FB)
      IF(.NOT.BI) X=.5*(A+B)
      FX=F(X)-Y
      IF(FX.EQ..0) RETURN
      IF(FX*FA.GT..0) GOTO 2

      B=X
      FB=FX
      GOTO 1

    2 A=X
      FA=FX
      GOTO 1

    3 X=A
      RETURN

    4 X=B
      RETURN

    5 X=(A+B)*.5
      RETURN

      END
