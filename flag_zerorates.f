      SUBROUTINE FLAG_ZERORATES(NFIRNA, NLANA, RATCO, NRANK, 
     >                          IMAXPOP, EN, POPMIN, ZERO_RATES)
C******************************************************************
C***  Check if all Rate Coefficients in one column are non-zero
C***  (otherwise: the matrix is singular!)
C***  and store logical flag ZERO_RATES for later use
C***  Note: This subroutine is called once per each element
C***  Called from:
C***  STEAL -  LINPOP - COMA 
***   STEAL -  POPZERO - NLTEPOP 
C******************************************************************

      DIMENSION RATCO(NRANK,NRANK), EN(NRANK)
      LOGICAL ZERO_RATES(2)

C***  Threshold value with hysteresis: 
C***    only if the estimated popnumber is lower, the level can be flagged
ccc      POPLIMIT = 0.01 * POPMIN ! version before  4-Mar-2015
      POPLIMIT = MAX(1.E-20 * POPMIN, 1.E-45)

C***  check colum J

C***  test modification, wrh 29-Aug-2008 16:23:38:
C***  Only the hightest levels can be switched off!
C***  I.e., for a level to be switched off, 
C***  the next-higher level must be already switched off 
C***  For this criterion, the loop over levels must run backwards:
      DO J = NLANA, NFIRNA, -1

C***     Never flag the IMAXPOP level!
C***     Never flag a level which was strongly populated last iteration
           IF (J .EQ. IMAXPOP .OR. EN(J) .GT. 1000*POPMIN) THEN
              ZERO_RATES(J) = .FALSE.
              CYCLE
           ENDIF

C***     Flag if diagonal element too small, else divison below will fail!
         IF (EN(J) .LT. POPMIN .AND. ABS(RATCO(J,J)) .LT. 1.E-200) THEN
            ZERO_RATES(J) = .TRUE.
            CYCLE
         ENDIF

C***     Sum up the gains of level J
         GAINS = .0
         DO I = NFIRNA, NLANA
            IF (I .EQ. J) CYCLE
            IF (EN(I) .GT. 1.1*POPMIN) 
     >         GAINS = GAINS + EN(I) * RATCO(I,J)
         ENDDO

C***     expected popnumber of level J
         ENJ = - GAINS / RATCO(J,J)
C***     mysteriously usefull??
         ENJ = ABS(ENJ)

C***     flag level, depending on predicted and on last popnumber
C***        the second condition provides some hysteresis 
         ZERO_RATES(J) = ENJ .LT. POPLIMIT
     >       .OR. (EN(J) .LE. 1.1*POPMIN .AND. ENJ .LE. 100.*POPLIMIT)

ccc version before  4-Mar-2015
ccc     >          .OR. (EN(J) .LE. 1.1*POPMIN .AND. ENJ .LE. 1.1*POPMIN)

C***     This version: Only the highest levels can be flagged
         IF (J .EQ. NLANA) CYCLE
         IF (.NOT. ZERO_RATES(J+1) ) ZERO_RATES(J)=.FALSE. 

      ENDDO

      RETURN
      END
