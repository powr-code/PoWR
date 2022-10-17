      SUBROUTINE INHIBIT (POPNUM, N, ND, NCHARG, RNE, 
     $                   NATOM, ABXYZ, NFIRST, NLAST, POPMIN)
C*******************************************************************************
C***  THE POPULATION NUMBERS ARE MODIFIED BY THIS SUBROUTINE,
C***     IN ORDER TO AVOID NUMBERS WHICH ARE TOO SMALL OR NEGATIVE
C***  CALLED FROM: STEAL, EXTRAP, MODIFY
C!!!  NOTE : NORMALLY RENORM IS A LOGICAL, BUT THE CRAY PREPROZESSOR (fpp) 
C!!!         HAS A BUG. NOW IRENORM IS USED WITH 
C!!!         RENORM = FALSE --> IRENORM = 0
C!!!         RENORM = TRUE  --> IRENORM = 1
C*******************************************************************************
 
      DIMENSION POPNUM(ND,N)
      DIMENSION NCHARG(N), RNE(ND)
      DIMENSION ABXYZ(NATOM), NFIRST(NATOM), NLAST(NATOM)
C!!!      LOGICAL RENORM

C**********************************************************************
C***  SET POPMIN : SMALLER POPNUMBERS ARE SET TO THIS VALUE
C***     New on 29-Sep-1999 13:51:23 (wrh)
C***     POPMIN is now set by the following CARDS option: 
C***     POPMIN = ...
C***     The Default value (1.E-100) is set in DECSTE
C**********************************************************************

C**********************************************************************
C***  INHIBIT NEGATIVE (AND/OR SMALL) POP.NUMBERS
C**********************************************************************
      NSMALL=0
      NEGWARN=0

      DO 19 L=1,ND
C!!!      RENORM = .FALSE.
      IRENORM = 0

      DO 16 NA=1,NATOM

      DO 16 J=NFIRST(NA),NLAST(NA)
      POPLJ=POPNUM(L,J)

C***  INHIBIT SMALL POP.NUMBERS 
      IF (POPLJ .LT. POPMIN .AND. POPLJ .GE. 0.0) THEN
C!!!            RENORM = .TRUE.
            IRENORM = 1
            NSMALL = NSMALL + 1
            POPNUM(L,J) = POPMIN
            ENDIF
C***  INHIBIT NEGATIVE POP. NUMBERS
      IF (POPLJ .LT. .0) THEN
C!!!            RENORM = .TRUE.
            IRENORM = 1
            NEGWARN = NEGWARN + 1
            POPNUM(L,J) = POPMIN
            ENDIF
   16 CONTINUE
 
C***  RENORMALIZATION OF THE POPULATION NUMBERS
C!!!      IF (RENORM) THEN
      IF (IRENORM .EQ. 1) THEN
         DO 15 NA=1,NATOM
         NFIRNA=NFIRST(NA)
         NLANA=NLAST(NA)
         SUM=0.0

         DO 18 J=NFIRNA,NLANA
   18    SUM = SUM + POPNUM(L,J)

         SUM = SUM / ABXYZ(NA)

         DO 15 J=NFIRNA,NLANA
         POPNUM(L,J) = POPNUM(L,J) / SUM
   15    CONTINUE

C***     RE-ADJUST THE ELECTRON DENSITY
         RNEL=0.0
         DO 20 J=1,N
         RNEL = RNEL + NCHARG(J) * POPNUM(L,J)
   20    CONTINUE
         RNE(L)=RNEL
      ENDIF
   19 CONTINUE
 

C**********************************************************************
C***  PRINTOUT OF WARNINGS AND ERROR MESSAGES  ************************
C**********************************************************************
      IF (NSMALL .GT. 0) PRINT 27, NSMALL, POPMIN
   27 FORMAT (10X, I4, 
     >  ' WARNINGS: POP. NUMBERS .LT.', 1PE8.1, ' SET TO THIS VALUE')
      IF (NEGWARN .GT. 0) PRINT 17, NEGWARN, POPMIN
   17 FORMAT (10X,I4,' WARNINGS: NEGATIVE POP. NUMBERS ARE SET TO ',
     $        1PE8.1 )
  
      RETURN
      END
