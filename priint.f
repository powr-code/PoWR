      SUBROUTINE PRIINT (XJC,EDDI,R,ND,XLAMBDA,NF,LSTEP,JOBNUM,MODHEAD)
C***********************************************************************
C***  PRINT INTENSITIES
C***********************************************************************

      DIMENSION XJC(ND),XLAMBDA(NF),EDDI(3,ND),R(ND)
      CHARACTER MODHEAD*100
      CHARACTER*8 NAME
 
      PRINT 2,MODHEAD,JOBNUM
    2 FORMAT (1X,  A  ,20X,'JOB NO.',I3,
     $            //,10X,'MEAN INTENSITY',/,10X,14('-'),/,
     $ ' FREQUENCY     DEPTH       J-NUE          T-RAD   ',
     $ '         EDDINGTON     SPHERICITY      EDDINGTON',/,
     $ '   INDEX       INDEX     (ERG/CM+2)       (KELVIN)',
     $ '         FACTOR F      FACTOR  R.R.Q    FACTOR H/J',/)
 
      DO 4 K=1,NF
      WRITE (NAME, '(A3, I4, A1)') 'XJC', K, ' '
      CALL READMS (3,XJC,ND,NAME, IERR)
      IF (K <= 999) THEN
        WRITE (NAME, '(A4, I3, A1)') 'EDDI', K, ' '
      ELSE
        WRITE (NAME, '(A4, I4)') 'EDDI', K
      ENDIF
      CALL READMS (3,EDDI,3*ND,NAME, IERR)
      DO 1 L=1,ND
      IF (((L-1)/LSTEP)*LSTEP .NE. (L-1) .AND. L.NE.ND) GOTO 1
      TRAD=TRADFUN (XLAMBDA(K),XJC(L))
      IF (L.EQ.ND .OR. L.EQ.1) THEN
      PRINT 5,K,L,XJC(L),TRAD,EDDI(1,L),EDDI(2,L)*R(L)*R(L),EDDI(3,L)
    5 FORMAT (2I10,E15.3,F15.0,F15.3,F15.3,F15.3)
      ELSE
      PRINT 5,K,L,XJC(L),TRAD,EDDI(1,L),EDDI(2,L)*R(L)*R(L)
      ENDIF
    1 CONTINUE
    4 PRINT 6
    6 FORMAT (1X)
 
      RETURN
      END
