      SUBROUTINE PRIMO (ND,N,RNE,LEVEL,POPNUM,JOBNUM,MODHEAD,LSPOP,
     $                  IFRO,ITO)
C***********************************************************************
C***********************************************************************

      DIMENSION RNE(ND),POPNUM(ND,N)
      CHARACTER MODHEAD*100,LEVEL(N)*10

      PRINT 1,MODHEAD,JOBNUM,IFRO,ITO
    1 FORMAT (1X,  A,  20X,'JOB NO.',I3,
     $ //,20X,'RELATIVE NON-LTE POPULATION NUMBERS MODIFIED FROM',
     $ ' POINT',I3,' TO POINT',I3,/,20X,79('-'))
      J1=1
    4 J2=MIN0(N,J1+9)
      PRINT 2, (LEVEL(J),J=J1,J2)
    2 FORMAT (//,'  L EL.DENS.',10(2X,A10))
      IF (N.LT.J2) PRINT 11
   11 FORMAT (1X)
      PRINT 11
      DO 3 L=1,ND
      IF(((L-1)/LSPOP)*LSPOP.NE.(L-1) .AND. L.NE.ND) GOTO 3
      PRINT 9,L,RNE(L),(ALOG10(POPNUM(L,J)),J=J1,J2)
    9 FORMAT (I3,F7.3,2X,10F12.2)
    3 CONTINUE
      IF (J2.EQ.N) RETURN
      J1=J1+10
      GOTO 4
      END
