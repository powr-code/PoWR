      SUBROUTINE PRIMINT (XJC,ND,XLAMBDA,NF,LSINT,JOBNUM,MODHEAD)
C***********************************************************************
C***  PRINT MEAN INTENSITIES
C***********************************************************************

      DIMENSION XJC(ND,NF),XLAMBDA(NF)
      LOGICAL KFIRST
      CHARACTER MODHEAD*100

      PRINT 2,MODHEAD,JOBNUM
    2 FORMAT (1X,  A  ,20X,'JOB NO.',I3,
     $            //,10X,'MEAN INTENSITY',/,10X,14('-'),/,
     $ ' FREQUENCY     WAVELENGTH    DEPTH       J-NUE',
     $ '          T-RAD   ',  /,
     $ '   INDEX       (ANGSTROEM)   INDEX     (ERG/CM+2)',
     $ '          (KELVIN)',/)
      DO 4 K=1,NF
      KFIRST=.TRUE.
      DO 1 L=1,ND
      IF (((L-1)/LSINT)*LSINT .NE. (L-1) .AND. L.NE.ND) GOTO 1
      XLAM=XLAMBDA(K)
      XJCLK=XJC(L,K)
      TRAD=TRADFUN (XLAM,XJCLK)
      IF (KFIRST) THEN
          PRINT 5, K,XLAM,L,XJCLK,TRAD
    5     FORMAT (I10,1X,F15.2,I8,E15.3,F15.0)
          KFIRST=.FALSE.
      ELSE
          PRINT 15, L,XJCLK,TRAD
   15     FORMAT (26X,I8,E15.3,F15.0)
      ENDIF
    1 CONTINUE
    4 PRINT 6
    6 FORMAT (1X)

      RETURN
      END
