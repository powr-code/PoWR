      SUBROUTINE PRIOPA (XLAM,K,ND,LSTEP,R,
     $ OPA,ETA,THOMSON,IWARN,MAINPRO,MAINLEV,JOBNUM,MODHEAD)
C***********************************************************************
C***  PRINTOUT OF THE CONTINUUM OPACITIES ETC.
C***********************************************************************

      DIMENSION OPA(ND),ETA(ND),THOMSON(ND),R(ND),IWARN(ND)
      CHARACTER*10 MAINPRO(ND),MAINLEV(ND)
      CHARACTER MODHEAD*100
 
      IF (IHELP.EQ.5HIHELP) GOTO 1
      IHELP=5HIHELP
      PRINT 2,MODHEAD,JOBNUM
    2 FORMAT (1X,A,20X,'JOB NO.',I3,//,10X,
     $ 'CONTINUOUS OPACITY, EMISSIVITY AND SOURCE FUNCTION',
     $ /,10X,50('-'),//,
     $ ' FREQUENCY  DEPTH      OPACITY  THOMSON   OPTICAL    R(TAU=1)  '
     $ 'MAIN CONTRIBUTION      LASER    EMISSIVITY   SOURCE FUNCTION',/,
     $ '   INDEX    INDEX   (PER RSTAR) FRACTION   DEPTH               '
     $ 'PROCESS     LEVEL      WARNING    (...)        TRAD/KELVIN',/)
 
    1 PRINT 3
    3 FORMAT (1H )
      TAU=.0
      RTAU1=.0
 
      DO 6 L=1,ND
      IF (L.EQ.ND) GOTO 7
      TAUOLD=TAU
      TAU=TAU+0.5*(OPA(L)+OPA(L+1))*(R(L)-R(L+1))
      IF( TAUOLD.GE.1. .OR. TAU.LT.1. )  GOTO 7
      Q=(1.-TAUOLD)/(TAU-TAUOLD)
      RTAU1=(1.-Q)*R(L)+Q*R(L+1)
    7 IF(((L-1)/LSTEP)*LSTEP.NE.(L-1) .AND. L.NE.ND) GOTO 6
      SDENOM = OPA(L) * (1.-THOMSON(L))
C***  source function as radiation temperature; negative if S < 0
      IF (SDENOM == .0) THEN
        TRAD= .0
      ELSE
        S=ETA(L)/OPA(L)/(1.-THOMSON(L))
        TRAD = SIGN(TRADFUN (XLAM,ABS(S)), S)
      ENDIF
      IF (L.EQ.ND) THEN
            PRINT 9,K,L,OPA(L),THOMSON(L),TAU,RTAU1
     $                  ,MAINPRO(L),MAINLEV(L),IWARN(L),ETA(L),TRAD
            ELSE
            PRINT 5,K,L,OPA(L),THOMSON(L),TAUOLD
     $                  ,MAINPRO(L),MAINLEV(L),IWARN(L),ETA(L),TRAD
            ENDIF
    6 CONTINUE
      RETURN
 
    5 FORMAT (I6,I10,1P,E15.3,0P,F7.3,2X,1P,E10.2, 10X ,3X,
     $                          A10,5X,A10,5X,A1,1P,E11.3,0P,F13.0)
    9 FORMAT (I6,I10,1P,E15.3,0P,F7.3,2X,1P,E10.2,0P,F10.3,3X,
     $                          A10,5X,A10,5X,A1,1P,E11.3,0P,F13.0)
 
      END
