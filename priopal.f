      SUBROUTINE PRIOPAL(KARTE,XLAM,ND,OPA,OPAL,ETA,ETAL,R,JOBNUM,LSOPA,
     $          MODHEAD)
C***********************************************************************
C***  PRINTOUT OF THE LINE OPACITIES ETC.
C***********************************************************************

      DIMENSION OPA(ND),OPAL(ND),ETA(ND),ETAL(ND),R(ND)
      CHARACTER MODHEAD*100, KARTE*80
C***  WPI = SQRT(PI)
      DATA WPI /1.772454/
      DATA IHELP / 0 /
 
C***  Print header only once at the beginning
      
      IF (IHELP .EQ. 0) GOTO 1
      IHELP=1
      PRINT 2,MODHEAD,JOBNUM
    2 FORMAT (1H1,1X,  A  , 4X,'AFTER JOB NO.',I3,//,10X,
     $ 'LINE OPACITY, EMISSIVITY AND SOURCE FUNCTION',
     $ /,10X,44('-'),//,
     $ ' LINE       DEPTH LINE OPACITY LINE/CONT.   TOTAL OPT.DEPTH  ',
     $ 'R (TAU=1)  LINE EMISS.   LINE SOURCE F.    TOTAL SOURCE F.',/,
     $ '            INDEX  (PER RSTAR)               (LINE CENTER)   ',
     $ '                          TRAD/KELVIN       TRAD/KELVIN   ',/)

    1 PRINT 3
    3 FORMAT (1H )
      TAU=.0
      RTAU1=.0
      DO 6 L=1,ND
      OPALC=OPAL(L)/WPI
      ETALC=ETAL(L)/WPI
      OPATOT=OPA(L)+OPALC
      ETATOT=ETA(L)+ETALC
      IF (L.EQ.ND) GOTO 7
      TAUOLD=TAU
      OPATOTP=OPA(L+1)+OPAL(L+1)/WPI
      TAU=TAU+0.5*(OPATOT+OPATOTP)*(R(L)-R(L+1))
      IF( TAUOLD.GE.1. .OR. TAU.LT.1. )  GOTO 7
      Q=(1.-TAUOLD)/(TAU-TAUOLD)
      RTAU1=(1.-Q)*R(L)+Q*R(L+1)
    7 IF(((L-1)/LSOPA)*LSOPA.NE.(L-1) .AND. L.NE.ND) GOTO 6
      S=ETALC/OPALC
      TRADLIN=TRADFUN (XLAM,S)
      S=ETATOT/OPATOT
      TRADTOT=TRADFUN (XLAM,S)
      IF (L.EQ.ND) GOTO 8
      PRINT 5,      L,OPALC,OPALC/OPA(L),TAUOLD,   ETALC,TRADLIN,TRADTOT
    5 FORMAT ( 10X, I6,1P,3E14.2,    11X,    E14.3,0P,2F16.0)
    6 CONTINUE
    8 PRINT 9,KARTE,L,OPALC,OPALC/OPA(L),TAU,RTAU1,ETALC,TRADLIN,TRADTOT
    9 FORMAT (1X,A9,I6,1P,3E14.2,0P,F11.3,1P,E14.3,0P,2F16.0)
      RETURN
      END
