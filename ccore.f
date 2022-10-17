      SUBROUTINE CCORE (WCHARM,NF,GAMMAC,DELTAC,IPRICC,MODHEAD,JOBNUM,
     $      SCOLD,RADIUS,XLAMBDA,ND,T,RNE,POP1,ENTOT,RSTAR,
     $      OPA,ETA,THOMSON,IWARN,MAINPRO,MAINLEV,NOM,
     $      NDIM,N,LEVEL,NCHARG,WEIGHT,ELEVEL,EION,EINST,SIGMAKI,
     $      MAXATOM,SIGMATHK,SEXPOK,EDGEK,KODAT,XDATA,
     $      KONTNUP,KONTLOW,LASTKON, FILLFAC, WJC, OPC, WJCMIN,
     >      BPLOCC, LPLOCC, KPLOCC, KANAL)

 
C*******************************************************************************
C***  - CALLED FROM: SUBROUTINE LINPOP
C***  - DETERMINES THE FREQUENCY WEIGHTS WCHARM REPRESENTING THE APPROXIMATE
C***       LAMBDA-OPERATORS AT EACH CONTINUUM FREQUENCY POINT
C***  - NOTE: THIS SUBROUTINE ALSO PROVIDES "SCOLD" = THE CONTINUUM
C***          SOURCE FUNCTION CALCULATED FROM OLD POP.NUMBERS 
C***  Clumping: CCORE should be called with ENTOT = clump density
C***  OPA and ETA are scaled back with FILLFAC for evaluating WCHARM
C***  This does not affect SCOLD, as the scaling factor cancels out
C*******************************************************************************
 
      DIMENSION WCHARM(ND,NF),SCOLD(NF,ND), WJC(ND,NF)
      DIMENSION RADIUS(ND),OPA(ND),ETA(ND),THOMSON(ND),XLAMBDA(NF)
      CHARACTER*8 VERSION, OPC
      LOGICAL BPLOCC

C***  VALID VERSIONS: 'LOCAL', 'GLOBAL', 'OAB'
      VERSION = 'LOCAL'
      IF (DELTAC .LT. .0) VERSION='OAB'
 
      IF (OPC(1:7) .EQ. 'DIAGTAU') THEN
        VERSION = 'DIAGTAU'
      ELSE IF (OPC(1:4) .EQ. 'DIAG') THEN
        VERSION = 'DIAG'
      ENDIF
      write (0,*) 'VERSION=', VERSION

C***  LOOP OVER ALL CONTINUUM FREQUENCIES  *************************************
      DO 4 K=1,NF
      XLAM=XLAMBDA(K)
      CALL COOP (XLAM,ND,T,RNE,POP1,ENTOT,RSTAR,
     $           OPA,ETA,THOMSON,IWARN,MAINPRO,MAINLEV,NOM,KODAT,
     $           NDIM,N,MAXATOM,LEVEL,NCHARG,WEIGHT,ELEVEL,EION,EINST,
     $           DUMMY,DUMMY,DUMMY,
     $           DUMMY,DUMMY,DUMMY,
     $           SIGMATHK,SEXPOK,EDGEK,
     $           K,NF,SIGMAKI,RADIUS,
     $           KONTNUP,KONTLOW,LASTKON,XDATA)                  
 
      DO L=1, ND
        OPA(L) = OPA(L) * FILLFAC
        ETA(L) = ETA(L) * FILLFAC
      ENDDO

C***  BRANCH FOR  GAMMAC .EQ. .0  -----------------
      IF (GAMMAC .EQ. .0) THEN
         DO 11 L=1,ND
         WCHARM(L,K)=.0
   11    CONTINUE
         GOTO 10
         ENDIF
C***  ---------------------------------------------

C******************************
C***  'GLOBAL' - VERSION:   ***
C******************************

      IF (VERSION .EQ. 'GLOBAL') THEN
C***  FIRST: WCHARM IS FILLED WITH THE OPTICAL DEPTH TAU-NUE
C***     = RADIAL OPTICAL DEPTH TOWARDS THE NEXT BOUNDARY
      WCHARM(1,K)=.0
      DO 5 L=2,ND
    5 WCHARM(L,K)=WCHARM(L-1,K) +
     +                      .5*(RADIUS(L-1)-RADIUS(L))*(OPA(L-1)+OPA(L))

      TAUMAX=WCHARM(ND,K)
      DO 6 L=1,ND
    6 WCHARM(L,K)=AMIN1(WCHARM(L,K), TAUMAX-WCHARM(L,K))
      

C****************************** 
C***  'LOCAL' - VERSION:    ***
C******************************

      ELSE IF (VERSION .EQ. 'LOCAL') THEN
      DO 7 L=1,ND
C***  TAU=LOCAL OPACITY * TYPICAL DISTANCE
      DR=AMIN1(RADIUS(1)-RADIUS(L),RADIUS(L)-RADIUS(ND))
      WCHARM(L,K)=OPA(L)*AMIN1(DR,1.)
    7 CONTINUE


C**************************
C***  'OAB' - VERSION   ***
C**************************

      ELSE IF (VERSION .EQ. 'OAB') THEN
C***  BOUNDARIES:
      WCHARM ( 1,K) = .0
      WCHARM (ND,K) = .0
C***  INNER POINTS:
      DO 3 L=2, ND-1
      DTP = 0.5*(OPA(L+1)+OPA(L)) * (RADIUS(L)-RADIUS(L+1)) / GAMMAC
      DTM = 0.5*(OPA(L-1)+OPA(L)) * (RADIUS(L-1)-RADIUS(L)) / GAMMAC
      DT = 0.5 * (DTM + DTP)
      BP = (1. - EXP(-DTP)) / (DT*DTP)
      BM = (1. - EXP(-DTM)) / (DT*DTM)
      WCHARM(L,K) = 1. / (1. + BP + BM)
      IF (WCHARM(L,K) .LT. WJCMIN) WCHARM(L,K)=.0
    3 CONTINUE


C************************************
C***  DIAG or DIAGTAU - VERSION   ***
C************************************
C***  For Gamma damping, the WCHARM-like factor WJC is first transformed in 
C***    a TAU-like quantity, then reduced and transformed back
      ELSE IF (VERSION .EQ. 'DIAG' .OR. VERSION .EQ. 'DIAGTAU') THEN
C***    BOUNDARIES: No amplification
        WCHARM ( 1,K) = .0
        WCHARM (ND,K) = .0
C***    INNER POINTS:
        DO L=2, ND-1
          TAU = -ALOG(1.-WJC(L,K))
          WCHARM(L,K) = 1.-EXP(-TAU/GAMMAC)
          IF (WCHARM(L,K) .LT. WJCMIN) WCHARM(L,K)=.0
        ENDDO

C*************************
C***  INVALID VERSION  ***
C*************************
      ELSE
      CALL REMARK ('VERSION NOT DECLARED')
      STOP 'ERROR'
      ENDIF

C***  VERSIONS 'GLOBAL' AND 'LOCAL': 
C***    GENERATE THE WEIGHT FACTORS FROM THE OPTICAL DEPTHS
      IF (VERSION .EQ. 'GLOBAL' .OR. VERSION .EQ. 'LOCAL') THEN
      DO 2 L=1,ND
      IF (WCHARM(L,K) .GT. .0) WCHARM(L,K)=WCHARM(L,K)**DELTAC
      IF (WCHARM(L,K) .GT. GAMMAC*0.1 ) THEN
            WCHARM(L,K)=1.-EXP(-WCHARM(L,K)/GAMMAC)
            ELSE
            WCHARM(L,K)=.0
            ENDIF
    2 CONTINUE
      ENDIF

 
C***  SOURCE FUNCTION WITH OLD POP.NUMBERS (WITHOUT THOMSON OPACITY)
   10 CONTINUE
      DO 1 L=1,ND
C***  LASER SECURITY
      OPAL=OPA(L)*(1.-THOMSON(L))
      IF (OPAL .GT. .0) THEN
            SCOLD(K,L)=ETA(L)/OPAL
            ELSE
            WCHARM(L,K)=.0
            SCOLD(K,L)=.0
            ENDIF
    1 CONTINUE
 
    4 CONTINUE
C***  ENDLOOP  *****************************************************************
 
      IF (IPRICC .EQ. 1)
     $      CALL PRICC (ND,NF,WCHARM,GAMMAC,DELTAC,MODHEAD,JOBNUM)
 
      write (0,*) 'CCORE: BPLOCC=', BPLOCC
      IF (BPLOCC) THEN
        CALL PLOCC (LPLOCC, KPLOCC, 
     >              ND, NF, WCHARM, GAMMAC, DELTAC, MODHEAD, JOBNUM, 
     >              OPC, KANAL)
      ENDIF

      RETURN
      END
