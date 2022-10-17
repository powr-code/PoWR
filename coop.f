      SUBROUTINE COOP (XLAM,ND,T,RNE,POPNUM,ENTOT,RSTAR,
     $                 OPA,ETA,THOMSON,IWARN,MAINPRO,MAINLEV,NOM,KODAT,
     $                 NDIM,N,MAXATOM,LEVEL,NCHARG,WEIGHT,ELEVEL,EION,
     $                 EINST,ALPHA,SEXPO,
     $                 ADDCON1, ADDCON2, ADDCON3, 
     $                 IGAUNT,SIGMATHK,SEXPOK,EDGEK,K,NF,SIGMAKI,RADIUS,
     $                 KONTNUP,KONTLOW,LASTKON,XDATA)
C***********************************************************************
C***  NON-LTE CONTINUOUS OPACITY AT GIVEN FREQUENCY FOR ALL DEPTH POINTS
C***  CALLED FROM VARIOUS PLACES: 
C***                  WRCONT
C***                  COMO
C***                  FORMAL
C***                  FORMAL - FORMCMF
C***                  STEAL - LINPOP - CCORE
C***                  STEAL - LINPOP - OPAROSS
C***                  WRCONT - DIFDTDR - OPAROSS
C***                  COMO - DIFDTDR - OPAROSS
C***                  COLI - DIFDTDR - OPAROSS
C***                  STEAL - TAUSCAL (2x) - OPAROSS
C***                  STEAL - ENSURETAUMAX - TAUSCAL (2x) - OPAROSS
C***  SIMILAR SUBROUTINES: CMFCOOP, COOPFRQ, DCOOP
C***  This version (23-Mar-2007) assumes that KODAT positions
C***  (i.e. KODATIND) give the atomic number (NCORECHARGE)
C***********************************************************************
 
C***  MAXIMUM ION CHARGE WHICH MAY OCCUR
      PARAMETER ( MAXION = 27 )
C***  MAXIMUM X-RAY DATA
      PARAMETER ( MAXXDAT = 10)
C***  Dimension of the core-charge data locally provided here
      PARAMETER (MAXATOMDIM = 26)

      DIMENSION XDATA(MAXXDAT)
      DIMENSION NCHARG(N),WEIGHT(N),ELEVEL(N),EION(N)
      DIMENSION NOM(N)
      DIMENSION KODAT(MAXATOM)
      DIMENSION SIGMATHK(MAXATOM,MAXATOM),SEXPOK(MAXATOM,MAXATOM)
      DIMENSION EDGEK(MAXATOM,MAXATOM)
      DIMENSION EINST(NDIM,NDIM)
      DIMENSION POPNUM(ND,N)
      DIMENSION OPA(ND),ETA(ND),THOMSON(ND)
      CHARACTER(8), DIMENSION(ND) :: IWARN
      DIMENSION RADIUS(ND)
      DIMENSION T(ND),RNE(ND),ENTOT(ND)
      DIMENSION KONTNUP(LASTKON),KONTLOW(LASTKON)
      DIMENSION SIGMAKI(NF,LASTKON)
      DIMENSION GFF(0:MAXION)
      DIMENSION KODATIND(MAXATOMDIM)
      CHARACTER*10 LEVEL(N),MAINPRO(ND),MAINLEV(ND)
      LOGICAL XRAYS, KSHELL

C***  Output of laser warnings for bound-free transitions
      DATA NWARN /0/ ! no warning has been issued yet
      SAVE NWARN

C***  C1 = H * C / K    ( CM * KELVIN )
      DATA C1 / 1.4388 /
C***  C2 = 2 * H * C    ( G * CM**3 / S**2 )
      DATA C2 / 3.9724E-16 /
C***  SIGMAE = ELCTRON SCATTERING CROSS SECTION  ( CM**2 )
      DATA SIGMAE / 0.6652E-24 /
C***  C3 = RECIPROCAL STATISTICAL WEIGHT OF FREE ELECTRON
      DATA C3 / 2.07E-16 /
C***  CFF = COEFFICIENT FOR FREE-FREE CROSS SECTION ( ALLEN PAGE 100, CM**5 )
      DATA CFF / 1.370E-23 /

      W  = 1.E8/XLAM
      W3 = W*W*W
 
C***  K-SHELL ABSORPTION CROSS-SECTIONS PROVIDED BY DATOM FILE ?
      KSHELL = .FALSE.
      DO NA=1, MAXATOM
        DO ISTATE=1, MAXATOM
         IF (SIGMATHK(NA,ISTATE) .NE. .0) THEN
            KSHELL = .TRUE.
            EXIT
         ENDIF
        ENDDO
      ENDDO

C***  X-RAY EMISSION SWITCHED ON BY CARD OPTION ?
      XRAYS = .FALSE.
      IF (XDATA(1) .NE. 0.) XRAYS = .TRUE.

C***  ONLY NEEDED FOR K-SHELL OR XRAY BRANCH
      IF (XRAYS .OR. KSHELL) THEN 
C***     Establish KODAT index for each used element
C***     First find number of used elements
         NATOMMAX = NOM(N)
         IF (MAXATOM .GT. MAXATOMDIM) THEN
            WRITE (0,*) '*** ERROR: MAXATOMDIM TOO SMALL'
            STOP 'ERROR IN COOP'
         ENDIF

C***     Now find for each NA the corresponding KODAT index
         DO NA=1, NATOMMAX
            KODATIND(NA) = 0
            DO J = 1, MAXATOM
               IF (NA .EQ. KODAT(J)) KODATIND(NA) = J 
            ENDDO
            IF (KODATIND(NA) .EQ. 0) THEN
               WRITE (0,*) '*** ERROR: ELEMENT NOT FOUND'
               STOP 'ERROR IN COOP'
            ENDIF
            IF (KODATIND(NA) .GT. MAXATOMDIM) THEN
               WRITE (0,*) '*** ERROR: NCORECHARGE NOT FOUND'
               STOP 'ERROR IN COOP'
            ENDIF
         ENDDO
      ENDIF

C***  PARAMETER FOR X-RAY SOURCE, AND SOME PREPARATIONS
      IF (XRAYS) THEN
         XFILL = XDATA(1)
         XRAYT = XDATA(2)
         XMINR = XDATA(3)
         DIFFEMEXP = XDATA(4)
         EXPFACXRAY = EXP(-C1*W/XRAYT)
         PRESIGXRAY = CFF / W3 / SQRT(XRAYT)
         IF (XDATA(5) .NE. 0.) THEN
            XFILL2 = XDATA(5)
            XRAYT2 = XDATA(6)
            XMINR2 = XDATA(7)
            EXPFACXRAY2 = EXP(-C1*W/XRAYT2)
            PRESIGXRAY2 = CFF / W3 / SQRT(XRAYT2)
         ENDIF
C***    Calculate number of free electrons, assuming full ionization
C***    Because of number conservation, this can be done for any depth point
         RNEXRAY = .0
         DO J = 1, N
          RNEXRAY = RNEXRAY + POPNUM(1,J) * KODATIND(NOM(J)) 
         ENDDO
      ENDIF
 
C***  LOOP OVER ALL DEPTH POINTS  --------------------------------------
      DO 1 L=1,ND
      NBFLASER = 0
   55 CONTINUE
      OPAMAX=.0
      OPAL=.0
      ETAL=.0
      IWARN(L)='        '
      TL=T(L)
      ROOTTL=SQRT(TL)
      T32=TL*ROOTTL
 
C***  BOUND-FREE  ******************************************************
C***  I = LOW      J = UP
      DO 5 KON=1,LASTKON
      J=KONTNUP(KON)
      I=KONTLOW(KON)
      EDGE=ELEVEL(J)+EION(I)-ELEVEL(I)
      IF (W .LT. EDGE) GOTO 5
 
C***  CALCULATE SIGMA, THE FREQUENCY-DEPENDENT CROSS SECTION
C***  IF ( K .GT. 0 ) IT IS ASSUMED THAT THE BOUND-FREE CROSS SECTIONS SIGMA
C***  HAVE BEEN ALREADY CALCULATED ( ARRAY SIGMAKI )
      IF (K .GT. 0) THEN
            SIGMA=SIGMAKI(K,KON)
            ELSE
            SIGMATH=EINST(I,J)*1.E-18
            CALL PHOTOCS (SIGMA,SIGMATH,EDGE,W,ALPHA,SEXPO,
     >                    ADDCON1, ADDCON2, ADDCON3, 
     >                    IGAUNT,KON)
            ENDIF
 
C***  RECIPROCAL STATISTICAL WEIGHT OF FREE ELECTRON
      WE=C3*RNE(L)*ENTOT(L)/T32
      G=WEIGHT(I)/WEIGHT(J)*WE*EXP(C1*(EDGE-W)/TL)
      EMINDU=G*POPNUM(L,J)*SIGMA

C***  Set emissivities zero if both levels are equal (=POPMIN)
      IF (POPNUM(L,I) .EQ. POPNUM(L,J)) EMINDU = .0

      SUM=POPNUM(L,I)*SIGMA-EMINDU
C***  LASER WARNING IF STIMULATED EMISSION EXCEEDS ABSORPTION IN THIS TRANSITION
C***  IF THIS BF opacity is negative, it is not added to the sum
      IF (SUM.LT. .0) THEN
         IWARN(L)='*       '
         IF (NWARN .EQ. 0)  THEN
            WRITE (0, 90) L, LEVEL(I), LEVEL(J)
   90       FORMAT ('*** WARNING FROM Subr. COOP: ',
     >       'LASERING BOUND-FREE CONTINUA SUPPRESSED',
     >       /,'*** THIS OCCURED FOR THE FIRST TIME AT DEPTH INDEX', I7, 
     >       /,'*** BETWEEN LEVELS ', A10, ' AND ', A10 )
             NWARN = NWARN + 1
         ENDIF
      ELSE
         OPAL=OPAL+SUM
         ETAL=ETAL+EMINDU
      ENDIF
      IF(SUM .LT. OPAMAX) GOTO 5
        OPAMAX=SUM
        MAINPRO(L)='BOUND-FREE'
        MAINLEV(L)=LEVEL(I)
    5 CONTINUE

C***  K-SHELL IONISATION  **********************************************
      IF (KSHELL) THEN      
C***     LOOP OVER ALL LEVELS 
         LASTNOMJ = -1
         LASTISTATE = -1
         DO 6 J=1,N
            NOMJ = NOM(J)
            ISTATE = NCHARG(J) + 1
C***        ARE THERE K-SHELL-DATA FOR CURRENT ELEMENT?
            IF (SIGMATHK(NOMJ,ISTATE) .EQ. 0.) GOTO 6

C***        IS RADIATION HARDER THAN K-SHELL EDGE ? 
            IF (W .LT. EDGEK(NOMJ,ISTATE)) GOTO 6

C***        K-SHELL IONIZATION NEEDS IONS WITH AT LEAST 3 ELECTRONS LEFT
            IF (KODATIND(NOMJ) - NCHARG(J) .LT. 3) THEN  
               WRITE (0,*) 'UNEXPECTED INCONSISTENCY WITH K-SHELL DATA'
               STOP 'ERROR in Subr. COOP'
            ENDIF

            IF (LASTNOMJ .NE. NOMJ .OR. LASTISTATE .NE. ISTATE) THEN
               CALL KSIGMA (SIGMAK, SIGMATHK(NOMJ,ISTATE), 
     >                   EDGEK(NOMJ,ISTATE), W, SEXPOK(NOMJ,ISTATE))
               LASTNOMJ = NOMJ
               LASTISTATE = ISTATE
            ENDIF

            SUM = POPNUM(L,J) * SIGMAK
            OPAL = OPAL + SUM
            IF (SUM .GT. OPAMAX) THEN
               OPAMAX=SUM
               MAINPRO(L) = 'K-SHELL'
               MAINLEV(L) = LEVEL(J)(:2)
               WRITE (MAINLEV(L)(4:5), '(I2)') ISTATE 
            ENDIF
    6    CONTINUE
      ENDIF

C***  X-RAY EMISSION -------------------------------------------------
      IF (XRAYS) THEN
C***     X-RAY SOURCE: FREE-FREE BREMSSTRAHLUNG 
         IF ((XFILL .GT. 0.) .AND. (RADIUS(L) .GE. XMINR)) THEN
            DO I=1, N
               NCHARI = KODATIND(NOM(I)) 
               SIGMAFF= PRESIGXRAY * FLOAT(NCHARI * NCHARI)
               SUM = RNEXRAY * ENTOT(L) * POPNUM(L,I) * SIGMAFF * XFILL
C***           For the differential emission measure option, 
C***           the exponential factor is replaced by a function 
C***           that depends on the DEM exponent - see documentation
               IF (DIFFEMEXP .EQ. .0) THEN
                  FDEM = EXPFACXRAY
               ELSE IF (DIFFEMEXP .EQ. 1.5) THEN
                  XXMAX = C1*W/XRAYT 
                  XXMIN = C1*W/TL
                  FDEM  = (EXPFACXRAY - EXP(-XXMIN)) / XXMAX
C*                normalization
                  FDEM = FDEM * 0.5 / ((XRAYT/TL)**0.5 - 1.)
               ELSE IF (DIFFEMEXP .EQ. 2.5) THEN
                  XXMAX = C1*W/XRAYT 
                  XXMIN = C1*W/TL
                  FDEM  = EXPFACXRAY - EXP(-XXMIN)
     >                  + XXMAX * EXPFACXRAY - XXMIN * EXP(-XXMIN)
                  FDEM = FDEM / (XXMAX*XXMAX)
C*                normalization
                  FDEM = FDEM * 1.5 / ((XRAYT/TL)**1.5 - 1.)
               ELSE
                  WRITE (0,*) '*** XRAY DIFFERENTIAL EMISSION MEASURE:' 
                  WRITE (0,*) '*** INVALID EXPONENT:', DIFFEMEXP
                  STOP '*** FATAL ERROR IN SUBR. COOP'
               ENDIF
               EMINDU = SUM * FDEM
               SUM = SUM - EMINDU
               OPAL = OPAL + SUM
               ETAL = ETAL + EMINDU 
               IF (SUM .GE. OPAMAX) THEN
                  OPAMAX=SUM
                  MAINPRO(L)='XRAYSOURCE'
                  MAINLEV(L)=LEVEL(I)
               ENDIF
            ENDDO
         ENDIF
         IF ((XFILL2 .GT. 0.) .AND. (RADIUS(L) .GE. XMINR2)) THEN
            DO I=1, N
               NCHARI = KODATIND(NOM(I)) 
               SIGMAFF= PRESIGXRAY2 * FLOAT(NCHARI * NCHARI)
               SUM = RNEXRAY * ENTOT(L) * POPNUM(L,I) * SIGMAFF *XFILL2
               EMINDU = SUM * EXPFACXRAY2
               SUM = SUM - EMINDU
               OPAL = OPAL + SUM
               ETAL = ETAL + EMINDU 
               IF (SUM .GE. OPAMAX) THEN
                  OPAMAX=SUM
                  MAINPRO(L)='XRAYSOURCE'
                  MAINLEV(L)=LEVEL(I)
               ENDIF
            ENDDO
         ENDIF
      ENDIF
 
C***  FREE-FREE  *******************************************************
C***  PRECALCULATE FREE-FREE GAUNT FACTORS FOR THE DIFFERENT ION CHARGES
      GFF(0)=.0
      DO 10 ION=1, MAXION
      CALL GAUNTFF (GIII,ION,XLAM,TL)
      GFF(ION)=GIII*FLOAT(ION*ION)
   10 CONTINUE

C***  PRECALCULATE SIGMAFF, LEAVING OUT THE FACTOR NCHARGE*NCHARGE
      PRESIG=CFF/W3/ROOTTL
      EXPFAC=EXP(-C1*W/TL)
      DO 3 I=1,N
      NCHARI=NCHARG(I)
      IF (NCHARI .GT. MAXION) THEN
         WRITE (0,*) '*** ERROR COOP: MAXION TOO SMALL'
         STOP '*** ERROR IN COOP'
      ENDIF
      SIGMAFF=PRESIG*GFF(NCHARI)
      SUM=RNE(L)*ENTOT(L)*POPNUM(L,I)*SIGMAFF
      EMINDU=SUM*EXPFAC
      SUM=SUM-EMINDU
      OPAL=OPAL+SUM
      ETAL=ETAL+EMINDU
      IF (SUM.LT.OPAMAX) GOTO 3
      OPAMAX=SUM
      MAINPRO(L)='FREE-FREE'
      MAINLEV(L)=LEVEL(I)
    3 CONTINUE

C***  If total true continuum opacity is negative: 
C***  re-do the whole calculation, but skip lasering bound-free transitions 
C***  Note: in the version with strict suppression (wrh  4-Apr-2019)
C***        this condition should never be met, and is therefore
C***        commented 
cc      IF (OPAL .LE. .0 .AND. NBFLASER .EQ. 0) THEN
cc         NBFLASER = 1 
cc         GOTO 55
cc      ENDIF

C***  ATENTION NOTICE: LINE IS IMPORTANT FOR ALL ETAL MODIFIED BEFORE 
      ETAL=ETAL*C2*W3
 
C***  THOMSON SCATTERING ***********************************************
      SUM=RNE(L)*SIGMAE
      IF (SUM.LT.OPAMAX) GOTO 4
      MAINPRO(L)='THOMSON'
      MAINLEV(L)='ELECTRON'
    4 OPAL=OPAL+SUM
C***  THOMSON = RELATIVE FRACTION FROM THE TOTAL OPACITY
      THOMSON(L)=SUM/OPAL
 
      OPA(L)=OPAL*ENTOT(L)*RSTAR
      ETA(L)=ETAL*ENTOT(L)*RSTAR
    1 CONTINUE
C***  ENDLOOP  ---------------------------------------------------------
 
      RETURN
      END
