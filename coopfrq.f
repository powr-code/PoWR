      SUBROUTINE COOPFRQ (NF,OPAC,ETAC,XLAMBDA,EXPFAC,SIGMAKI,N,NCHARG,
     $                   WEIGHT,ELEVEL,EION,NFEDGE,EN,NOM,RSTAR,ENTOTL,
     $                   RNEL,TL,SIGMAFF,MAXION,RL,XDATA,
     $                   SIGMATHK,SEXPOK,EDGEK,KODAT,MAXATOM,
     $                   KONTNUP,KONTLOW,LASTKON,OPATHOM)
C***********************************************************************
C***  NON-LTE CONTINUOUS OPACITY AT CURRENT DEPTH FOR ALL FREQUENCIES
C***  NOTE: ONLY TRUE OPACITY, WITHOUT THOMSON SCATTERING TERM.
C***  THOMSON OPACITY IS PREPARED AS AN EXTRA VARIABLE: OPATHOM
C***  This version (23-Mar-2007) assumes that KODAT positions
C***  (i.e. KODATIND) give the atomic number (NCORECHARGE)
C***  Called from: STEAL --> LINPOP --> COMA
C***          and: WRSTART --> GREY --> OPAGREY
C***********************************************************************
 
      DIMENSION NOM(N)
      DIMENSION SIGMATHK(MAXATOM,MAXATOM),SEXPOK(MAXATOM,MAXATOM)
      DIMENSION EDGEK(MAXATOM,MAXATOM)
      DIMENSION XDATA(10)
      DIMENSION KODAT(MAXATOM)
      DIMENSION NCHARG(N),WEIGHT(N),ELEVEL(N),EION(N),EN(N)
      DIMENSION KONTNUP(LASTKON),KONTLOW(LASTKON),NFEDGE(LASTKON)
      DIMENSION OPAC(NF),ETAC(NF),XLAMBDA(NF),EXPFAC(NF)
      DIMENSION SIGMAKI(NF,LASTKON),SIGMAFF(NF,0:MAXION)
      LOGICAL XRAYS, KSHELL
C***  Dimension of the core-charge data locally provided here
      PARAMETER ( MAXATOMDIM = 26)
      DIMENSION KODATIND(MAXATOMDIM)

C***  Output of laser warnings for bound-free transitions
      DATA NWARN /0/ ! no warning has been issued yet
      SAVE NWARN
 
C***  C1 = H * C / K    ( CM * ANGSTROEM )
      DATA C1 / 1.4388 /
C***  C2 = 2 * H * C    ( G * CM**3 / S**2 )
      DATA C2 / 3.9724E-16 /
C***  C3 = RECIPROCAL STATISTICAL WEIGHT OF FREE ELECTRON
      DATA C3 / 2.07E-16 /
C***  CFF = COEFFICIENT FOR FREE-FREE CROSS SECTION (ALLEN PAGE 100)
      DATA CFF / 1.370E-23 /
C***  SIGMAE = ELCTRON SCATTERING CROSS SECTION  ( CM**2 )
      DATA SIGMAE / 0.6652E-24 /


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

C***  PARAMETER FOR X-RAY SOURCE
      XRAYS = .FALSE.
      IF (XDATA(1) .NE. 0.) THEN
        XRAYS = .TRUE.
        XFILL = XDATA(1)
        XRAYT = XDATA(2)
        XMINR = XDATA(3)
        DIFFEMEXP = XDATA(4)
        IF (XDATA(5) .NE. 0.) THEN
           XFILL2 = XDATA(5)
           XRAYT2 = XDATA(6)
           XMINR2 = XDATA(7)
        ENDIF
      ENDIF

      T32=TL*SQRT(TL)
 
ccC***  Safety against negative opacities from lasering bound-free trans.
ccC***  Note: In contrast to the other opacity routines, in COOPFRQ
ccC***   all lasering b-f continue are suppressed *at all frequencies* 
ccC**    if a negative "true" opacity is encountered at any frequency 
cc      NBFLASER = 0
cc   55 CONTINUE

      DO 10 K=1,NF
      OPAC(K)=.0
      ETAC(K)=.0
   10 CONTINUE
 
C***  BOUND-FREE  ******************************************************

C***  LOOP OVER ALL CONTINUUM TRANSITIONS
      DO 5 KON=1,LASTKON
      NUP=KONTNUP(KON)
      LOW=KONTLOW(KON) 
      EDGE=ELEVEL(NUP)+EION(LOW)-ELEVEL(LOW)
      IF (C1*EDGE/TL .GT. 700.) THEN
        WE = 0.
      ELSE
        EXPEDGE=EXP(C1*EDGE/TL)
        WE=C3*RNEL*ENTOTL/T32 *WEIGHT(LOW)/WEIGHT(NUP)*EXPEDGE
      ENDIF
      NFLOW=NFEDGE(KON)

C***  LOOP OVER ALL CONTINUUM FREQUENCY POINTS WITH XLAMBDA(K) < EDGE
      DO 11 K=1,NFLOW
      SIGMA=SIGMAKI(K,KON)
C***  RECIPROCAL STATISTICAL WEIGHT OF FREE ELECTRON
      G=WE*EXPFAC(K)
      EMINDU=G * SIGMA * EN(NUP)

C***  Set emissivities zero if both levels are equal (=POPMIN)
      IF (EN(LOW) .EQ. EN(NUP)) EMINDU = .0

      SUM=    EN(LOW)*SIGMA-EMINDU

C***  LASER bf continua are skipped!  3-Feb-2016
cc      IF (SUM .LT. .0 .AND. NBFLASER .EQ. 1) THEN
cc         IF (NWARN .EQ. 0)   WRITE (0, 90) 
cc   90    FORMAT ('*** WARNING FROM Subr. COOPFRQ: ',
cc     >   'LASERING BOUND-FREE CONTINUA SUPPRESSED')
cc         NWARN = NWARN + 1
cc      ELSE
      IF (SUM .GT. .0) THEN
         OPAC(K)=OPAC(K)+SUM
         ETAC(K)=ETAC(K)+EMINDU
      ENDIF
   11 CONTINUE
 
    5 CONTINUE

C***  ONLY NEEDED FOR K-SHELL OR XRAY BRANCH
      IF (XRAYS .OR. KSHELL) THEN 
C***     Establish KODAT index for each used element
C***     First find number of used elements
         NATOMMAX = NOM(N)
         IF (MAXATOM .GT. MAXATOMDIM) THEN
            WRITE (0,*) '*** ERROR: MAXATOMDIM TOO SMALL'
            STOP 'ERROR IN COOPFRQ'
         ENDIF
C***     Now find for each NA the corresponding KODAT index
C***     in order to know the core charge
         DO NA=1, NATOMMAX
            KODATIND(NA) = 0
            DO II = 1, MAXATOM
               IF (NA .EQ. KODAT(II)) KODATIND(NA) = II
            ENDDO
            IF (KODATIND(NA) .EQ. 0) THEN
               WRITE (0,*) '*** ERROR: ELEMENT NOT FOUND'
               STOP 'ERROR IN COOPFRQ'
            ENDIF
            IF (KODATIND(NA) .GT. MAXATOMDIM) THEN
               WRITE (0,*) '*** ERROR: NCORECHARGE NOT FOUND'
               STOP 'ERROR IN COOPFRQ'
            ENDIF
            IF (KODATIND(NA) .GT. MAXATOMDIM) THEN
               WRITE (0,*) '*** ERROR: NCORECHARGE NOT FOUND'
               STOP 'ERROR IN COOPFRQ'
            ENDIF
         ENDDO
      ENDIF

C***  K-SHELL IONISATION  **********************************************
      IF (KSHELL) THEN
C***  LOOP OVER ALL LEVELS 
      DO 6 J=1,N
         NOMJ = NOM(J)
         ISTATE = NCHARG(J) + 1
C***     ARE THERE K-SHELL-DATA FOR CURRENT ELEMENT?
         IF (SIGMATHK(NOMJ,ISTATE) .EQ. 0.) GOTO 6

C***     K-SHELL IONIZATION NEEDS IONS WITH AT LEAST 3 ELECTRONS LEFT
         IF (KODATIND(NOMJ) - NCHARG(J) .LT. 3) THEN  
            WRITE (0,*) 'UNEXPECTED INCONSISTENCY WITH K-SHELL DATA'
            STOP 'ERROR in Subr. COOPFRQ'
         ENDIF

         WK = 1.E8 / EDGEK(NOMJ,ISTATE)
         DO K=1,NF
C***        IS RADIATION TOO SOFT FOR K-SHELL-IONISATION? 
            IF (XLAMBDA(K) .GT. WK) EXIT
            W = 1.E8 / XLAMBDA(K)
            CALL KSIGMA (SIGMAK, SIGMATHK(NOMJ,ISTATE), 
     >                   EDGEK(NOMJ,ISTATE), W, SEXPOK(NOMJ,ISTATE))
  
            SUM = EN(J) * SIGMAK
            OPAC(K) = OPAC(K) + SUM
         ENDDO
    6 CONTINUE
      ENDIF

C***  X-RAY SOURCE: FREE-FREE BREMSSTRAHLUNG HeIII ASSUMED
      IF (XRAYS) THEN
      IF ((XFILL .GT. 0.) .AND. (RL .GE. XMINR)) THEN
        DO 15,K=1,NF
          W = 1.E8/XLAMBDA(K)
          W3 = W * W * W
          EXPFACXRAY = EXP(-C1*W/XRAYT)
          SUM    = 2. * ENTOTL * 4. * CFF / W3 / SQRT(XRAYT)
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
                  STOP '*** FATAL ERROR IN SUBR. COOPFRQ'
               ENDIF
               EMINDU = SUM * FDEM
C*             In the DEM branch, the opacity is set to zero 
               IF (DIFFEMEXP .EQ. .0) THEN
                 OPAX   = (SUM-EMINDU) * XFILL
                 OPAC(K)=OPAC(K)+OPAX
               ENDIF
          ETAX   = EMINDU * XFILL
          ETAC(K)=ETAC(K)+ETAX
   15   CONTINUE
      ENDIF
      IF ((XFILL2 .GT. 0.) .AND. (RL .GE. XMINR2)) THEN
        DO 16,K=1,NF
          W = 1.E8/XLAMBDA(K)
          W3 = W * W * W
          EFACTX = EXP(-C1*W/XRAYT2)
          SUM    = 2. * ENTOTL * 4. * CFF / W3 / SQRT(XRAYT2)
          EMINDU = SUM * EFACTX
          OPAX   = (SUM-EMINDU) * XFILL2
          ETAX   = EMINDU * XFILL2
          OPAC(K)=OPAC(K)+OPAX
          ETAC(K)=ETAC(K)+ETAX
 16     CONTINUE
      ENDIF
      ENDIF


C***  FREE-FREE  *******************************************************
C***  Loop 3 is most time consuming in WRSTART (53 percent
C***  Tested on 30-Jan-1997 19:52:38, Lars
      DO 30 K=1,NF
      W=1.E8/XLAMBDA(K)
      W3=W*W*W
      DO 3 I=1,N
      SUM=RNEL*ENTOTL*EN(I)*SIGMAFF(K,NCHARG(I))
      EMINDU=SUM*EXPFAC(K)
      SUM=SUM-EMINDU
      OPAC(K)=OPAC(K)+SUM
      ETAC(K)=ETAC(K)+EMINDU
    3 CONTINUE
      ETAC(K)=ETAC(K)*C2*W3*ENTOTL*RSTAR
      OPAC(K)=OPAC(K)*ENTOTL*RSTAR
   30 CONTINUE

ccC***  If total true continuum opacity is negative at at least one freq.:
ccC***  re-do the whole calculation, but skip lasering bound-free trans.
cc      OPACMIN = OPAC(1) 
cc      DO K=2, NF
cc         IF (OPAC(K) .LT. OPACMIN) OPACMIN = OPAC(K) 
cc      ENDDO
cc      IF (OPACMIN .LE. .0 .AND. NBFLASER .EQ. 0) THEN
cc         NBFLASER = 1
cc         GOTO 55
cc      ENDIF

C***  THOMSON SCATTERING OPACITY (FREQUENCY INDEPENDENT! ) ************
      OPATHOM = RNEL * SIGMAE * ENTOTL * RSTAR
 
      RETURN
      END
