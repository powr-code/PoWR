      SUBROUTINE CMFCOOP (XLAM,ND,T,RNE,POPNUM,ENTOT,RSTAR, LEVEL,
     >       OPA,ETA,THOMSON,NDIM,N,NCHARG,WEIGHT,ELEVEL,EION,
     >       EINST,ALPHA,SEXPO, ADDCON1, ADDCON2, ADDCON3, 
     >       IGAUNT,SIGMA1I,KONTLOW,KONTNUP,LASTKON,NATOM,KONTHLP, 
     >       DENSCON,BPLOT,BPLOT2,IPLOT,K,KCL,KCU,KCDELTA,
     >       OPACL,OPACU,ETACL,ETACU,XLAM0LN,ALN, MAXXDAT, XDATA, 
     >       SIGMATHK,SEXPOK,EDGEK, MAXATOM, NOM, KODAT, RADIUS)
C***********************************************************************
C***  NON-LTE CONTINUOUS OPACITY AT GIVEN FREQUENCY FOR ALL DEPTH POINTS
C***  OPTIMIZED VERSION OF SUBR. COOP, CALLED FROM MAIN PROGRAM CMF
C***  This version (23-Mar-2007) assumes that KODAT positions
C***  (i.e. KODATIND) give the atomic number (NCORECHARGE)
C***********************************************************************
 
C***  MAXIMUM ION CHARGE WHICH MAY OCCUR
      PARAMETER ( MAXION = 27 )

      DIMENSION NCHARG(N),WEIGHT(N),ELEVEL(N),EION(N)
      DIMENSION NOM(N)
      DIMENSION KODAT(MAXATOM)
      DIMENSION SIGMATHK(MAXATOM,MAXATOM),SEXPOK(MAXATOM,MAXATOM)
      DIMENSION EDGEK(MAXATOM,MAXATOM)
      DIMENSION KONTLOW(LASTKON),KONTNUP(LASTKON),KONTHLP(LASTKON)
      DIMENSION EINST(NDIM,NDIM)
      DIMENSION POPNUM(ND,N)
      DIMENSION OPA(ND),ETA(ND),THOMSON(ND), RADIUS(ND)
      DIMENSION T(ND),RNE(ND),ENTOT(ND)
      DIMENSION SIGMA1I(LASTKON)
      DIMENSION GFF(0:MAXION)
      DIMENSION OPACL(ND),OPACU(ND)
      DIMENSION ETACL(ND),ETACU(ND)
      DIMENSION XDATA(MAXXDAT)
      CHARACTER(10), DIMENSION(N) :: LEVEL
      LOGICAL XRAYS, KSHELL

      LOGICAL BPLOT,BPLOT2,BFCALC

C***  Output of laser warnings for bound-free transitions
      DATA NWARN /0/ ! no warning has been issued yet
      SAVE NWARN

C***  tiefenabh. clumping nach goetz
      DIMENSION DENSCON(ND)

C***  Dimension of the core-charge data locally provided here
      PARAMETER (MAXATOMDIM = 26)
      DIMENSION KODATIND(MAXATOMDIM)
      
      REAL :: XLAM, W, W3
 
C***  C1 = H * C / K    ( CM * ANGSTROEM )
      DATA C1 / 1.4388 /
C***  C2 = 2 * H * C    ( G * CM**3 / S**2 )
      DATA C2 / 3.9724E-16 /
C***  SIGMAE = ELCTRON SCATTERING CROSS SECTION  ( CM**2 )
      DATA SIGMAE / 0.6652E-24 /
C***  C3 = RECIPROCAL STATISTICAL WEIGHT OF FREE ELECTRON
      DATA C3 / 2.07E-16 /
C***  CFF = COEFFICIENT FOR FREE-FREE CROSS SECTION ( ALLEN PAGE 100 )
      DATA CFF / 1.370E-23 /
 
      W=1.E8/XLAM
      W3=W*W*W

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
               STOP 'ERROR IN CMFCOOP'
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

C***  This version has been optimized by Goetz in oder to avoid 
C***  the calculation of cont. opacities at each fine-frequency point. 
C***  Instead, the opacity is calculated only for K indices which are 
C***  spaced wider than by one counter.    
C***  Default is not to calculate BF-Opacity
      BFCALC = .FALSE.

C***  New Interval if K is outside the actual Interval
      IF (K .GT. KCU) THEN
C         WRITE (39,'(I6,1X,F3.0)') KCL, -3.
         KCL = KCU
C         WRITE (39,'(I6,1X,F3.0)') KCL, -3.
         DO L=1, ND
            OPACL(L) = OPACU(L)
            ETACL(L) = ETACU(L)
         ENDDO
         KCU = KCU+KCDELTA
         IF (BPLOT) WRITE (39,'(I6,1X,F3.0)') KCU, -3.
         BFCALC = .TRUE.
      ENDIF
C***  New Start of Interpolation Procedure if (KCL = KCU)
      IF (KCL .EQ. KCU) THEN
         BFCALC = .TRUE.
         IF (BPLOT) WRITE (39,'(I6,1X,F3.0)') KCU, -3.
      ENDIF

C***  PRECALCULATE SIGMA1I, THE PHOTOIONIZATION CROSS SECTION AT XLAMKCU
C***  FOR ALL CONTINUA WITH  (EDGE .LE. WU)
      IF (BFCALC) THEN
         XLAMKCU = EXP(XLAM0LN + FLOAT(KCU)*ALN)
         WU=1.E8/XLAMKCU
         W3U=WU*WU*WU
         KHELP=0
         DO 17 KON=1,LASTKON
            NUP=KONTNUP(KON)
            LOW=KONTLOW(KON)
            EDGE=ELEVEL(NUP)+EION(LOW)-ELEVEL(LOW)
ccc         in der folgenden Zeile W --> WU verbessert (?) - wrh  6-Dec-2006 
            IF (WU .LT. EDGE) GOTO 17
            KHELP=KHELP+1
            KONTHLP(KHELP)=KON
            SIGMATH=EINST(LOW,NUP)*1.E-18
            CALL PHOTOCS (SIGMA1I(KON),SIGMATH,EDGE,WU,ALPHA,SEXPO,
     >                    ADDCON1, ADDCON2, ADDCON3, 
     >                    IGAUNT,KON)
 17      CONTINUE
      ENDIF

C***  LOOP OVER ALL DEPTH POINTS  --------------------------------------
      DO 1 L=1,ND
      NBFLASER = 0
   55 CONTINUE
      OPAL=.0
      ETAL=.0
      TL=T(L)
      ROOTTL=SQRT(TL)
      T32=TL*ROOTTL

      ENEL = RNE(L) * ENTOT(L) * DENSCON(L)  
 
C***  BOUND-FREE  ******************************************************
C***  I = LOW      J = UP
      IF (BFCALC) THEN
         OPACU(L) = 0.
         ETACU(L) = 0.
         DO 5 KH=1,KHELP
            KONT=KONTHLP(KH)
            J=KONTNUP(KONT)
            I=KONTLOW(KONT)
            EDGE=ELEVEL(J)+EION(I)-ELEVEL(I)
 
C***  RECIPROCAL STATISTICAL WEIGHT OF FREE ELECTRON
            WE=C3 * ENEL / T32
            G=WEIGHT(I)/WEIGHT(J)*WE*EXP(C1*(EDGE-WU)/TL)
            SIGMA=SIGMA1I(KONT)
            EMINDU=G*POPNUM(L,J)*SIGMA
C***        Set emissivities zero if both levels are equal (=POPMIN)
            IF (POPNUM(L,I) .EQ. POPNUM(L,J)) EMINDU = .0

            SUM=POPNUM(L,I)*SIGMA - EMINDU
            IF (BPLOT2) THEN
               IF (L .EQ. IPLOT) THEN
                  SUMBF = SUM
                  EMINDUBF = EMINDU
               ENDIF
            ENDIF


ccC***  LASER WARNING IF STIMULATED EMISSION EXCEEDS ABSORPTION 
ccC***    IN THIS TRANSITION
ccC***  IF TOTAL CONT. OPA WAS < 0 AT FIRST TRIAL. THIS BF TRANSITION IS SKIPPED
cc      IF (SUM .LT. .0 .AND. NBFLASER .EQ. 1) THEN

C***    Note: the above version, which makes use of NBLASER. has been
C***          replaced by the more radical version that any lasering b-f 
C***          transition is immediately disregarded - wrh  4-Apr-2019
        IF (SUM .LT. .0) THEN 
          IF (NWARN .EQ. 0)   WRITE (0, 90) L, LEVEL(I), LEVEL(J)
   90     FORMAT ('*** WARNING FROM Subr. CMFCOOP: ',
     >    'LASERING BOUND-FREE CONTINUA SUPPRESSED',
     >    /,'*** THIS OCCURED FOR THE FIRST TIME AT DEPTH INDEX', I7,
     >    /,'*** BETWEEN LEVELS ', A10, ' AND ', A10 )
         NWARN = NWARN + 1
        ELSE
           OPACU(L)=OPACU(L)+SUM
           ETACU(L)=ETACU(L)+EMINDU
        ENDIF

 5       CONTINUE


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
            OPACU(L) = OPACU(L) + SUM
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
C***  The following statement implies that the filling factor refers 
C***      to material of clump density. The reasons are given in a WR-Memo 
C***      -- Inserted by wrh  2-Aug-2006 19:32:10
               SUM = SUM * DENSCON(L)
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
C*             In the DEM branch, the opacity is set to zero
               IF (DIFFEMEXP .EQ. .0)
     >               OPACU(L) = OPACU(L) + (SUM-EMINDU)
               ETACU(L) = ETACU(L) + EMINDU 
            ENDDO
         ENDIF
C***     X-RAY SOURCE: FREE-FREE BREMSSTRAHLUNG 2
         IF ((XFILL2 .GT. 0.) .AND. (RADIUS(L) .GE. XMINR2)) THEN
            DO I=1, N
               NCHARI = KODATIND(NOM(I))
               SIGMAFF= PRESIGXRAY2 * FLOAT(NCHARI * NCHARI)
               SUM = RNEXRAY * ENTOT(L) * POPNUM(L,I) * SIGMAFF *XFILL2
C***  The following statement implies that the filling factor refers 
C***      to material of clump density. The reasons are given in a WR-Memo 
C***      -- Inserted by wrh  2-Aug-2006 19:32:10
               SUM = SUM * DENSCON(L)
               EMINDU = SUM * EXPFACXRAY2
               SUM = SUM - EMINDU
               OPACU(L) = OPACU(L) + SUM
               ETACU(L) = ETACU(L) + EMINDU 
            ENDDO
         ENDIF
      ENDIF
 
C***  FREE-FREE  *******************************************************
C***  PRECALCULATE FREE-FREE GAUNT FACTORS FOR THE DIFFERENT ION CHARGES
         GFF(0)=.0
         DO 10 ION=1,MAXION
            CALL GAUNTFF (GIII,ION,XLAMKCU,TL)
            GFF(ION)=GIII*FLOAT(ION*ION)
 10      CONTINUE

C***  PRECALCULATE SIGMAFF, LEAVING OUT THE FACTOR NCHARGE*NCHARGE
         PRESIG =CFF/W3U/ROOTTL
         EXPFAC=EXP(-C1*WU/TL)
         DO 3 I=1,N
            SIGMAFF=PRESIG*GFF(NCHARG(I))
            SUM = ENEL * POPNUM(L,I) * SIGMAFF
            EMINDU=SUM*EXPFAC
            SUM=SUM-EMINDU
            IF (BPLOT2) THEN
               IF (I .EQ. N) THEN
                  IF (L .EQ. IPLOT) THEN
                     SUMFF = SUM
                     EMINDUFF = EMINDU
                     SIGMAFFFF = SIGMAFF
                     ENELFF = ENEL
                     POPNUMFF = POPNUM(L,I)
                     GFFFF = GFF(NCHARG(I))
                  ENDIF
               ENDIF
            ENDIF
            OPACU(L)=OPACU(L)+SUM
            ETACU(L)=ETACU(L)+EMINDU
 3       CONTINUE

      ENDIF

C***  If total true continuum opacity is negative:
C***  re-do the whole calculation, but skip lasering bound-free transitions
C***  Note: in the version with strict suppression (wrh  4-Apr-2019)
C***        this condition should never be met, and is therefore
C***        commented
cc      IF (OPACU(L) .LE. .0 .AND. NBFLASER .EQ. 0) THEN
cc         NBFLASER = 1
cc         GOTO 55
cc      ENDIF

C***  Interpolation of the BF+FF-Opacity
      IF (K .EQ. KCU) THEN
         OPAL = OPAL+ OPACU(L)
         ETAL = ETAL+ ETACU(L)
C***    New Start of Interpolation Procedure if (KCL = KCU)
         IF (KCL .EQ. KCU) THEN
            OPACL(L) = OPACU(L)
            ETACL(L) = ETACU(L)
         ENDIF
      ELSEIF (K .GE. KCL) THEN
         FL = FLOAT(KCU-K)/FLOAT(KCU-KCL)
         FU = 1.-FL         
         OPAL = OPAL + FL*OPACL(L) + FU*OPACU(L)
         ETAL = ETAL + FL*ETACL(L) + FU*ETACU(L)
      ELSE
         STOP 'CMFCOOP: FATAL ERROR'
      ENDIF

      ETAL=ETAL*C2*W3
 
C***  THOMSON SCATTERING ***********************************************
C***  Clump density and filling factor cancel out for Thomson scattering
      SUM=RNE(L)*SIGMAE
      if (l .eq. iplot) then
        sumt = sum
      endif
      OPAL=OPAL+SUM
C***  THOMSON = RELATIVE FRACTION FROM THE TOTAL OPACITY
      THOMSON(L)=SUM/OPAL
 
      OPA(L)=OPAL*ENTOT(L)*RSTAR
      ETA(L)=ETAL*ENTOT(L)*RSTAR
    1 CONTINUE
C***  ENDLOOP  ---------------------------------------------------------
 
      IF (BPLOT2) THEN
        WRITE (92,'(I8,32(E15.8,1X))') K, XLAM,
     >        SUMBF, EMINDUBF, SUMFF, EMINDUFF, SUMT, 
     >        SIGMAFFFF, ALOG10(ENELFF), POPNUMFF, GFFFF, PRESIG
      ENDIF

      RETURN
      END
