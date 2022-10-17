      SUBROUTINE ZONEINT (LRED, LBLUE, EMINT, CEMINT, TAUSUM, TAUSUMC,
     >                   PJPJ, ZFINE, KINDEX,
     >                   OPAFINE, OPAFC, OPALFIN, ETAFINE, ETAFC,
     >                   ETALFIN, SFINE, CSFINE, 
     >                   RRAY, ZRAY, XCMF, 
     >                   OPARAY, OPALRAY, ETARAY, ETACRAY, ETALRAY,
     >                   MAXXN, LTOT, DELXLAP, NBLINE, NLOBLFI, NLOBLLA, 
     >                   NDADDIM, IVERSION, LINPRO, AVOIGT,
     >                   TAU, TAUC, DTAU, DTAUC, WTAU, WTAUC,
     >                 XCMFFINE, POROLENGTHFINE, TAUMAX, XMAX, XMAXLIN,
     >                 DXMAX, BIRONLINES, OPAFE, ETAFE, NDDIM, NFLDIM, 
     >                 XCMFBLUE, XCMFRED, DXCMF, CORE, POROLENGTHRAY,
     >                 PHITAB, NFDIMPHITAB, NLDIMPHITAB, IPOINTERPHITAB, 
     >                 RADIUS, ND, MAXLAP, DD_VDOPDU_RAY, NATOM, 
     >                 IND_ELLINE,  
     >                 DD_VDOPDU_FINE_NORMFAC,
     >                 DD_VDOPDU_FINE_SQRD, DD_VDOPDU_FINE,
     >                 GRIEMPAR, KODAT, VDOP, ZINTER, NMOD, MAXMOD)

C***********************************************************************
C***  CALLED FROM: OBSFRAM
C***  FORMAL INTEGRATION ACROSS THE SCATTERING ZONE
C***  IVERSION = 0: INTEGRATION IN TAU - improved version 2-SEP-1998
C***  IVERSION = 1: INTEGRATION IN Z
C***********************************************************************

      IMPLICIT NONE

      REAL, DIMENSION(NDADDIM, NATOM), INTENT(IN) :: DD_VDOPDU_RAY
      REAL, DIMENSION(MAXXN, NATOM) :: DD_VDOPDU_FINE_SQRD, 
     >                                 DD_VDOPDU_FINE_NORMFAC, DD_VDOPDU_FINE
      INTEGER, DIMENSION (MAXLAP) :: IND_ELLINE
      REAL OPALRAY(NDADDIM,NBLINE), ETALRAY(NDADDIM,NBLINE)
      REAL OPALFIN(MAXXN,NBLINE), ETALFIN(MAXXN,NBLINE)
      REAL OPARAY(LTOT), ETARAY(LTOT), ETACRAY(LTOT)
      REAL RRAY(LTOT), ZRAY(LTOT), XCMF(LTOT), DELXLAP(NBLINE)  
      INTEGER IPOINTERPHITAB(NBLINE)
      REAL, DIMENSION(MAXLAP,NDDIM,MAXMOD) :: AVOIGT, GRIEMPAR
      REAL OPAFINE(MAXXN),ETAFINE(MAXXN),SFINE(MAXXN),ZFINE(MAXXN)
      REAL OPAFC(MAXXN), ETAFC(MAXXN), CSFINE(MAXXN)
      REAL TAU(MAXXN),DTAU(MAXXN),WTAU(MAXXN)
      REAL TAUC(MAXXN),DTAUC(MAXXN),WTAUC(MAXXN), XCMFFINE(MAXXN)
      REAL OPAFE(NDDIM,NFLDIM,MAXMOD), ETAFE(NDDIM,NFLDIM,MAXMOD)
      CHARACTER(8), DIMENSION(NBLINE) :: LINPRO 
      REAL, DIMENSION(NBLINE) :: XMAXLIN
      REAL POROLENGTHRAY(LTOT), POROLENGTHFINE(MAXXN)
      INTEGER, DIMENSION(NATOM) :: KODAT
      REAL ZINTER(2)      
      REAL, DIMENSION
     > (-NFDIMPHITAB:NFDIMPHITAB, NDDIM, NLDIMPHITAB,MAXMOD) :: PHITAB

      LOGICAL LASER, BIRONLINES, CORE, bHYROGEN, DEBUG

C***  next declarations needed only with IMPLICIT NONE
      REAL P, Q, WPIINV, XI, EPS, EPS2, DXFINE, DXMAX, DZFINE, RFINE
      REAL XCMFBLUE, XKCMF, XCMFRED, PJPJ, TEMP, DXCMF, WIF_P, WIF_M
      REAL ETAFEL, OPAFEL, OPAFEI, OPAFELM, ETAFELM,  ETAFEI
      REAL EMINT, CEMINT, RFINERIP, PHI, XI_VDOP, CORRFAC_MACROCLUMP
      REAL TAUSUM, TAUSUMC, WB, EXPTAUP, PR, QR, EXPTAU, EXPTAUM 
      REAL EXPTAUC, EXPTAUCP, EXPTAUCM, DZ, W, WAC, WBC
      REAL TAUMAX, VDOP, XMAX
      REAL RADIUS(NDDIM)
      REAL STARKHOLTSMARK, STARKPROF, SDOT, STARKVOIGT, WA

      INTEGER L, NFINE, LRIP, LMAX, LTOT, LRED, LBLUE, NPOINT, IPOINT 
      INTEGER I, NA, NATOM, KCMF_M, KCMF_P, NLOC, NLOBLFI, NLOBLLA
      INTEGER LL, IMOD, NMOD, IVERSION, LAST, KINDEX, ND, MAXXN, NBLINE
      INTEGER MAXMOD, NDADDIM, NDDIM, NFLDIM, NFDIMPHITAB, NLDIMPHITAB
      INTEGER MAXLAP

C***  WPIINV = 1. / SQRT(PI)
      DATA WPIINV / 0.564189583549 /
C***  The following two parameters are to make the new TAUVERSION 
C***     laser-resistant. 
      DATA EPS  / 1.E-10 /
      DATA EPS2 /  0.5   /
      SAVE EPS, EPS2

C***  DEBUG option; slows down the code!
      DEBUG = .FALSE. 

      LRIP = 1   ! initialization for depth interpolation (STARKPROF)
    
C***  LMAX is required because OPAFE, ETAFE was not mirrored at symmetry plane
      IF (CORE) THEN
         LMAX = LTOT
      ELSE
         LMAX = (LTOT+1) / 2
      ENDIF


C***  Speed up the formal: ZFINE only established at KINDEX=1
C***  All other interpolations are done (for each KINDEX) subsequently
      IF (KINDEX .EQ. 1) THEN

         DXFINE = DXMAX
      
C**      Loop on sub-intervals [L-1,L] of XCMF rough array.
C**      Refines every sub interval whose frequency-range is bigger than DXMAX 

C**      NFINE  - total number of inserted points in fine array
         NFINE = 0
      
      DO L=LRED+1,LBLUE
C***    NPOINT = number of fine points in current sub-interval [L-1,L] 
C***              including ZRAY(L-1), but excluding ZRAY(L)
        
C***    DXFINE defines the maximum Doppler shift that is allowed over an
C***    integration step along the ray. DXFINE must be small enough to resolve
C***    a Gaussian profile (0.3 Doppler units by default). When the *local*
C***    Doppler width is larger (e.g. by depth-dependent microturbulence or 
C***    Doppler width), this resolution can be relaxed according to the 
C***    finest resolution needed, as derived from the array DD_VDOPDU_RAY.
        DXFINE = DXMAX * MINVAL(DD_VDOPDU_RAY(L,1:NATOM))
        NPOINT = MAX(1, NINT(ABS(XCMF(L-1)-XCMF(L))/DXFINE))
    
        IF (NFINE+NPOINT .GE. MAXXN) THEN
            WRITE (0,'(A,I5,A,I5)') 
     >       'DIMENSION MAXXN=', MAXXN, ' INSUFFICIENT', MAXXN 
            STOP 'ERROR IN ZONEINT'
        ENDIF

        DZFINE = (ZRAY(L-1) - ZRAY(L)) / FLOAT(NPOINT)
        DO IPOINT=1, NPOINT
           ZFINE(NFINE+IPOINT) = ZRAY(L-1) - (IPOINT-1) * DZFINE
        ENDDO        
        NFINE = NFINE + NPOINT
      
      ENDDO
C***  End of Loop over the coarse intervals [L,L+1]
C***  Add endpont of last interval
      NFINE = NFINE+1
      ZFINE(NFINE) = ZRAY(LBLUE)

      ENDIF 
C**************************************************************

c      if (nfine .le. 1) write (0,*) '!!!! nfine .le 1 !!!!' 

c      DO i=1, nfine
c         write (0,*) 'i, zfine =', i, zfine(i)
c      enddo
c
c      stop 'test'

ccccccccccccccccccc begin steinbruch ccccccccccccccccccccccc      DXFINE = DXMAX
      
C** Interpolates the other arrays (ETA, OPA...) 

      L=2 ! meaning that ZRAY(L-1) > ZFINE(I) > ZRAY(L)
      DO I=1, NFINE
         IF (ZFINE(I) .LT. ZRAY(L)) L = L + 1
         CALL SPLINPO_FAST
     >        (XI, ZFINE(I), XCMF, ZRAY, LTOT, L, DEBUG)
            XCMFFINE(I) = XI
c            write (0,'(A,I3,F10.5,I3,F10.5,I5)') 'i, zfine, L, xcmf =', 
c     >           i, zfine(i), L, xi, KINDEX

C***     INTERPOLATION WEIGHT, LINEARLY IN RADIUS R
         RFINE=SQRT(PJPJ+ZFINE(I)*ZFINE(I))
         IF (RRAY(L) .NE. RRAY(L-1)) THEN
           P=(RFINE-RRAY(L-1))/(RRAY(L)-RRAY(L-1))
         ELSE
           P = 0.5
         ENDIF
         Q=1.-P

         POROLENGTHFINE(I) = 
     >      Q * POROLENGTHRAY(L-1) + P * POROLENGTHRAY(L)
          
C***     Interpolation of depth dependent VDOP arrays (linearly in radius)
         DO NA=1, NATOM
           DD_VDOPDU_FINE(I, NA) = 
     >         Q * DD_VDOPDU_RAY(L-1, NA) + P * DD_VDOPDU_RAY(L, NA)
           DD_VDOPDU_FINE_SQRD(I,NA) = 
     >         DD_VDOPDU_FINE(I,NA) * DD_VDOPDU_FINE(I,NA)
C***  This array contains the normalization factors for the gaussians
           DD_VDOPDU_FINE_NORMFAC(I,NA) = WPIINV / DD_VDOPDU_FINE(I,NA)
          ENDDO
         
C***      INTERPOLATION WEIGHTS MODIFIED: LINEAR INTERPOLATION IN R**2
          Q=Q*RRAY(L-1)*RRAY(L-1)/(RFINE*RFINE)
          P=P*RRAY(L  )*RRAY(L  )/(RFINE*RFINE)

          TEMP = Q * OPARAY(L-1) + P * OPARAY(L)
          OPAFINE(I) = TEMP
          OPAFC  (I) = TEMP

          ETAFINE(I) = Q * ETARAY (L-1) + P * ETARAY (L)
          ETAFC  (I) = Q * ETACRAY(L-1) + P * ETACRAY(L)

          
C***  Add Iron Opacities
          IF (BIRONLINES) THEN
C***        Find CMF-Frequency Index of the current fine frequency XI
C***        and interpolate OPAFE, ETAFE for that frequency
C***        at depth points (L-1) and (L)

C***        Frequency Indices and Interpolation Weights
            XKCMF  = (XCMFBLUE -  XI) / DXCMF + 1. 
            IF (XCMFBLUE .LT. XI) THEN
                WRITE (0,*) '*** INTERNAL ERROR DETECTED IN ZONEINT:'
                WRITE (0,*) '*** XCMFBLUE .LT. XI'
                WRITE (0,*) '*** XCMFBLUE=', XCMFBLUE
                WRITE (0,*) '*** XI      =', XI
                WRITE (0,*) '*** Frequency range in FORMCMF too small??'
                STOP        '*** FATAL ERROR ***'
            ENDIF  
            IF (XCMFRED .GT. XI) THEN            
                WRITE (0,*) '*** INTERNAL ERROR DETECTED IN ZONEINT:'
                WRITE (0,*) '*** XCMFRED .GT. XI'
                WRITE (0,*) '*** XCMFRED=', XCMFRED
                WRITE (0,*) '*** XI      =', XI
                WRITE (0,*) '*** Frequency range in FORMCMF too small??'
                STOP        '*** FATAL ERROR ***'
            ENDIF   
            KCMF_M = INT(XKCMF)
            KCMF_P = KCMF_M + 1
            WIF_P  = XKCMF - FLOAT(KCMF_M)
            WIF_M  = 1. - WIF_P

C***        Opacities at XKCMF at (L)
            IF (L .LE. LMAX) THEN
                LL = L
            ELSE
                LL = 2 * LMAX - L
            ENDIF

C***     Check if current point lies in second-model domain
            IMOD=1
            IF (NMOD .EQ. 2) THEN
               IF ((ZFINE(I)-ZINTER(1))*(ZFINE(I)-ZINTER(2)) .LT. .0)
     >            IMOD=2
            ENDIF

            OPAFEL  = WIF_M * OPAFE(LL,KCMF_M,IMOD) 
     >              + WIF_P * OPAFE(LL,KCMF_P,IMOD)
            ETAFEL  = WIF_M * ETAFE(LL,KCMF_M,IMOD) 
     >              + WIF_P * ETAFE(LL,KCMF_P,IMOD)
C***     Opacities at XKCMF at (L-1)
            IF (L-1 .LE. LMAX) THEN
                LL = L-1
            ELSE
                LL = 2 * LMAX - L-1
            ENDIF
            OPAFELM = WIF_M * OPAFE(LL,KCMF_M,IMOD) 
     >              + WIF_P * OPAFE(LL,KCMF_P,IMOD)
            ETAFELM = WIF_M * ETAFE(LL,KCMF_M,IMOD) 
     >              + WIF_P * ETAFE(LL,KCMF_P,IMOD)
            
C***        Final Radius-Interpolation
            OPAFEI = Q * OPAFELM + P * OPAFEL
            ETAFEI = Q * ETAFELM + P * ETAFEL
C***        Adding up Iron to the total Opacity/Emissivity
            OPAFINE(I) = OPAFINE(I) + OPAFEI
            ETAFINE(I) = ETAFINE(I) + ETAFEI
         ENDIF
C***     -- end of iron branch ---

         DO NLOC=NLOBLFI,NLOBLLA
            OPALFIN(I,NLOC) = Q * OPALRAY(L-1,NLOC) +
     >                          P * OPALRAY(L,NLOC)
            ETALFIN(I,NLOC) = Q * ETALRAY(L-1,NLOC) +
     >                          P * ETALRAY(L,NLOC)
         ENDDO
      ENDDO
C***  End of Loop over all fine depth points

C***  Macroclumping correction for the continuum opacity
      DO I = 1, NFINE
       CALL MACROCLUMP (OPAFC(I), POROLENGTHFINE(I), CORRFAC_MACROCLUMP)
       OPAFC(I) = OPAFC(I) * CORRFAC_MACROCLUMP
       ETAFC(I) = ETAFC(I) * CORRFAC_MACROCLUMP
      ENDDO
      
C******************************************************************
C***  Z version
C***  NOTE: This version has no problems with LASER effects.
C***        The weight function exp(-tau) is NOT incorporated in the 
C***        quadrature weights. This implies a loss of accuracy.
C******************************************************************
      IF (IVERSION .EQ. 1) THEN
C***  Opacities and optical depths ...
      DO I=1, NFINE
C***     This setting of RFINE assures that L-interpolation weigts
C***     have not been calculated yet (for Subr. STARKPROF)
         RFINERIP = .0

C***     Check if current point lies in second-model domain
         IMOD=1
         IF (NMOD .EQ. 2) THEN
            IF ((ZFINE(I)-ZINTER(1))*(ZFINE(I)-ZINTER(2)) .LT. .0) 
     >         IMOD=2
         ENDIF

C***     ADD LINE PLUS CONTINUUM OPACITIES
         DO NLOC=NLOBLFI,NLOBLLA
            NA = IND_ELLINE(NLOC) !NA ist the element index corresponding to the line
            XI = XCMFFINE(I) - DELXLAP(NLOC)
C***        Max. bandwidth holds for all profile types!
            IF (ABS(XI) .GT. XMAXLIN(NLOC)) CYCLE
            IF (LINPRO(NLOC) .EQ. '        ') THEN
               PHI = EXP(-XI*XI / DD_VDOPDU_FINE_SQRD(I,NA))
               PHI = PHI * DD_VDOPDU_FINE_NORMFAC(I,NA)
               OPAFINE(I)=OPAFINE(I)+PHI*OPALFIN(I,NLOC)
               ETAFINE(I)=ETAFINE(I)+PHI*ETALFIN(I,NLOC)

            ELSE IF (LINPRO(NLOC) .EQ. 'VOIGT   '
     >          .OR. LINPRO(NLOC) .EQ. 'BRD-HeI '
     >          .OR. LINPRO(NLOC) .EQ. 'Q-STARK ') THEN
C***           XI Argument of VOIGTH needs to account for VDOP
               XI_VDOP = XI/DD_VDOPDU_FINE(I,NA)
C***           STARKVOIGT interpolates appropriate value for AV 
C***             as function of radius and for primary and second model, 
C***             then calls the Voigt function VOIGTH 
               PHI = STARKVOIGT (XI_VDOP, AVOIGT(1,1,IMOD), NLOC, MAXLAP, 
     >               ZFINE(I), LRIP, RFINERIP, PR, QR, RADIUS, ND, 
     >               PJPJ)
C***           Profile is normalized appropriately
               PHI = PHI / DD_VDOPDU_FINE(I,NA)
               OPAFINE(I) = OPAFINE(I) + PHI * OPALFIN(I,NLOC)
               ETAFINE(I) = ETAFINE(I) + PHI * ETALFIN(I,NLOC)
            ELSE IF (LINPRO(NLOC) == 'L-STARK ') THEN
C***           Stark broadening: depth-dependent HOLTSMARK PROFILE
               bHYROGEN = (NA == KODAT(1))
C***           STARKHOLTSMARK interpolates appropriate value for GRIEMPAR
               PHI = STARKHOLTSMARK (XI, VDOP, DD_VDOPDU_FINE(I,NA),
     >                               GRIEMPAR(1,1,IMOD), NLOC, MAXLAP, 
     >                               ZFINE(I), LRIP, RFINERIP, PR, QR,
     >                               RADIUS, ND, PJPJ, bHYROGEN)
               OPAFINE(I) = OPAFINE(I) + PHI * OPALFIN(I,NLOC)
               ETAFINE(I) = ETAFINE(I) + PHI * ETALFIN(I,NLOC)
            ELSE IF (LINPRO(NLOC)(:3) .EQ. 'BRD') THEN
C***              PRESSURE BROADENING: Hydrogen and Helium / tabulated
C*** If depth-dependent VDOP active, data is tabulated appropriately in STARKBROAD 
               PHI = STARKPROF (XI, IPOINTERPHITAB(NLOC), ZFINE(I),
     >                  LRIP, RFINERIP, PR, QR,
     >                  PHITAB(-NFDIMPHITAB,1,1,IMOD),
     >                  NFDIMPHITAB, NLDIMPHITAB, RADIUS, ND, NDDIM,
     >                  PJPJ, DXMAX)
               OPAFINE(I) = OPAFINE(I) + PHI*OPALFIN(I,NLOC)
               ETAFINE(I) = ETAFINE(I) + PHI*ETALFIN(I,NLOC)
            ELSE
               WRITE (0, '(2A)') 
     >          '*** UNDEFINED LINE PROFILE TYPE: ', LINPRO(NLOC)
               STOP '*** FATAL ERROR IN ZONEINT'
            ENDIF
         ENDDO
 
C***  Macroclumping correction for the total opacity
      CALL MACROCLUMP (OPAFINE(I), POROLENGTHFINE(I), 
     >                  CORRFAC_MACROCLUMP)
      OPAFINE(I) = OPAFINE(I) * CORRFAC_MACROCLUMP
      ETAFINE(I) = ETAFINE(I) * CORRFAC_MACROCLUMP

C***  ESTABLISH DTAU(I) = OPTICAL DEPTH INCREMENT BETWEEN I-1 AND I
C***         AND TAU(I) = OPTICAL DEPTH SCALE
         IF (I .EQ. 1) THEN
            DTAU (1) = 0.
            DTAUC(1) = 0.
            TAU  (1) = TAUSUM
            TAUC (1) = TAUSUMC
         ELSE
            DTAU (I) = 0.5 * 
     >              (OPAFINE(I-1)+OPAFINE(I)) * (ZFINE(I-1)-ZFINE(I))  
            DTAUC(I) = 0.5 * 
     >              (OPAFC  (I-1)+OPAFC  (I)) * (ZFINE(I-1)-ZFINE(I))

            TAU (I) = TAU (I-1) + DTAU (I)
            TAUC(I) = TAUC(I-1) + DTAUC(I)
         ENDIF
      ENDDO   !end of fine frequency loop

      TAUSUM  = TAU (NFINE)
      TAUSUMC = TAUC(NFINE)


C***      S means: EMISSIVITY * EXP(-TAU)
          DO 14 I=1, NFINE
            IF (TAU(I) .LT. 700.) THEN
               SFINE(I) = ETAFINE(I) * EXP(-TAU (I))
            ELSE
              SFINE(I) = 0.
            ENDIF
            IF (TAUC(I) .LT. 700.) THEN
              CSFINE(I) = ETAFC  (I) * EXP(-TAUC(I))
            ELSE
              CSFINE(I) = 0.
            ENDIF
   14     CONTINUE

          WTAU (1) = 0.5 * (ZFINE(1) - ZFINE(2))
          WTAUC(1) = WTAU(1)
          DO 7 I = 2, NFINE-1
            WTAU (I) = 0.5 * (ZFINE(I-1) - ZFINE(I+1))
            WTAUC(I) = WTAU(I)
    7     CONTINUE
          WTAU (NFINE) = 0.5 * (ZFINE(NFINE-1) - ZFINE(NFINE))
          WTAUC(NFINE) = WTAU(NFINE)

C***  INTEGRATION SUM, USING CRAY VECTOR FUNCTION SDOT (SCALAR PRODUCT)
      CEMINT = CEMINT + SDOT(NFINE,CSFINE,1,WTAUC,1)
       EMINT =  EMINT + SDOT(NFINE, SFINE,1,WTAU ,1)



C******************************************************************
C***  TAU version (NEW! 2-SEP-1998, wrh)
C***        This version combines the advantages of both the (old) TAU 
C***        and the Z version, as it is accurate AND laser-resistant. 
C***        The weight function exp(-tau) is incorporated in the 
C***        quadrature weights for better accuracy (equivalent to the 
C***        TAU version) 
C***        In LASER intervals, the normal Z integration is restored. 
C***        This is still possible, because the divison by opa (S=ETA/OPA)
C***        is not performed in advance, but ETA is kept and 1/OPA is 
C***        implied in the weights. For this switching of the integretion 
C***        methods, it is necessary to apply the endpoint-terms exp(-tau)
C***        in the integration weights each time, although they cancel out 
C***        in case of subsequent tau-version intervals
C******************************************************************
      ELSE IF (IVERSION .EQ. 0) THEN
         EPS = 1.E-10

C***     LINE + CONTINUUM ***************************************
 
C***  ESTABLISH DTAU(I) = OPTICAL DEPTH INCREMENT BETWEEN I-1 AND I
C***         AND TAU(I) = OPTICAL DEPTH SCALE
C***     AND CALCULATE INTEGRATION WEIGHTS EXP(-TAU) DTAU
C***     THE LOOP RUNS FOR I=0 FIRST IN ORDER TO PREPARE LINE OPACITIES AT I=1

         TAU(1) = TAUSUM
         EXPTAU = EXP(-TAU(1))
         DO 5 I = 0, NFINE-1
C***        This setting of RFINE assures that L-interpolation weigts 
C***        have not been calculated yet (for Subr. STARKPROF)
            RFINERIP = .0

C***        Check if current point lies in second-model domain
            IMOD=1
            IF (NMOD .EQ. 2) THEN
               IF ((ZFINE(I+1)-ZINTER(1))*(ZFINE(I+1)-ZINTER(2)) .LT. .0) 
     >            IMOD=2
            ENDIF

C***        ADD LINE PLUS CONTINUUM OPACITIES
            DO NLOC=NLOBLFI,NLOBLLA
               NA = IND_ELLINE(NLOC) !NA ist the element index corresponding to the line
               XI = XCMFFINE(I+1) - DELXLAP(NLOC)
C***           Max. bandwidth holds for all profile types!
cc               IF (ABS(XI) .GT. XMAX) CYCLE
               IF (ABS(XI) .GT. XMAXLIN(NLOC)) CYCLE
               IF (LINPRO(NLOC) .EQ. '        ') THEN
C***              DOPPLER BROADENING ONLY: GAUSS PROFILE
                  PHI = EXP(-XI*XI / DD_VDOPDU_FINE_SQRD(I+1,NA))
                  PHI = PHI * DD_VDOPDU_FINE_NORMFAC(I+1,NA)
                  OPAFINE(I+1) = OPAFINE(I+1) 
     >                         + PHI * OPALFIN(I+1,NLOC)
                  ETAFINE(I+1) = ETAFINE(I+1) 
     >                         + PHI * ETALFIN(I+1,NLOC)

               ELSE IF (LINPRO(NLOC) .EQ. 'VOIGT   '
     >             .OR. LINPRO(NLOC) .EQ. 'BRD-HeI '
     >             .OR. LINPRO(NLOC) .EQ. 'Q-STARK ') THEN
                  XI_VDOP = XI/DD_VDOPDU_FINE(I+1,NA)
                  PHI = STARKVOIGT (XI_VDOP, AVOIGT(1,1,IMOD), NLOC, MAXLAP, 
     >                  ZFINE(I+1), LRIP, RFINERIP, PR, QR, RADIUS, ND, 
     >                  PJPJ)
C***              Profile is normalized appropriately
                  PHI = PHI / DD_VDOPDU_FINE(I+1,NA)
                  OPAFINE(I+1) = OPAFINE(I+1) + PHI * OPALFIN(I+1,NLOC)
                  ETAFINE(I+1) = ETAFINE(I+1) + PHI * ETALFIN(I+1,NLOC)

               ELSE IF (LINPRO(NLOC) == 'L-STARK ') THEN
C***              Stark broadening: depth-dependent HOLTSMARK PROFILE
                  bHYROGEN = (NA == KODAT(1))
C***              STARKHOLTSMARK interpolates appropriate value for GRIEMPAR
                  PHI = STARKHOLTSMARK (XI, VDOP, DD_VDOPDU_FINE(I+1,NA),
     >                              GRIEMPAR(1,1,IMOD), NLOC, MAXLAP, 
     >                              ZFINE(I+1), LRIP, RFINERIP, PR, QR,
     >                              RADIUS, ND, PJPJ, bHYROGEN)
                  OPAFINE(I+1) = OPAFINE(I+1) + PHI * OPALFIN(I+1,NLOC)
                  ETAFINE(I+1) = ETAFINE(I+1) + PHI * ETALFIN(I+1,NLOC)
               ELSE IF (LINPRO(NLOC)(:3) .EQ. 'BRD') THEN
C***              PRESSURE BROADENING: H I and He II / tabulated
C*** If depth-dependent VDOP active, data is tabulated appropriately in STARKBROAD 
                  PHI = STARKPROF (XI, IPOINTERPHITAB(NLOC), ZFINE(I+1), 
     >                  LRIP, RFINERIP, PR, QR, 
     >                  PHITAB(-NFDIMPHITAB,1,1,IMOD),
     >                  NFDIMPHITAB, NLDIMPHITAB, RADIUS, ND, NDDIM,
     >                  PJPJ, DXMAX)
                  OPAFINE(I+1) = OPAFINE(I+1) + PHI*OPALFIN(I+1,NLOC)
                  ETAFINE(I+1) = ETAFINE(I+1) + PHI*ETALFIN(I+1,NLOC)
               ELSE
                  WRITE (0, '(2A)') 
     >             '*** UNDEFINED LINE PROFILE TYPE: ', LINPRO(NLOC)
                  STOP '*** FATAL ERROR IN ZONEINT'
               ENDIF
            ENDDO

C***  Macroclumping correction for the total opacity
           CALL MACROCLUMP (OPAFINE(I+1), POROLENGTHFINE(I+1), 
     >                  CORRFAC_MACROCLUMP)
           OPAFINE(I+1) = OPAFINE(I+1) * CORRFAC_MACROCLUMP
           ETAFINE(I+1) = ETAFINE(I+1) * CORRFAC_MACROCLUMP
         
            IF (I .EQ. 0) GOTO 5

            DTAU (I+1) = 0.5 * 
     >              (OPAFINE(I)+OPAFINE(I+1)) * (ZFINE(I)-ZFINE(I+1))
            TAU (I+1) = TAU (I) + DTAU (I+1)


C***        WB COMES FROM THE INTERVAL I-1, I (ZERO for first interval!)
            IF (I .EQ. 1) THEN
               WB = .0
               EXPTAUP  = EXP(-TAU (2))
            ELSE
               EXPTAUM  = EXPTAU
               EXPTAU   = EXPTAUP
               EXPTAUP  = EXP(-TAU(I+1))
C***           LASER condition for that interval (LINE)
               LASER = ABS(DTAU (I)) .LT. EPS  .OR. 
     >              OPAFINE(I-1) .LT. EPS2 * OPAFC(I-1)    .OR.
     >              OPAFINE(I  ) .LT. EPS2 * OPAFC(I)
               IF (LASER) THEN 
                  DZ = ZFINE(I-1) - ZFINE(I)
                  WB =  EXPTAU  * DZ * 0.5
               ELSE            
                  WB = (-EXPTAU - (EXPTAU  - EXPTAUM ) / DTAU (I)) 
     >                 / OPAFINE(I)
               ENDIF
            ENDIF

C***        WA COMES FROM THE INTERVAL I, I+1
C***        LASER condition for that interval (LINE)
            LASER = ABS(DTAU (I+1)) .LT. EPS  .OR. 
     >           OPAFINE(I+1) .LE. EPS2 * OPAFC(I+1)    .OR.
     >           OPAFINE(I  ) .LE. EPS2 * OPAFC(I  )
            IF (LASER) THEN
               DZ = ZFINE(I) - ZFINE(I+1)
               WA  = EXPTAU  * DZ * 0.5
            ELSE
               WA  = (EXPTAU + (EXPTAUP - EXPTAU ) / DTAU (I+1)) 
     >               /  OPAFINE(I)
            ENDIF
          
            W = WA  + WB
            EMINT = EMINT + ETAFINE(I) * W 

            IF (TAU(I+1) .GT. TAUMAX) THEN
               LAST = I + 1
               GOTO 51
            ENDIF 

    5    CONTINUE

C***     Final Interval N-1, N (or LAST-1, LAST, if TAUMAX was reached)

         LAST = NFINE
   51    CONTINUE


C***     LASER condition for that interval (LINE)
         LASER = ABS(DTAU (LAST)) .LT. EPS  .OR. 
     >           OPAFINE(LAST-1) .LT. EPS2 * OPAFC(LAST-1)  .OR.
     >           OPAFINE(LAST  ) .LT. EPS2 * OPAFC(LAST  )
         IF (LASER) THEN 
            DZ = ZFINE(LAST-1) - ZFINE(LAST)
            W = EXPTAUP * DZ * 0.5
         ELSE
            W = (-EXPTAUP - (EXPTAUP-EXPTAU)/DTAU(LAST) ) 
     >          / OPAFINE(LAST) 
         ENDIF

         EMINT = EMINT + ETAFINE(LAST) * W 
         TAUSUM  = TAU(LAST)


C***     ONLY CONTINUUM  *******************************************

C***  ESTABLISH DTAUC(I) = OPTICAL DEPTH INCREMENT BETWEEN I-1 AND I
C***         AND TAUC(I) = OPTICAL DEPTH SCALE
C***         Note: This calculation must hurry one step ahead
C***     AND CALCULATE INTEGRATION WEIGHTS EXP(-TAU) DTAU


         TAUC (1) = TAUSUMC
         EXPTAUC  = EXP(-TAUC(1))

         DO 55 I = 1, NFINE-1

            DTAUC(I+1) = 0.5 * 
     >              (OPAFC(I)+OPAFC(I+1)) * (ZFINE(I)-ZFINE(I+1))
            TAUC(I+1) = TAUC(I) + DTAUC(I+1)

C***        WB COMES FROM THE INTERVAL I-1, I (ZERO for the first interval!)
            IF (I .EQ. 1) THEN
               WBC = .0
               EXPTAUCP = EXP(-TAUC(2))
            ELSE
               EXPTAUCM = EXPTAUC
               EXPTAUC  = EXPTAUCP
               EXPTAUCP = EXP(-TAUC(I+1))
C***           LASER condition for that interval (Continuum)
               LASER = ABS(DTAUC (I)) .LT. EPS  .OR. 
     >              OPAFC(I-1) .LT. .0  .OR.   
     >              OPAFC(I  ) .LT. .0    
               IF (LASER) THEN 
                  DZ = ZFINE(I-1) - ZFINE(I)
                  WBC = EXPTAUC * DZ * 0.5
               ELSE            
                  WBC = (-EXPTAUC - (EXPTAUC - EXPTAUCM) / DTAUC(I)) 
     >                  / OPAFC  (I)
               ENDIF
            ENDIF

C***        WA COMES FROM THE INTERVAL I, I+1
C***        LASER condition for that interval (Continuum)
            LASER = ABS(DTAUC (I+1)) .LT. EPS  .OR. 
     >           OPAFC(I+1) .LT. 0.  .OR.   
     >           OPAFC(I  ) .LT. 0.    
            IF (LASER) THEN
               DZ = ZFINE(I) - ZFINE(I+1)
               WAC = EXPTAUC * DZ * 0.5
            ELSE
               WAC = (EXPTAUC + (EXPTAUCP - EXPTAUC) / DTAUC(I+1)) 
     >               / OPAFC  (I)
            ENDIF
          
            W = WAC + WBC
            CEMINT = CEMINT + ETAFC(I) * W 

         IF (TAUC(I+1) .GT. TAUMAX) THEN
            LAST = I + 1
            GOTO 52
         ENDIF 

   55    CONTINUE

C***     Final Interval N-1, N (or LAST, if TAUMAX is reached)

         LAST = NFINE

   52    CONTINUE

C***     LASER condition for that interval (Continuum)
         LASER = ABS(DTAUC (LAST)) .LT. EPS  .OR. 
     >           OPAFC(LAST-1) .LE. 0.  .OR.   
     >           OPAFC(LAST  ) .LE. 0.    

         IF (LASER) THEN 
            DZ = ZFINE(LAST-1) - ZFINE(LAST)
            W = EXPTAUCP * DZ * 0.5
         ELSE
            W = (-EXPTAUCP - (EXPTAUCP - EXPTAUC) / DTAUC(LAST)) 
     >          / OPAFC(LAST) 
         ENDIF
         CEMINT = CEMINT + ETAFC(LAST) * W 
         TAUSUMC = TAUC(LAST)

C***********************************************************************
      ELSE
          STOP 'IVERSION INVALID IN ZONEINT'
      ENDIF

      RETURN
      END
