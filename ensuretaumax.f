      SUBROUTINE ENSURETAUMAX(bENSURETAUXMAX, !determines whether taumax iteration is allowed or not
     >                        HYSTACC,
     >                        IHSSTATUS,
     >    TAUMAX,       !TAU should have this value at inner boundary
     >    TAUACC,       !Accuracy level for TAUMAX (read from CARDS file)
     >    bTauStrict,   !if .FALSE. TAU scale is only changed if TAU(ND) is too low
     >    VMIN,         !Minimum velocity (can be changed in this routine)
     >    VELO,         !Velocity field vector
     >    GRADI,        !Velocity gradient vector
     >    RADIUS,       !radius values (in Rstar) for depth points
     >    XMDOT,        !log Mdot in M_sun/yr
     >    GEDDRAD,
     >    RHO,          !mass density per depth point
     >    TAUROSScont,  !Rosseland optical depth vector (continuum only)
     >    RMAX,         !outer boundary of calculations (VFINAL is reached here)
     >    RCONorg,      !old connection point between hydrostatic and beta law domain (updated here)
     >    T,            !temperature at each depth point  (TNEW in STEAL main prog.)
     >    TEFF,         !effective temperature Tstar at Rstar
     >    TOLD,         !old temperature one iteration ago  (for interpolation)
     >    TOLD2,        !old temperature two iterations ago (for interpolation)
     >    TOLD3,        !old temperature three iterations ago (for interpolation)
     >    CKONVER,      !char vector indicating if depth point is converged
     >    bTauMaxSafe,  !if true taumax is only adjusted if all depth points are conv.
     >    LASTTAU,      !Number of STEAL jobs since last TAUMAX iteration
     >    GLOG,         !log g
     >    GEFFLOG,      !log geff (only > 0 if geff or Eddington Gamma predefined in the CARDS file)
     >    bFixGEFF,     !if true, GEFF was fixed by the CARDS input
     >    XMSTAR,       !stellar mass in solar masses
     >    RSTAR,        !stellar radius at TAUMAX (in cm)
     >    XMU,          !relative mass (in AMU?) (called XMASS in INITEL)
     >    ENTOT,        !Total number of electrons per depth point
     >    RNE,          !relative electron density
     >    ND,           !Number of depth points
     >    NDDIM,        !Maximum number of depth points
     >    NP,           !Number of impact parameters
     >    NPDIM,        !Maximum number of impact parameters
     >    P,            !Impact parameter grid
     >    Z,            !z values (via z=sqrt(r^2-p^2)
     >    VTURB,        !turbulence velocity in km/s
     >    ARAD,         !radiative acceleration (from COLI)
     >    APRESS,       !pressure acceleration (from STEAL->INITFCORR)
     >    AMECH,        !mechanical acceleration (inerta, from STEAL->INITFCORR)
     >    AGRAV,        !gravitational acceleration (from STEAL->INITFCORR)
     >    XJC,
     >    bThinWind,    !Thin Wind CARDS option set? (TRUE/FALSE)
     >    ThinCard,     !Thin-Karte aus dem CARDS file
     >    RadiusGridParameters,     !Lines from CARDS for radius grid settings
     >    DENSCON,      !clumping factor per depth point (from CLUMP_STRUCT)
     >    FILLFAC,      !filling factor per depth point (= 1./DENSCON)
     >    DENSCON_FIX,  !fixed clumping factor value from CARDS
     >    DENSCON_LINE, !DENSCON-Line from CARDS
     >    NFIRST,
     >    NLAST,
     >    NATOM,        !Number of different atoms in DATOM
     >    ABXYZ,        !Abundances by number
     >    ATMASS,       !
     >    GEDD,         !Eddington Gamma (only read if fixed)
     >    bUNLU,        !True if temperature corrections are allowed
     >    HTOTL,
     >    HTOTCMF0,
     >    FQLIMIT,
     >    bFULLHYDROSTAT,
     >    bGAMMARADMEAN,
     >    TAUROSS,
     >    CORMAX,       !max corrections in this steal job
     >    TauCorLimits, !correction limits from CARDS for tau adjustment
     >    INCRIT,       !specifies RGRID criterion for each depth point
     >    VCRIT,        !specifies velocity criterion for each depth point
     >    SCRATCH,      !Scratcharray for popnum interpolation stuff
     >    ReduceTauCorrections,    !CARDS option for tau reducing has been set
     >    GEDDreduce,   !reduce GEDD change by mixing with old value
     >    QIONMEAN,
     >    GAMMARADMEAN,
     >    GRSTATIC,     !Gammarad(L) values from last TAUMAX iteration
     >    bGEddFix,     !True if GEdd has been fixed via CARDS (see DECSTAR)
     >    bTauUpdated,  !returns if tau iteration was done or not
     >    JOBNUM,
     >    POPMIN,
     >    XJCorg,
     >    NCOLIP,
     >    NEWWRC,
     >    NTOLD,
     >    DEPARTNDorg,
C     
C         !Parameters after here are only needed to call TAUSCAL
     >    NDIM,MAXKONT,POPNUM,POP1,POP2,    !5     
     >    POP3,N,EN,LEVEL,NCHARG,WEIGHT,    !11
     >    ELEVEL,EION,EINST,ALPHA,SEXPO,    !16
     >    ADDCON1, ADDCON2, ADDCON3,        !19
     >    IGAUNT,NOM,NF,                    !22
     >    XLAMBDA,FWEIGHT, TAUTHOM,         !25
     >    MAXATOM,SIGMATHK,SEXPOK,EDGEK,    !29
     >    KODAT,KONTNUP,KONTLOW,LASTKON     !33
     >    )
C***********************************************************************
C***  THIS SUBROUTINE adjusts the velocity field such that the innermost
C***   (=maximum) value of the Rosseland continuum optical depth has the
C***   same way as set in the CARDS file. The exact way of how the 
C***   velocity is changed depends on several options.
C    
C     TAUMAX-Line: FIX or MIN are needed to enable the routine at all
C     HYDROSTATIC INTEGRATION - If this line is set, the inner part is
C                               calculated such that the hydrostatic
C                               equation is fulfilled 
C                               (g_eff is calculated via a_thom only)
C     HYDROSTATIC INTEGRATION FULL - Same as above but with a_rad 
C                                    instead of a_thom
C    
C     There are several minor options which can change significantly how
C     much this routine changes the current stratification. A bad set of
C     parameters can even ruin a model which would otherwise run fine. The
C     details of all CARDS options are given in the PoWR manual.
C
C     This subroutine should only be called if TAUMAX has been specified
C     in the CARDS file and the HYDRO option (HD solution) is not used
C***********************************************************************

      IMPLICIT NONE
      INCLUDE 'interfacebib.inc'

      INTEGER, INTENT(IN) :: ND, NDDIM, NP, NPDIM, NATOM, JOBNUM,
     >                       N, NDIM, NF, MAXKONT, MAXATOM, LASTKON,
     >                       NCOLIP, NEWWRC, NTOLD
      REAL, INTENT(IN) :: TAUMAX, TAUACC, TEFF, RSTAR, VTURB, 
     >                    CORMAX, ReduceTauCorrections, GEDDreduce,
     >                    XMDOT, POPMIN
      REAL, DIMENSION(2), INTENT(IN) :: TauCorLimits
      REAL, INTENT(INOUT) :: GEDD, GLOG, GEFFLOG, RMAX, RCONorg, 
     >                       VMIN, XMSTAR, GEDDRAD, HYSTACC
      REAL, DIMENSION(ND), INTENT(IN) :: TAUROSS             !full tau from COLI (not changed here)
      REAL, DIMENSION(ND), INTENT(INOUT) :: RADIUS, VELO, GRADI, ENTOT,
     >                                      T, TAUTHOM, TAUROSScont,
     >                                      RHO, XMU, RNE, GRSTATIC,
     >                                      TOLD, TOLD2, TOLD3
      REAL, DIMENSION(ND-1), INTENT(IN) :: ARAD, APRESS, AMECH, AGRAV, 
     >                                     HTOTL, HTOTCMF0
      REAL, DIMENSION(ND, NF), INTENT(INOUT) :: XJC
      INTEGER, INTENT(INOUT) :: LASTTAU, IHSSTATUS
      CHARACTER(1), DIMENSION(ND) :: CKONVER

C***  Info on depth point criteria from GEOMESH>RGRID and 
C***   the velocity calculation criterion
      CHARACTER(8), DIMENSION(ND), INTENT(INOUT) :: INCRIT, VCRIT     

      LOGICAL, INTENT(IN) :: bENSURETAUXMAX, bThinWind, bTauStrict, 
     >                       bGEddFix, bUNLU, bTauMaxSafe
      LOGICAL, INTENT(OUT) :: bTauUpdated

      INTEGER :: NDVALUE, NC
      REAL, DIMENSION(1) :: OLDRADI, OLDVELO

      !Variables from Parameters that are only needed for TAUSCAL (and POPNUM interpolation)
      INTEGER, DIMENSION(N), INTENT(IN) :: NCHARG, NOM
      INTEGER, DIMENSION(NATOM), INTENT(IN) :: NFIRST, NLAST
      INTEGER, DIMENSION(MAXKONT), INTENT(IN) :: 
     >                IGAUNT, KONTNUP, KONTLOW
      INTEGER, DIMENSION(MAXATOM), INTENT(IN) :: KODAT
      REAL, DIMENSION(N), INTENT(IN) :: WEIGHT, ELEVEL, EION
      REAL, DIMENSION(N), INTENT(INOUT) :: EN, DEPARTNDorg
      REAL, DIMENSION(N, N), INTENT(IN) :: EINST
      REAL, DIMENSION(ND, N), INTENT(INOUT) :: POPNUM, POP1, POP2, POP3,
     >                                         SCRATCH
      REAL, DIMENSION(MAXKONT), INTENT(IN) :: 
     >                ALPHA, SEXPO, ADDCON1, ADDCON2, ADDCON3
      REAL, DIMENSION(NATOM), INTENT(IN) :: ABXYZ, ATMASS
      REAL, DIMENSION(MAXATOM, MAXATOM), INTENT(IN) :: 
     >                SIGMATHK, SEXPOK, EDGEK
      REAL, DIMENSION(NF), INTENT(IN) :: FWEIGHT, XLAMBDA
      CHARACTER(10), DIMENSION(N), INTENT(IN) :: LEVEL

      REAL, DIMENSION(ND, NF), INTENT(INOUT) :: XJCorg

      !Internal variables
      INTEGER :: I, J, K, L, ITTAU, Nct, LastDAMP, IT_failsafe, LL,
     >           NA, IPAR, NPAR, LCON, LfM, LTMCM, NHYSTBAD, NTAUIV,
     >           LCHECK, IERR
      INTEGER, PARAMETER :: MAXMASSIT = 100   !Maximum number of mass iterations
      INTEGER, PARAMETER :: MAXITTAU  = 30    !Maximum number of iterations
      REAL, PARAMETER :: GEDDreduceMIN = 0.001  !Minimum GEDDreduce value
C***  The full Eddington factor may not be higher than GAMMARADmax 
C***    to have some safety distance to the Eddington limit, 
C***    which helps keeping the density stratification stable. 
C***      A maximum value of 0.9 is also used in TLUSTY 
      REAL, PARAMETER :: GAMMARADmax = 0.9
      REAL, DIMENSION(MAXITTAU) :: VMINHIST   !history of vmin steps
      REAL, DIMENSION(MAXITTAU) :: TAUMAXHIST !history of TAUROSScont(ND)

C***  Attention! The following block of vectors is loccaly dimensioned
C***             ("HEAP variables" since dimension ND is formal parameter!)
      REAL, DIMENSION(ND) :: RADIUSorg, VELOorg, TAUSCALorg, RI,
     >                       OPAROSSorg, Torg, T1old, T2old, T3old,
     >                       GEDDTHOM, GEFFL, GAMMARAD, RHELP,
     >                       ENTOTorg, DR, ARADHELP, RNEorg
      REAL, DIMENSION(ND-1) :: ARADmod
          
      REAL :: ATMEAN, FM, STEPDAMP, RSTARSU, RNEL, XMSTARG,
     >        VFINAL, VMINCAND, BETA, BETA2, BETA2FRACTION, TAUCURMAX,
     >        VPAR1, VPAR2, HSCALE, VPAR1_2, VPAR2_2, FQLIMIT,
     >        RCON, VMINOLD, VMINOLD2, CORLOG, VCON, ARADL, FQ,
     >        TAUCHECK, TROLD, TROLD2, TAU23, R23, GEDDL,
     >        VMIN_over, VMIN_under, VMIN_min, VMIN_max, DARADDR,
     >        DAGRAVDR, RL, XLSTAR, XLSTARS, Qhopf, Qorg, 
     >        P1, P2, RHOINT, APR, GFR, VPAR1cand, VPAR2cand,
     >        GAMMARADMEAN, VL, TL, TscaleFac, RINT, TAUINT,
     >        facMass, DIFF1, DIFF2, ENTOTL, AGRAVL, MECHRATIOLOG,
     >        HSCALEold, TauGammaEffect, QIONMEAN, fGAMMACOR, XMG,
     >        GAMMARADMEANcalc, GEDDreducenow, RCRIT, ENEL,
     >        GLOGorg, XMSTARorg, GEFFLOGorg, HYSTRATIO, GEDDreduceL,
     >        DIFFTAUMAX_J, DIFFTAUMAX_MIN, 
     >        QVorg, QVnew, QRHOnew, QRHOorg, Textrap, TAUTHICK
      REAL, DIMENSION(3) :: TMCORRMAX

C***  Vector, true if HYSTACC fulfilled per depth point 
      LOGICAL, DIMENSION(ND) :: HYSTGOOD
      
      CHARACTER(8), DIMENSION(ND) :: VELOCRITERION
      CHARACTER(40) :: ACTPAR
      CHARACTER(4) :: TMENFCRIT, CTFMT
      CHARACTER(8) :: TAUFLAB

      LOGICAL bNewVelo, bHydroStat, bNoRGrid, bNoDetails, 
     >        bTauInterpolation, bUseBisection, bTauCycle,
     >        bReduceDone, bThinImprove, bFixGEFF, bHScaleOnly,
     >        bFULLHYDROSTAT, bGAMMARADMEAN, bPrintChanges, 
     >        bForceTAUMAX, bPrintHYST, bTauInterval, 
     >        bTMTHOM, bTMNoOS, bHYSTloose, bUseENTOT, bRCON,
     >        bARADCHECK, bOldThinGrid, bVDAMP, bWRCnext
      
      REAL, EXTERNAL :: WRVEL   !velocity field function 

C***  Tiefenabhaengiges Clumping nach Goetz Graefener
      REAL, DIMENSION(ND), INTENT(INOUT) :: DENSCON, FILLFAC       
      REAL, INTENT(IN) :: DENSCON_FIX
      CHARACTER(80), INTENT(IN) :: DENSCON_LINE, ThinCard

C***  !contains all RADIUS-GRID CARDS (for subr. RGRID)
      CHARACTER(80), DIMENSION(3) :: RadiusGridParameters           

      REAL, DIMENSION(NP), INTENT(INOUT) :: P
      REAL, DIMENSION(ND*NP), INTENT(INOUT) :: Z

C***  The following variables are just declared to call PRIMOD
      LOGICAL :: OLDTEMP, BTWOT, TTABLE, FALSEdummy
      REAL :: TFAC, TEFFdummy, ZEROdummy, ZEROdummy2
      INTEGER :: JOBNOLD, JOBNOLD2
      CHARACTER(100) :: MODHEAD, MODOLD, MODOLD2

C***  Constants
      REAL, PARAMETER :: AMU = 1.66E-24        !Atomic mass unit (gramm)     
      REAL, PARAMETER :: PI4 = 12.5663706144   !PI4 = 4*PI
      REAL, PARAMETER :: STEBOL = 5.6705E-5    !STEFAN-BOLTZMANN CONSTANT (cgs)
      REAL, PARAMETER :: XLSUN = 3.85E33       !Solar Luminosity (CGS-Units)
      REAL, PARAMETER :: RSUN = 6.96E10        !Solar Radius (cm)
      REAL, PARAMETER :: BOLTZ = 1.38E-16      !BOLTZMANN CONSTANT (ERG/DEG)
      REAL, PARAMETER :: GCONST = 6.670E-8  !GRAVITATION CONSTANT (CGS UNITS)
      REAL, PARAMETER :: XMSUN = 1.989E33   !XMSUN = Solar Mass (g)
      REAL, PARAMETER :: RGAS = BOLTZ / AMU !Gas Constant (CGS) = BOLTZ / AMU

      !File and channel handles (=KANAL)
      INTEGER, PARAMETER :: hOUT = 6        !write to wruniqX.out (stdout)
      INTEGER, PARAMETER :: hCPR = 0        !write to wruniqX.cpr (stderr)

      !COMMON-Block fuer Geschwindigkeitsfeldparameter
      COMMON /VELPAR/ VFINAL, VMINCAND, BETA, 
     >                VPAR1cand, VPAR2cand, RCON, HSCALE, 
     >                BETA2, BETA2FRACTION, VPAR1_2, VPAR2_2

      !This common block is needed for WRVEL
      COMMON /COMRADI/ OLDRADI
      COMMON /COMVELO/ bNewVelo, NDVALUE, OLDVELO

C***  COMMON /COMTEFF/  TRANSFERS THE EFF. TEMPERATURE TO: PRIMOD
      COMMON /COMTEFF/ TEFFdummy,ZEROdummy,ZEROdummy2,FALSEdummy

      TEFFdummy = TEFF
      ZEROdummy = 0.
      ZEROdummy2 = 0.
      FALSEdummy = .FALSE.
      bForceTAUMAX = .FALSE.
      
C***  interpolate POPNUM, GEFFL over ENTOT (if true), or RADIUS (if false)
      bUseENTOT = .FALSE.      

      bPrintHYST = .FALSE.     !default: do not print hydrostatic situation
      bPrintChanges = .FALSE.  !default: do not print changes to V and T
      
C***  if true, TAUMAX iteration is performed in only every NTAUIV-th job
      bTauInterval = .FALSE.    
      NTAUIV = 5

C***  if true, overshooting is prevented/damped by further GAMMARAD reducing
C***   Changed to FALSE by wrh  4-Apr-2019, as suggested by Andreas 
      bTMNoOS = .FALSE.          
      
C***  if true, iteration is done based on TAUTHOM instead of TAUROSScont
      bTMTHOM = .FALSE.         

C***  if true, ensuretaumax is not enforced if TAUMAX is good, but HYST not
      bHYSTloose = .FALSE.      

C***  if true, damping applies to the whole v(r) changes, otherwise toVMIN only
      bVDAMP = .TRUE.           

C***  Determine whether the old grid is already based on a VELTHIN call
C***  (i.e. HYDROSTATIC INTEGRATION card was used)
C***  In this case, the VCRIT vector should be filled and at
C***  least one point (i.e. the innermost) should be marked
C***  as being from a former hydrostatic (HS) integration
      bOldThinGrid = (VCRIT(ND)(1:2) == 'HS')

C***  Determine if next job should be WRCONT
      bWRCnext = (NCOLIP >= NEWWRC)
      IF (.NOT. bWRCnext) THEN
        DO L=1, ND
C***      If any of the Scharmer iterations is not converged,
C***      there will also be a WRCONT afterwards.
C***      Therefore we can also do a stratification update
C***      (This can be switched off by the SAFE option on the TAUMAX card)
          IF (CKONVER(L) /= 'C') THEN
            bWRCnext = .TRUE.
            EXIT
          ENDIF
        ENDDO
      ENDIF
      
C***  Ensure useful format for TAU output, even for very large TAUMAX      
      IF (TAUMAX > 9999.) THEN
        CTFMT = 'G9.5'
      ELSE      
        CTFMT = 'F8.4'
      ENDIF

           
      
      NDVALUE = ND
      RSTARSU = RSTAR / RSUN

C***  NC = number of core rays
      NC = NP - ND

C***  ATMEAN = mean atomic weight (without electrons)
      ATMEAN = .0
      DO NA=1, NATOM
        ATMEAN = ATMEAN + ABXYZ(NA)*ATMASS(NA)
      ENDDO
C***  XMU = mean particle mass (atoms and electrons)
      DO L=1, ND
        XMU(L) = ATMEAN / (1. + RNE(L))
      ENDDO

C***  RHO = mass density from current ENTOT
      DO L=1, ND
        RHO(L) = ENTOT(L) * AMU * ATMEAN
      ENDDO

C***  RI = interstice radius grid; DR = interval (L, L+1)
      DO L=1, ND-1
        RI(L) = 0.5 * ( RADIUS(L) + RADIUS(L+1) )
        DR(L) = RADIUS(L) - RADIUS(L+1)                
      ENDDO

C***  Mass flux per square centimeter at R* 
C***  The constant "3.02" converts M_sun/yr/(4pi Rsun^2) gram/cm^2/s 
      FM = 10.**(XMDOT+3.02) / RSTARSU**2 

      RCON = RCONorg
      IF (VMIN < 0.) THEN
          VMIN = VMINCAND
      ENDIF
      
      !This calculates the continuum tau-rosseland vector
      CALL TAUSCAL(RSTAR, ND, RADIUS, RNE, ENTOT,
     >              T, POPNUM, NDIM, N, EN, LEVEL, NCHARG, WEIGHT,
     >              ELEVEL, EION, EINST, ALPHA, SEXPO,
     >              ADDCON1, ADDCON2, ADDCON3, 
     >              IGAUNT, NOM, NF,
     >              XLAMBDA, FWEIGHT,
     >              TAUTHOM, TAUROSScont,
     >              MAXATOM, SIGMATHK, SEXPOK, EDGEK, KODAT,
     >              KONTNUP, KONTLOW, LASTKON, 
     >              DENSCON, FILLFAC
     >     )
     
C***  Consistency check: Is VFINAL kept? Stop otherwise!
C***  (The terminal velocity can only be changed via WRSTART or
C***   the HYDROSOLVE subroutine.)
      IF (ABS(VFINAL - VELO(1)) > 1.) THEN
        WRITE (hCPR,'(A,A)') '*** ERROR: VFINAL from CARDS file ',
     >          'does not match with the current terminal velocity!'
        WRITE (hCPR,'(A)') ' Please run a new WRSTART if ' //
     >          'you want to force a new terminal velocity!'
        STOP 'FATAL ERROR in ENSURETAUMAX'        
      ENDIF
      

C-  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -    
C***         start of ENSURETAUMAX criteria check
C-  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -    
     
C***  Accuracy check of hydrostatic equation in hydrostatic part      
C***  For each inner point, the fulfillment of the hydrostatic
C***  integration is checked. If the equation is violated by more then
C***  what is specified in the HYSTACC constant (see header), the 
C***  depth point is considered as "bad"
C***  If the HYDROSTATIC INTEGRATION card is set with the FULL option,
C***  one or more "bad" points will lead to a velocity update, regardless
C***  of whether TAUMAX is accurate or not.
      HYSTGOOD = .TRUE.
      NHYSTBAD = 0
      DO L=1, ND-1        
        IF (AMECH(L) > .0) THEN
           MECHRATIOLOG = LOG10(AMECH(L)/AGRAV(L))
        ELSE
           MECHRATIOLOG = -3.
        ENDIF
C***    Only points inside of RCON are considered
C***    Furthermore such points where AMECH plays a non-negligible 
C***      role and therefore the situation is not really hydrostatic 
C***      are also not considered  (MECHRATIOLOG check)
        IF (RI(L) < RCON .AND. MECHRATIOLOG < -2.) THEN
          HYSTRATIO = MIN(GAMMARADmax, ARAD(L)/AGRAV(L)) 
     >               + APRESS(L)/AGRAV(L) - 1.
          IF (ABS(HYSTRATIO) > HYSTACC) THEN
            HYSTGOOD(L) = .FALSE.
            IF (bFULLHYDROSTAT) THEN
              ! Enforce TAUMAX iteration even with good TAUMAX if 
              ! hydrostatic equation is not fulfilled.
              ! (can be switched off with the "LOOSE" option)
              bForceTAUMAX = .NOT. bHYSTloose
            ENDIF
            NHYSTBAD = NHYSTBAD + 1
          ENDIF
          IF (bPrintHYST) THEN
            WRITE (hCPR,'(A, I3, 3X, F7.4, 2X, L1, 1(3X, F12.5))') 
     >        'HyS: ', L, HYSTRATIO, HYSTGOOD(L), VELO(L)
          ENDIF
        ENDIF
      ENDDO      
      
      IF (NHYSTBAD == 0 .AND. IHSSTATUS >= 0) THEN
C***    Convergence criterion for HYDROSTATIC INTEGRATION:
C***    If all the hydrostatic eq. is matched with the required
C***    accuracy this convergence criterion is fulfilled
        IHSSTATUS = 1
      ENDIF

C***  Setup proper output labeling
C***  if (hardcoded) switch bTMTHOM is .TRUE., the TAUMAX value
C***  specified in the CARDS file refers to the THOMSON tau scale.
C***  Otherwise it refers to the full Rosseland continuum tau. (default)
      IF (bTMTHOM) THEN
        TAUCURMAX = TAUTHOM(ND)
        TAUFLAB = ' TAUTHOM'
      ELSE
        TAUCURMAX = TAUROSScont(ND)
        TAUFLAB = ' TAUcont'
      ENDIF
      
C***  Interval option for velocity update
C***  Instead of updating the velocity field as soon as one of the other
C***  criteria enforces this, it is also possible to restrict the run of
C***  this routine to every N-th STEAL job.
C***  To switch this on, the bTauInterval logical must be .TRUE.
      IF (.NOT. bWRCnext) THEN
            WRITE (hCPR,FMT='(A,F6.1,A,I2,A,F10.5)') 
     >      'TAUMAX update postponed; current values: ' //
     >        'TAUMAX=', TAUMAX,
     >        '  ' // TAUFLAB // '(',ND ,')=', TAUCURMAX
            bTauUpdated = .FALSE.
         RETURN
      ENDIF
           
C***  Check for CORRLIMIT option on the TAUMAX card           
C***  If set, the CORRLIMIT value is compared to CORMAX from the 
C***  Scharmer iteration and the velocity field update is skipped if 
C***  the corrections are higher than the specified limit.
C***  Note that CORRLIMIT is interpreted as a LOGARITHMIC value
      IF (bENSURETAUXMAX .AND. TauCorLimits(1) > -99.) THEN
        IF (CORMAX > 1.E-100) THEN
          CORLOG = LOG10(CORMAX)
        ELSE
          CORLOG = -100.
        ENDIF
        IF (CORLOG > TauCorLimits(1)) THEN
          WRITE (hCPR, FMT='(A, A, F6.2, A, F6.2, A)') 
     >      'ENSURETAUMAX skipped due to high corrections:',
     >      ' LOG COR= ', CORLOG, '  (Limit: ', TauCorLimits(1), ')'
          WRITE (hCPR,FMT='(A,F10.5,A,I2,A,F10.5)') 
     >        ' current diagnosis: TAUMAX=', TAUMAX,
     >        '  ' // TAUFLAB // '(',ND ,')=', TAUCURMAX
          bTauUpdated = .FALSE.
          RETURN
        ENDIF
      ENDIF
      
C***  Check for SAFE option on TAUMAX CARD
C***  If the "SAFE" is set, the velocity field update is only performed
C***  if all depth points during the Scharmer iteration are converged
      IF (bTauMaxSafe) THEN
        DO L=1, ND
          IF (CKONVER(L) /= 'C') THEN
            WRITE (hCPR, FMT='(A)') 
     >        'ENSURETAUMAX skipped due non-converged depth points.'
            WRITE (hCPR,FMT='(A,F10.5,A,I2,A,F10.5)') 
     >        ' current diagnosis: TAUMAX=', TAUMAX,
     >        '  ' // TAUFLAB // '(',ND ,')=', TAUCURMAX
            bTauUpdated = .FALSE.
            RETURN
          ENDIF
        ENDDO
      ENDIF

C***  Flux ratio check: (turned off by default)
C***  If the FQL option has been set on the TAUMAX card, the update 
C***  of the velocity stratification is only performed if the flux
C***  consistency is better than the specified limit
C***  (Example: TAUMAX FIX EPS=0.01 FQL=0.05 )
      IF (FQLIMIT > 0.) THEN
        DO L=2, ND-2      !Aussenpunkte weglassen
          !Flusskonstanz pruefen
          FQ = HTOTL(L) / HTOTCMF0(L)
C          WRITE (hCPR,*) 'FQ(L) ', L , FQ
          IF (FQ > (1. + FQLIMIT) .OR. FQ < (1. - FQLIMIT)) THEN
            WRITE (hCPR, FMT='(A)') 
     >        'ENSURETAUMAX skipped because flux conservation too poor'
            WRITE (hCPR,FMT='(A,F10.5,A,I2,A,F10.5)') 
     >        ' current diagnosis: TAUMAX=', TAUMAX,
     >        '  ' // TAUFLAB // '(',ND ,')=', TAUCURMAX
            bTauUpdated = .FALSE.
            RETURN
          ENDIF
        ENDDO
      ENDIF
      
C***  The original ARAD vector is not touched, but to allow
C***  that modifications on the values entering the calculations
C***  can be performed, the original vector is copied into
C***  ARADmod. The rest of this routine then uses ARADmod instead of ARAD
      DO L=1, ND
        ARADmod(L) = ARAD(L)
      ENDDO
      
C***  Interpolate over negative ARAD (and negative FLUX) points
      bARADCHECK = .FALSE.
      LL=0
      LCHECK=0
      DO L=ND-1, 1, -1
        IF (.NOT. bFULLHYDROSTAT) EXIT
        IF (bGAMMARADMEAN) EXIT
        bARADCHECK = .TRUE.
        TAUINT = 0.5 * (TAUROSS(L) + TAUROSS(L+1))
        LCHECK = LCHECK + 1
        IF (ARAD(L) > 0. .AND. HTOTL(L) > 0.) THEN
C***      Store "good" ARAD values in HELP arrays for interpolation  
          LL = LL + 1
          ARADHELP(LL) = ARAD(L)
          RHELP(LL) = RI(L)
C***      Note that HELP arrays are inside out in radius direction
        ENDIF
      ENDDO      

C***  Interpolation
C***  for boundaries: take last "good" value
C***  There must always be at least two good points
C***  (if LL == LCHECK all ARAD values are fine and nothing needs to be done!)
      IF (LL < LCHECK .AND. LL >= 2) THEN
        WRITE (hCPR, FMT='(A)') 
     >      'WARNING: ENSURETAUMAX uses interpolated a_rad.'
        DO L=1, ND-1
          IF (RI(L) < RHELP(1)) THEN
            ARADmod(L) = ARADHELP(1)
          ELSEIF (RI(L) > RHELP(LL)) THEN
            ARADmod(L) = ARADHELP(LL)
          ELSE 
            CALL SPLINPOX(ARADmod(L),RI(L),ARADHELP,RHELP,LL)
          ENDIF
        ENDDO               
      ELSEIF (bARADCHECK .AND. LL < 2) THEN
        WRITE (hCPR, FMT='(A)') 
     >   'WARNING: BAD ARAD => ENSURETAUMAX SKIP FORCED'        
          WRITE (hCPR,FMT='(A,F10.5,A,I2,A,F10.5)') 
     >        ' current diagnosis: TAUMAX=', TAUMAX,
     >        '  ' // TAUFLAB // '(',ND ,')=', TAUCURMAX
        bTauUpdated = .FALSE.
        RETURN       
      ENDIF
      
C***  Check the accuracy range for TAUMAX     
      TAUCHECK = TAUMAX - TAUCURMAX
      IF (bTauStrict) THEN
C***     STRICT option (implied by FIX) demands TAUMAX to be approached
C***     NOT STRICT allows TAUROSScont(ND) > TAUMAX
C***     This can be switched by writing 
C***          TAUMAX=x.x MIN     instead of   TAUMAX=x.x FIX  
C***     in the CARDS file
         TAUCHECK = ABS(TAUCHECK)
      ENDIF
      
      IF (.NOT. bENSURETAUXMAX) THEN
C***    TAUMAX cards line does neither include FIX, nor MIN 
C***     => no update of the velocity field
        bTauUpdated = .FALSE.
        WRITE (hCPR,FMT=69) TAUMAX, TAUFLAB, ND, TAUCURMAX
   69       FORMAT ('ENSURETAUMAX: not enabled - diagnosis only',/,
     >          'TAUMAX=', F10.5,'  ', A12, '(',I2 ,')=', F10.5)        
        RETURN
      ELSEIF (TAUCHECK <= TAUACC .AND. .NOT. bForceTAUMAX) THEN
C***    TAUMAX is still in specified accuracy range
C***      and no update is forced from the hydrostatic integration
C***    => EXIT subroutine
        bTauUpdated = .FALSE.
        WRITE (hCPR,FMT=70) TAUMAX, TAUFLAB, ND, TAUCURMAX
   70       FORMAT ('ENSURETAUMAX: no corrections needed',/,
     >          'TAUMAX=', F10.5,'  ', A12, '(',I2 ,')=', F10.5)        
        RETURN
      ELSE
C***    Velocity field will be updated 
C***    => print current TAUSCAL values and update criterion
        TMENFCRIT = 'TAU '
        IF (bForceTAUMAX .AND. NHYSTBAD > 0) THEN
          WRITE (hCPR, FMT='(A, I3, A)') 'Hydrostatic equation not '
     >      // ' fulfilled for ', NHYSTBAD, ' depth points.'
          TMENFCRIT = 'HYST'
        ENDIF
        WRITE (hCPR, FMT='(A,F10.5,A,I2,A,F10.5,2A)')
     >     'ENSURETAUMAX enforced: TAUMAX=', TAUMAX,
     >     '  ' // TAUFLAB // '(', ND, ')=', TAUCURMAX,
     >     '   Criterion=', TMENFCRIT
      ENDIF
      
C-  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -    
C***         end of ENSURETAUMAX criteria check
C-  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -    

C***  Save old Radius, velocity, and TAUROSScont for later POPNUM and 
C***  Temperature interpolation onto a changed grid
      DO L=1, ND
        RADIUSorg(L) = RADIUS(L)
        VELOorg(L) = VELO(L)
        TAUSCALorg(L) = TAUROSScont(L)
        ENTOTorg(L) = ENTOT(L)
        Torg(L) = T(L)
        RNEorg(L) = RNE(L)
      ENDDO
      DO K=1, NF
       DO L=1, ND
          XJCorg(L,K) = XJC(L,K)
        ENDDO
      ENDDO
      XMSTARorg = XMSTAR
      GLOGorg = GLOG
      GEFFLOGorg = GEFFLOG
      GAMMARADMEANcalc = GAMMARADMEAN
C***  Mis-use of the scratch array which is dimensioned in STEAL as
C***  (NDIM+2)*(NDIM+2)
      CALL BACKUPPOPNUM(ND, N, POPNUM, SCRATCH)      

C***  Calculate Eddington Gamma and g_eff
      XLSTAR = PI4 * STEBOL * RSTAR*RSTAR * TEFF*TEFF*TEFF*TEFF
      XLSTARS = XLSTAR / XLSUN

C***  If the previous new solution of the HYDROSTATIC equation 
C***   resulted in a change of the sign of TAUMAX - TAUCURMAX, 
C***   the mixing fraction of the new ARAD to the old one is reduced
C***   and the solution is repeated ("overshooting")           
      GEDDreducenow = GEDDreduce
  22  CONTINUE   !ENTRY TO JUMP BACK IF automatic GEDD reduce is enforced
  
C***  Restore all values after a possible
C***   GAMMARAD-reducing caused by "overshooting"
      XMSTAR = XMSTARorg
      GLOG = GLOGorg
      GEFFLOG = GEFFLOGorg
      GAMMARADMEAN = GAMMARADMEANcalc
      DO L=1, ND
        RADIUS(L) = RADIUSorg(L)
        VELO(L) = VELOorg(L)
        IF (bTMTHOM) THEN
          TAUTHOM(L) = TAUSCALorg(L)
        ELSE
          TAUROSScont(L) = TAUSCALorg(L)
        ENDIF
        ENTOT(L) = ENTOTorg(L)
        T(L) = Torg(L)
        RNE(L) = RNEorg(L)
        XMU(L) = ATMEAN / (1. + RNE(L))
      ENDDO
      
C***  If reduce option is enabled, only this fraction of the new
C      GAMMARAD is used (GEDDRAD contains value from last step)
C     Note: This works only if the MEAN option is set!
      IF (GEDDreducenow /= 1. .AND. GEDDRAD > 0.) THEN
        GAMMARADMEAN = GEDDreducenow * GAMMARADMEAN 
     >                    + (1. - GEDDreducenow) * GEDDRAD
      ENDIF

C***  Limit GAMMARADMEAN to a maxmimum value (otherwise situation can get unstable)
      GAMMARADMEAN = MIN(GAMMARADMEAN, GAMMARADmax)         

      fGAMMACOR = 1.          

      
C***  Update of  g_eff or g_grav plus M_star
C***  Note: Mass update can provide a correction factor to GAMMARAD
C***  GAMMARADMEAN and QIONMEAN were calculated in INITFCORR
      IF (.NOT. bFixGEFF) THEN
        IF (.NOT. bGEddFix) THEN
          IF (bFULLHYDROSTAT) THEN
            GEDD = GAMMARADMEAN
          ELSE
            GEDD = 10.**(-4.51) * QIONMEAN * XLSTARS / XMSTAR
            IF (GEDD >= 1.) THEN
              WRITE (hCPR,'(A,A)') '*** Star exceeds Eddington limit!',
     >          'Stellar mass might be too low for this Luminosity.'
              WRITE(hCPR,'(A,F7.2)') '  Suggestion: Set MSTAR > ', 
     >                   10**(-4.51) * QIONMEAN * XLSTARS
              STOP 'FATAL ERROR in ENSURETAUMAX'
            ENDIF
          ENDIF
        ENDIF
        GEFFLOG = LOG10( (10.**GLOG) * (1. - GEDD) ) 
      ELSEIF (bGEddFix) THEN
        !GEFF and GEDD fixed by wrstart
        GLOG = LOG10(10**GEFFLOG / (1. - GEDD))
        XMSTARG = 10.**GLOG * RSTAR * RSTAR / GCONST
        XMSTAR = XMSTARG / XMSUN
      ELSE
        !GEFF was fixed by wrstart => Calculate GEDD and new mass
        CALL CALCMASSFROMGEFF(bFixGEFF, bFULLHYDROSTAT,
     >                        XMSTAR, GLOG, GEFFLOG, 
     >                        ARAD, APRESS, AGRAV, RADIUS, ND,
     >                        RSTAR, RCON, RI, TAUROSS,
     >                        GAMMARADMEAN, GEDD,
     >                        LOG10(XLSTARS), QIONMEAN, fGAMMACOR)
        IF (bFULLHYDROSTAT) THEN
          GEDD = GAMMARADMEAN
        ENDIF
      ENDIF      

      
C***  calculate depth-dependend GEDD now and store it in the
C***  GEFFL vector (g_eff per depth point). This is the vector 
C***  which is later passed to the VELTHIN routine where the 
C***  hydrostatic integration is performed.
C***  Note: Back interpolation from ARAD, AGRAV on interstices
C***        to ARADL, AGRAVL on full depth points
      DO L=1, ND
        IF (bFULLHYDROSTAT) THEN
C***      HYDROSTATIC INTEGRATION with FULL option is set
          IF (.NOT. bGEddFix) THEN
C***        EDDINGTON-GAMMA is not fixed via (direct or indirect) CARDS input          
            IF (bGAMMARADMEAN) THEN
C***          MEAN option has been set on the HYDROSTATIC INTEGRATION card            
              GAMMARAD(L) = GAMMARADMEAN
            ELSE
              IF (RADIUS(L) > RI(1)) THEN
                ARADL = ARADmod(1)
                AGRAVL = AGRAV(1)
              ELSEIF (RADIUS(L) < RI(ND-1)) THEN
                CALL SPLINPOX(ARADL, RADIUS(ND-1), 
     >                        ARADmod, RI, ND-1, DFDX=DARADDR)
                CALL SPLINPOX(AGRAVL, RADIUS(ND-1), 
     >                        AGRAV, RI, ND-1, DFDX=DAGRAVDR)
                ARADL = ARADmod(ND-1) - DARADDR * DR(ND-1)
                AGRAVL = AGRAV(ND-1) - DAGRAVDR * DR(ND-1)
              ELSE
                CALL SPLINPOX(ARADL, RADIUS(L), ARADmod, RI, ND-1)
                AGRAVL = GCONST * XMSTAR * XMSUN / 
     >                   (RADIUS(L)*RSTAR)**2
              ENDIF
              
C***          GAMMARAD calculation is corrected by a factor fGAMMACOR  
C***          which differs from one if the MASS has been updated 
C***          (This happens if LOG GEFF is fixed)
              GAMMARAD(L) = ARADL / AGRAVL * fGAMMACOR

C***          Depth-dependend GAMMARAD reducing (i.e. mixing with older value)              
C***          GRSTATIC(L) - vector with GAMMARAD from LAST TAUMAX iteration
C***          GEDDRAD     - GAMMARADMEAN from LAST TAUMAX iteration
              IF (GEDDreducenow /= 1. .AND. GRSTATIC(L) > 0.) THEN
                GEDDreduceL = GEDDreducenow
                IF ( .NOT. HYSTGOOD(L)) THEN
                  !Less strict reducing for depth points where the hydrostatic 
                  !  equation is not yet fulfilled
                  GEDDreduceL = MIN(1., GEDDreducenow * 2.)
                ENDIF
                GAMMARAD(L) = GAMMARAD(L) * GEDDreduceL
     >             + (1. - GEDDreduceL) * GRSTATIC(L)
              ELSEIF (GEDDreducenow /= 1. .AND. GEDDRAD > 0.) THEN
                !After the model start we might have only GEDDRAD (= old GAMMARADMEAN)
                GAMMARAD(L) = GAMMARAD(L) * GEDDreducenow
     >             + (1. - GEDDreducenow) * GEDDRAD
              ENDIF
              
C***          Make sure that eventual GAMMARAD vector does not exceed the GAMMARADmax limit
              GAMMARAD(L) = MIN(GAMMARAD(L), GAMMARADmax)              
            ENDIF
          ELSE
C***        EDDINGTON-GAMMA is fixed (directly or implicit)
C***        => we must take this value          
            GAMMARAD(L) = GEDD
          ENDIF
          
C***      Prepare finale GEFFL vector (for VELTHIN call)
          GEFFL(L) = (10**GLOG) * (1. - GAMMARAD(L))
          IF (GEFFL(L) < 0.) THEN
C***        avoid any values which could lead to a non-monotonic velocity law
            GEFFL(L) = 0.
          ENDIF
        ELSE 
C***      HYDROSTATIC INTEGRATION is given without FULL option        
          IF (bGAMMARADMEAN) THEN
C***        MEAN option has been set on the HYDROSTATIC INTEGRATION card            
            GEFFL(L) = 10**(GEFFLOG)
          ELSE
            GEFFL(L) = (10**GLOG) 
     >        - RNE(L)/(ATMEAN*QIONMEAN) * ( (10**GLOG) - (10**GEFFLOG) ) 
          ENDIF
        ENDIF
      ENDDO
      
      
C***  Initialize Iteration values
      ITTAU = 0         ! TAUMAX iteration counter
      TROLD = 0.        ! Initvalue for last-iteration TAUMAX
      VMINOLD = 0.      ! Value from last iteration
      VMINOLD2 = 0.     ! VMINOLD-Value from last iteration
      
      bHydroStat = .TRUE.      !default value: hydrostatic domain is achieved

      !prepare bisection    
      bUseBisection = .FALSE.   !FALSE = bisection method is currently not in use
      VMIN_over = -1.0
      VMIN_under = -1.0
      VMIN_min = VMIN
      VMIN_max = VMIN
      IT_failsafe = 0
      LastDAMP = 2
      STEPDAMP = 1.0

C***  if true, the velocity corrections have been damped      
      bReduceDone = .FALSE.

      WRITE (hCPR,FMT=72) GEDD, GLOG, GEFFLOG
   72   FORMAT ('Current Eddington Gamma: ', F6.3,
     >          '   log g_grav = ', F10.5, 
     >          '   log g_eff = ', F10.5)
     
      IF (RCON > RADIUS(1)) THEN
        !Give a warning as we might overwrite a hydro solution here
        WRITE (hCPR,*) 'ENSURETAUMAX WARNING> Old velocity might not '
     >    // 'be a beta-type law'
        RCON = RADIUS(1)    !this is needed to avoid SPLINPOX errors
        bRCON= .FALSE.
      ELSE
        bRCON= .TRUE.
      ENDIF

C***  Major update Mar 2016: In contrast to earlier versions,
C***  the whole iteration is now always performed on the old grid
C***  and only afterwards a new grid is established if this gives
C***  a better representation. Therefore major parts after this
C***  line were re-written and could be streamlined!

C***  Note that VMIN is not initialized in this routine since 
C***  this is already done by REMOST and DECSTE which ensure
C***  that the first VMIN is the current VELO(ND)
      
C     ------- Main TAU iteration loop ----------------------
      mainitloop: DO 
          ITTAU = ITTAU + 1
          bTauCycle = .FALSE.  

          !RESET VMIN if outside of boundaries 0...VFINAL
          IF (VMIN < 0.) VMIN = 1.E-4
          IF (VMIN > VFINAL) VMIN = 0.9 * VFINAL
        
          IF (ITTAU == MAXITTAU) THEN
C***        TAUMAX iteration did not converge 
C***        detect from the history the VMIN value which gave 
C***        the closest approximation to TAUMAX
            DIFFTAUMAX_MIN = ABS(TAUMAXHIST(1)-TAUMAX) 
            IT_failsafe = 1
            DO J=2, ITTAU-1
              DIFFTAUMAX_J = ABS(TAUMAXHIST(J)-TAUMAX)
              IF (DIFFTAUMAX_J < DIFFTAUMAX_MIN) THEN
                DIFFTAUMAX_MIN = DIFFTAUMAX_J
                IT_failsafe = J
              ENDIF  
            ENDDO
            WRITE(hCPR,*) 'TAUMAX iteration not converged: ',
     >                     'best iteration used, no reduction'
            VMIN = VMINHIST(IT_failsafe)
C***        Damping the correction may be counter-productive, when 
C***        the relation between VMIN and TAUMAX was non-monotonic!
C***        Therefore:
            bReduceDone = .TRUE.
          ENDIF

C***      Transfer VMIN to VELPAR common block          
          VMINCAND = VMIN

C***      VMIN history and min/max values          
          VMIN_min = MIN(VMIN, VMIN_min)
          VMIN_max = MAX(VMIN, VMIN_max)
          VMINHIST(ITTAU) = VMIN    
                                                       
C***      Initialize VELPAR common block parameters 
C***      VPAR1, VPAR2, RCON, HSCALE
C***      (for 2beta-laws also VPAR1_2 and VPAR2_2)
C***      for HYDROSTATIC INTEGRATION only the scale height H_0 is calculated
          bHScaleOnly = bThinWind
          CALL INITVEL (RMAX, TEFF, GEFFLOG, RSTAR, XMU(ND),
     >                  VTURB, bHScaleOnly, bHydroStat) 

C***      bThinWind = HYDROSTATIC INTEGRATION card has been found            
          IF (bThinWind) THEN
C***        The VELTHIN routine solves the HYDROSTATIC EQUATION
C***        and connects the result to a BETA law in the outer part.
C***        The default connection criterion a continouus gradient
C***        between both velocity laws
C***        (can be changed to a fraction of the sonic speed via CARDS option)
            CALL VELTHIN (T, RADIUS, VELO, ND, RSTAR, RMAX,
     >                    GEFFL, XMU, VTURB, ThinCard, VELOCRITERION)
     
          ELSE
C***        Simple exponential law with constant scale height 
C***        (i.e. barometric formula) in the inner part:
C***        This is directly calculated in the WRVEL function
C***        The necessary parameters were calculated in INITVEL and are
C***        transfered via the VELPAR common block     
            bNewVelo = .TRUE.
            DO L=1, ND
              VELO(L) = WRVEL(RADIUS(L))
C***          Fill velocity criterion info vector
              IF (RADIUS(L) <= RCON) THEN
                VELOCRITERION(L) = 'STATIC  '
              ELSE 
                VELOCRITERION(L) = ''
              ENDIF
            ENDDO
                
          ENDIF

C***      Adjust the density vectors after the velocity update
          DO L=1, ND, 1
            RL = RADIUS(L)
            RHO(L) = FM / (VELO(L) * 1.E5 * RL * RL) 
            ENTOT(L) = RHO(L) / (AMU*ATMEAN)
            
C***        Furthermore, the informative VCRIT model vector 
C***        is updated here: 
            SELECTCASE (VELOCRITERION(L))
              CASE ('HYSTINT ')
                VCRIT(L) = 'HS      '
              CASE ('STATIC  ')
                VCRIT(L) = 'ST      '
              CASE ('BETA    ')
                VCRIT(L) = 'B       '
              CASE ('2BETA   ')
                VCRIT(L) = '2B      '
              CASE ('SQRT    ')
                VCRIT(L) = 'R       '
              CASE DEFAULT
                IF (RADIUS(L) > RCON) THEN
                  IF (BETA2FRACTION > 0.) THEN
                    VCRIT(L) = '2B      '
                  ELSE
                    VCRIT(L) = 'B       '
                  ENDIF
                ENDIF
            ENDSELECT
          ENDDO

C***      The new velocity CAN lead to a different depth-dependend clumping
C***      As we want to calculate TAUROSS afterwards with the new DENSCON,
C***      it is not possible to account here for effects in DENSCON which 
C***      stem from an update of the TAU scale. This can only be accounted 
C***      for in the final CLUMP_STRUCT call at the end of this routine
C***      and thus still lead to a small Inconsistency which has to be
C***      resolved in the next JOB with a TAUMAX iteration.
          CALL CLUMP_STRUCT (DENSCON, FILLFAC, ND, DENSCON_FIX, 
     >                 VELO, TAUSCALorg, DENSCON_LINE, RADIUS, T, XMU)

C***      Calculation of the new TAU scale (continuum opacities only)
          CALL TAUSCAL(RSTAR, ND, RADIUS, RNE, ENTOT,
     >              T, POPNUM, NDIM, N, EN, LEVEL, NCHARG, WEIGHT,
     >              ELEVEL, EION, EINST, ALPHA, SEXPO,
     >              ADDCON1, ADDCON2, ADDCON3, 
     >              IGAUNT, NOM, NF,
     >              XLAMBDA, FWEIGHT,
     >              TAUTHOM, TAUROSScont,
     >              MAXATOM, SIGMATHK, SEXPOK, EDGEK, KODAT,
     >              KONTNUP, KONTLOW, LASTKON, 
     >              DENSCON, FILLFAC
     >      )
          IF (bTMTHOM) THEN
            TAUCURMAX = TAUTHOM(ND)
          ELSE
            TAUCURMAX = TAUROSScont(ND)
          ENDIF
C***      check if current VMIN yields higher or lower TAUMAX than wanted
          IF (TAUCURMAX > TAUMAX) THEN
             VMIN_OVER = VMIN
          ELSE
             VMIN_UNDER = VMIN
          ENDIF 
 
          !Automatic GEDD reducing for overshooting situations
          IF ((ITTAU == 1) .AND. bFULLHYDROSTAT .AND. (.NOT. bFixGEFF)
     >      .AND. bTMNoOS
     >      .AND.  (GEDDreducenow > GEDDreduceMIN) .AND.
     >      ((TAUCURMAX-TAUMAX)*(TAUSCALorg(ND)-TAUMAX) < 0.)
     >       ) THEN
            WRITE (hCPR,'(A,2(' // CTFMT // ',A))') 
     >             'Reduding of GAMMARAD changes '
     >          // 'enforced due to Tau overshooting (', 
     >           TAUSCALorg(ND), ' vs. ', TAUCURMAX, ')'
            GEDDreducenow = GEDDreducenow / 2. 
            WRITE (hCPR,'(A,F7.4)') 'New GAMMARAD is accounted with '
     >         // ' a factor of GEDDreduce = ', GEDDreducenow
            GOTO 22
          ENDIF
 
          IF (bTMTHOM) THEN
            TAUMAXHIST(ITTAU) = TAUTHOM(ND)
          ELSE
            TAUMAXHIST(ITTAU) = TAUROSScont(ND)     !save taumax history
          ENDIF          

  
          !- - - - - - - - Pruefung des berechneten TAUMAX-Wertes - - - - - - - - -
          IF (TAUMAX > .0) THEN
            !TAUMAX adjustment (for velocity field) starts here

            IF ((IT_failsafe > 0) .AND. bReduceDone) THEN              
              !The failsafe iteration is of course the last iteration
              WRITE (hCPR,FMT='(A,I3,A,' // CTFMT // ',A,G12.4)') 
     >          '*** TAUMAX ITERATION: FAILSAFE=', IT_failsafe,
     >           TAUFLAB // '=', TAUCURMAX, '  VMIN=', VMIN 
               EXIT
            ELSEIF (ReduceTauCorrections /= 1. .AND. bReduceDone) THEN              
              !If correction damping has been done this was the last iteration
              WRITE (hCPR,FMT='(A,F4.2,A,' // CTFMT // ',A,G12.4,A)') 
     >          '*** TAUMAX ITERATION: REDUCE=', ReduceTauCorrections,
     >           TAUFLAB // '=', TAUCURMAX, '   VMIN=', VMIN, 
     >          ' damped'
               EXIT
            ELSE
              WRITE (hCPR,FMT='(A,I3,A,' // CTFMT // ',A,G12.4)') 
     >          ' *** TAUMAX ITERATION: ITTAU=', ITTAU, 
     >           TAUFLAB // '=', TAUCURMAX, '   VMIN=', VMIN 
            ENDIF

            !store old values to track current trend
            VMINOLD2 = VMINOLD
            VMINOLD = VMIN
            TROLD2 = TROLD
            TROLD = TAUCURMAX           

            IF (TAUCURMAX < 0.) THEN
              !total failure of last step => interpolation between the two steps before
              VMIN = 0.5 * (VMINHIST(ITTAU - 2) + VMINHIST(ITTAU - 1) )
              CYCLE
            ENDIF            

C***        If not yet converged:
            IF (ABS(TAUCURMAX-TAUMAX) > TAUACC) THEN
C***           Check if bisection is possible               
               IF (VMIN_OVER > .0 .AND. VMIN_UNDER > .0) THEN
                  VMIN = 0.5 * (VMIN_OVER + VMIN_UNDER)
                  BTAUCYCLE = .TRUE.
               ENDIF 
            ENDIF

C***          TAUROSScont TOO LARGE: New vmin is at v(TAUMAX) using spline interpolation
              IF (bTauStrict .AND. (.NOT. bTauCycle) .AND.
     >            (TAUCURMAX > TAUMAX+TAUACC)) THEN
                VMIN = VMIN*(TAUCURMAX/TAUMAX)**(.5*STEPDAMP)
c                !SPLINPO offen makes a too large step inside, try other things first
c                IF ((VMINOLD2 > 0.) .AND. (TROLD2 > 0.) .AND. 
c     >                           (TROLD2 < (TAUMAX-TAUACC))) THEN
c                  !Only last step was too large, try half step
c                  VMIN = (VMINOLD2 + VMINOLD) / 2.
c                ELSEIF (bTMTHOM) THEN
c                  ! standard method: 
c                  ! calculate new v_min by using v(tau) interpolation
c                  CALL SPLINPOX(VMIN, TAUMAX, VELO, TAUTHOM, ND)
c                ELSE
c                  ! standard method: 
c                  ! calculate new v_min by using v(tau) interpolation
c                  CALL SPLINPOX(VMIN, TAUMAX, VELO, TAUROSScont, ND)
c                ENDIF
                bTauCycle = .TRUE.
              ENDIF

C***          TAUROSScont TOO SMALL: New vmin by scaling with sqrt(tau(ND)/TAUMAX)
              !STEPDAMP can damp this adjustment by using higher roots
              IF ((.NOT. bTauCycle) 
     >             .AND. (TAUCURMAX < TAUMAX-TAUACC)) THEN
                VMIN = VMIN*(TAUCURMAX/TAUMAX)**(.5*STEPDAMP)
                bTauCycle = .TRUE.                
              ENDIF

            IF (bTauCycle .AND. (ITTAU < MAXITTAU)) THEN
              !check if steps should be made smaller            
              IF ((ITTAU - LASTDAMP) > 8) THEN
                Nct = 0
                DO J=LastDAMP, ITTAU
                  !check how often the trend for TAUROSScont(ND) has changed
                  IF (((TAUMAXHIST(J) - TAUMAX)
     >              * (TAUMAXHIST(J-1) - TAUMAX)) < 0. ) THEN
                    Nct = Nct + 1
                  ENDIF
                ENDDO
                IF (Nct > 4) THEN
                  !make smaller (damped) steps if there have been
                  ! more than four changes in the last eight iterations
                  STEPDAMP = STEPDAMP / 2.
                  LastDAMP = ITTAU
                ENDIF
              ENDIF
              CYCLE
            ENDIF

          ENDIF

          IF ((ReduceTauCorrections /= 1.) .AND. 
     >        (.NOT. bReduceDone)) THEN
            IF (bVDAMP) THEN
C***        Velocity field damping (can harm TAU value)     
C***        currently spoils VELOCRITERION information (can this be fixed?)
              
C***          Mix new velocity field with (1 - Reduce) of original field
C***          and adjust density (but keep grid spacing)
              DO L=1, ND
                VELO(L) = ReduceTauCorrections * VELO(L)
     >                      + (1. - ReduceTauCorrections) * VELOorg(L)
                RHO(L) = FM / (VELO(L) * 1.E5 * RADIUS(L)**2) 
                ENTOT(L) = RHO(L) / (AMU*ATMEAN)
              ENDDO
              
C***          Update VMIN, mix RCON with old value
              VMIN = VELO(ND)
              IF (bRCON) 
     >           RCON = ReduceTauCorrections * RCON
     >                        + (1. - ReduceTauCorrections) * RCONorg
              
C***          Any velocity update might affect the clumping stratification              
              CALL CLUMP_STRUCT (DENSCON, FILLFAC, ND, DENSCON_FIX, 
     >                 VELO, TAUROSScont, DENSCON_LINE, RADIUS, T, XMU)     
     
              !Berechnung der neuen TAU-Skala (continuum opacities only)
              CALL TAUSCAL(RSTAR, ND, RADIUS, RNE, ENTOT,
     >                T, POPNUM, NDIM, N, EN, LEVEL, NCHARG, WEIGHT,
     >                ELEVEL, EION, EINST, ALPHA, SEXPO,
     >                ADDCON1, ADDCON2, ADDCON3, 
     >                IGAUNT, NOM, NF,
     >                XLAMBDA, FWEIGHT,
     >                TAUTHOM, TAUROSScont,
     >                MAXATOM, SIGMATHK, SEXPOK, EDGEK, KODAT,
     >                KONTNUP, KONTLOW, LASTKON, 
     >                DENSCON, FILLFAC
     >        )
              IF (bTMTHOM) THEN
                TAUCURMAX = TAUTHOM(ND)
              ELSE
                TAUCURMAX = TAUROSScont(ND)
              ENDIF
              
              WRITE (hCPR,FMT='(A,F4.2,A,' // CTFMT // ',A,G12.4,A)') 
     >          '*** TAUMAX ITERATION: REDUCE=', ReduceTauCorrections,
     >           TAUFLAB // '=', TAUCURMAX, '   VMIN=', VMIN,
     >          ' mixed'              
              
              EXIT !exit mainitloop in case of VELO mixing
            ELSE     
C***        Classical damping method (VMIN + new iteration)
              !Reduce total vmin corrections if CARDS option has been set (default: by 50%)
              ! Note: This does NOT help to achieve the specified TAUMAX value, but
              !       damps the total changes in order to damp possible following
              !       corrections in the STEAL-COLI-iteration
              VMIN = ReduceTauCorrections * VMINHIST(ITTAU) +
     >                 ( 1. - ReduceTauCorrections ) * VMINHIST(1)
              bReduceDone = .TRUE.
              CYCLE            
            ENDIF
          ENDIF

          EXIT !exit TAU loop if no cycling criterion was met

      ENDDO mainitloop
C     ------- Main TAU iteration loop ends here ----------------------


C**** ---------------------------------------------------------------
C****  Check if a new grid is significantly better
C**** ---------------------------------------------------------------


      RCONorg = RCON    !Write back updated connection point

C***  Update WRVEL interpolation vectors with solution of the TAU iteration
      DO L=1, ND
        OLDRADI(L) = RADIUS(L)
        OLDVELO(L) = VELO(L)
      ENDDO
C***  To ensure pure interpolation in WRVEL (inside GEOMESH->RGRID)
C***  the following logical MUST be .FALSE. and RCON must be > RMAX
C***  (Both values are transfered via COMMON blocks.)
      bNewVelo = .FALSE.
      RCON = 1.1 * RADIUS(1)
      
C***  Calculate "quality" of the old radius grid
      QVorg = 0.
      QRHOorg = .0
      DO L=1, ND-1
        QVorg = MAX(QVorg, VELO(L) - VELO(L+1))
        QRHOorg = MAX(QRHOorg, LOG10(RHO(L+1)/RHO(L)))
      ENDDO

C***  To test for a new grid, bNoRGrid needs to be .FALSE.      
      bNoRGrid = .FALSE.    
      
      gridloop: DO
      
C***    Creation of radius grid (unless bNoRGrid == .TRUE.)  plus
C***    impact parameter grid, i.e. P(NPDIM) and Z(NDDIM,NPDIM) arrays
C***    updates the INCRIT vector which prints the depth point criteria

        CALL GEOMESH (RADIUS,INCRIT,P,Z,ND,NDDIM,NP,NPDIM,RMAX, 
     >                RadiusGridParameters, bNoRGrid, NC)         
            
C***          setup velocity and density stratification for new radius grid
        DO L=1, ND, 1
          IF (.NOT. bNoRGrid) THEN
            CALL SPLINPOX(VELO(L), RADIUS(L), OLDVELO, OLDRADI, ND)
          ENDIF
          RHO(L) = FM / (VELO(L) * 1.E5 * RADIUS(L) * RADIUS(L)) 
          ENTOT(L) = RHO(L) / (AMU*ATMEAN)
        ENDDO            

C***    Calculate "quality" of the new radius grid
        IF (.NOT. bNoRGrid) THEN
          QVnew = 0.
          QRHOnew = .0
          DO L=1, ND-1
            QVnew = MAX(QVnew, VELO(L) - VELO(L+1))
            QRHOnew = MAX(QRHOnew, LOG10(RHO(L+1)/RHO(L)))
          ENDDO

C***      If the new grid does not improve neither the density
C***      spacing, nor the velocity spacing, the old grid is
C***      kept (requires recalculation of GEOMESH and ENTOT)
          IF (QRHOnew > 0.9*QRHOorg .AND. QVnew > 0.9*QVorg) THEN
            DO L=1, ND
              RADIUS(L) = OLDRADI(L)
              VELO(L) = OLDVELO(L)
            ENDDO
            bNoRGrid = .TRUE.
            CYCLE gridloop
          ENDIF

        ENDIF
        
        EXIT gridloop
      
      ENDDO gridloop
      
C**** ---------------------------------------------------------------
C****  End of grid check and (potential) update
C**** ---------------------------------------------------------------

      RCON = RCONorg    !Restore updated connection point
            
C***  Backup old Temperatures for interpolation
      T1old = TOLD
      T2old = TOLD2
      T3old = TOLD3
      
C***  The status of bNoRGrid AFTER the gridloop can now be used
C***  to determine whether we have a new grid (i.e. bNoRGrid == .FALSE.)
C***  or the old grid is still used. Only in the first case, radius
C***  interpolations have to be performed
      IF (.NOT. bNoRGrid) THEN
            
C***    Interpolation of temperature on the new radius grid
        CALL INTERPOLATETEMP(T, Torg, RADIUS, RADIUSorg, ND)

C***    Interpolation of Popnumbers on new radius grid
        CALL INTERPOLATEPOPNUM(POPNUM, SCRATCH, POPMIN,
     >                         RADIUS, RADIUSorg, ENTOT, ENTOTorg,
     >                         N, ND, ABXYZ, NFIRST, NLAST, NATOM, 
     >                         bUseENTOT)

C***    re-calculate RNE(L) and XMU(L) from new POPNUMs
        DO L=1, ND
          RNEL=0.0
          DO J=1,N
            RNEL = RNEL + NCHARG(J) * POPNUM(L,J)
          ENDDO
          RNE(L)=RNEL
          XMU(L) = ATMEAN / (1. + RNE(L))
        ENDDO

C***    --- Interpolation of old temperatures on Radius-Grid ---
        CALL INTERPOLATETEMP(TOLD , T1old, RADIUS, RADIUSorg, ND)
        CALL INTERPOLATETEMP(TOLD2, T2old, RADIUS, RADIUSorg, ND)
        CALL INTERPOLATETEMP(TOLD3, T3old, RADIUS, RADIUSorg, ND)
        
        WRITE (hCPR,*) '*** ENSURETAUMAX: New radius grid established'
        
      ENDIF


      
C***  Additional tau interpolation from old to new scale
C***  in the optical thick part (defined by TAU > TAUTHICK)
      TAUTHICK = 0.1

C***  Calculate departure coefficients for the innermost
C***  depth point of the old stratification
      ENEL = RNEorg(ND) * ENTOTorg(ND)
      CALL LTEPOP (N, EN, Torg(ND), ENEL, 
     >             WEIGHT, NCHARG, EION, ELEVEL, NOM,
     >             ABXYZ, NFIRST, NLAST, NATOM) 
      DO J=1, N
        DEPARTNDorg(J) = SCRATCH(ND,J)/EN(J)
      ENDDO
    
      DO L=1, ND 
        IF (TAUROSScont(L) < TAUTHICK) CYCLE
        IF (TAUROSScont(L) <= TAUSCALorg(ND)) THEN
          CALL SPLINPOX(T(L), TAUROSScont(L), Torg, TAUSCALorg, ND)
          CALL SPLINPOX(TOLD(L), TAUROSScont(L), T1old, TAUSCALorg, ND)
          CALL SPLINPOX(TOLD2(L), TAUROSScont(L), T2old, TAUSCALorg, ND)
          CALL SPLINPOX(TOLD3(L), TAUROSScont(L), T3old, TAUSCALorg, ND)
          DO J=1, ND
            CALL SPLINPOX(POPNUM(L,J), TAUROSScont(L), 
     >                                SCRATCH(1,J), TAUSCALorg, ND)
          ENDDO
        ELSE
          Textrap= ((TAUROSScont(L)+0.67)/(TAUSCALorg(ND)+0.67))**(0.25)
          T(L)     = Torg(ND)  * Textrap
          TOLD(L)  = T1old(ND) * Textrap
          TOLD2(L) = T2old(ND) * Textrap
          TOLD3(L) = T3old(ND) * Textrap

          ENEL = RNE(L) * ENTOT(L)
C***      LTEPOP returns the LTE population numbers
C***        for depth point L to the vector EN  (here not ENLTE!)
          CALL LTEPOP (N, EN, T(L), ENEL, 
     >                 WEIGHT, NCHARG, EION, ELEVEL, NOM,
     >                 ABXYZ, NFIRST, NLAST, NATOM) 
          DO J=1, N
C***        Apply old departure coefficient to LTE popnums
            POPNUM(L,J) = EN(J) * DEPARTNDorg(J)
          ENDDO
        ENDIF        
        RNEL=0.0
        DO J=1, N
          RNEL = RNEL + NCHARG(J) * POPNUM(L,J)
        ENDDO
        RNE(L)=RNEL
        XMU(L) = ATMEAN / (1. + RNE(L))
      ENDDO
            
            
C***  Interpolation of older popnumbers on the new radius grid
C***  and interpolation on TAU in the optical thick part
      CALL REGRIDOLDPOP(POP1, T1old, SCRATCH, ND, N, T,
     >                  POPMIN, TAUTHICK, bNoRGrid,
     >                  RADIUS, RADIUSorg, ENTOT, ENTOTorg,
     >                  RNE, RNEorg, ABXYZ, NFIRST, NLAST, 
     >                  NATOM, WEIGHT, NCHARG, EION, ELEVEL, NOM, 
     >                  TAUROSScont, TAUSCALorg,
     >                  EN, DEPARTNDorg, bUseENTOT)

      IF (NTOLD > 1) THEN      
        CALL REGRIDOLDPOP(POP2, T2old, SCRATCH, ND, N, T,
     >                    POPMIN, TAUTHICK, bNoRGrid,
     >                    RADIUS, RADIUSorg, ENTOT, ENTOTorg,
     >                    RNE, RNEorg, ABXYZ, NFIRST, NLAST, 
     >                    NATOM, WEIGHT, NCHARG, EION, ELEVEL, NOM, 
     >                    TAUROSScont, TAUSCALorg,
     >                    EN, DEPARTNDorg, bUseENTOT)
      ENDIF
      IF (NTOLD > 2) THEN      
        CALL REGRIDOLDPOP(POP3, T3old, SCRATCH, ND, N, T,
     >                    POPMIN, TAUTHICK, bNoRGrid,
     >                    RADIUS, RADIUSorg, ENTOT, ENTOTorg,
     >                    RNE, RNEorg, ABXYZ, NFIRST, NLAST, 
     >                    NATOM, WEIGHT, NCHARG, EION, ELEVEL, NOM, 
     >                    TAUROSScont, TAUSCALorg,
     >                    EN, DEPARTNDorg, bUseENTOT)
      ENDIF

C     ------- end of POPNUM interpolation ----------------------

            
C***  Final calculation of the velocity gradient
      CALL GRADIFF(ND, VELO, GRADI, RADIUS)

C***  Final calculation of the TAU scale (continuum opacities only)
C***  now including all grid updates
      CALL TAUSCAL(RSTAR, ND, RADIUS, RNE, ENTOT,
     >             T, POPNUM, NDIM, N, EN, LEVEL, NCHARG, WEIGHT,
     >             ELEVEL, EION, EINST, ALPHA, SEXPO,
     >             ADDCON1, ADDCON2, ADDCON3, 
     >             IGAUNT, NOM, NF,
     >             XLAMBDA, FWEIGHT,
     >             TAUTHOM, TAUROSScont,
     >             MAXATOM, SIGMATHK, SEXPOK, EDGEK, KODAT,
     >             KONTNUP, KONTLOW, LASTKON, 
     >             DENSCON, FILLFAC)
            
C***  Final calculation of the DENSCON and FILLFAC vectors,
C***  now including all damping and grid update effects
      CALL CLUMP_STRUCT (DENSCON, FILLFAC, ND, DENSCON_FIX, 
     >                   VELO, TAUROSScont, DENSCON_LINE, RADIUS, T, XMU)     
      
C***  Calculate Tau2/3 (always on Rosseland continuum scale)
      TAU23=0.666666666666
      IF (TAUROSScont(ND) < TAU23) THEN
         R23=1.
      ELSE
         CALL LIPO (R23,TAU23,RADIUS,TAUROSScont,ND)
      ENDIF
      

      
C***  Update GRSTATIC vector for the MODEL file:
C***  This vector stores GAMMARAD (=ARAD/AGRAV) in the MODEL and will be 
C***  used in the next call of ENSURETAUMAX to allow for the potential
C***  mixing of the current GAMMARAD with the one from the last run.
C***  (This is a technical requirement for the REDUCE option
C***    on the HYDROSTATIC INTEGRATION card)
      DO L=1, ND
        IF (bFULLHYDROSTAT) THEN
          IF (RADIUS(L) > RADIUSorg(1)) THEN
            GRSTATIC(L) = GAMMARAD(1)
          ELSEIF (RADIUS(L) < RADIUSorg(ND)) THEN
            GRSTATIC(L) = GAMMARAD(ND)
          ELSE
            CALL SPLINPOX(GRSTATIC(L),RADIUS(L),GAMMARAD,RADIUSorg,ND)
          ENDIF
        ELSE
          !If Thomson is used, just fill with mean electron Gamma to
          ! avoid any problems if FULL might be switched on later
          GRSTATIC(L) = GEDD
        ENDIF
      ENDDO      
      
C***  Check if final velocity stratification is monotonic
      monocheck: DO I=ND-1, 1, -1
        IF (VELO(I) <= VELO(I+1)) THEN
          WRITE (hCPR,*) ' WARNING: velocity is not monotonic'
          DO J=1, ND
            WRITE (hCPR,FMT='(A,I2,A,F12.6,A,A)') 
     >       "VELO(",J,")=",VELO(J),'  ',VELOCRITERION(J)
          ENDDO
          EXIT monocheck
        ENDIF
      ENDDO monocheck

C***  Determine and print maximum corrections
C***  (Debug option: Print all corrections via bPrintChanges = .TRUE.)
      TMCORRMAX = 0.
      LTMCM = 0.
      DO L=1, ND
        IF (RADIUS(L) > RADIUSorg(1)) THEN
          VL = VELOorg(1)
          TL = Torg(1)
          ENTOTL = ENTOTorg(1)
        ELSE
          CALL SPLINPOX(VL, RADIUS(L), VELOorg, RADIUSorg, ND)
          CALL SPLINPOX(TL, RADIUS(L), Torg, RADIUSorg, ND)
          CALL SPLINPOX(ENTOTL, RADIUS(L), ENTOTorg, RADIUSorg, ND)
        ENDIF
        IF (ABS(VELO(L)/VL-1.) > TMCORRMAX(1)) THEN
          LTMCM = L
        ENDIF
        TMCORRMAX(1) = MAX(TMCORRMAX(1), ABS(VELO(L)/VL-1.))
        TMCORRMAX(2) = MAX(TMCORRMAX(2), ABS(T(L)/TL-1.))
        TMCORRMAX(3) = MAX(TMCORRMAX(3), ABS(ENTOT(L)/ENTOTL-1.))
        IF (bPrintChanges) THEN
          WRITE(hCPR,'(A,I2,6(3X,G15.8))') 'corr V and T: ', 
     >      L, ABS(VELO(L)/VELOorg(L)-1.), ABS(VELO(L)/VL-1.), 
     >      ABS(T(L)/Torg(L)-1.), ABS(T(L)/TL-1.),
     >      ABS(ENTOT(L)/ENTOTorg(L)-1.), ABS(ENTOT(L)/ENTOTL-1.)
        ENDIF
      ENDDO
      WRITE(hCPR,'(A, I3, 2(F15.8))') 
     >  'Maximum corrections from ENSURETAUMAX at L to V & ENTOT: ', 
     >       LTMCM , TMCORRMAX(1), TMCORRMAX(3)

C***  PRINTOUT OF VARIOUS MODEL SPECIFICATIONS
      OLDTEMP = .FALSE.
      TTABLE = .FALSE.
      BTWOT = .FALSE.
      MODHEAD=' '
      MODOLD = ' ' 
      MODOLD2 = ' ' 
      JOBNOLD = 0
      JOBNOLD2 = 0 
      TFAC = 1.0

      bNoDetails = .TRUE.  !do not print depth-dependend values
      CALL PRIMOD (ND,RADIUS,INCRIT,ENTOT,T,VELO,GRADI,NP,
     >             OLDTEMP, MODHEAD, JOBNUM, 
     >             MODOLD, JOBNOLD, TTABLE, TAUROSScont,
     >             R23, TEFF, bThinWind, ITTAU, MAXITTAU, RCON, 
     >             BTWOT, MODOLD2, JOBNOLD2, TFAC,
     >             BETA, VPAR1cand, VPAR2cand, DENSCON, BETA2, 
     >             BETA2FRACTION, HSCALE, bNoDetails,.FALSE.,VTURB)

      GEDDRAD = GAMMARADMEAN        !store new GAMMARADMEAN in GEDDRAD (saved in MODEL)
      bTauUpdated = .TRUE.
      LASTTAU = 0.

      VMIN = VELO(ND)   !ensure that new VMIN value is really the lowest velocity value
      VMINCAND = VMIN      

      RETURN
      END
