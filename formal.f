C***  MAIN PROGRAM FORMAL  ****************************************************
      SUBROUTINE FORMAL
 
C*******************************************************************************
C***  FORMAL INTEGRAL IN THE OBSERVERS FRAME, YIELDING EMERGENT FLUX PROFILES
C*******************************************************************************
 
      IMPLICIT NONE
      character onechar*1

C***  IRON: ADD GENERIC ION TO MAXATOM
      INTEGER, PARAMETER :: MAXATOM =          26 

      INTEGER, PARAMETER :: MAXLAP  =           15000
      INTEGER, PARAMETER :: MAXSUBL =            5500
      INTEGER, PARAMETER :: MAXIND  = 22000 + MAXSUBL
      INTEGER, PARAMETER :: MAXFEIND  =       1500 
C***  NDIM = Number-of-original-levels + MAXSUBLevels (without 2*)
      INTEGER, PARAMETER :: NDIM    =  1560 + MAXSUBL
      INTEGER, PARAMETER :: NFLDIM  =          300000   
      INTEGER, PARAMETER :: NFODIM  =          250000   
      INTEGER, PARAMETER :: MAXMOD  =               2
      INTEGER, PARAMETER :: NFDIM   =    2*NDIM + 400 
      INTEGER, PARAMETER :: MAXKONT =         NFDIM/2 
      INTEGER, PARAMETER :: MAXKODR =            NDIM 
      INTEGER, PARAMETER :: NDDIM   =             100 
      INTEGER, PARAMETER :: NDADDIM =    2*NDDIM + 10 
      INTEGER, PARAMETER :: NPDIM   =             130 
      INTEGER, PARAMETER :: MAXXN   =            4000
      INTEGER, PARAMETER :: NREDMAX =            6000
      INTEGER, PARAMETER :: MAXSTRI =              20 
      INTEGER, PARAMETER :: MAXXDAT =              10 
      INTEGER, PARAMETER :: NPHIMAX =             500 

C***  MAXIMUM ION CHARGE WHICH MAY OCCUR (SEE ALSO SUBR. GAUNTFF)
      INTEGER, PARAMETER :: MAXION = 27 
      
C***  IRON: COMMON BLOCK FOR IRON-SPECIFIC DATA
C***  include "dimblock"
      INTEGER, PARAMETER :: INDEXMAX = 1E7, NFEREADMAX = 3E5    !std
C      INTEGER, PARAMETER :: INDEXMAX = 4E7, NFEREADMAX = 5E5    !vd20
C      INTEGER, PARAMETER :: INDEXMAX = 1E8, NFEREADMAX = 6E5    !xxl
      
C***  ARRAYS FOR TREATMENT OF LINE OVERLAPS:
      REAL, DIMENSION(MAXLAP) :: XLAMLAP, DELXLAP 
      INTEGER, DIMENSION(MAXLAP) :: INDLAP, IPOINTERPHITAB
      REAL, DIMENSION(NDDIM,MAXLAP,MAXMOD) :: OPAL, ETAL 
      REAL, DIMENSION(NDDIM,MAXLAP) :: OPALRAY, ETALRAY
      REAL, DIMENSION(MAXXN,MAXLAP) :: OPALFIN, ETALFIN
      REAL, DIMENSION(MAXLAP,NDDIM,MAXMOD) :: AVOIGT, GRIEMPAR

C***  ARRAYS FOR MULTIPLET HANDLING:
      INTEGER, DIMENSION(MAXSUBL) :: NSUBLOW, NSUBNUP

C***  HANDLING OF DIELECTRONIC RECOMBINATION / AUTOIONIZATION (SUBR. DATOM)
      INTEGER, PARAMETER :: MAXAUTO = 2850
      INTEGER, DIMENSION(MAXAUTO) :: LOWAUTO, IONAUTO, KRUDAUT
      REAL, DIMENSION(MAXAUTO) :: WAUTO, EAUTO, AAUTO
      INTEGER, DIMENSION(MAXKODR) :: KODRNUP, KODRLOW

C***  VECTORS FOR USE IN SUBR. ZONEINT
      REAL, DIMENSION(MAXXN) :: TAU, DTAU, WTAU, TAUC, DTAUC, WTAUC,
     >                          XCMFFINE, POROLENGTHFINE

      REAL, DIMENSION(NDIM) :: WEIGHT, ELEVEL, EION, ENLTE
      INTEGER, DIMENSION(NDIM) :: NCHARG, MAINQN, NOM, IONGRND
      REAL, DIMENSION(NDIM, NDIM) :: EINST
      REAL, DIMENSION(4, NDIM) :: ALTESUM
      REAL, DIMENSION(4, MAXIND) :: COCO
      REAL, DIMENSION(MAXATOM) :: ATMASS, STAGE
      INTEGER, DIMENSION(MAXATOM) :: KODAT, NFIRST, NLAST, KODATIND
      REAL, DIMENSION(NFDIM, MAXMOD) :: XLAMBDA
      REAL, DIMENSION(NDDIM, MAXMOD) :: OPA, ETA
      REAL, DIMENSION(NDDIM) :: ETANOTH, THOMSON, DELW, ADELW,
     >                          XJCIND, OPAFEFT, ETAFEFT
      INTEGER, DIMENSION(NDDIM) :: IWARN
      REAL, DIMENSION(NDDIM):: VEC_SECMOD
      REAL, DIMENSION(NDDIM,MAXMOD) :: RADIUS, ENTOT, T, RNE,
     >                                  VELO, VDU, GRADI, TAURCONT
      REAL, DIMENSION(NDDIM) :: RADIUS_MERGED 
      REAL, DIMENSION(NDDIM,MAXMOD) :: VDU_ORIG, 
     >               T_ORIG, RNE_ORIG, ENTOT_ORIG
      REAL, DIMENSION(MAXMOD) :: RCON, XMDOT
      REAL, DIMENSION(MAXATOM,MAXMOD) :: ABXYZ
      REAL, DIMENSION(NDDIM,NFDIM,MAXMOD) :: XJC
      REAL, DIMENSION(NDDIM,NPDIM) :: U, UK
      REAL, DIMENSION(3,NDDIM) :: EDDI ! For dimensioning, used in FORMCMF
      REAL, DIMENSION(NDADDIM) :: OPARAY, ETARAY, ETACRAY, 
     >                                    ZRAY, XCMF, RRAY
      REAL, DIMENSION(MAXXN) :: OPAFINE, ETAFINE, OPAFC, ETAFC,
     >                          SFINE, CSFINE, ZFINE
      REAL, DIMENSION(NFODIM) :: PROFILE, CPROFILE, DLAM, PROFILEC
      REAL, DIMENSION(MAXION) :: TFEEXC
      REAL, DIMENSION(NPDIM) :: A,WE, WX
      REAL, DIMENSION(NPDIM,NPDIM) :: BX
      LOGICAL, DIMENSION(NDIM,NDDIM,MAXMOD) :: ZERO_RATES

C***  ATTENTION: B AND C MUST BE LOCATED SUBSEQUENTLY IN THE MEMORY]
C***  @TODO: GET RID OF THIS CONDITION ASAP!!!
      REAL, DIMENSION(NPDIM,NPDIM) :: B
      REAL, DIMENSION(NPDIM) :: C

      
      REAL, DIMENSION(NPDIM,MAXMOD) :: PGRID
      REAL, DIMENSION(NPDIM) :: PGRID_ORIG, PGRID_MERGED, EMINT_P
      INTEGER NP_ORIG, NDNP_ORIG
      REAL, DIMENSION(NDDIM,NPDIM,MAXMOD) :: ZGRID, ZGRID_ORIG
      REAL, DIMENSION(NDDIM,NPDIM) :: ZGRID_MERGED

      REAL, DIMENSION(NDDIM,NDIM,MAXMOD) :: POPNUM, POPNUM_ORIG, 
     >                                       POPLTE
      REAL, DIMENSION(MAXKONT) :: ALPHA, SEXPO, 
     >                            ADDCON1, ADDCON2, ADDCON3
      INTEGER, DIMENSION(MAXKONT) :: KONTNUP, KONTLOW 
      CHARACTER*8 IGAUNT(MAXKONT), KEYCBF(MAXKONT)
      INTEGER, DIMENSION(MAXIND) :: INDNUP, INDLOW, MULTIIND


C***  ARRAYS FOR CMF (ELECTRON SCATTERING REDISTRIBUTION OPTION)
      REAL, DIMENSION(NDDIM) :: TA, TB, TC, UB, GA, H, QQ, SCMF,
     >                          VA, VB, PP, W0
      REAL, DIMENSION(NDDIM,NPDIM) :: V
      REAL, DIMENSION(NDDIM,NFLDIM) :: XJNUE, ETANCK, THOMCK,
     >                                 OPAK, ETAK
      REAL, DIMENSION(NDDIM,NFLDIM,MAXMOD) :: ETACK, OPACK, ETACCK
      REAL, DIMENSION(NFLDIM,MAXMOD) :: BCORE, DBDR
      REAL, DIMENSION(NFLDIM) :: XPLOT, YPLOT, YSCRATCH

C***  ARRAYS TO HANDLE DEPTH-DEPENDENT REDISTRIBUTION INTEGRAL
      REAL, DIMENSION(NDDIM) :: VDUEL, WREDI0
      REAL, DIMENSION(NREDMAX,NDDIM) :: WREDI

C***  X-RAY DATA
      REAL, DIMENSION(MAXXDAT) :: XDATA
      REAL, DIMENSION(MAXATOM,MAXATOM) :: SIGMATHK, SEXPOK, EDGEK

C***  Density Contrast
      REAL, DIMENSION(NDDIM,MAXMOD) :: DENSCON, FILLFAC 
      REAL, DIMENSION(NDDIM) :: ENTOTDENS, POROLENGTH
      REAL, DIMENSION(NDADDIM) :: POROLENGTHRAY

C***  Line profiles for pressure broadening
C***  NLDIMPHITAB = max. number of H + He lines in one spectral range 
C***  NFDIMPHITAB = max. number of frequency points in one half-profile 
      INTEGER, PARAMETER :: NLDIMPHITAB = 100
      INTEGER, PARAMETER :: NFDIMPHITAB = 2000
      REAL, DIMENSION 
     > (-NFDIMPHITAB:NFDIMPHITAB, NDDIM, NLDIMPHITAB,MAXMOD) :: PHITAB
      REAL, DIMENSION(2*NFDIMPHITAB+1) :: PHISCRATCH
      REAL, DIMENSION(MAXLAP) :: XMAXLIN
      CHARACTER(132), DIMENSION(0:MAXSTRI) :: STRING1
      CHARACTER(100), DIMENSION(MAXMOD) :: MODHEAD
      CHARACTER(100) :: MODHEADT
      CHARACTER(10), DIMENSION(NDDIM) :: MAINPRO, MAINLEV
      CHARACTER(10), DIMENSION(NDIM) :: LEVEL
      CHARACTER(10), DIMENSION(MAXATOM) :: ELEMENT
      CHARACTER(4), DIMENSION(MAXIND) :: KEYCBB
      CHARACTER(2) :: EL_MIN
      CHARACTER(2), DIMENSION(MAXATOM) :: SYMBOL
      CHARACTER(8), DIMENSION(MAXLAP) :: LINPRO
      CHARACTER(10) :: XUNIT
      CHARACTER(20) :: FREQIN 
      CHARACTER(80) :: MANIPOP_OPTIONS, DD_VDOP_LINE, NOWIND_LINE, 
     >                 VDOPPLOT_LINE, TRANSDWLLINE, RCOROTLINE, 
     >                 MACROCLUMPLINE
      CHARACTER(256) :: KARTE, PATH_VCSSB, PATH_LEMKE_DAT

      LOGICAL :: PLOT, FIN, REDIS, IDENT
      LOGICAL :: BROAD, bBIGBANDLIMIT
      LOGICAL :: ABSWAV, BCONT, BNOCONT
      LOGICAL :: CORE, LINELIST, BCALIBRATED  
      LOGICAL :: BDD_VDOP, BPLOTVDOP, BPLOTVMIC, BMICROTURB
      INTEGER :: INDCUT, ND_MERGED, NP_MERGED

C     abarnisk 04062004: BCALIBRATED ob profil kalibriert werden soll
      CHARACTER(8) :: FUNIT
C***  IND_ORIGLEV(LEVEL) = index of "mother level" of multiplets 
C***                       (identical for normal lines)      
      INTEGER, DIMENSION(NDIM) :: IND_ORIGLEV

C***  For rotation    
      INTEGER, DIMENSION(NPDIM) :: NPHI
      REAL, DIMENSION(NPHIMAX,NPDIM) :: PHIWEIGHT, PHIARR    
 
C***  Variables needed for SECONDMODEL_DEFINE
      CHARACTER SECONDMODEL_LINE*1000, SECONDMODEL_PATH*400
      LOGICAL SECONDMODEL_CHANGED, IGNORE_DIFF
      REAL, DIMENSION(2) :: SECMOD_RRANGE

      REAL, DIMENSION(NPHIMAX) :: PHI_VEC

C***  VARIABLES READ BY FORMOSA
      INTEGER, DIMENSION(MAXMOD) :: ND, NP, NF, JOBNUM
      REAL,    DIMENSION(MAXMOD) :: RSTAR, VDOP_MODEL, TEFF

C***  FOR depth-dependent VDOP:
      REAL :: VMICFRAC_DEFAULT, VDOP, VDOP_FIRSTMOD

C*** vectors in km/s, Doppler units, and squared Doppler units are prepared
C     to avoid recalculation in long loops.
      REAL, DIMENSION(NDDIM,MAXMOD) :: 
     >                       DD_VMIC, DD_VMICDU, DD_VMICDU_ORIG
      REAL, DIMENSION (NDDIM, MAXATOM, MAXMOD) :: 
     >     DD_VDOP, DD_VDOPDU, DD_VDOP_ORIG, DD_VDOPDU_ORIG
      REAL, DIMENSION (MAXXN, MAXATOM) :: DD_VDOPDU_FINE_NORMFAC,
     >                  DD_VDOPDU_FINE_SQRD, DD_VDOPDU_FINE
      REAL, DIMENSION(NDADDIM, MAXATOM) :: DD_VDOPDU_RAY

C*** Variables for limb-darkening plot
      REAL, DIMENSION(NFODIM) :: WEIGHT_LIMBDARK
      CHARACTER(132) :: LIMB_LINE
      LOGICAL :: BLIMB

      REAL :: DXLAM, DX_WINDROT, WEIGHTSUM

      REAL, DIMENSION(NPDIM) :: SUMFILT, SUMINT
      
      INTEGER, DIMENSION(MAXLAP) :: IND_ELLINE
      
C***  FROM PREPRAY -> OBSFRAM
      REAL BCOREL, DBDRL

C***  Intersection points with secondmodel domain
      REAL, DIMENSION(2, NPDIM, NPHIMAX) :: ZINTER

C***  To Store Feautrier Matrices ST(94*95,89)
      REAL, DIMENSION((NPDIM+1)*NPDIM,NDDIM) :: ST
      LOGICAL :: BELIFI

C***  SPZ is for the Input Option SET_POP_ZERO
C***  - can be sed to suppress specific lines 
      INTEGER, PARAMETER :: MAXSPZ = 10
      CHARACTER(10), DIMENSION(MAXSPZ) :: SPZ1, SPZ2

C***  Iron data variables      
      REAL, DIMENSION(NFEREADMAX) :: FEDUMMY
      INTEGER, DIMENSION(MAXFEIND) :: INDRB, INDRF, IFRBSTA, IFRBEND,
     >                                IFELOW, IFENUP, INDFEACT
      REAL, DIMENSION(MAXFEIND) :: SIGMAACT, SIGMAACTUL, SIGNU3ACT,
     >                             SIGMAINT, SIGMAINTUL
      REAL, DIMENSION(INDEXMAX) :: SIGMAFE
ccc     SIGMAFEUL, SIGNU3FE unused - from Andreas
ccc      REAL, DIMENSION(INDEXMAX) :: SIGMAFE, SIGMAFEUL, SIGNU3FE
      INTEGER, DIMENSION(NDIM) :: INDRBS, IFRBSSTA, IFRBSEND, 
     >                            INDFESACT, IFES
      REAL, DIMENSION(NDDIM,NFLDIM,MAXMOD) :: OPAFE, ETAFE
      REAL, DIMENSION(NDDIM,MAXFEIND) :: OPAFEI, ETAFEI
      LOGICAL, DIMENSION(NDDIM) :: bFELASER
      LOGICAL :: BFECHECK, BFEMODEL, BIRONLINES, BFEWARNING,
     >           BAIRWAVELENGTHSET, BAIRWAVELENGTH, 
     >           bDDOPAFORMCMF, bDDFECONVOL, BVSINI_AT_RCOROT

      REAL :: VSINI, TAUMAX, XMAXMIN, DXMAX, VDOPFE,
     >        DXFE, XLAM0FE, VMAX, OSMIN, RMAX, XMAX, YMAX,
     >        XMAXBROAD, VSIDU, RANGEBLUE, VMAXDU,
     >        RANGERED, RCOROT, XLAM, XLAMLN,
     >        DISP, BWESEX, ALN, XRANGERED, XRANGEBLUE, XLAP,
     >        FREMAX, FREMIN, DUMMY, XCMFRED, XCMFBLUE, DXCMF, FNUEC,
     >        XLAMREF, DXOBS, FINCRI, BEGLAM, XOBS0, XO, PJPJ, 
     >        PWEIGHT, XLAM2, DISMO, FACLOG, FAC, CEMINT, EMINT, 
     >        POPMIN, ATMEAN, XMUL, AMACHL, TAUMINBROAD
      INTEGER LSOPA, LSDWL, LSPRO, IFIRST, IERR, IDUMMY,
     >        LPHISTA, LPHIEND, IMOD, JPFIRST, JPLAST, KANAL1, LINE,
     >        LASTIND, LASTKDR, LASTFE, KWORDS, K,
     >        L, N, I, J, NMOD, NATOM, NAUTO, LASTKON, NA, NDN, NFL,
     >        NSTRING, IDTERM, NLPHITAB, IVERSION, NPHIROT_DEFAULT,
     >        NBLINE, NMULTI, NBL, IND, LBREF, IRANGEBLUE, IRANGERED,
     >        LOW, NUP, INDREF, IWARNJ0, JP, NFOBS, IFOBR,
     >        LPHI, IRAY, LTOT, NSTRI, LASTSELF, IVDOPSTATUS,
     >        JPFIRST_ORIG, JPLAST_ORIG, LPHISTA_ORIG, LPHIEND_ORIG, 
     >        NBFIRST, NBLAST 

      INTEGER, EXTERNAL :: IDX, ISRCHFLE, ISRCHFLT, ISRCHFGT, ISRCHFGE
     
C***  Operating system:
      CHARACTER(8) :: OPSYS
      COMMON / COMOS / OPSYS

      CHARACTER(10) :: TIM1, TIM2

C***  Physical constants      
      REAL, PARAMETER :: CLIGHT = 2.99792458E5      ! IN KM/SECOND
      REAL, PARAMETER :: CLIGHT2 = 2.99792458E18    ! in Angstroem/sec
      REAL, PARAMETER :: PARSEC = 3.08561E18    !PARSEC IN CM
      REAL, PARAMETER :: CONMV = 48.64 !Calibration CONSTANT FOR F-NUE IN MV
      REAL, PARAMETER :: PI = 3.14159 
      REAL, PARAMETER :: WPIINV = 0.564189583549       

      !File and channel handles (=KANAL)
      INTEGER, PARAMETER :: hOUT   = 6    !write to wruniqX.out (stdout)
      INTEGER, PARAMETER :: hCPR   = 0    !write to wruniqX.cpr (stderr)

C***  Initialize values (only DATA statements here!)
      DATA LSOPA,LSDWL,LSPRO,IFIRST /-1, -1, -1, 1/
      DATA VSINI / 0. /
      DATA BVSINI_AT_RCOROT / .TRUE. /
C***  X-UNITS OF PLOT IN ANGSTROEM (ALTERNATIVELY: MICROMETER)
      DATA XUNIT / 'ANGSTROEM ' /

C***  Link data to identify program version
      CHARACTER(30) :: LINK_DATE
      CHARACTER(10) :: LINK_USER
      CHARACTER(60) :: LINK_HOST
      COMMON / COM_LINKINFO / LINK_DATE, LINK_USER, LINK_HOST

C***  For progress-bar
      INTEGER NRAYDONE, NPHISUM, IMPATIENCE
C********* END OF DECLARATIONS *****************************************

C***  Write Link Data (Program Version) to CPR file
      WRITE (hCPR,'(2A)') '>>> FORMAL started: Program Version from '
     >                 ,LINK_DATE
      WRITE (hCPR,'(4A)') '>>> created by '
     >                 , LINK_USER(:IDX(LINK_USER))
     >     ,' at host ', LINK_HOST(:IDX(LINK_HOST))

      CALL INSTALL
ccc not in gfortran      CALL TIME(TIM1)

C***  TAUMAX = MAXIMUM OPTICAL DEPTH WHERE INTEGRATION IS TRUNCATED
C***  XMAX   = COMOVING-FRAME BANDWIDTH OF THE SCATTERING ZONE
C***  DXMAX  = MAXIMUM CMF FREQUENCY STEP IN THE SCATTERING ZONE
C***  TAUMINBROAD = MININUM TAU VALUE REQUIRED FOR EXTENDING THE 
C***                BROADENING CALCULATION RANGE IN SUBR. BANDWIDTH
C***  DEFAULT VALUES - MAY BE CHANGED BY INPUT OPTIONS!
      TAUMAX = 7.0
      XMAXMIN = 3.5
      DXMAX  = 0.3
      TAUMINBROAD = 0.01       !was 0.1 before Jan 2017 
      MANIPOP_OPTIONS = ' '      
      NMOD = 1

C***  Initialize "STRING COMMENT"
      FREQIN = ' '

C***  Initialize NOWIND option
      NOWIND_LINE = 'NONE'

C***  Initialize LIMBDARKENIG option
      LIMB_LINE = 'OFF'

C***  Initialize MACROCLUMP option
      MACROCLUMPLINE = 'NONE'

C***  Initialize porosity length
      DO L=1, NDDIM
         POROLENGTH(L) = .0
      ENDDO

C***  Initialize BELIFI; To Store Feautrier Matrices in Memory is now default
      BELIFI = .FALSE.

C***  OPEN CARDS-FILE
      OPEN (UNIT=2, FILE='CARDS', STATUS='UNKNOWN')

      CALL DATOM (NDIM, N, LEVEL, NCHARG, WEIGHT, 
     >            ELEVEL, EION, MAINQN,
     >            EINST, ALPHA, SEXPO,
     >            ADDCON1, ADDCON2, ADDCON3, 
     >            IGAUNT, COCO, KEYCBB, ALTESUM,
     >            INDNUP, INDLOW, LASTIND, MAXIND, MAXATOM, NATOM,
     >            ELEMENT, SYMBOL, NOM, KODAT, ATMASS, STAGE,
     >            SIGMATHK, SEXPOK, EDGEK, NFIRST,
     >            NLAST, NAUTO, MAXAUTO, LOWAUTO, WAUTO, EAUTO, AAUTO,
     >            IONAUTO, KRUDAUT, KONTNUP, KONTLOW, LASTKON, MAXKONT,
     >            IONGRND, KODRNUP, KODRLOW, LASTKDR, MAXKODR, KEYCBF,
C***  IRON: ADDITIONAL PARAMETERS FOR IRON-GROUP LINE BLANKETING
     >            'FORMAL', INDEXMAX, NFEREADMAX, MAXFEIND,
     >             LASTFE, SIGMAFE, INDRB, INDRF,
     >             IFENUP, IFELOW, IFRBSTA, IFRBEND, FEDUMMY,
     >             VDOPFE, DXFE, XLAM0FE, SIGMAINT, BFEMODEL)
      BIRONLINES = BFEMODEL    ! Default as found in DATOM
    
C***  KODATIND gives the Core Charge Number of an element
C***  This is needed to identify H and He for pressure broadening
      DO NA=1, NATOM
         KODATIND(NA) = 0
         DO J = 1, MAXATOM
            IF (NA .EQ. KODAT(J)) KODATIND(NA) = J
         ENDDO
         IF (KODATIND(NA) .EQ. 0) THEN
            WRITE (0,*) '*** ERROR: ELEMENT NOT FOUND'
            STOP 'ERROR IN FORMAL'
         ENDIF
         IF (KODATIND(NA) .GT. MAXATOM) THEN
            WRITE (0,*) '*** ERROR: NCORECHARGE NOT FOUND'
            STOP 'ERROR IN FORMAL'
         ENDIF
      ENDDO

C***  READING OF THE MODEL FILE
      IMOD=1
      CALL FORMOSA (ND(IMOD), RADIUS(1,IMOD), NP(IMOD), 
     >          PGRID(1,IMOD), ZGRID(1,1,IMOD),
     >          ENTOT(1,IMOD), RNE(1,IMOD), ABXYZ(1,IMOD), NATOM, 
     >          T(1,IMOD), VELO(1,IMOD), NF(IMOD),
     >          XLAMBDA(1,IMOD), GRADI(1,IMOD),
     >          POPNUM_ORIG(1,1,IMOD), RSTAR(IMOD), VDOP_MODEL(IMOD), 
     >          JOBNUM(IMOD), N, NDDIM, NPDIM, NFDIM, 
     >          MODHEAD(IMOD), TEFF(IMOD),
     >          MAXXDAT, XDATA, XJC(1,1,IMOD), IMOD, 
     >          DENSCON(1,IMOD), FILLFAC(1,IMOD), TAURCONT(1,IMOD),
     >          ZERO_RATES(1,1,IMOD), RCON(IMOD), NDIM, XMDOT(IMOD) )

      VMAX = VELO(1,1)
      RMAX = RADIUS(1,1)
C***  Default (if not specified otherwise in FORMAL_CARDS):
      VDOP = VDOP_MODEL(1)

      NDN = ND(IMOD) * N
      DO I=1, NDN
         POPNUM(I,1,IMOD) = POPNUM_ORIG(I,1,IMOD)
      ENDDO

      CALL POPMIN_NULLING (ZERO_RATES(1,1,IMOD), POPNUM(1,1,IMOD), 
     >                     ND(IMOD), N)

C***  The no. of core-intersecting impact parameters might be modified 
C***  in subr. ROTATION_PREP; 

C***  The radius-grid might be changed by subr. SECONDMODEL_PREP
C***  as well as bu subr. NOWIND
C***  therefore, P and Z are cloned and restored after one RANGE 
C***  has been completed 
      NP_ORIG = NP(1) 
      DO JP=1, NP_ORIG
         PGRID_ORIG(JP) = PGRID(JP,1)
      ENDDO

      ZGRID_ORIG = ZGRID

C***  Printout of Model Parameters

      WRITE (*,'(2A)') 'MODHEAD=', MODHEAD(IMOD)
      CALL PRI_PAR (TEFF(IMOD), RSTAR(IMOD), 
     >         VELO(1,IMOD), DENSCON(1,IMOD), XMDOT(IMOD))

C***  Open file for plots of the wind-rotation geometric grid
      OPEN (65, FILE='windgrid.plot', STATUS='UNKNOWN')

C***  OPEN FILE 1 = 'PLOT' FOR DIRECT PLOT TRANSFER
      KANAL1 = 1
      CALL JSYMSET ('G1','TRANSFER')
      CALL REMARK ('PLOT DATA TO BE ROUTED')
      OPEN (KANAL1, FILE='PLOT', STATUS='UNKNOWN')

C***  MASS STORAGE ON FILE 7 FOR THE FEAUTRIER MATRICES (CONT.)
      CALL OPENMS (7, IDUMMY, IDUMMY, 0, IERR)
 
C***  DEFAULT OPTIONS
      BROAD = .FALSE.
      IDENT = .TRUE.
      OSMIN = 0.05
      REDIS = .TRUE.
      BCONT = .FALSE.
      BWESEX = 1.
      FIN=.FALSE.
      DISP= .0
      BDD_VDOP = .FALSE.
      RCOROTLINE = 'DEFAULT'
      VMICFRAC_DEFAULT = 0.05
      IVERSION=0
      ABSWAV=.TRUE.
      LINELIST = .FALSE.
      BNOCONT = .FALSE.
      BCALIBRATED = .FALSE.
      FUNIT = 'FLAM10PC'
      BAIRWAVELENGTHSET=.FALSE.
      BAIRWAVELENGTH=.FALSE.
      BPLOTVDOP = .FALSE.
      BMICROTURB = .FALSE.
      BLIMB = .FALSE.
      DD_VDOP_LINE = ''
      VDOPPLOT_LINE = ''
      PATH_VCSSB     = 'default'
      PATH_LEMKE_DAT = 'default'
      bBIGBANDLIMIT = .TRUE.
      bDDOPAFORMCMF = .TRUE.
      bDDFECONVOL = .TRUE.
      IVDOPSTATUS = 1
      IGNORE_DIFF = .FALSE. 

C***  LOOP FOR EVERY SPECTRAL RANGE TO BE SYNTHESIZED    ---------------------
    1 CONTINUE

C***  DECFORM reads from FORMAL_CARDS all options till a new
C***  spectral range is opened with a LINE, BLEND or RANGE option
      CALL DECFORM (KARTE, LSOPA, RANGERED, RANGEBLUE,
     >              VSINI, DISP, REDIS, BWESEX,
     >              LSPRO, LSDWL, FIN, VDOP, IDTERM, 
     >              JPFIRST_ORIG, JPLAST_ORIG,
     >              IVERSION, TAUMAX, XMAXMIN, DXMAX, PATH_VCSSB,
     >              IDENT, OSMIN, BROAD, MODHEAD(1), JOBNUM(1),
     >              STRING1, MAXSTRI, NSTRING, ABSWAV, BCONT, FUNIT,
     >              LINELIST, 
     >              BIRONLINES, BNOCONT, SPZ1, SPZ2, MANIPOP_OPTIONS, 
     >              XUNIT, BCALIBRATED, FREQIN, MACROCLUMPLINE, 
     >              PATH_LEMKE_DAT, BAIRWAVELENGTHSET, BAIRWAVELENGTH, 
     >              MAXSPZ, TRANSDWLLINE, RCOROTLINE, 
     >              DD_VDOP_LINE, BPLOTVDOP, BMICROTURB,  
     >              LIMB_LINE, VDOPPLOT_LINE, bBIGBANDLIMIT, 
     >              NOWIND_LINE, TAUMINBROAD, 
     >              bDDOPAFORMCMF, bDDFECONVOL, 
     >              LPHISTA_ORIG, LPHIEND_ORIG,
     >              BVSINI_AT_RCOROT, DX_WINDROT, SECONDMODEL_LINE)

C***  If no more ranges requested, terminate the program
      IF (FIN) GOTO 20

      XMAX = XMAXMIN   
      
      IF (.NOT. BFEMODEL) BIRONLINES = .FALSE.

      WRITE (0,'(/,A)') '----------- Starting with range '//FREQIN

C***  If JPFIRST and/or JPLAST were specified: 
C***     make sure that theu fall into the allowed range
      JPFIRST = MAX(1,       JPFIRST_ORIG)
      JPFIRST = MIN(JPFIRST, NP(1)-1)
      JPLAST  = MIN(NP(1)-1, JPLAST_ORIG)
      JPLAST  = MAX(JPLAST, JPFIRST)

      LPHISTA = 1
      LPHIEND = 1


C***  If iron lines are disabled by the option NO-IRONLINES :
C***  WARNING issued here to the cpr-File, and by Subr. PRIPRO to formal.out 
      IF (.NOT. BIRONLINES  .AND. BFEMODEL) then
        BFEWARNING = .TRUE.
        WRITE(0,'(A)') 'WARNING: MODEL CONTAINS IRON LINES THAT ARE '
     >              // 'SUPPRESSED BY THE OPTION "NO-IRON LINES"'
      ELSE 
        BFEWARNING = .FALSE.
      ENDIF    

C***  Preparation of depth-dependent ("DD") VDOP. 
C***  DD is activatd by specifying VMIC 
C***  If DD is not activated: 
C***     VDOP is constant, taken from MODEL file, or
C***     overwritten by FORMAL_CARDS option VDOP 
C***  If DD is activated: VDOP refers to the smallest Doppler-broadening
C***     of any line at any depth for ensuring sufficient resolution

C***  VDOP_STRUCT must be REPEATED FOR THE SECOND MODE (different T)
      IMOD = 1
      CALL VDOP_STRUCT (BDD_VDOP, DD_VDOP_LINE, 
     >   DD_VDOP(1,1,IMOD), VDOP, VELO(1,IMOD), 
     >   T(1,IMOD), ND(IMOD), NDDIM, NATOM, MAXATOM,  
     >   DD_VDOPDU(1,1,IMOD), 
     >   VMICFRAC_DEFAULT, ATMASS, XMAX, XMAXMIN,
     >   SYMBOL, VDOPFE, BMICROTURB, BIRONLINES,
     >   DD_VMIC(1,IMOD), TAURCONT, RADIUS(1,IMOD), EL_MIN,  
     >   DD_VMICDU(1,IMOD), bDDFECONVOL, IVDOPSTATUS)

C***  Note for Multi-Model mode: all dimensionless frequencies
C***  refer to VDOP of model 1
      VSIDU  = VSINI / VDOP
      VMAXDU = VMAX / VDOP

C*******************************************************************
C***  Preparation in case of SECONDMODEL
C*******************************************************************
      IF (SECONDMODEL_LINE .NE. '')
     >  CALL SECONDMODEL_DEFINE (SECONDMODEL_LINE, NMOD, 
     >        SECONDMODEL_PATH, SECONDMODEL_CHANGED, 
     >        RMAX, SECMOD_RRANGE, IGNORE_DIFF)

C***  Error stop if DWL-Plot requested but multi-model active 
      IF (LSDWL .GT. 0 .AND. NMOD .GT. 1) THEN
        WRITE (0,*) 'DWL-PLOT NOT POSSIBLE FOR MORE THAN ONE MODEL'
        STOP 'ERROR detected in main program FORMAL'
      ENDIF

C***  Macroclumping 
C***  Note: in case of SECOND MODEL, the same POROLENGTH vector will be
C***        used there
      IF (MACROCLUMPLINE(:10) .EQ. 'MACROCLUMP') THEN
         CALL PREPMACROCLUMP (MACROCLUMPLINE, DENSCON, VELO, RADIUS,
     >            TAURCONT, ND, POROLENGTH)
      ENDIF

      IF (SECONDMODEL_CHANGED .AND. NMOD .EQ. 2) THEN
        
        CALL COPY_SECONDMODEL 
     >        (SECONDMODEL_PATH, IGNORE_DIFF, BIRONLINES)
        WRITE (0,*) 'SECOND MODEL COPIED'
        IMOD=2
        CALL FORMOSA (ND(IMOD), RADIUS(1,IMOD), NP(IMOD), 
     >          PGRID(1,IMOD), ZGRID(1,1,IMOD),
     >          ENTOT(1,IMOD), RNE(1,IMOD), ABXYZ(1,IMOD), NATOM, 
     >          T(1,IMOD), VELO(1,IMOD), NF(IMOD),
     >          XLAMBDA(1,IMOD), GRADI(1,IMOD),
     >          POPNUM_ORIG(1,1,IMOD), RSTAR(IMOD), VDOP_MODEL(IMOD), 
     >          JOBNUM(IMOD), N, NDDIM, NPDIM, NFDIM, 
     >          MODHEAD(IMOD), TEFF(IMOD),
     >          MAXXDAT, XDATA, XJC(1,1,IMOD), IMOD,
     >          DENSCON(1,IMOD), FILLFAC(1,IMOD), TAURCONT(1,IMOD),
     >          ZERO_RATES(1,1,IMOD), RCON(IMOD), NDIM, XMDOT(IMOD) )

        VMAX = MAX (VELO(1,1),VELO(1,2))

        NDN = ND(IMOD) * N
        DO I=1, NDN
           POPNUM(I,1,IMOD) = POPNUM_ORIG(I,1,IMOD)
        ENDDO

        CALL POPMIN_NULLING (ZERO_RATES(1,1,IMOD), POPNUM(1,1,IMOD), 
     >                       ND(IMOD), N)

C***    Printout of Model Parameters (second model)
        WRITE (*,'(2A)') 'SECOND MODEL READ FROM ', 
     >      SECONDMODEL_PATH(:IDX(SECONDMODEL_PATH))
        WRITE (*,'(2A)') 'MODHEAD=', MODHEAD(IMOD)
        CALL PRI_PAR (TEFF(IMOD), RSTAR(IMOD), 
     >         VELO(1,IMOD), DENSCON(1,IMOD), XMDOT(IMOD))

C***    Preparation of depth-dependent ("DD") VDOP. 
C***    This is repeated here for the SECOND MODEL (different T(r)!)
        VDOP_FIRSTMOD = VDOP
        CALL VDOP_STRUCT (BDD_VDOP, DD_VDOP_LINE, 
     >   DD_VDOP(1,1,IMOD), VDOP, VELO(1,IMOD), 
     >   T(1,IMOD), ND(IMOD), NDDIM, NATOM, MAXATOM,  
     >   DD_VDOPDU(1,1,IMOD), 
     >   VMICFRAC_DEFAULT, ATMASS, XMAX, XMAXMIN,
     >   SYMBOL, VDOPFE, BMICROTURB, BIRONLINES,
     >   DD_VMIC(1,IMOD), TAURCONT, RADIUS(1,IMOD), EL_MIN, 
     >   DD_VMICDU(1,IMOD), bDDFECONVOL, IVDOPSTATUS)

C***     VDOP should not be overwritten by SECOND MODEL
         VDOP = VDOP_FIRSTMOD 
      ENDIF
C***  End of preparations in case of SECONDMODEL *******************

C***  Manipulation of Popnumbers, in order to simulate the emission of a hot 
C***    component (wrh,  5-Feb-2001)
      IF (MANIPOP_OPTIONS .NE. ' ') THEN
         IF (NMOD .EQ. 1) THEN 
            CALL MANIPOP (ENLTE, WEIGHT, NCHARG, EION, ELEVEL, NOM, 
     >              ABXYZ, NFIRST, NLAST, NATOM, POPNUM, RNE, ENTOT, 
     >              N, ND, MANIPOP_OPTIONS, DENSCON, T, level)
         ELSE
            WRITE (0,*) '*** MANIPOP and SECONDMODEL options ',
     >               'cannot be combined!'
            STOP ' *** ABORT in FORMAL ***'
         ENDIF
      ENDIF

C***  Manipulation of Popnumbers: Zero setting for specified levels
      DO IMOD=1, NMOD
         CALL SET_POP_ZERO (LEVEL, NDIM, N, SPZ1, SPZ2, POPNUM(1,1,IMOD),
     >                 ND(IMOD), NMOD, MAXSPZ)
      ENDDO

C***  INTRODUCING DIMENSIONLESS VELOCITY UNITS
      DO IMOD=1, NMOD
        DO L=1, ND(IMOD)
          VDU (L,IMOD)=VELO(L,IMOD) / VDOP
        ENDDO
      ENDDO

C*******************************************************************
C***  Reading of FORMAL-CARDS is now continued for the spectral range
C***  (current single LINE, or till -BLEND closes the range)
C***
C***  For the SECOND-MODEL mode, the atomic date in the FORMAL_CARDS
C***  are for both models; the DATOM data must be identical anyhow.
C***  Subr. MULTIPLE - MULTISPLI adds additional SUBLEVELS with relative  
C***  LTE population numbers to the POPNUM array; this requires T(L) 
C***  and thus an internal loop over IMOD. 
C***  Moreover PREFORM might encounter a DRTRANSIT option which adds 
C***  auto-ionizing levels, requiring T, ENTOT and RNE and 
C***  also an internal loop IMOD=1, NMOD 
C*******************************************************************

      CALL PREFORM (KARTE, N, ELEVEL, LINE, INDLOW, INDNUP, LASTIND,
     >                CLIGHT, VDOP, INDLAP, XLAMLAP, DELXLAP, ALN, 
     >                XLAM, NBLINE, MAXLAP, MAXIND, MAXATOM, 
     >                LEVEL, WEIGHT, EINST, NDIM, POPNUM,
     >                T, ND, NOM, NCHARG, EION, ENTOT, RNE,
     >                MAXSUBL, NSUBLOW, NSUBNUP, BROAD, 
     >                LINPRO, AVOIGT, NMOD, NDDIM, MAXMOD, DENSCON, 
     >                MAINQN, MULTIIND, 
     >                NMULTI, DD_VDOP, NATOM,
     >                IND_ORIGLEV)

C***  The reference wavelenght XLAM is the wavelength of the first line
C***  specified in the current BLEND range. If the BLEND is void of
C***  lines, XLAM is taken as the middle of the range. 
      IF (NBLINE == 0) THEN
        XLAM = (RANGERED + RANGEBLUE) / 2.
      ENDIF

C***  X-Scale: Frequency in harmonic Doppler units referring to XLAM
      XLAMLN = ALOG (XLAM)

C***  Restrict wavelength range if RANGE option given,
C***  remove lines outside range
      IF (RANGERED .NE. RANGEBLUE) THEN
         XRANGERED  = ALOG (RANGERED  / XLAM) / ALN
         XRANGEBLUE = ALOG (RANGEBLUE / XLAM) / ALN

         IRANGERED = ISRCHFLE(NBLINE,XLAMLAP,1,RANGERED)
         NBLINE = NBLINE+1-IRANGERED
         DO I=1, NBLINE
            INDLAP(I) =  INDLAP(I+IRANGERED-1)
            XLAMLAP(I) = XLAMLAP(I+IRANGERED-1)
            DELXLAP(I) = DELXLAP(I+IRANGERED-1)
            LINPRO(I) =  LINPRO(I+IRANGERED-1)
            DO IMOD=1, NMOD
               DO L=1, ND(IMOD)
                  AVOIGT(I,L,IMOD) =  AVOIGT(I+IRANGERED-1,L,IMOD)
               ENDDO
            ENDDO 
         ENDDO
         IF (RANGEBLUE .GT. XLAMLAP(NBLINE)) THEN 
            IRANGEBLUE = ISRCHFLE(NBLINE,XLAMLAP,1,RANGEBLUE)
            NBLINE = IRANGEBLUE - 1
         ENDIF
         IF (NBLINE .LT. 1) THEN
            WRITE (0,*) '*** WARNING: RANGE CONTAINS NO LINES'
         ENDIF
      ELSE
         XRANGERED  = .0
         XRANGEBLUE = .0
      ENDIF

C***  ONLY LINELIST TO BE GENERATED -- NO CALCULATION!
      IF (LINELIST) THEN
        WRITE (*,*) ' ******** LISTONLY OPTION, NO CALCULATION ! *****'
C*      If RANGE option is used together with LISTONLY,
C*      the definition of FREMIN, FREMAX would be missing
        IF (RANGERED .NE. RANGEBLUE) THEN
           FREMIN = XRANGERED
           FREMAX = XRANGEBLUE
        ENDIF
        GOTO 10
      ENDIF

C***  To each line will be assigned the CMF bandwidth that must be
C***  covered; this will be done in Subr. STARKBROAD -> BANDWIDTH
C***  and stored in the vector XMAXLIN
C***  XMAXBOAD is initialzed here and might be incremented in BANDWIDTH
      XMAXBROAD = XMAX

C***  LOOP OVER DETECTED AND ALL BLENDING LINES  -----------------------
      DO 25 NBL=1,NBLINE
        IND=INDLAP(NBL)
        LOW=INDLOW(IND)
        NUP=INDNUP(IND)

C***    IND_ELLINE(LEVEL) gives corresponding element index = NA.         
        IND_ELLINE(NBL) = NOM(LOW)
        NA = NOM(LOW)
        XLAP=XLAMLAP(NBL)
C***    Initialize bandwidth of all lines to XMAX for the case that
C***    STARKBROAD is not called (no broadening requested)
        XMAXLIN(NBL) = XMAX
        DO IMOD=1, NMOD
          CALL LIOP (EINST(NUP,LOW), WEIGHT(LOW), WEIGHT(NUP), LOW,NUP,
     >               ND(IMOD), XLAP, ENTOT(1,IMOD), POPNUM(1,1,IMOD),
     >               RSTAR(IMOD), OPAL(1,NBL,IMOD), ETAL(1,NBL,IMOD),
     >               VDOP)
        ENDDO


C***    OPTION: PRINT LINE OPACITIES ***************************************
        IF (LSOPA.GT.0) THEN
C***      Not allowed for Multi-Model Model !
C***      Not accounting for clumping!
          IF (NMOD .GT. 1) THEN
            WRITE (0,*) 'ERROR : PRINTING OF LINE OPACITIES NOT ',
     >                  'POSSIBLE,'
            WRITE (0,*) '        IF MORE THAN ONE MODEL IS SPECIFIED'
            STOP 'ERROR IN SUBR. FORMAL'
          ENDIF

C***    CALL COOP AND BACKJC ONLY FOR THE USE OF OPA, ETA, THOMSON, XJCIND
C***    IN THE ROUTINE PRIOPAL
          CALL COOP (XLAP, ND, T(1,1), RNE(1,1),
     >               POPNUM(1,1,1), ENTOT(1,1), RSTAR(1),
     >               OPA(1,1), ETANOTH, THOMSON, IWARN, MAINPRO, 
     >               MAINLEV, NOM, KODAT,
     >               NDIM, N, MAXATOM, LEVEL, NCHARG, WEIGHT, 
     >               ELEVEL, EION, EINST,
     >               ALPHA, SEXPO,
     >               ADDCON1, ADDCON2, ADDCON3,
     >               IGAUNT, SIGMATHK, SEXPOK, EDGEK, 
     >               0, NF(1), DUMMY,
     >               RADIUS, KONTNUP, KONTLOW, LASTKON,XDATA)
          CALL BACKJC (XJC(1,1,1), ND, NF(1), 
     >                 XJCIND, XLAMBDA(1,1), 
     >                 XLAP, RADIUS)
          DO 117 L=1, ND(1)
            ETA(L,1) = ETANOTH(L) + 
     >                    OPA(L,1) * THOMSON(L) * XJCIND(L)
  117     CONTINUE

          CALL PRIOPAL(KARTE, XLAP, ND, 
     >                   OPA(1,1), OPAL(1,NBL,1), 
     >                   ETA(1,1), ETAL(1,NBL,1),
     >                   RADIUS, JOBNUM(1), 
     >                   LSOPA, MODHEAD(1))
        ENDIF 
C***    END OF PRINT-OPACITY-BLOCK *************************************

   25 CONTINUE
C***  END LOOP OVER ALL LINES IN BLEND -------------------------

C***********************************************************************
C***  FORMCMF: Co-moving frame calculation of the frequeny
C***           redistribution by electron scattering
C***********************************************************************

      DO IMOD=1, NMOD
         IF (IMOD .EQ. 2) WRITE (0,'(A)') 
     >      'NEXT: FORMCMF for SECOND MODEL'

         VMAX = AMAX1(VMAX, VELO(1,IMOD))
         CALL FORMCMF (TEFF(IMOD), U, UK, NDDIM, NFLDIM,
     >        ZGRID(1,1,IMOD), OPA(1,IMOD), ETA(1,IMOD), ETANOTH, 
     >        ETACK(1,1,IMOD), ETACCK(1,1,IMOD), 
     >        OPACK(1,1,IMOD), ETANCK,
     >        XCMFRED, XCMFBLUE, DXCMF, XMAX, RMAX, VMAXDU,
     >        THOMSON, THOMCK, OPAL(1,1,IMOD),
     >        ETAL(1,1,IMOD),  
     >        RADIUS(1,IMOD), ND(IMOD), NP(IMOD), PGRID(1,IMOD),
     >        TA, TB, TC, UB, GA, H, QQ, SCMF, V, VA, VB, PP,
     >        BCORE(1,IMOD), DBDR(1,IMOD), VDU(1,IMOD), 
     >        GRADI(1,IMOD), VDOP,
     >        W0, NBLINE, DELXLAP, XRANGERED, XRANGEBLUE,
     >        XJNUE, OPAK, ETAK,
     >        XPLOT, YPLOT, XJC(1,1,IMOD), NFDIM, XJCIND,
     >        LINE, MODHEAD(IMOD), JOBNUM(IMOD), NREDMAX, WREDI,
     >        WREDI0, VDUEL, BWESEX,
     >        XLAMBDA(1,IMOD), XLAM, T(1,IMOD), 
     >        RNE(1,IMOD), POPNUM(1,1,IMOD), 
     >        ENTOT(1,IMOD), RSTAR(IMOD), MAINPRO,
     >        MAINLEV, NOM, KODAT, NDIM, N, 
     >        MAXATOM, LEVEL, NCHARG,
     >        WEIGHT, ELEVEL, EION, EINST, ALPHA, SEXPO, ADDCON1,
     >        ADDCON2, ADDCON3, IGAUNT, 
     >        SIGMATHK, SEXPOK, EDGEK, NF(IMOD),
     >        KONTNUP, KONTLOW, LASTKON, XDATA, REDIS, LBREF, INDREF,
     >        NPDIM, A, B, C, WE, BX, WX, EDDI, IWARNJ0,
     >        IWARN, FNUEC, LSDWL, ST, BELIFI, 
     >        DENSCON(1,IMOD), FILLFAC(1,IMOD), ENTOTDENS, 
     >        bBIGBANDLIMIT,
C *** IRON: Additional Parameter used in SUBROUTINE FECHECK called from FORMCMF
     >        INDRB, IFRBSTA, IFRBEND, LASTFE,
     >        CLIGHT, VDOPFE, DXFE, XLAM0FE,
     >        INDFEACT, BFECHECK,
C *** IRON: Additional Parameter used in SUBROUTINE CMFFEOP called from FORMCMF
     >        SIGMAFE, OPAFE(1,1,IMOD), ETAFE(1,1,IMOD), 
     >        IFENUP, IFELOW, BIRONLINES,
     >        SIGMAACT, OPAFEI, ETAFEI, NFL, BAIRWAVELENGTHSET, 
     >        BAIRWAVELENGTH, VSIDU,
     >        DD_VDOPDU(1,1,IMOD), NATOM, YSCRATCH,
     >        MAXLAP, IND_ELLINE, bDDOPAFORMCMF, bDDFECONVOL)   

      ENDDO
      WRITE (0,*) 'FORMCMF finished for '//FREQIN
C***********************************************************************

C***  A second model, if involved, might have a different radius grid. 
C***  All relevant quantities are now interpolated to the grid of the 
C***  main model. 
C***  VDU and DD_VDOPDU are not only specific for the current RANGE
C***  and therefore saved before being scaled, and will be restored 
C***  when the current range is done. 
      IF (NMOD .EQ. 2 .OR. NOWIND_LINE .NE. 'NONE') THEN
         VDU_ORIG       = VDU
         DD_VDOP_ORIG   = DD_VDOP
         DD_VDOPDU_ORIG = DD_VDOPDU
         DD_VMICDU_ORIG = DD_VMICDU
         T_ORIG         = T
         ENTOT_ORIG     = ENTOT
         RNE_ORIG       = RNE
      ENDIF

      IF (NMOD .EQ. 2) THEN
         CALL MERGE_RGRID (RADIUS, NDDIM, MAXMOD, ND, RADIUS_MERGED,
     >                     ND_MERGED, PGRID, NP, NPDIM, PGRID_MERGED, 
     >                     NP_MERGED, ZGRID_MERGED, SECMOD_RRANGE)
C***     Note: Number of impact parameters NP might have been enhanced 
C***           by subr. MERGE_RGRID
         JPFIRST = MAX(1,           JPFIRST_ORIG)
         JPFIRST = MIN(JPFIRST,     NP_MERGED-1)
         JPLAST  = MIN(NP_MERGED-1, JPLAST_ORIG)
         JPLAST  = MAX(JPLAST, JPFIRST)

      ELSE
C***  "Normal" case without SECOND MODEL 
         ND_MERGED     = ND(1)
         NP_MERGED     = NP(1)
         RADIUS_MERGED = RADIUS(1:ND(1), 1)
         PGRID_MERGED  = PGRID(1:NP(1), 1)
         ZGRID_MERGED  = ZGRID(1:NDDIM, 1:NPDIM,1)
      ENDIF

C***  Default (no WINDROT, no SECONDMODEL)
      DO JP = 1, NP_MERGED
            NPHI(JP) = 1
            PHIWEIGHT(1,JP) = 1.
      ENDDO


C***  NOWIND option: omit (part of) the wind in the formal integral;
C***     this subroutine modifies RADIUS_MERGED etc. 
      IF (NOWIND_LINE .NE. 'NONE') THEN 
         CALL NOWIND (NOWIND_LINE, RCON, NATOM, ATMASS, ABXYZ,
     >                ND, RNE, T, RADIUS, VELO, VDOP, TAURCONT, 
     >                ND_MERGED, RADIUS_MERGED, PGRID_MERGED, 
     >                NP_MERGED, ZGRID_MERGED, IERR)
         IF (IERR .GT. 0) GOTO 2
C***     Note: Number of impact parameters NP might have been enhanced 
C***           by subr. MERGE_RGRID
         JPFIRST = MAX(1,           JPFIRST_ORIG)
         JPFIRST = MIN(JPFIRST,     NP_MERGED-1)
         JPLAST  = MIN(NP_MERGED-1, JPLAST_ORIG)
         JPLAST  = MAX(JPLAST, JPFIRST)
      ENDIF

C*******************************************************************
C***  Preparation in case of wind rotation
C*******************************************************************
      IF (VSINI .GT. 0.) THEN
         WRITE(0,'(A,F7.1)') 'Wind Rotation: VSINI [km/s] =', VSINI
C***     ROTATION_PREP evaluates RCOROT, inserts more core-intersecting 
C***     points in PGRID_MERGED, re-calculates ZGRID_MERGED, 
C***     and establishes / modifies the array of azimuthal angles PHIARR 
        
         CALL ROTATION_PREP (RCOROTLINE, RADIUS, VELO, TAURCONT, 
     >                    ND, RCOROT, NPHI, VSIDU, XMAX, 
     >                    ND_MERGED, RADIUS_MERGED,
     >                    NP_MERGED, PGRID_MERGED, ZGRID_MERGED,
     >                    NDDIM, NPDIM, NPHIMAX, PHIWEIGHT, PHIARR, 
     >                    BVSINI_AT_RCOROT, DX_WINDROT)        

C***  If JPFIRST and/or JPLAST were specified: 
C***     make sure that theu fall into the allowed range
C***     Note: Number of impact parameters NP_MERGED might have been enhanced 
C***           by subr. ROTATION_PREP
         JPFIRST = MAX(1,       JPFIRST_ORIG)
         JPFIRST = MIN(JPFIRST, NP_MERGED-1)
         JPLAST  = MIN(NP_MERGED-1, JPLAST_ORIG)
         JPLAST  = MAX(JPLAST, JPFIRST)

         CALL PLOT_WINDROT_GRID (PGRID_MERGED, NPDIM, JPFIRST, JPLAST, 
     >      LPHISTA_ORIG, LPHIEND_ORIG, NPHI, 
     >      PHIARR, NPHIMAX, XPLOT, YPLOT)

          IF (LPHISTA_ORIG .GT. 0 .OR. LPHIEND_ORIG .LT. 999) 
     >       WRITE (0,'(A,/,A, 2I4)') 
     >       '**** TEST RUN WITH RESTRICTED RANGE OF ANGLE INTEGRAL:', 
     >       '**** LPHISTA, LPHIEND =',  LPHISTA_ORIG, LPHIEND_ORIG

      ENDIF

      IF (NMOD .EQ. 2 .OR. NOWIND_LINE .NE. 'NONE') THEN
         CALL RESCALE_SECMOD (NDDIM, NFLDIM, MAXLAP, MAXMOD, MAXATOM,
     >       ND, NATOM, NFL, NBLINE, RADIUS, RADIUS_MERGED, ND_MERGED,
     >       VDU, OPAL, ETAL, ETACK, ETACCK, OPACK, OPAFE, ETAFE,
     >       DD_VDOPDU, DD_VDOP, DD_VMICDU, T, ENTOT, RNE, VEC_SECMOD, 
     >       RSTAR, NMOD)
      ENDIF

      IF (NMOD .EQ. 2) THEN
         CALL SECONDMODEL_PREP (ZINTER, NPHI, PGRID_MERGED, NP_MERGED, 
     >        NPDIM, NPHIMAX, PHIARR, PHIWEIGHT, PHI_VEC, 
     >        SECONDMODEL_LINE, 
     >        JPFIRST, JPLAST, LPHISTA_ORIG, LPHIEND_ORIG)

C***  Test output in case of a single ray
         IF (JPFIRST == JPLAST .AND. LPHISTA_ORIG == LPHIEND_ORIG) THEN
            WRITE (0,'(A,I3,A,F8.3,A,I3,A,F8.3)') 
     >      '*** Single Ray at P(', JPFIRST, ')=', 
     >      PGRID_MERGED(JPFIRST), 
     >      ' and PHI(', LPHISTA_ORIG, ')=',
     >      PHIARR(LPHISTA_ORIG,JPFIRST)*180./PI
            WRITE (0,'(A,2F8.3)') 'Intersection points with second model:'
     >       // ' Z1, Z2 =', ZINTER(1,JPFIRST,LPHISTA_ORIG), 
     >                       ZINTER(2,JPFIRST,LPHISTA_ORIG)
         ENDIF
      ENDIF

C*************************************************************
C***  Line broadening preparation (if requested)
C************************************************************************
      IF (BROAD) THEN 
         DO IMOD=1, NMOD
C***        initialize counter for line-broadening profile table
            NLPHITAB = 0
            DO NBL=1, NBLINE
               IND=INDLAP(NBL)
               LOW=INDLOW(IND)
               NUP=INDNUP(IND)
               NA = NOM(LOW)
               CALL STARKBROAD (KODATIND, NOM, 
     >          PHITAB(-NFDIMPHITAB,1,1,IMOD), 
     >          NFDIMPHITAB, NLDIMPHITAB, NLPHITAB,
     >          ND_MERGED, RADIUS_MERGED, OPA, OPAL(1,NBL,1), LINPRO(NBL), 
     >          AVOIGT, NBL, MAXLAP, IPOINTERPHITAB(NBL),
     >          XMAX, XMAXBROAD, XMAXLIN(NBL), NDDIM, DXMAX, PATH_VCSSB, 
     >          LOW, NUP, PHISCRATCH, XLAMLAP(NBL), ALN, T(1,IMOD), 
     >          ENTOT(1,IMOD), RNE(1,IMOD), LEVEL, MAINQN, 
     >          NCHARG, POPNUM(1,1,IMOD), N, VDOP, ELEVEL, EION, 
     >          PATH_LEMKE_DAT, DD_VDOP(1,NA,IMOD), 
     >          DD_VDOPDU(1,NA,IMOD), 
     >          IND_ORIGLEV, BDD_VDOP, DD_VMICDU(1,IMOD), 
     >          GRIEMPAR, TAUMAX, TAUMINBROAD, IMOD, MAXMOD, EINST, NDIM)

C***          Ensure XMAX covers maximum broadening  
              XMAX = MAX(XMAX, XMAXBROAD)

            ENDDO
         ENDDO
       ENDIF
C************************************************************************

ccc   **** test output *****************************
      IF (.false.) THEN
      do imod=1, nmod
         write (0,*) 'GRIEMPAR of MODEL', imod
         do nbl=1, nbline
            do l=1, nd(imod)
               if (griempar(nbl,l,imod) .gt. .0) then
                  onechar = '+'
               elseif (griempar(nbl,l,imod) .lt. .0) then
                  onechar = '-'
               else
                  onechar = '0'
               endif
               write (0,'(a1,$)') onechar
            enddo
            write (0,'(1x,I2,2x,A)') nbl, linpro(nbl)
         enddo
      enddo
      ENDIF
ccc   **** end of test output *****************************

C***  DEFINING THE OBSERVER'S FRAME FREQUENCY RANGE (FREMIN, FREMAX) 
C***  Note: this range must be smaller than covered by FORMCMF
      IF (RANGERED .NE. RANGEBLUE) THEN
         FREMIN = XRANGERED
         FREMAX = XRANGEBLUE
         RANGERED  = .0
         RANGEBLUE = .0
      ELSE
         FREMIN = XCMFRED  + 1.1*VMAXDU
         FREMAX = XCMFBLUE - 1.1*VMAXDU
      ENDIF

      WRITE (0,'(A,F11.4,A,F11.4,A)')
     > 'XCMFBLUE=',XCMFBLUE," (",XLAM * EXP( XCMFBLUE * ALN )," Ang)"
      WRITE (0,'(A,F11.4,A,F11.4,A)')
     > 'FREMAX  =',FREMAX," (" ,XLAM * EXP( FREMAX * ALN )," Ang)"
      WRITE (0,'(A,F11.4,A,F11.4,A)')
     > 'FREMIN  =',FREMIN," (",XLAM * EXP( FREMIN * ALN )," Ang)"
      WRITE (0,'(A,F11.4,A,F11.4,A)')
     > 'XCMFRED =',XCMFRED," (",XLAM * EXP( XCMFRED * ALN )," Ang)"

      IF (LSDWL.GT.0) THEN
          CALL PRIDWL (VELO(1,1), GRADI(1,1), VDOP,
     >                 POPNUM(1,1,1), ENTOT(1,1),
     >                 RADIUS, RSTAR(1), OPAL(1,1,1),
     >                 EINST, WEIGHT, FNUEC, XLAP, 
     >                 LOW, NUP, ND, NDIM, DELW,
     >                 LSDWL, OPA(1,1), THOMSON, ADELW,
     >                 JOBNUM(1), MODHEAD(1))
        GOTO 501
      ENDIF
           
C***  Prepare LOOP FOR EVERY OBSERVER'S FRAME-FREQUENCY  -----------------------
C***
C***  Define DXOBS = frequency spacing in Doppler units
      DXOBS = DXMAX
      NFOBS = IFIX((FREMAX-FREMIN)/DXOBS) + 2
      IF (NFOBS .GT. NFODIM) THEN
            WRITE (0,*) '********************************************'
            WRITE (0,*) 'ERROR: NFODIM TOO SMALL, REQUIRED:', NFOBS
            WRITE (0,*) '********************************************'
            WRITE (*,*) '********************************************'
            WRITE (*,*) 'ERROR: NFODIM TOO SMALL, REQUIRED:', NFOBS
            WRITE (*,*) '********************************************'
          STOP 'FATAL ERROR detected in FORMAL'
          ENDIF

C**  DEFINING ZERO-POINT AND INCREMENT OF THE OBSERVER'S FRAME FREQUENCY
      XOBS0 = FREMAX + DXOBS

      WRITE (0,'(1X,A,I6)')
     >      'Number of Observers-frame frequencies : ', NFOBS

      DO K = 1, NFOBS
          PROFILE(K)=0.
         CPROFILE(K)=0.
         PROFILEC(K)=0.
      ENDDO

C***  Decode LIMBDARKENING options and 
C***  prepare limb-darkening weights if requested
      IF (LIMB_LINE .NE. 'OFF') 
     >   CALL LIMBDARK_PREP (LIMB_LINE, XOBS0, DXOBS, NFOBS, XLAM, ALN, 
     >                      WEIGHT_LIMBDARK, BLIMB)

C***  Prepare total number of rays (for progress-bar)
      NPHISUM = 0
      DO JP= JPFIRST, JPLAST
         NPHISUM = NPHISUM + NPHI(JP)
      ENDDO
      IMPATIENCE = MAX (1, NPHISUM * NFOBS / 50)

      WRITE (0,'(A,I4,I4,I4)') 'JPFIRST, JPLAST, NP:', 
     >                          JPFIRST, JPLAST, NP_MERGED
      WRITE (0,'(A,I6)') 'Total number of azimuth-angle points: ',
     >               NPHISUM
      WRITE (0,'(A,/,A)') 'Progress Bar for this range:', 
     > '0% 10% 20%  30%  40%  50%  60%  70%  80%  90% 100%'

      NRAYDONE = 0


C***  LOOP OVER IMPACT PARAMETERS *********************************
      DO 14 JP=JPFIRST, JPLAST
         PJPJ = PGRID_MERGED(JP) * PGRID_MERGED(JP)
C***     IRAY = starting position of depth vector with impact parameter JP
C***            within the z-array if used as 1-dimensional field 
         IRAY=ND_MERGED * (JP-1) + 1

C***     for limb-darkening function
         EMINT_P(JP) = .0

C***     Prepare Loop over azimuth angles 
C***      - if no wind rotation and no second-model: only 1 point) 
C***      - if LPHISTA or LPHIEND are specified, make sure that 
C***        that they fall into the allowed range at current JP
          IF  (LPHISTA_ORIG .GT. 0 .OR. LPHIEND .LT. 999) THEN
              LPHISTA = MAX(1,       LPHISTA_ORIG)
              LPHISTA = MIN(LPHISTA, NPHI(JP))
              LPHIEND = MIN(NPHI(JP),LPHIEND_ORIG)
              LPHIEND = MAX(LPHIEND, LPHISTA)
          ELSE
              LPHISTA = 1
              LPHIEND = NPHI(JP)
          ENDIF

C***     Normalization of PHIWEIGHT to unity
C***       i.e.phi integral always returns the average over the P-ring
C***      (or part of that ring, if LPHISTA or LPHIEND are restricted)

         WEIGHTSUM = .0
         DO LPHI=1, NPHI(JP)
             WEIGHTSUM = WEIGHTSUM + PHIWEIGHT(LPHI,JP)
         ENDDO

         DO LPHI=1, NPHI(JP)
             PHIWEIGHT(LPHI,JP) = PHIWEIGHT(LPHI,JP) / WEIGHTSUM
         ENDDO

C***     Loop over azimuth angles ***********************************
         DO 21 LPHI = LPHISTA, LPHIEND

C***         Loop over obsframe-frequencies  ************************
             DO 13 K=1, NFOBS 

             XO=XOBS0 - K*DXOBS
             XLAMREF = EXP (XLAMLN + XO * ALN)
             DLAM(K)= XLAMREF - XLAM

             CALL PREPRAY (ZGRID_MERGED(IRAY,1), COS(PHIARR(LPHI,JP)),
     >                PGRID_MERGED, ND_MERGED, NDDIM, NP_MERGED, JP,XO, 
     >                LTOT, PWEIGHT, 
     >                CORE, VDU, 
     >                RADIUS_MERGED, OPAL, ETAL, RRAY, OPARAY, ETARAY, 
     >                ETACRAY, OPALRAY, ETALRAY, NFLDIM, ZRAY, XCMF,
     >                NDADDIM, NBLINE, MAXLAP, REDIS, ETACK, ETACCK, OPACK, 
     >                     XCMFRED, XCMFBLUE, DXCMF, BCORE, BCOREL, 
     >                     DBDR, DBDRL, K, POROLENGTH, POROLENGTHRAY,
     >                     RCOROT, VSIDU, DD_VDOPDU_RAY, DD_VDOPDU, NATOM, 
     >                     NBFIRST, NBLAST, 
     >                     MAXATOM, BVSINI_AT_RCOROT, NMOD, MAXMOD, 
     >                     ZINTER(1,JP,LPHI), XMAX, DELXLAP)
C***         If no p-integral is calculated (test case), omit
C***         integration weight p*dp
             IF (JPFIRST .EQ. JPLAST) PWEIGHT = 1.

             CALL OBSFRAM (LTOT, CORE, XMAX, XMAXLIN, 
     >               EMINT, CEMINT, BCOREL, DBDRL,
     >               TAUMAX, PJPJ, ZFINE,  
     >               OPAFINE, OPAFC, OPALFIN, ETAFINE,
     >               ETAFC, ETALFIN, SFINE, CSFINE, RRAY, 
     >               OPARAY, OPALRAY, ETARAY, ETACRAY, 
     >               ETALRAY, ZRAY, XCMF, MAXXN,
     >               NDADDIM, NDDIM, DELXLAP, NBLINE,
     >               IVERSION, LINPRO, AVOIGT,
     >               TAU, TAUC, DTAU, DTAUC, WTAU, WTAUC,
     >               XCMFFINE, POROLENGTHFINE, MAXLAP, DXMAX, 
     >               BIRONLINES, OPAFE, ETAFE, NFLDIM, 
     >               XCMFBLUE, XCMFRED, DXCMF, POROLENGTHRAY,
     >               PHITAB, NFDIMPHITAB, NLDIMPHITAB, IPOINTERPHITAB,
     >               RADIUS_MERGED, ND_MERGED, 
     >               DD_VDOPDU_RAY, NATOM, NBFIRST, NBLAST,
     >               IND_ELLINE, DD_VDOPDU_FINE_NORMFAC,
     >               DD_VDOPDU_FINE_SQRD, 
     >               DD_VDOPDU_FINE,
     >               GRIEMPAR, KODAT, VDOP, 
     >               INDCUT, ZINTER(1,JP,LPHI), NMOD, MAXMOD, K)

           PROFILEC(K)= PROFILEC(K) +  EMINT*PHIWEIGHT(LPHI,JP)*PWEIGHT
           CPROFILE(K)= CPROFILE(K) + CEMINT*PHIWEIGHT(LPHI,JP)*PWEIGHT
           IF (BLIMB) EMINT_P(JP) = EMINT_P(JP) + 
     >                EMINT*PHIWEIGHT(LPHI,JP)*WEIGHT_LIMBDARK(K)

C***       Progress bar
           NRAYDONE = NRAYDONE + 1
           IF ( (NRAYDONE / IMPATIENCE) * IMPATIENCE .EQ. NRAYDONE) THEN
              WRITE (0,'(A,$)') 'X'
           ENDIF

   13      CONTINUE
C***       End of frequency loop ********************************************

   21     CONTINUE
C***      End of phi-loop ***********************************************

   14 CONTINUE
C***  End of p-loop *************************************************

      WRITE (0,'(/,A)') 'Range done: ' // FREQIN 
C*************************************************

C***  Normalization: 
C***   PROFILEC = calibrated flux
C***   PROFILE  = normalized flux
C***  CPROFILE  = contimuum  flux
      DO K=1, NFOBS
         PROFILE(K) = PROFILEC(K) / CPROFILE(K)
      ENDDO
 
C***  Transformation of the wavelength scale to air, if appropriate
      IF (BAIRWAVELENGTH) THEN
        DO K=1, NFOBS
          XO=XOBS0 - K*DXOBS
          XLAMREF = EXP (XLAMLN + XO * ALN)
          XLAM2   = XLAMREF * XLAMREF
          DLAM(K) = DLAM(K) - XLAMREF*(2.735182E-4 + 131.4182
     >                 / XLAM2 + 2.76249E8 / (XLAM2*XLAM2))
        ENDDO
      ENDIF

C***  OUTPUT FOR THE DETECTED LINE: 
C***  WRITE BUFFERED STRING ARRAY WITH INDIVIDUAL COMMENTS 
C***  (OPTION 'STRING') TO FILE 'PLOT' (KANAL = 1)
      DO 11 NSTRI=0, MIN0(NSTRING,MAXSTRI)
        WRITE (KANAL1,'(A)') STRING1(NSTRI)(:IDX(STRING1(NSTRI)))
        IF (NSTRI .GT. 0) 
     >   WRITE (*     ,'(A)') STRING1(NSTRI)(:IDX(STRING1(NSTRI)))
   11 CONTINUE
      LOW=INDLOW(LINE)
      NUP=INDNUP(LINE)
      IF (LSDWL .LE. 1) THEN
        CALL TRAPLO (1, PROFILE, DLAM, NFOBS,
     >               LINE, FREQIN, MODHEAD(1), JOBNUM(1),
     >               DISP, N, NBLINE, IDENT, OSMIN,
     >               XLAM, XLAMLAP, INDLAP, INDLOW, INDNUP, LEVEL,
     >               ELEVEL, WEIGHT, EINST, NDIM, ABSWAV, 
     >               .FALSE., FUNIT, NMOD, 
     >               LPHISTA, LPHIEND, NPHI,  
     >               XUNIT, .FALSE., BAIRWAVELENGTH, 
     >               RSTAR(1), KODATIND, MAXATOM, NOM)
      ENDIF
      
      IF (BCONT) THEN 
           WRITE(0,'(A,A)') 'Plotting Continuum flux'
           WRITE (KANAL1,'(A)') 'PLOT: CONTINUUM '//FREQIN
           WRITE (KANAL1,'(A)') '* CONTINUUM '//FREQIN
          CALL TRAPLO (1, CPROFILE, DLAM, NFOBS, LINE,
     >          FREQIN, MODHEAD(1), JOBNUM(1), DISP, N, NBLINE, 
     >          IDENT, OSMIN,
     >          XLAM, XLAMLAP, INDLAP, INDLOW, INDNUP, LEVEL,
     >          ELEVEL, WEIGHT, EINST, NDIM, ABSWAV, 
     >          BCONT, FUNIT, NMOD, 
     >          LPHISTA, LPHIEND, NPHI, 
     >          XUNIT, BCALIBRATED, BAIRWAVELENGTH,
     >          RSTAR(1), KODATIND, MAXATOM, NOM)
      ENDIF
C***    Plotten des kompletten, nichtrektifizierten spektrums
      IF (BCALIBRATED) THEN
          WRITE(0,'(A,A)') 'Plotting flux-calibrated spectrum'
          WRITE (KANAL1,'(A)') 'PLOT: CALIBRATED '//FREQIN
          WRITE (KANAL1,'(A)') '* CALIBRATED '//FREQIN
          CALL TRAPLO (1, PROFILEC, DLAM, NFOBS, LINE,
     >         FREQIN, MODHEAD(1), JOBNUM(1),
     >         DISP, N, NBLINE, 
     >         IDENT, OSMIN,
     >         XLAM, XLAMLAP, INDLAP, INDLOW, INDNUP, LEVEL,
     >         ELEVEL, WEIGHT, EINST, NDIM, ABSWAV, 
     >         .TRUE., FUNIT, NMOD, 
     >         LPHISTA, LPHIEND, NPHI, 
     >         XUNIT, BCALIBRATED, BAIRWAVELENGTH,
     >         RSTAR(1), KODATIND, MAXATOM, NOM)
      ENDIF
      

  501 CONTINUE
      IF (LSDWL .GT. 1)
     > CALL TRADWL (1, DELW, ADELW, VELO(1,1), ENTOT(1,1), ND,
     >              LINE, MODHEAD(1), JOBNUM(1),
     >              LEVEL(NUP), XLAM, XPLOT, TAURCONT, RADIUS,
     >              TRANSDWLLINE)
   10 CONTINUE
 
C***  PRINTOUT OF PROFILE TABLE ONLY WHEN NOT "DWL" MODE ACTIVE
      IF (LSDWL.LE.1) THEN
         CALL PRIPRO (XLAM, VDOP, NFOBS, PROFILE,
     >         XOBS0, DXOBS, JOBNUM(1),
     >         VSINI, MODHEAD(1), DLAM, LSPRO, IFIRST,
     >         NPHI, LPHISTA, LPHIEND,
     >         DXMAX, TAUMAX, XMAX, TAUMINBROAD, JPFIRST, JPLAST, 
     >         PGRID_MERGED,  
     >         PWEIGHT, ELEVEL, NDIM, INDNUP, INDLOW,
     >         LASTIND, INDLAP, MAXLAP, FREMAX, FREMIN,
     >         NBLINE, RMAX, XLAMLAP,
     >         WEIGHT, LEVEL, EINST, VMAX,
     >         IVERSION, LINE, REDIS,
     >         BWESEX, BROAD, LINPRO, AVOIGT, BNOCONT, NMOD, 
     >         DENSCON, FILLFAC, BFEWARNING, SPZ1,
     >         SPZ2,MANIPOP_OPTIONS,ND, 
     >         MULTIIND, NMULTI, BAIRWAVELENGTH, MAXSPZ, RCOROT,
     >         BDD_VDOP, DD_VDOP_LINE, BMICROTURB, EL_MIN, 
     >         IVDOPSTATUS, BVSINI_AT_RCOROT)
      ENDIF

      IF (BPLOTVDOP) THEN
        CALL PLOTVDOP (RADIUS, TAURCONT, VELO, DD_VDOP, 
     >                 ND, NATOM, BPLOTVDOP, SYMBOL, BMICROTURB,
     >                 DD_VMIC, DD_VDOP_LINE, VDOPPLOT_LINE, MODHEAD(1),
     >                 BIRONLINES)
      ENDIF

      IF (BLIMB) 
     >   CALL LIMBDARK_OUTPUT (KANAL1, LIMB_LINE, NFOBS, EMINT_P,
     >                            NP_MERGED, PGRID_MERGED)

C***  Entry in case of error which enforced skipping the  current range
    2 CONTINUE

C***  The no. of core-intersecting impact parameters might have been modified 
C***  in subr. ROTATION_PREP; therefore, P and Z were cloned and 
C***  are now restored for the potential next RANGE

      NP(1) = NP_ORIG 
      DO JP=1, NP_ORIG
         PGRID(JP,1) = PGRID_ORIG(JP)
      ENDDO

      ZGRID = ZGRID_ORIG

C***  In case of a second model, several vectors of both models have been
C***  saved before being rescaled to the merged radius grid 
C***  by subr. RESCALE_SECMOD, and must now be restored for 
C***  a possible further range
      IF (NMOD .EQ. 2 .OR. NOWIND_LINE .NE. 'NONE') THEN
         VDU       = VDU_ORIG
         DD_VDOP   = DD_VDOP_ORIG
         DD_VDOPDU = DD_VDOPDU_ORIG
         DD_VMICDU = DD_VMICDU_ORIG
         T         = T_ORIG
         ENTOT     = ENTOT_ORIG
         RNE       = RNE_ORIG
      ENDIF


      GOTO 1
C***  ENDLOOP   (RANGE COMPLETED - continue with next RANGE) -----------

C***********************************************************************
      
 20   CLOSE (1)
 
      CALL CLOSMS (7, IERR)

C***  CLOSE THE CARDS-FILE  
      CLOSE (2)

C***  CLOSE THE FILE for plots of the wind-rotation geometry grid 
      CLOSE (65)

      CALL JSYMSET ('G0', '0')

      CALL STAMP (OPSYS, 'FORMAL', TIM1)

      STOP 'O.K.'
      END
 

