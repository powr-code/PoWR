C***  MAIN PROGRAM STEAL  ******************************************************
      SUBROUTINE STEAL
C*******************************************************************************
C***  STATISTICAL EQUATIONS WITH APPROXIMATE LAMBDA-OPERATORS
C*******************************************************************************

      IMPLICIT NONE

C***  SET ARRAY DIMENSION PARAMETERS
C***  IRON: ADD GENERIC ION TO MAXATOM
      INTEGER, PARAMETER :: MAXATOM =          26 
C*** !!! ATTENTION !! USE of COMMON BLOCK in LINPOP !!!
C*** NDIM and NFDIM are also set in LINPOP 
      INTEGER, PARAMETER :: NDIM    =        1560 
      INTEGER, PARAMETER :: NFDIM   = 2*NDIM + 400 
      INTEGER, PARAMETER :: MAXAUTO =        2850 
C*** !!! ATTENTION !! USE of COMMON BLOCK in LINPOP !!!
C*** MAXIND is also set in LINPOP 
      INTEGER, PARAMETER :: MAXIND  =       20000 
      INTEGER, PARAMETER :: MAXINDE = MAXAUTO+MAXIND 
      INTEGER, PARAMETER :: MAXNDR  =     MAXAUTO 
      INTEGER, PARAMETER :: NDDIM   =          89 
      INTEGER, PARAMETER :: NPDIM   =          94 
      INTEGER, PARAMETER :: NFLDIM  =          40 
      INTEGER, PARAMETER :: MAXHIST =        4000 
      INTEGER, PARAMETER :: MAXXDAT =          11 

C***  To store the Source funktion on the COLI-Fine grid
      INTEGER, PARAMETER :: MAXFINE =      160000 

C***  MAXIMUM ION CHARGE WHICH MAY OCCUR (SEE ALSO SUBR. GAUNTFF)
      INTEGER, PARAMETER :: MAXION = 27 

      INTEGER, PARAMETER :: NDIMP2  = NDIM + 2 
      INTEGER, PARAMETER :: MAXKONT = NDIM 
      INTEGER, PARAMETER :: MAXKODR = NDIM 

      INTEGER, PARAMETER ::MAXLAP = 40

C***  NUMBER OF LEVELS WHICH CAN BE PLOTTED BY PLOTPOP
      INTEGER, PARAMETER ::MAXSETS = 25
C***  NUMBER OF PLOT OPTIONS WHICH ARE STORED FOR DEFERRED EXECUTION
C***  (JNUE, JLINE PLOTS), AND MAXIMUM LENGTH OF PLOT VECTORS
      INTEGER, PARAMETER :: MAXPLOTOPT = 100 
      INTEGER, PARAMETER :: MAXPLOTN   = 19000

C***  NUMBER OF ENTRYS STORED IN THE GAMMA HISTORY
      INTEGER, PARAMETER :: MAXGAHIST = 100

C***  IRON: COMMON BLOCK ARRAY DIMENSIONS
C***  include "dimblock"
      INTEGER, PARAMETER :: INDEXMAX = 1E7, NFEREADMAX = 3E5    !std
C      INTEGER, PARAMETER :: INDEXMAX = 4E7, NFEREADMAX = 5E5     !vd20
C      INTEGER, PARAMETER :: INDEXMAX = 1E8, NFEREADMAX = 6E5     !xxl

C***  Fine Integration in STEAL - SETXJFINE
      INTEGER, PARAMETER :: IFF_MAX =   80000      !std
C      INTEGER, PARAMETER :: IFF_MAX =  200000      !vd20
C      INTEGER, PARAMETER :: IFF_MAX =  300000      !xxl

      INTEGER, PARAMETER :: TINYINT = SELECTED_INT_KIND(2)

C***  IFF_MAX_MS must be 1/8 * IFF_MAX !!!
C      INTEGER, PARAMETER :: IFF_MAX_MS =   IFF_MAX / 8 


C***  HANDLING OF DIELECTRONIC RECOMBINATION / AUTOIONIZATION (SUBR. DATOM)
      INTEGER, DIMENSION(MAXAUTO) :: LOWAUTO, IONAUTO
      REAL, DIMENSION(MAXAUTO) :: WAUTO, EAUTO, AAUTO, KRUDAUT
      REAL, DIMENSION(NDIM) :: DRRATEN, RDIEL, RAUTO, DRJLW, 
     >                         DRJLWE, DRLJW
      REAL, DIMENSION(MAXKODR) :: KODRNUP, KODRLOW

      REAL, DIMENSION(NPDIM) :: P
      REAL, DIMENSION(NDDIM,NPDIM) :: Z

      CHARACTER(MAXHIST*8) :: MODHIST
      REAL, DIMENSION(MAXXDAT) :: XDATA
      REAL, DIMENSION(NDDIM) :: GRADI, TAUTHOM, TAUROSS, TAUROSScont,
     >                          TNEW, TOLD, ENTOT, RNE, RADIUS, ETA,
     >                          OPA, THOMSON, VELO, OLDVELO, OLDRADI,
     >                          OPALROSS, XJTOTL, XKTOTL, XNTOTL,
     >                          HTOTL, HTOTCMF0, HTOTM, HTOTG, HTOTOBS,
     >                          AGRAV, AMECH, ARAD, APRESS, DTKUBAT,
     >                          ENTOTDENS, FTCOLI, TOLD2, TOLD3,
     >                          GRSTATIC, VMACH, FLUXERR, CORRS
      INTEGER, DIMENSION(NDDIM) :: ITNE, IWARN
      REAL, DIMENSION(NDDIM, NDIM) :: POPNUM, DEPART, POPLTE,
     >                                POP1, POP2, POP3
      REAL, DIMENSION(NDIM) :: WEIGHT, ELEVEL, EION, ENLTE, ENOLD, 
     >                         DEPARTNDorg
      INTEGER, DIMENSION(NDIM) :: NCHARG, MAINQN, NOM, IONGRND
      REAL, DIMENSION(NDIM, NDIM) :: CRATE, DCRATEDT, RRATE, EINST
      LOGICAL, DIMENSION(NDIM, NDIM) :: BCOLLIDONE
      REAL, DIMENSION(4, NDIM) :: ALTESUM
      REAL, DIMENSION(NFLDIM) :: PHI, PWEIGHT
      REAL, DIMENSION(MAXNDR) :: ELEVDR, WEIGHTDR
      INTEGER, DIMENSION(MAXNDR) :: NCHARGDR, NOMDR
      INTEGER, DIMENSION(MAXAUTO) :: INDNUPDR
      REAL, DIMENSION(MAXATOM) :: ABXYZ, ATMASS, STAGE
      INTEGER, DIMENSION(MAXATOM) :: KODAT, NFIRST, NLAST, IMAXPOP
      REAL, DIMENSION(MAXATOM,MAXATOM) :: SIGMATHK, SEXPOK, EDGEK
      REAL, DIMENSION(MAXKONT) :: SIGMA1I, ALPHA, SEXPO,
     >                            ADDCON1, ADDCON2, ADDCON3
      REAL, DIMENSION(NFDIM) :: FWEIGHT, EMFLUX, EMCOLI 
      REAL, DIMENSION(NDDIM,NFDIM) :: XJC, WCHARM, XJCorg
      REAL, DIMENSION(NDDIM,MAXINDE) :: XJL
      REAL, DIMENSION(NFDIM,NDDIM) :: SCOLD
      REAL, DIMENSION(NFDIM,0:MAXION) :: SIGMAFF
      REAL, DIMENSION(NFDIM,MAXKONT) :: SIGMAKI
      INTEGER, DIMENSION(MAXKONT) :: KONTNUP, KONTLOW, NFEDGE
      CHARACTER*8 IGAUNT(MAXKONT), KEYCBF(MAXKONT)
C     SPLIT: Groessere SCRATCH-Dimensionierung noetig
C OLD> DIMENSION SCRATCH(2*NDIMP2), ATEST(NDIMP2,NDIMP2)
      REAL, DIMENSION(NDIMP2*NDIMP2) :: SCRATCH
      REAL, DIMENSION(NDIMP2) :: BRS, BRSP, BRY, EN, V1, V2, DTEST
      REAL, DIMENSION(NDIMP2,NDIMP2) :: RATCO, DM, AOFF, ATEST, BTEST
      REAL, DIMENSION(MAXPLOTN) :: XPLOT, YPLOT
      LOGICAL :: BPRIUNLU, BPGAHIST, BAG,
     >           BGAMMACFIX, BGAMMALFIX, BGAMMARFIX, BGAMMADFIX,
     >           BPGAHISTE, BSMOOTHTC, BRUDZERO, BNEWTONRESET

      REAL, DIMENSION(MAXINDE) :: XJLAPP
      INTEGER, DIMENSION(10) :: NGAMC, NGAMR, NGAML, NGAMD
      REAL, DIMENSION(10) :: AGAMC, AGAMR, AGAML, AGAMD
      REAL, DIMENSION(2) :: TauCorLimits
      LOGICAL, DIMENSION(MAXAUTO) :: DRXJL
      LOGICAL, DIMENSION(MAXIND) :: LINE
      LOGICAL :: NOTEMP, TPLOT, TPLTAU, CONVERG, MOREJOBS,
     >           NOEXTRAP, STHLP, NODATOM, OLDSTART, NOCON, COMPO,
     >           DRNORUD, NOPOP, BAUTO, BAUTO_ABORT, BINBOX,
     >           BPRICORR, BCOREX, VPLOT, BUNLU, BHTOTERR, BITCONT,
     >           BGFIN, BRDIVFAILED, BTRACE_POPNUM, bHLPHIST
      CHARACTER(100) :: MODHEAD
      CHARACTER(255) :: HISTENTRY
      CHARACTER(80) :: COMMAND
      CHARACTER(20) :: ACTPAR
      CHARACTER(16) :: BUFFER16
      CHARACTER(10), DIMENSION(NDIM) :: LEVEL
      CHARACTER(10), DIMENSION(NDDIM) :: MAINPRO, MAINLEV
      CHARACTER(10), DIMENSION(MAXATOM) :: ELEMENT
      CHARACTER(10), DIMENSION(MAXSETS, NDIM) :: LEVELPL, LEVELPLDEP
      CHARACTER(10), DIMENSION(MAXNDR) :: IONDR, LEVELDR
      CHARACTER(10) :: PRILEVRA
      CHARACTER(8), DIMENSION(NDDIM) :: INCRIT, VCRIT
      CHARACTER(8), DIMENSION(NFDIM) :: KEY
      CHARACTER(8) :: NAME, BUFFER8, JOBSTR      
      CHARACTER(4), DIMENSION(MAXIND) :: KEYCBB
      CHARACTER(4) :: KODE
      CHARACTER(2) :: WRTYPE
      CHARACTER(2), DIMENSION(MAXATOM) :: SYMBOL
      CHARACTER(1), DIMENSION(NDDIM) :: CKONVER
      CHARACTER(256) :: TPLOTOPT
      CHARACTER(256), DIMENSION(MAXPLOTOPT) :: PLOTOPT
      REAL :: RMAX, GEDD, RCSAVE, DENSCON_FIX, CORRLAST
      CHARACTER(1200) :: UNLUTECLINE
      CHARACTER(120) :: DENSCON_LINE,
     >                 ThinCard, HydroCard, AlphaCard, MacroCard
      CHARACTER(8) :: OPC, GEFFKEY
      LOGICAL :: BFEWING, bTDIFFUS, BPLOCC, BTALTER
      REAL, DIMENSION(MAXIND) :: XLAMBLUE, XLAMRED, DEXFAC
      INTEGER, DIMENSION(MAXLAP,MAXIND) :: IBLENDS
      REAL, DIMENSION(NFLDIM,MAXLAP) :: PHIL
      REAL, DIMENSION(MAXLAP) :: BETA
      INTEGER, DIMENSION(2*NDDIM+8) :: IADR16
      REAL, DIMENSION(NDDIM) :: DTLOCAL
      REAL, DIMENSION(NDIM) :: DTINT, DTRMAX
      LOGICAL, DIMENSION(MAXKONT) :: NRB_CONT
      LOGICAL, DIMENSION(NDIM,NDDIM) :: ZERO_RATES
C***  New CARDS options added by ansander
      LOGICAL :: bENSURETAUMAX, bHYDROSOLVE, bBLOCKINVERSION,
     >           bFixGEFF, bTauStrict, bTauUpdated, bThinWind, 
     >           bGEddFix, bHydroDone, bNewVelo, bSUMMARY, bNoARAD,
     >           bSaveTauCont, BOOLdummy, bUpdateMass, bUseTWOPNT,
     >           bLateTau, bModHeadUpdate, bHydroHelp, bRELSCHARMERACC,
     >           bOLDVELO, bTauMaxSafe, bFULLHYDROSTAT, bGAMMARADMEAN,
     >           bBackupNow, bForceColiPP, bUCPP, bFeTCORR, bNoFeTCORR
      REAL :: VMIN, RCON, TAUACC, TAUMAX, ReduceTauCorrections
      CHARACTER(80), DIMENSION(3) :: RadiusGridParameters           !contains all RADIUS-GRID CARDS (for subr. RGRID)      

      REAL :: TEFF, GLOG, GEFFLOG, XMSTAR, GF, RSTAR, VDOP, XMAX,
     >        VDOPFE, DXFE, XLAM0FE, TOTOUT, TOTIN, OPARND, VMAX,
     >        EDDIHOUTJMEAN, HTOTOUTMINUS, EPSILON, REDUCE, GEDDreduce,
     >        BRRESET, SMALLPOP, POPMIN, POPRANG, DELTAC, COREX,
     >        XLAM_FINE_START, XLAM_FINE_END, WJCMIN, FQLIMIT,
     >        WORKRATIO, MG, GAMMAC, Rcritical,
     >        DUNLU_LOC, DUNLU_INT, DUNLU_RMAX, DUNLU_TB, TBTAU, TAUINT,
     >        GAMMAR, GAMMAL, GAMMAD, TNDCORR, DFEINDR,
     >        CORMAX, RTCMAX, RELIT, XLOGL, RTRANS, XMDOT, GEDDRAD,
     >        TEFFdummy, ZEROdummy, ZEROdummy2, GFLSAV, VTURB,
     >        VPAR1, VPAR2, HSCALE, BETA2, BETA2FRACTION, CLUMP_SEP,
     >        VPAR1_2, VPAR2_2, VFINAL, ATMEAN, OPALINE_SCALE,
     >        GAMMARADMEAN, QIONMEAN, DTDRIN, HTOTND, HNDCORFAC, 
     >        FLUXEPS, HYSTACC

      REAL, PARAMETER :: AMU = 1.6605E-24   ! Atomic mass unit in gramm

      INTEGER :: N, L, LAST, JOBNUM, NCOLIP, IDX, KANAL, NFL, LASTIND,
     >           I, J, NATOM, NAUTO, LASTKON, LASTKDR, LASTFE, LAINDHI, 
     >           ND, NP, NF, NF2, IFF_N_MS, LSRAT, JOBMAX, K, KV, LST,
     >           NDDUMMY, IHIST, IFRRA, ITORA, IPRICC, IPRILC, LSEXPO,
     >           LSTAU, LSPOP, IFLUX, IDAT, NPLOT, NEWWRC, NOLAP, ITBR,
     >           NATOUT, NUP, LOW, IPRINTZERORATES, NPLOTDEP, NPLOTOPT,
     >           LPLOCC, NSCHAR, KPLOCC, NLINE, IPLOT_XJLAPP, IPLOTOPT,
     >           IPLOT_XJCAPP, LPLOT_XJCAPP, NITER_PLOT_JAPP, INDDR,
     >           NLASER, IND, IERR, NBRCANC, NKONVER, MAXFEACT, iAMT,
     >           IWARN_NEG_XJCAPP, IWARN_NEG_XJLAPP, ITSUM, IDUMMY, 
     >           NEXTHYDRO, LASTHYDRO, MFORM, iBLOCKINVERSION, iZRType,
     >           LASTBACKUP, LASTTAU, NBACKUP, NTOLD, NKONV_THRESH,
     >           IHSSTATUS


C***  Hydro Stuff
      REAL, DIMENSION(NDDIM-1) :: ATHOM, ACONT
      REAL, DIMENSION(NDDIM) :: RHO, XMU, ALPHAF

C***  Unsoeld-Lucy terms for TEMPEQ
      REAL, DIMENSION(NDDIM) :: OPASMEAN, QFJMEAN, OPAJMEAN,
     >                          OPAPMEAN, QOPAHMEAN, OPASMEANTC, 
     >                          OPAJMEANTC, OPALAMBDAMEAN, SMEAN

      REAL, DIMENSION(MAXIND) :: XRED, XBLUE, SLOLD, DETAL, OPAL,
     >                           SLNEW, DOPAL, SCOLIND, SCNEIND,
     >                           OPACIND, XRED0, XBLUE0, OPALOLD,
     >                           XLAMZERO
      INTEGER, DIMENSION(MAXIND) :: INDNUP, INDLOW, NBLENDS
      REAL, DIMENSION(4,MAXIND) :: COCO

C*** !!! Der folgende COMMON-BLOCK steht auch in LINPOP !!!!!
      COMMON /COMIND/  INDNUP,INDLOW,XRED,XBLUE, SLOLD,
     > DETAL,OPAL,SLNEW,DOPAL,SCOLIND,SCNEIND,OPACIND,
     > XRED0,XBLUE0,COCO,OPALOLD,XLAMZERO,NBLENDS

      REAL, DIMENSION(NFDIM) :: XLAMBDA, XLAMBDA2, XJCAPP, OPAC,
     >                          SCNEW, DOPA, DETA, ETAC, EXPFAC,
     >                          XJCLP1, OPAC1, XKC, XKC2, 
     >                          XJCAPPNEW, FWTEST
      REAL, DIMENSION(NFDIM, 4) :: XJC_PLOTDATA_L

C*** !!! Der folgende COMMON-BLOCK steht auch in LINPOP !!!!!
      COMMON /COMNF/  XLAMBDA, XLAMBDA2, XJCAPP, OPAC, SCNEW, DOPA, 
     > DETA, ETAC, EXPFAC, XJCLP1, OPAC1, XKC, XKC2, XJCAPPNEW, 
     > XJC_PLOTDATA_L, FWTEST

      COMMON /COMRADI/ OLDRADI
      COMMON /COMVELO/ bNewVelo, NDDUMMY, OLDVELO
      
C***  COMMON /COMTEFF/  TRANSFERS THE EFF. TEMPERATURE FOR SUBR. PRIMOD
      COMMON /COMTEFF/ TEFFdummy,ZEROdummy,ZEROdummy2,BOOLdummy

C***  IRON: COMMON BLOCK FOR IRON-SPECIFIC DATA 
C*** !!! Der folgende COMMON-BLOCK steht auch in LINPOP !!!!!
      INTEGER, PARAMETER :: MAXFEIND  =       2500 
      INTEGER, DIMENSION(MAXFEIND) :: INDRB, INDRF, IFRBSTA, IFRBEND,
     >                                IFENUP, IFELOW, INDFEACT
      REAL, DIMENSION(MAXFEIND) :: SIGMAACT, SIGMAINT, FERATLU, 
     >                             FERATUL, FERATLU0, FERATUL0
      COMMON /IRON/ INDRB, INDRF, IFRBSTA, IFRBEND, IFENUP, IFELOW,
     >              INDFEACT, SIGMAACT, SIGMAINT, FERATLU, FERATUL,
     >              FERATLU0, FERATUL0

      REAL, DIMENSION(NFEREADMAX) :: FEDUMMY
      REAL, DIMENSION(NDDIM) :: OPAFE, ETAFE
      REAL, DIMENSION(NDDIM, MAXIND) :: WFELOW, WFENUP
      REAL, DIMENSION(INDEXMAX) :: SIGMAFE

      LOGICAL :: BFECHECK, BFEMODEL, BNUEFE 
      LOGICAL, DIMENSION(MAXIND) :: RUDLINE

C***  Depth dependent Clumping
      REAL, DIMENSION(NDDIM) :: DENSCON, FILLFAC

C***  Special Variables defined on the COLI-Fine grid
      REAL, DIMENSION(MAXFINE) :: SFINE_OLD, SFINE_NEW
      REAL, DIMENSION(MAXIND) :: XLAMMIN, XLAMMAX
      REAL, DIMENSION(MAXKONT) :: KONTHLP
      INTEGER, DIMENSION(MAXIND) :: LINEINDEX
      REAL, DIMENSION(MAXIND) :: XLAMSOR, XKMIN, XKMAX, XKMID, XKRED
      REAL, DIMENSION(MAXLAP) :: XKRED_CORE, XKBLUE_CORE, ETAL, 
     >                           PWEIGHTCL, WS
      INTEGER, DIMENSION(MAXLAP) :: NUPACT, LOWACT, LIND, LINDS
      LOGICAL, DIMENSION(MAXIND) :: BLASERL
      REAL, DIMENSION(MAXINDE) :: XJLAPPNEW
      CHARACTER(80) :: MESSAGE, MESSAGE2
C***  Plot of Approximate Radiation Fields
      LOGICAL :: BXJLAPPNEW, BXJCAPPNEW, BNEWOPER, 
     >           BPLOTAPP, BXJLAPPCORE
      REAL, DIMENSION(NDDIM,5) :: XJL_PLOTDATA
      REAL, DIMENSION(NDDIM,4) :: XJC_PLOTDATA_I

C***  Diagonal Weights for Continua
      REAL, DIMENSION(NDDIM,NFDIM) :: WJC

C***  Variables to handle Fine Integration in Steal
C***  FF_INFO : Information on the Frequency grid
C***         (1) : XLAM0
C***         (2) : ALN
C***         (3) : XLAM_FINE_START
C***         (4) : XLAM_FINE_END
C***               (3) + (4) as demanded by the Input Cards
C***         (5) : Start index of Integration
C***         (6) : End index of Integration
C***         (7) : IFF_N Total number of Fine Frequencies
C***         (8) : XLAM_FINE_START
C***         (9) : XLAM_FINE_END
C***               (8) + (9) wavelength of first and last fine point
C***        (10) : KSPACE
C***  BFF_ACT : True, if actual K is in Fine Interval
C***  IFF_DK stores in a compressed way the Delta-K values
C***    IFF_DK can have values between -100 and 100
      REAL, DIMENSION(10) :: FF_INFO
C      INTEGER, DIMENSION(IFF_MAX) :: IFF_DK, IFF_WCHARM
      INTEGER (KIND=TINYINT), DIMENSION(IFF_MAX) :: IFF_DK
      INTEGER (KIND=TINYINT), DIMENSION(IFF_MAX) :: IFF_WCHARM
      REAL, DIMENSION(IFF_MAX) :: WCHARM_FINE
      
C***  Gamma History
C***  JOBNUM, PROGRAM(0,1), 4x(DEPTH,LEVEL,CORR), MAXCORR, 
C***  DEPTH(T), MAXCORR(T), NEWTON:AVER,MAX,CANCELLATIONS,NR. NOT CONV.
      REAL, DIMENSION(26,MAXGAHIST) :: GAHIST

C***  Auto Gamma
      REAL, DIMENSION(7) :: AG

      CHARACTER(10) :: TIM1

      INTEGER :: ITWARN, ITMAX
      COMMON / COMITWA / ITWARN, ITMAX
 
C***  Operating system:
      CHARACTER(8) :: OPSYS
      COMMON / COMOS / OPSYS

C***  Link data to identify program version
      CHARACTER(30) :: LINK_DATE 
      CHARACTER(10) :: LINK_USER 
      CHARACTER(60) :: LINK_HOST
      COMMON / COM_LINKINFO / LINK_DATE, LINK_USER, LINK_HOST

C***  File and channel handles (=KANAL)
      INTEGER, PARAMETER :: hOUT = 6        !write to wruniqX.out (stdout)
      INTEGER, PARAMETER :: hCPR = 0        !write to wruniqX.cpr (stderr)
      INTEGER, PARAMETER :: hMODEL = 3      !write to MODEL file
      INTEGER, PARAMETER :: hHIST = 21      !write to MODHIST file

      INCLUDE 'interfacebib.inc'

C***  Initialization of variables; otherwise, optimized compilation 
C***    might not allocate memory for them unless compiled with "save"
      XLOGL = -99.
      CORMAX = -99.
      HTOTND = 0.
      CONVERG = .FALSE.
      FLUXERR = .0
      BGFIN = .FALSE.
      IHSSTATUS = -1
C***  Set Counters for Warnings
      IWARN_NEG_XJCAPP = 0
      IWARN_NEG_XJLAPP = 0


      bHLPHIST = .FALSE.
      bTauUpdated = .FALSE. !Default setting: Tau was not updated (can be overwritten bei ENSURETAUMAX() )
      bHydroDone = .FALSE.  !Default setting: Hydro routine was not done (overwritten by HYDROSOLVE() )
      bHydroHelp = .FALSE.
      bSaveTauCont = .FALSE.
      bBackupNow = .FALSE.

C***  Write Link Data (Program Version) tp CPR file
      WRITE(hCPR,'(2A)') '>>> STEAL started: Program Version from ', 
     >      LINK_DATE
      WRITE(hCPR,'(4A)') '>>> created by ', LINK_USER(:IDX(LINK_USER)),
     >      ' at host ', LINK_HOST(:IDX(LINK_HOST))

      CALL INSTALL
                                 
      IF (OPSYS .EQ. 'CRAY' .OR. OPSYS .EQ. 'SGI') THEN
        CALL CLOCK(TIM1)
      ELSE
        CALL TIME(TIM1)
      ENDIF

      !Put dummy stuff into COMMON block data to ensure storage reservation
      NDDUMMY = NDDIM
      bNewVelo = .FALSE.
      DO L=1, NDDIM
        OLDRADI(L) = 1337.
        OLDVELO(L) = 1338.
      ENDDO
      WFELOW = 0.
      WFENUP = 0.

C***  DEFINE CHANNEL 2 FOR PLOTS
      KANAL = 2
      OPEN (KANAL, FILE='PLOT', STATUS='UNKNOWN')

C***************************************************************************
C***  XMAX = LINE BANDWIDTH IN DOPPLER UNITS, USED FOR:                  ***   
C***   A) MAX. LINE CORE WIDTH;    B) OVERLAP CRITERION.                 ***
C***  IF THE OVERLAP OPTION IS ACTIVE, XMAX SHOULD BE CHOSEN             ***
C***  CONSISTENTLY WITH PROGRAM CMF                                      ***
      XMAX = 4.5
C***  NFL = NUMBER OF LINE FREQUENCIES (USED FOR LINE-CORE INTEGRATION)  ***
      NFL = 31
C*************************************************************************** 

C***  Default: negative DRDRIN indicates that DTDR at inner boundary has not yet been calculated
      DTDRIN = -99.

      BPRICORR = .TRUE.
      bForceColiPP = .FALSE.

C***  READING THE ATOMIC DATA FROM FILE DATOM
      CALL       DATOM (NDIM,N,LEVEL,NCHARG,WEIGHT,ELEVEL,EION,MAINQN,
     $                  EINST,ALPHA,SEXPO,
     $                  ADDCON1, ADDCON2, ADDCON3, 
     $                  IGAUNT,COCO,KEYCBB,ALTESUM,
     $                  INDNUP,INDLOW,LASTIND,MAXIND,MAXATOM,NATOM,
     $                  ELEMENT,SYMBOL,NOM,KODAT,ATMASS,STAGE,
     $                  SIGMATHK,SEXPOK,EDGEK,NFIRST,
     $                  NLAST,NAUTO,MAXAUTO,LOWAUTO,WAUTO,EAUTO,AAUTO,
     $                  IONAUTO,KRUDAUT,KONTNUP,KONTLOW,LASTKON,MAXKONT,
     $                  IONGRND,KODRNUP,KODRLOW,LASTKDR,MAXKODR,KEYCBF,
C***  IRON: ADDITIONAL PARAMETERS FOR IRON-GROUP LINE BLANKETING
     >                  'STEAL', INDEXMAX, NFEREADMAX, MAXFEIND,
     >                  LASTFE, SIGMAFE, INDRB, INDRF,
     >                  IFENUP, IFELOW, IFRBSTA, IFRBEND, FEDUMMY, 
     >                  VDOPFE, DXFE, XLAM0FE, SIGMAINT, BFEMODEL)
 
C***  READING OF THE MODEL FILE ----------------------------------------
C***  IRON: READIN XJF'S
      LAINDHI = LASTIND + NAUTO
      CALL   REMOST (ND,NDDIM,NP,NPDIM,P,Z,
     >               ENTOT,TOLD,NF,NFDIM,XLAMBDA,FWEIGHT,RADIUS,
     >               RSTAR,RCON,VELO,GRADI,VDOP,INCRIT,VCRIT,TOTOUT,
     >               TOTIN,EMFLUX, MAXAUTO,LASTIND,INDNUP,INDLOW,KEY,
     >               RNE,MODHEAD,JOBNUM,XJC,XJL,MODHIST,MAXHIST,LAST,
     >               EINST,NDIM,ABXYZ,NATOM, TEFF, GLOG, GEFFLOG,
     >               bFixGEFF, GEDD, bGEddFix, ENTOTDENS, XMDOT, 
     >               Rcritical, GEDDRAD, XMSTAR, NAUTO, DRXJL, KRUDAUT, 
     >               MAXXDAT, XDATA, HTOTL, BHTOTERR, HTOTND, DTDRIN,
     >               DENSCON, FILLFAC, 
     >               LASTFE, ARAD, bNoARAD, WFELOW, WFENUP, NCOLIP,
     >               FTCOLI, EMCOLI, GAHIST, 
     >               MAXGAHIST, MAXIND, MAXINDE, MAXKONT, LASTKON,
     >               XJTOTL, XKTOTL, XNTOTL, WJC, GF, OPARND, 
     >               OPASMEAN, OPASMEANTC, QFJMEAN, OPAJMEAN, SMEAN,
     >               OPAJMEANTC, OPAPMEAN, QOPAHMEAN,
     >               EDDIHOUTJMEAN, XLAMBDA2, NF2, HTOTOUTMINUS, 
     >               FF_INFO, IFF_DK, IFF_MAX, IFF_N_MS, 
     >               OPALAMBDAMEAN, TAUROSS, TAUROSScont, VMIN, N, 
     >               POPNUM, SIGMAKI, ATMASS, ACONT, ATHOM, OPALROSS,
     >               RHO, XMU, ALPHAF, NEXTHYDRO, LASTHYDRO, VTURB,
     >               GRSTATIC, LASTBACKUP, LASTTAU, GEFFKEY, ZERO_RATES)      
     
      WRITE(hCPR,'(A,I7)') '>>> This is job number ', JOBNUM
      IF (TAUROSScont(ND) < 0.) THEN
        bSaveTauCont = .TRUE.    !save TAUROSScont into MODEL if calculated this time
      ENDIF
      IF (NEXTHYDRO == -13) THEN
        bHydroHelp = .TRUE.
      ENDIF

C***  DECODING INPUT DATA ******************************************
      CALL DECSTE (LSRAT, LSPOP, JOBMAX, EPSILON, REDUCE, IHIST,        
     >     IFRRA, ITORA, IPRICC, IPRILC, LSEXPO, LSTAU, IFLUX,          
     >     IDAT, LEVELPL, MAXSETS, NPLOT, NDIM, NEWWRC, LASTIND, 
     >     NGAMC,NGAMR,NGAML,NGAMD,AGAMC,AGAMR,AGAML,AGAMD,DELTAC,
     >     LINE, MAXIND, TPLOT, TPLOTOPT, JOBNUM, NOEXTRAP, MODHIST,
     >     STHLP,NODATOM,NSCHAR,NOLAP,ITBR,ITMAX,OLDSTART,COMPO,
     >     ELEMENT,SYMBOL,NATOM,NATOUT,BRRESET, 
     >     DRNORUD,NOPOP, BAUTO, BAUTO_ABORT, SMALLPOP, POPMIN, 
     >     BINBOX, POPRANG, 
     >     COREX, BCOREX, VPLOT, BUNLU, UNLUTECLINE, BPRIUNLU,  
     >     ND, IPRINTZERORATES, PRILEVRA,
     >     LEVELPLDEP, NPLOTDEP,  
     >     BITCONT, bTDIFFUS, BPGAHIST, AG, BAG, BPGAHISTE, 
     >     PLOTOPT, MAXPLOTOPT, NPLOTOPT, OPC, BTALTER, 
     >     BPLOCC, LPLOCC, KPLOCC, BRUDZERO, 
     >     NLINE, LINEINDEX, 
     >     BXJLAPPNEW, BXJCAPPNEW, BNEWOPER, IPLOT_XJLAPP, IPLOT_XJCAPP, 
     >     LPLOT_XJCAPP, NITER_PLOT_JAPP, BPLOTAPP, BNEWTONRESET, 
     >     XLAM_FINE_START, XLAM_FINE_END, BTRACE_POPNUM, 
     >     BGAMMACFIX, BGAMMALFIX, BGAMMARFIX, BGAMMADFIX, 
     >     BNUEFE, BXJLAPPCORE, WJCMIN, NKONV_THRESH,
     >     iBLOCKINVERSION, RSTAR, RMAX, TAUMAX, bENSURETAUMAX, TAUACC,
     >     bTauStrict, ReduceTauCorrections, TauCorLimits, FQLIMIT,
     >     RadiusGridParameters, VMIN, bThinWind, ThinCard, bSUMMARY,
     >     bHYDROSOLVE, bLateTau, HydroCard, AlphaCard, 
     >     DENSCON_FIX, DENSCON_LINE, MFORM, XMSTAR, GLOG, bUpdateMass,
     >     WRTYPE, MODHEAD, bModHeadUpdate, bOLDVELO, bTauMaxSafe,
     >     bFULLHYDROSTAT, bGAMMARADMEAN, GEDDreduce, iZRType, iAMT,
     >     CLUMP_SEP, MacroCard, OPALINE_SCALE, NBACKUP, GEFFKEY, bUCPP,
     >     FLUXEPS, HYSTACC, IHSSTATUS, bNoFeTCORR, bUseTWOPNT, 
     >     bRELSCHARMERACC)

      bBLOCKINVERSION = (iBLOCKINVERSION /= 0)  !temporary line until full ion split is implemented

C***  COLI++ after this STEAL job forced by CARDS option     
      IF (bUCPP) THEN
        bForceColiPP = .TRUE.
      ENDIF
      
cccc  Only for historical reasons, the variable NOTEMP still exists 
ccc   as a parameter for some subroutines. It is always: 
      NOTEMP=.TRUE.

      !Avoid temperature corrections right after hydro step
      IF (BUNLU .AND. ((LASTHYDRO == 1) 
     >                  .OR. (bLateTau .AND. LASTHYDRO == 2))) THEN
        BUNLU = .FALSE.
      ENDIF
      
C***  TCORR ALTERNATE implementation (from Goetz, slightly modified)
      IF (BUNLU .AND. BTALTER) THEN
        IF (MODULO(NCOLIP,2) /= 0) THEN
          BUNLU = .FALSE.
          WRITE (hCPR,*) 'NO T-CORRECTIONS (ALTERNATING)'
        ENDIF
      ENDIF

C***  The option RUDLINES_ZERORADRATES sets all f=0 for all rudimental
C***  lines. This might be usefull in order to achieve optimum consistency
C***  between radiation transfer and statistical equations
      IF (BRUDZERO) THEN
         DO LOW=1, N
            DO NUP=1,N
               IF (EINST(LOW,NUP) .EQ. -2.) EINST(NUP,LOW)=.0
            ENDDO
         ENDDO
      ENDIF

      IF (JOBNUM .GT. 3 .AND. JOBNUM .LT. 11 .AND. BITCONT) THEN
        STHLP = .TRUE.
        bHLPHIST = .TRUE.
        CALL JSYMSET ('G1','WRCONT')
      ENDIF
      IF (JOBNUM .GT. 3 .AND. JOBNUM .LT. 30 .AND. BITCONT) THEN
        NEWWRC = 1
      ENDIF

C***  Move Entries in GAMMA HISTORY
      IF (.NOT. STHLP) THEN
        DO I=MAXGAHIST-1, 1, -1
          DO J=1, 26
            GAHIST(J,I+1) = GAHIST(J,I)
          ENDDO
        ENDDO
      ENDIF

C***  This Job is STEAL
      GAHIST(2,1) = 0.
      
C***  Get last Correction
      IF (GAHIST(2,2) /= 0.) THEN
        LST = 3
      ELSE
        LST = 2
      ENDIF
      IF (GAHIST(19,LST) > 0.) THEN
        CORRLAST = GAHIST(19,LST)
      ELSE
        CORRLAST = 1.
      ENDIF
      CORRLAST = MIN(CORRLAST, 1.)

C***  Preparations for Fluxcorrection in TEMPCORR (UnLu-Method)
      CALL INITFCORR (TEFF, RSTAR, TOLD, RNE, bNoARAD,
     >                ABXYZ, ATMASS, ELEMENT, NATOM, WRTYPE,
     >                RADIUS, VELO, GRADI, VTURB, ENTOT, ND,
     >                HTOTM, HTOTG, HTOTOBS, HTOTL, ACONT, ATHOM,
     >                AGRAV, AMECH, ARAD, APRESS, WORKRATIO,
     >                CLUMP_SEP, OPALROSS, MacroCard, DENSCON,
     >                FILLFAC, OPALINE_SCALE,
     >                FTCOLI, XJTOTL, XKTOTL, XNTOTL, HTOTCMF0, 
     >                XMSTAR, MFORM, GLOG, ATMEAN, VMACH,
     >                TAUROSS, QIONMEAN, GAMMARADMEAN, RCON)

C***  Calculate ABS(DTDR) at inner boundary (e.g. for use in COLI)
C***  This works best if one full COLI was already done, since then the full
C***  opacity at the inner boundary is given. 
C***  If this is not the case, i.e. OPARND is not jet calculated, this
C***  value is approximated by the Rosseland continuum opacity
      IF (OPARND <= 0.) THEN      
        DO J=1,N
          EN(J)=POPNUM((J-1)*N+L,1)
        ENDDO
C***    Calculate OPAROSS (aka continuum opacity) with clump density
        CALL OPAROSS (OPARND, EN, TOLD(ND), RNE(ND), ENTOTDENS(ND),
     >                RSTAR, NDIM, N,
     >                LEVEL, NCHARG, WEIGHT, ELEVEL, EION, EINST,
     >                ALPHA, SEXPO, ADDCON1, ADDCON2, ADDCON3, 
     >                IGAUNT, NF, XLAMBDA, FWEIGHT, NOM,
     >                MAXATOM, SIGMATHK, SEXPOK, EDGEK, KODAT, 
     >                RADIUS(ND), KONTNUP, KONTLOW, LASTKON)
C***    Scale opacity down with filling factor
        OPARND = OPARND * FILLFAC(ND)
      ENDIF

C***  Start approximation for the temperature gradient at ND
C***  This should not be calculated before an OLD MODEL might have been adapted
C***  The first coli will calculate this gradient if still missing
      IF (DTDRIN < .0 .AND. JOBNUM .GT. 1) THEN
         DTDRIN = 0.1875 * OPARND * TEFF**4 / TOLD(ND)**3
      ENDIF

C***  UNSOELD-LUCY-PROCEDURE 
      bFeTCORR = .FALSE.
      IF (BUNLU) THEN
        CALL PREPKUBAT(ND, NF, N, NDIM, XLAMBDA, ENTOTDENS, RNE, 
     >                 POPNUM, ELEVEL, NCHARG, XJC, KONTNUP,
     >                 KONTLOW, SIGMAKI, FWEIGHT, WEIGHT, NOM, 
     >                 ABXYZ, NFIRST, NLAST, INDNUP, INDLOW,
     >                 KODAT, ENLTE, CRATE, DCRATEDT,
     >                 NATOM, EION, EINST, COCO, KEYCBB,
     >                 KEYCBF, IONGRND, ALTESUM, MAXIND,
     >                 SIGMATHK, SEXPOK, EDGEK, MAXATOM,
     >                 LASTKON, LASTIND, OPALROSS, TOLD, RADIUS,
     >                 RSTAR, OPALAMBDAMEAN, OPASMEANTC, 
     >                 DTKUBAT, BCOLLIDONE )

        IF (.NOT. bNoFeTCORR) THEN
C***      Added in Feb 2018 to wrh branch:        
C***      Modify Fe rates to account for current temperature changes       
C***       and switch off ALO for Fe transitions in this case
C***     (old behaviour can be forced via CARDS line "TCORR FERATES OFF")
          bFeTCORR = .TRUE.
          GAMMAD = 0.
          DO MG=1, 10
            AGAMD(MG) = 0.
          ENDDO
        ENDIF
        CALL TEMPCORR (TOLD, TNEW, ND, RADIUS, TEFF,
     >    HTOTL, XJTOTL, HTOTCMF0, FTCOLI, 
     >    DTLOCAL, DTINT, DTRMAX, FLUXERR,
     >    OPASMEANTC, QFJMEAN, OPAJMEANTC, OPAPMEAN, 
     >    QOPAHMEAN, EDDIHOUTJMEAN,
     >    HTOTOUTMINUS, HTOTND, OPARND, DTDRIN, UNLUTECLINE, 
     >    DUNLU_LOC, DUNLU_INT, DUNLU_RMAX, DUNLU_TB, TBTAU, TAUINT,
     >    OPALAMBDAMEAN, TAUROSS, ARAD, DTKUBAT, BTDIFFUS, HNDCORFAC,
     >    SMEAN, CORRS, FLUXEPS)
     
      ELSE
         DO L=1, ND
            TNEW(L) = TOLD(L)
         ENDDO
      ENDIF

      
      
      
      IF (DRNORUD .AND. (JOBNUM .GT. 1)) THEN
C***  TESTING THE XJL'S FOR DR-TRANSITIONS
C***  (DRXJL(INDDR)=.TRUE. AS DEFAULT FOR THOSE LINES WHICH ARE 
C***  SPECIFIED AS RUDIMENTALS IN THE DATOM-FILE)
         DO 10 INDDR = 1,NAUTO
            IF (.NOT. DRXJL(INDDR)) THEN
               WRITE (6,*) 'STEAL> Fatal Error: Radiation Field of ',
     >         'non-rudimental DR-Line not found in Model-File'
               STOP 'ERROR IN STEAL'
            ENDIF
   10    CONTINUE
      ENDIF
      IF (.NOT.DRNORUD) LAINDHI = LASTIND

      DO 1 MG=1,10
      IF (NGAMC(MG).LE.JOBNUM) GAMMAC=AGAMC(MG)
      IF (NGAMR(MG).LE.JOBNUM) GAMMAR=AGAMR(MG)
      IF (NGAML(MG).LE.JOBNUM) GAMMAL=AGAML(MG)
      IF (NGAMD(MG).LE.JOBNUM) GAMMAD=AGAMD(MG)
    1 CONTINUE


C***  HANDLING OF LASER LINES: INFORMATION TRANSFERRED FROM MAIN 
C***  PROGRAM "CMF" IN FILE "LASER_LINES"
      NLASER=0
      OPEN (51, FILE='LASER_LINES', STATUS='UNKNOWN')
   50 READ (51,'(A4,I6)',END=55) KODE,IND
      IF (KODE .EQ. 'LINE') THEN
        LINE(IND)=.FALSE.
        NLASER=NLASER+1
      ENDIF
      GOTO 50
   55 CLOSE (51)
      IF (NLASER .GT. 0) THEN
         WRITE (hCPR,56) NLASER
         WRITE (*,56) NLASER
   56    FORMAT (/,' STEAL> ZERO LINE CORE FOR ',I6,' LINES (LASER ',
     >           'CONDITION FROM >>>COLI<<<)')
      ENDIF

      !Ensure initial values for possibly unused old POPNUM arrays
      POP1 = 0.
      POP2 = 0.
      POP3 = 0.
      IF (STHLP .OR. OLDSTART .OR. (bLateTau .AND. LASTHYDRO == 1)) THEN
C***     READING THE POPNUMBERS OF THE LAST THREE ITERATIONS
         IERR = 1                                                        
         CALL READMS (hMODEL,POPNUM,ND*N,'POPNUM  ', IERR)
         IERR = 1                                                        
         IF (JOBNUM .GT. 5) THEN
           CALL READMS (hMODEL,POP1,ND*N,'POP1    ', IERR)
         ENDIF
         IF (IERR .EQ. -10) THEN
           BPRICORR = .FALSE.
         ENDIF
         IERR = 1                                                        
         IF (JOBNUM .GT. 8) THEN
           CALL READMS (hMODEL,POP2,ND*N,'POP2    ', IERR)
         ENDIF
         CALL READMS (hMODEL,POPLTE,ND*N,'POPLTE  ', IERR)
      ELSE
C***     READING THE LAST-ITERATION POPNUMBERS
         IERR = 1                                                        
         CALL READMS (hMODEL,POP1,ND*N,'POPNUM  ', IERR)
         IF (IERR .EQ. -10) THEN
           BPRICORR = .FALSE.
         ENDIF
         IERR = 1                                                        
         IF (JOBNUM .GT. 5) THEN
           CALL READMS (hMODEL,POP2,ND*N,'POP1    ', IERR) 
         ENDIF
      ENDIF
       
      IF (BUNLU) THEN
         GAHIST(26,1) = 2.
      ELSE
         GAHIST(26,1) = 0.
      ENDIF

C*** AUTO-GAMMA Artistics
      IF (BAG) 
     >  CALL ADJGAMMA(GAHIST, MAXGAHIST, AG, LASTHYDRO, BTALTER,
     >                GAMMAC, GAMMAL, GAMMAR, GAMMAD, BGFIN, GF, 
     >                BGAMMACFIX, BGAMMALFIX, BGAMMARFIX, BGAMMADFIX)

C***  IN CASE OF JOBNUM=1 (STARTJOB) OR OUTPUT REQUEST ONLY (STHLP),
C***  GAMMA'S ARE SET ZERO:
      IF ((JOBNUM .EQ. 1) .OR. STHLP) THEN
         GAMMAR= .0
         GAMMAL= .0
         GAMMAC= .0
         GAMMAD= .0
         WRITE (hCPR,*) 
     >       "STEAL: JOBNUM=1 or output-only -> All GAMMAs = 0"
      ENDIF
 
C***  SKIP ANY CALCULATIONS IN OUTPUT-ONLY MODE OR OLDSTART-INITIALIZATION
      IF (STHLP .OR. OLDSTART .OR. (bLateTau .AND. LASTHYDRO == 1)) THEN
        WRITE (hCPR,*) "This STEAL is a short one (skipping LINPOP)"
        IF (.NOT. STHLP) GOTO 20
      ENDIF

C***  IN CASE OF JOBNUM=1 (STARTJOB):
C***  CALCULATION OF POPNUMBERS WITH THE USUAL (LINEAR) RATE EQUATION
C***  AS START APPROXIMATION FOR SUBR. LINPOP
C***  Note: this branch is also invoked by PRINT RATES
C***  In this case, the wruniq chain is automatically stopped, and the
C***  MODEL file is not updated. 

      IF (JOBNUM .EQ. 1 .OR. LSRAT .NE. -1) THEN
cc      IF (JOBNUM .EQ. 1 .OR. (LSRAT .EQ. -1 .AND. GAMMAR.EQ.0. .AND.
cc     $    GAMMAL.EQ.0. .AND. GAMMAC.EQ.0. .AND. GAMMAD.EQ.0.)) THEN

          CALL POPZERO (TNEW, RNE, POPNUM, DEPART, ENTOTDENS, ITNE, 
     >                  NDIM, N, ENLTE,
     $                  WEIGHT,NCHARG,EION,ELEVEL,EN,EINST,LEVEL,
     $                  XLAMBDA,FWEIGHT,XJC,NF,XJL,IFRRA,ITORA,ALPHA,
     $                  SEXPO, ADDCON1, ADDCON2, ADDCON3, 
     $                  IGAUNT, MODHEAD, JOBNUM,
     $                  ND,LSRAT,CRATE,RRATE,RATCO,
     $                  SIGMAKI,ALTESUM,COCO,KEYCBB,NOM,NATOM,ABXYZ,
     $                  KODAT,NFIRST,NLAST,NATOUT,
     $                  NAUTO,MAXAUTO,LOWAUTO,WAUTO,EAUTO,AAUTO,
     $                  RDIEL,RAUTO,IONAUTO,IONGRND,
     $                  INDNUP,INDLOW,LASTIND,
     $                  KONTNUP,KONTLOW,LASTKON,KODRNUP,KODRLOW,LASTKDR,
     $                  KEYCBF,MAXATOM,SIGMATHK,SEXPOK,EDGEK,LAINDHI,
     $                  KRUDAUT, ZERO_RATES, POPMIN,
     >                  PRILEVRA, SCRATCH, SCRATCH(NDIM))

         write (hCPR,*) '***** this was subroutine POPZERO ***'
         IF (LSRAT .NE. -1) THEN
            WRITE (hCPR,*) 'STEAL run for printing RATES only' 
            WRITE (hCPR,*) 'Chain is now terminated' 
            WRITE (hCPR,*) 'MODEL file has not been updated' 
            JOBMAX = -1
ccc            CALL JSYMSET ('G1','STOP')
            GOTO 25
         ENDIF
      ENDIF
      IF (STHLP) GOTO 20

      !-------------------- Start of of LINPOP call -------------------------
cc      IF (.NOT. (LSRAT .EQ. -1  .AND. GAMMAR.EQ.0. .AND. GAMMAD.EQ.0.
cc     >          .AND. GAMMAL.EQ.0. .AND. GAMMAC.EQ.0. ) ) THEN

C***    RATE EQUATION WITH APPROXIMATE RADIATION TRANSFER (SCHARMER-METHOD)
C***    IN THIS BRANCH, PRIRAT MAY ONLY SHOW THE NETTO RATES
C***    CALCULATION OF NEW POPULATION NUMBERS, EL. DENSITY AND DEPARTURE COEFF.

C***    PREPARATION OF THE LINE-OVERLAP TREATMENT
C***    IRON: IRON-LINES ARE OMITTED (-->NBLENDS = 0)
        CALL OVERLAP (IBLENDS,MAXLAP,LASTIND,XLAMZERO,XLAMRED,XLAMBLUE,
     >                XMAX,EINST,NDIM,ELEVEL,INDNUP,INDLOW,NOLAP,VDOP,
     >                NBLENDS,LASTFE)
     
C***    Temperature correction at inner boundary if TDIFFUS is set in CARDS
        IF (bTDIFFUS) THEN
          CALL TDIFFUS(TNEW, RADIUS, ND, TEFF, OPARND, TNDCORR)
        ENDIF

C***    SOLUTION OF THE NON-LINEAR SET OF RATE EQUATIONS
        CALL  LINPOP (TNEW,RNE,ENTOT,ITNE,POPNUM,DEPART,POPLTE,POP1,
     $   N,ENLTE,WEIGHT,NCHARG,EION,ELEVEL,EN,EINST,LEVEL,
     $   FWEIGHT,XJC,NF,XJL,WCHARM,SCOLD,
     $   XJLAPP,
     $   DM, V1,V2,GAMMAC,GAMMAL,EPSILON,
     $   NOTEMP,DELTAC,GAMMAR,IPRICC,IPRILC,MODHEAD,JOBNUM,IFRRA,ITORA,
     $   RADIUS,RSTAR,OPA,ETA,THOMSON,IWARN,MAINPRO,MAINLEV,
     $   VELO,GRADI, VDOP,NFLDIM,PHI,PWEIGHT,
     $   LASTIND,LAINDHI,
     $   SIGMAKI,
     $   ND,LSRAT,CRATE,RRATE,RATCO,ALPHA,SEXPO, 
     $   ADDCON1, ADDCON2, ADDCON3, 
     $   IGAUNT, KEYCBB, NRB_CONT, ZERO_RATES, IPRINTZERORATES, POPMIN,
     $   LINE,ALTESUM,NFEDGE,NOM,NATOM,ABXYZ,KODAT,NFIRST,
     $   SCRATCH,ATEST,BTEST,DTEST,BRS,BRSP,BRY,XMAX,IBLENDS,MAXLAP,
     $   XLAMZERO,ITBR,NBRCANC,NLAST,SIGMAFF,MAXION,NFL,IONGRND,IONAUTO,
     $   NAUTO,MAXAUTO,LOWAUTO,WAUTO,EAUTO,AAUTO,DRRATEN,IADR16,RDIEL,
     $   RAUTO,DRJLW,DRJLWE,DRLJW,NSCHAR,MAXATOM,BETA,PHIL,NBLENDS,
     $   SIGMATHK,SEXPOK,EDGEK,XDATA,
     $   KONTNUP,KONTLOW,LASTKON,KODRNUP,KODRLOW,LASTKDR,KEYCBF,NATOUT,
     $   BRRESET,NOCON,ENOLD,TEFF, KRUDAUT, ELEMENT, 
     $   BAUTO, BAUTO_ABORT, CKONVER, SMALLPOP, BUNLU, 
     >   DENSCON, FILLFAC, ENTOTDENS,
     >   WFELOW, WFENUP, GAMMAD, NKONVER, WJC, 
     >   OPC, BPLOCC, LPLOCC, KPLOCC, KANAL, 
     >   LASTFE, WJCMIN, NKONV_THRESH,
C***  Fine-frequency grid quantitiew for SETXJFINE
     >   SFINE_OLD, SFINE_NEW, MAXFINE, KONTHLP, 
     >   SIGMA1I, NLINE, LINEINDEX, XLAMSOR, XLAMMIN, XLAMMAX,
     >   NUPACT, LOWACT, BLASERL, ETAL, LIND, LINDS,
C***  IRON
     >   VDOPFE, DXFE, XLAM0FE,
     >   MAXFEACT, BFECHECK, BFEWING,
     >   DFEINDR, SIGMAFE, OPAFE, ETAFE, 
     >   INDEXMAX, BFEMODEL, BNUEFE, 
C***
     >   NF2, NDDIM, XKMIN, XKMAX, XKMID, XKRED, 
     >   XKRED_CORE, XKBLUE_CORE, 
     >   BXJLAPPNEW, BXJCAPPNEW, BNEWOPER, BXJLAPPCORE, 
     >   XJL_PLOTDATA, XJC_PLOTDATA_I, 
     >   IPLOT_XJLAPP, IPLOT_XJCAPP, LPLOT_XJCAPP, NITER_PLOT_JAPP, 
     >   BPLOTAPP, PWEIGHTCL, WS, BRDIVFAILED, BNEWTONRESET, 
     >   IWARN_NEG_XJCAPP, IWARN_NEG_XJLAPP,
     >   XJLAPPNEW, 
     >   XLAM_FINE_START, XLAM_FINE_END, RUDLINE, IMAXPOP, 
C***  New Fine-spaced WCHARM handling
     >   IFF_MAX, IFF_MAX, FF_INFO, IFF_DK, IFF_WCHARM, WCHARM_FINE, 
     >   IFF_N_MS,
C***  Split Parameter (if blockwise DM inversion is used or not)
     >   bBLOCKINVERSION, AOFF, CORRLAST, bRELSCHARMERACC,
C***  Settings for ZERO_RATES, BROYDEN TWOPOINT, and AUTO_MODIFY
     >   iZRType, bUseTWOPNT, iAMT,
C***  Fe rate temperature correction
     >   CORRS, DEXFAC, bFeTCORR,
C***  For constant consistency check     
     >   NDIM, MAXIND, NFDIM, MAXFEIND)

cc      ENDIF
      !-------------------- End of LINPOP call -------------------------

C***  CYCLIC CHANGING OF THE ARRAY NAME FOR THE LAST THREE ITERATIONS
      !links: old  -- rechts: new
      CALL CHANGE ('POP3    '  ,'HELP    '  , 3)
      CALL CHANGE ('POP2    '  ,'POP3    '  , 3)
      CALL CHANGE ('POP1    '  ,'POP2    '  , 3)
      CALL CHANGE ('POPNUM  '  ,'POP1    '  , 3)
      CALL CHANGE ('HELP    '  ,'POPNUM  '  , 3)
      CALL CHANGE ('TOLD3   '  ,'THELP   '  , 3)
      CALL CHANGE ('TOLD2   '  ,'TOLD3   '  , 3)
      CALL CHANGE ('TOLD1   '  ,'TOLD2   '  , 3)
      CALL CHANGE ('T       '  ,'TOLD1   '  , 3)
      CALL CHANGE ('THELP   '  ,'T       '  , 3)
 
   20 CONTINUE
C***  REDUCED CORRECTIONS, IF OPTION IS SET
      IF (REDUCE.NE.1. .AND. JOBNUM.GT.1 .AND. .NOT. STHLP)
     $      CALL REDCOR (POPNUM,POP1,ND,N,RNE,NCHARG,REDUCE)
 
C***  PRINTOUT POP. NUMBERS (IF REQUESTED)
      IF (LSPOP.GT.0) CALL PRIPOP (LSPOP,WEIGHT,NCHARG,NOM,TNEW,BUNLU,
     $             ND,N,RNE,ITNE,LEVEL,POPNUM,DEPART,JOBNUM,MODHEAD,
     $             ABXYZ,NFIRST,NLAST,NATOM,NATOUT, SMALLPOP)
 
C***  PRINTOUT RELATIVE CORRECTIONS 
      IF ((JOBNUM .GT. 5 .OR. JOBNUM .GT. 1 .AND. .NOT. STHLP) .AND. 
     >     BPRICORR)
     $      CALL PRICORR (POPNUM, POP1, LEVEL, N, ND, MODHEAD, LSPOP, 
     $                   CORMAX, RTCMAX, JOBNUM, REDUCE, 
     >                   GAMMAC, GAMMAL, GAMMAR, GAMMAD, 
     $                   TOLD, TNEW, EPSILON, DELTAC, SMALLPOP, BUNLU, 
     $                   DUNLU_LOC, DUNLU_INT, DUNLU_RMAX, DUNLU_TB, 
     >                   bTDIFFUS,  TNDCORR, HNDCORFAC, GAHIST, MAXGAHIST, 
     >                   STHLP, IWARN_NEG_XJCAPP, IWARN_NEG_XJLAPP,
     >                   TBTAU, TAUINT)

C*** INHIBIT SMALL OR NEGATIVE POP. NUMBERS (WARNINGS ISSUED!)
      CALL INHIBIT (POPNUM, N, ND, NCHARG, RNE, 
     >              NATOM, ABXYZ, NFIRST, NLAST, POPMIN)

      IF (STHLP .OR. (bLateTau .AND. LASTHYDRO == 1)) GOTO 25

      IF (.NOT. (STHLP .OR. OLDSTART)) THEN
        ITSUM = 0
        ITMAX = 0
        DO L=1, ND
          ITSUM = ITSUM + ABS(ITNE(L))
          IF (ABS(ITNE(L)) .GT. ITMAX) ITMAX = ABS(ITNE(L))
        ENDDO
        RELIT = FLOAT(ITSUM) / FLOAT(ND)
        PRINT 36, RELIT
   36   FORMAT (' STEAL> AVERAGE NUMBER OF N.-R. ITERATIONS:', F6.2)
        IF (ITBR .LT. ITMAX) THEN
          IF (BNEWTONRESET) THEN
            PRINT 38, ITBR, NBRCANC, BRRESET
   38       FORMAT (' STEAL> BROYDEN UPDATE FROM ITERATION', I3,
     >          '  -  CANCELLED AT', I3, ' ITERATIONS',
     >          '  -  BRRESET = ',G15.1 ) 
          ELSE
            PRINT 338, ITBR, NBRCANC
  338       FORMAT (' STEAL> BROYDEN UPDATE FROM ITERATION', I3,
     >          '  -  CANCELLED AT', I3, ' ITERATIONS',
     >          '  -  NO NEWTON-RESETS' ) 
          ENDIF
        ENDIF
        GAHIST(22,1) = RELIT
        GAHIST(23,1) = ITMAX
        GAHIST(24,1) = FLOAT(NBRCANC)
        GAHIST(25,1) = FLOAT(NKONVER)
      ENDIF

C***  Information on Calculation of Approximate Radiation Fields
      IF (BXJCAPPNEW .OR. BXJLAPPNEW) THEN
        IF (BXJCAPPNEW .AND. BXJLAPPNEW) THEN
          WRITE (MESSAGE,'(A)') 'CONTINUA AND LINES'
        ELSEIF (BXJCAPPNEW) THEN
          WRITE (MESSAGE,'(A)') 'CONTINUA'
        ELSEIF (BXJLAPPNEW) THEN
          WRITE (MESSAGE,'(A)') 'LINES'
        ENDIF
        WRITE (MESSAGE2,'(A,F7.1,A,F7.1)')
     >    '    LAMBDA_START=', XLAM_FINE_START, 
     >    ' LAMBDA_END', XLAM_FINE_END
        MESSAGE = MESSAGE( :IDX(MESSAGE)) // MESSAGE2
        WRITE (*,'(A,A)') 
     >    ' STEAL> NEW APPROXIMATE RADIATION FIELDS CALCULATED FOR ',
     >    MESSAGE
      ENDIF

C***  Trace POPNUMs: Write POPNUMs to cpr-file
      IF (BTRACE_POPNUM) THEN
        WRITE (hCPR,'(A)')
     >    'TRACE POPNUMS: L, JOBNUM' 
        DO L=1, ND
          WRITE (hCPR,'(A,I2,1X,I7,1X,300(E8.2,1X))')
     >      'Trace Popnums ', L, JOBNUM, (POPNUM(L,I), I=1, N)
        ENDDO
      ENDIF


      !++++++++ OLD TAUMAX AND WRITMS LOCATION WAS HERE +++++++++


C********************************************************************
C***  PRINT / PLOT OUTPUT OPTIONS
C********************************************************************



C***  ENTRY FOR OUTPUT-ONLY (STEAL-HELP)
   25 CONTINUE

C***  Printout of Model Parameters
      ! (old small version, usually called at start only, use PRINTMODELSUMMARY for full model data)
      IF ((STHLP .OR. bSUMMARY) .AND. (.NOT. bNoARAD)) THEN
         CALL PRINTMODELSUMMARY (MODHEAD, ND, TEFF, TNEW, TAUROSS, 
     >                           RSTAR, XLOGL, XMDOT, RTRANS, VELO, 
     >                           VDOP, DENSCON, FILLFAC, GLOG, GEFFLOG,
     >                           GEDD, RADIUS, RCON, INCRIT, VCRIT,
     >                           VTURB, bThinWind, bHYDROSOLVE, 
     >                           WORKRATIO, RNE, ENTOT, XMU, XMSTAR, 
     >                           GRADI, NATOM, ABXYZ, ATMASS, ELEMENT,
     >                           SYMBOL, bFULLHYDROSTAT, 
     >                           GEDDRAD, ARAD, APRESS, AGRAV, 
     >                           GAMMARADMEAN, QIONMEAN, bFixGEFF)
      ELSEIF (bNoARAD) THEN
        WRITE (hCPR,'(A)') 'WARNING: *** MODEL SUMMARY NOT POSSIBLE: '
     >   // 'RADIATIVE ACCELERATIONS NOT YET ON MODEL FILE ***'
      ENDIF


C***  PLOT OF THE VELOCITY LAW
      IF (VPLOT) THEN
        CALL PLOTV (ND,RADIUS,VELO,MODHEAD,JOBNUM)
      ENDIF

      !++++++++++ CONVERGED PRINT OPTIONS WERE HERE +++++++

C***  PLOT OPTIONS **************************************
   30 CONTINUE

      !TEST plot => shift if successful
C      CALL PLOTRTAU1 (NF, XLAMBDA, K, ND, RADIUS, OPA, MAINPRO,
C     >                    MAINLEV, MODHEAD)
   
C***  COLLECTED PLOT OPTIONS
      DO IPLOTOPT=1, NPLOTOPT
       
C***     Issue comment line that can help for WRplot navigation
         WRITE (KANAL, '(A)') '* ' // 
     >       PLOTOPT(IPLOTOPT)(:IDX(PLOTOPT(IPLOTOPT)))

         CALL SARGV (PLOTOPT(IPLOTOPT),2,ACTPAR)
         
C***     **********
C***     *  JNUE  *
C***     **********
         IF (ACTPAR .EQ. 'JNUE') 
     >   CALL PLOTJNUE (KANAL, MODHEAD, JOBNUM, ND, NF, XLAMBDA,
     >     PLOTOPT(IPLOTOPT), MAXPLOTN, XPLOT, YPLOT, XJC )

C***     ***********
C***     *  JLINE  *
C***     ***********
         IF (ACTPAR .EQ. 'JLINE') 
     >   CALL PLOTJLINE (KANAL, MODHEAD, JOBNUM, ND, 
     >     PLOTOPT(IPLOTOPT), MAXPLOTN, XPLOT, YPLOT,
     >     XJL, LASTIND, INDLOW, INDNUP, ELEVEL, EINST, NDIM)

C***     *********
C***     *  ACC  *
C***     *********
         IF (ACTPAR .EQ. 'ACC') THEN
           CALL PLOTACC (PLOTOPT(IPLOTOPT), AGRAV, AMECH, ARAD, APRESS, 
     >                   ACONT, ATHOM, WORKRATIO, VELO, RADIUS, ND, 
     >                   ATMEAN, RNE, RCON, TNEW, TEFF, RSTAR, XMU,
     >                   XMSTAR, Rcritical, bFULLHYDROSTAT,
     >                   MODHEAD, JOBNUM, KANAL)
         ENDIF

C***     *********************
C***     *  EDDINGTON GAMMA  *
C***     *********************
         IF (ACTPAR == 'GAMMA') THEN
           CALL PLOTGAMMA (PLOTOPT(IPLOTOPT), TAUROSS,
     >                     AGRAV, AMECH, ARAD, APRESS, ACONT, ATHOM, 
     >                     WORKRATIO, VELO, RADIUS, ND, ATMEAN,
     >                     ENTOT, RNE, RCON, TNEW, TEFF, RSTAR, XMU,
     >                     XMSTAR, Rcritical, bFULLHYDROSTAT,
     >                     QIONMEAN, GAMMARADMEAN,
     >                     MODHEAD, JOBNUM, KANAL)
         ENDIF
         
C***     **********
C***     *  UNLU  *
C***     **********
         IF (ACTPAR == 'UNLU' .AND. BUNLU) THEN
           CALL PLOTUNLU (KANAL, PLOTOPT(IPLOTOPT), MODHEAD, JOBNUM, 
     >                    ND, TOLD, TNEW, DUNLU_LOC, DUNLU_TB, bTDIFFUS,
     >                    DTLOCAL, DTINT, DTRMAX, DTKUBAT)
         ENDIF

C***     **********
C***     *  HTOT  *
C***     **********
         IF (ACTPAR .EQ. 'HTOT') 
     >     CALL PLOTHSUM (HTOTL, HTOTM, HTOTG, HTOTOBS, HTOTCMF0,
     >                    ND, MODHEAD, JOBNUM, KANAL, TEFF, BHTOTERR,
     >                    FLUXEPS)

C***     **********
C***     *  FLUX  *
C***     **********
         IF (ACTPAR .EQ. 'FLUX') 
     >      CALL PLOTFLU (NF, XLAMBDA, EMFLUX, EMCOLI, MODHEAD, JOBNUM, 
     >                    KANAL, RSTAR, TOTOUT, PLOTOPT(IPLOTOPT))

C***     **********
C***     *  ALPHA *
C***     **********
         IF (ACTPAR == 'ALPHA')
     >     CALL PLOTALPHA(ND, RADIUS, ALPHAF, MODHEAD, JOBNUM, .FALSE.)
     
C***     ************
C***     *  SIGMAFE *
C***     ************
         IF (ACTPAR == 'SIGMAFE' .AND. JOBNUM >= 11) THEN
           CALL PLOTSIGMAFE(PLOTOPT(IPLOTOPT), MODHEAD, JOBNUM,
     >                      SIGMAFE, INDRB, MAXFEIND, 
     >                      LASTIND, LASTFE, IFRBSTA, IFRBEND, 
     >                      LEVEL, N, INDNUP, INDLOW,
     >                      INDEXMAX, NFEREADMAX,
     >                      VDOPFE, DXFE, XLAM0FE, .FALSE.)
         ENDIF
      ENDDO

C***  FURTHER PLOT OPTIONS, NOT YET COLLECTED (should be modernized!)

C***  POPNUMBER PLOT (IF REQUESTED)
      IF (NPLOT .GT. 0)
     $      CALL PLOTPOP (LEVELPL, NPLOT, N, ND, LEVEL, ENTOT, 
     $                    POPNUM, MODHEAD, JOBNUM, KANAL, MAXSETS,
     $                    BINBOX, POPRANG, NATOM, NFIRST, NLAST, NCHARG,
     >                    SYMBOL)

C***  DEPARTURE COEFFICIENT PLOT (IF REQUESTED)
      IF (NPLOTDEP .GT. 0)
     $      CALL PLOTDEP (LEVELPLDEP, NPLOTDEP, N, ND, LEVEL, ENTOT, 
     $                    DEPART, MODHEAD, JOBNUM, KANAL, MAXSETS,
     $                    BINBOX, POPNUM, POPLTE, NATOM, NFIRST, NLAST, 
     >                    NCHARG, SYMBOL)

C***  DIRECT TRANSFER OF TEMPERATURE STRATIFICATION PLOT (IF REQUESTED)
      IF (TPLOT) THEN
C         IF ((LSTAU .LE. 0) .AND. TPLTAU)
C     $       CALL TAUSCAL (RSTAR,ND,RADIUS,RNE,
C     $            ENTOT,TNEW,POPNUM,NDIM,N,EN,LEVEL,NCHARG,WEIGHT,
C     $            ELEVEL,EION,EINST,ALPHA,SEXPO,
C     $            ADDCON1, ADDCON2, ADDCON3, 
C     $            IGAUNT,NOM,NF,
C     $            XLAMBDA,FWEIGHT,TAUTHOM,TAUROSScont,
C     $            MAXATOM,SIGMATHK,SEXPOK,EDGEK,KODAT,
C     $            KONTNUP,KONTLOW,LASTKON, DENSCON, FILLFAC)
         CALL PLOTT (TPLOTOPT, ND, RADIUS, TAUROSS, TNEW,
     >               MODHEAD, JOBNUM, KANAL, BINBOX, UNLUTECLINE)
      ENDIF


C***  Plot of FS- and approximate radiation fields
      IF (BPLOTAPP) THEN
        CALL PLOTAPP(KANAL, ND, NF,
     >    XJL_PLOTDATA, XJC_PLOTDATA_I, XJC_PLOTDATA_L,
     >    IPLOT_XJLAPP, IPLOT_XJCAPP, LPLOT_XJCAPP, NITER_PLOT_JAPP)
      ENDIF

C********* End of PLOTs *******************

      IF (STHLP) GOTO 26

C     Possible point for new options that were previously only available in wrstart(?)

      NTOLD = 3
      IF ((bENSURETAUMAX .OR. bHYDROSOLVE) .AND. (JOBNUM > 1)) THEN
        CALL READMS (hMODEL,POP3, ND*N, 'POP3    ', IERR)
        CALL READMS (hMODEL,TOLD3,ND,   'TOLD3   ', IERR)
        IF (IERR == -10) NTOLD = 2
        CALL READMS (hMODEL,TOLD2,ND,   'TOLD2   ', IERR)
        CALL READMS (hMODEL,TOLD3,ND,   'TOLD3   ', IERR)
        IF (IERR == -10) NTOLD = 1
      ENDIF

      IF (JOBNUM > 1 ) THEN
          !ensures that the TAUMAX value given in the CARDS file is still reached at the inner boundary
          ! do not do this for the steal job that is run right after WRSTART (JOB = 1) because the
          ! adapter should be finished before this subroutine is called
          VMAX = VELO(1)
          CALL ENSURETAUMAX(bENSURETAUMAX, HYSTACC, IHSSTATUS,
     >              TAUMAX, TAUACC, bTauStrict,
     >              VMIN, VELO, GRADI, RADIUS, XMDOT, GEDDRAD,
     >              RHO, TAUROSScont, RMAX, RCON, TNEW, TEFF,
     >              TOLD, TOLD2, TOLD3, CKONVER, bTauMaxSafe,
     >              LASTTAU, GLOG, GEFFLOG, bFixGEFF, XMSTAR, RSTAR, 
     >              XMU, ENTOT, RNE, ND, NDDIM, NP, NPDIM, P, Z, VTURB,
     >              ARAD, APRESS, AMECH, AGRAV, XJC, bThinWind, 
     >              ThinCard, RadiusGridParameters, DENSCON, FILLFAC, 
     >              DENSCON_FIX, DENSCON_LINE, NFIRST, NLAST, NATOM, 
     >              ABXYZ, ATMASS, GEDD, BUNLU, HTOTL, HTOTCMF0, 
     >              FQLIMIT, bFULLHYDROSTAT, bGAMMARADMEAN, TAUROSS,
     >              CORMAX, TauCorLimits, INCRIT, VCRIT, SCRATCH,
     >              ReduceTauCorrections, GEDDreduce, QIONMEAN,
     >              GAMMARADMEAN, GRSTATIC, bGEddFix, bTauUpdated, 
     >              JOBNUM, POPMIN, XJCorg, NCOLIP, NEWWRC,
     >              NTOLD, DEPARTNDorg,
                    !The parameters hereafter are only needed for TAUSCAL
     >              NDIM,MAXKONT,POPNUM,POP1,POP2,POP3,N,EN,LEVEL,
     >              NCHARG,WEIGHT,ELEVEL,EION,EINST,ALPHA,SEXPO,
     >              ADDCON1, ADDCON2, ADDCON3, 
     >              IGAUNT,NOM,NF,
     >              XLAMBDA,FWEIGHT,TAUTHOM,
     >              MAXATOM,SIGMATHK,SEXPOK,EDGEK,KODAT,
     >              KONTNUP,KONTLOW,LASTKON
     >    )
          bSaveTauCont = .TRUE.
          IF (ABS(VMAX - VELO(1)) > 1.) THEN
            !Forced renewed EDDI file if VINF differs by more than 1 km/s
            bForceColiPP = .TRUE.
          ENDIF
      ENDIF


C***  UPDATING THE MODEL FILE
      CALL WRITMS (hMODEL,POPNUM,ND*N,'POPNUM  ',-1, IDUMMY, IERR)
      CALL WRITMS (hMODEL,RNE,ND,'RNE     ',-1, IDUMMY, IERR)
      CALL WRITMS (hMODEL,TNEW,ND,'T       ',-1, IDUMMY, IERR)
      IF (.NOT.((JOBNUM .EQ. 1) .OR. STHLP)) THEN
        CALL WRITMS (hMODEL,POPLTE,ND*N,'POPLTE  ',-1, IDUMMY, IERR)
        CALL WRITMS (hMODEL,SIGMAKI,NF*LASTKON,
     >                                'SIGMAKI ',-1, IDUMMY, IERR)
      ENDIF

      CALL WRITMS (hMODEL,ENTOT,     ND,  'ENTOT   ', -1, IDUMMY, IERR)
      DO L=1, ND
         RHO(L) = AMU * ATMEAN * ENTOT(L)
      ENDDO
      CALL WRITMS (hMODEL,RHO,       ND,  'RHO     ', -1, IDUMMY, IERR)
      CALL WRITMS (hMODEL,ZERO_RATES,N*ND,'ZERO_RAT', -1, IDUMMY, IERR)

      IF (bSaveTauCont .AND. (TAUROSScont(ND) > 0.)) THEN
        CALL WRITMS(hMODEL,TAUROSScont,ND,'TAURCONT',-1,IDUMMY,IERR)
      ENDIF

      IF (DTDRIN >= 0.) THEN
        CALL WRITMS (hMODEL, DTDRIN, 1, 'DTDRIN  ', -1, IDUMMY, IERR)
      ENDIF
      
      IF (bModHeadUpdate) THEN
        !Update MODHEAD entry if forced via CARDS
        CALL WRITMS (hMODEL,MODHEAD,13, 'MODHEAD ', -1, IDUMMY, IERR)
      ENDIF

      !Save updated TAUROSScont and density contrast
      IF (bTauUpdated .OR. bHydroDone) THEN
        !Save old popnums and temperatures because the radius grid has changed
        ! and these numbers have been interpolated on the new grid
        CALL WRITMS (hMODEL,POP1,ND*N, 'POP1    ',-1, IDUMMY, IERR)
        CALL WRITMS (hMODEL, TOLD, ND, 'TOLD1   ',-1, IDUMMY, IERR)
        IF (NTOLD > 1) THEN
          CALL WRITMS (hMODEL,POP2,ND*N, 'POP2    ',-1, IDUMMY, IERR)
          CALL WRITMS (hMODEL, TOLD2,ND, 'TOLD2   ',-1, IDUMMY, IERR)
        ENDIF
        IF (NTOLD > 2) THEN
          CALL WRITMS (hMODEL,POP3,ND*N, 'POP3    ',-1, IDUMMY, IERR)
          CALL WRITMS (hMODEL, TOLD3,ND, 'TOLD3   ',-1, IDUMMY, IERR)
        ENDIF
        !CALL WRITMS (hMODEL,TAUROSS,ND, 'TAUROSS ', -1, IDUMMY, IERR)
        CALL WRITMS (hMODEL,DENSCON,ND, 'DENSCON ', -1, IDUMMY,IERR)
        !Save updated grid
        CALL WRITMS (hMODEL,RADIUS,ND, 'R       ', -1, IDUMMY, IERR)
        CALL WRITMS (hMODEL,INCRIT,ND, 'INCRIT  ', -1, IDUMMY, IERR)
        CALL WRITMS (hMODEL,VCRIT,ND,  'VCRIT   ', -1, IDUMMY, IERR)
        CALL WRITMS (hMODEL,XMU  ,ND,  'XMU     ', -1, IDUMMY, IERR)
        CALL WRITMS (hMODEL,GEDDRAD,1, 'GEDDRAD ', -1, IDUMMY, IERR)
        CALL WRITMS (hMODEL,GRSTATIC,ND,'GRSTATIC',-1, IDUMMY, IERR)        
      ENDIF

      CALL WRITMS (hMODEL,P,NP,      'P       ', -1, IDUMMY, IERR)
      CALL WRITMS (hMODEL,Z,ND*NP,   'Z       ', -1, IDUMMY, IERR)
      
      !Save updated velocity field
      CALL WRITMS (hMODEL,VMIN,1 ,   'VMIN    ', -1, IDUMMY, IERR)
      IF (bTauUpdated .OR. bHydroDone) THEN
          CALL WRITMS (hMODEL,RCON,1 ,   'RCON    ', -1, IDUMMY, IERR)
      ENDIF
      IF (bHydroDone) THEN
        CALL WRITMS (hMODEL,RCSAVE,1 , 'RCSAVE  ', -1, IDUMMY, IERR)
      ENDIF
      CALL WRITMS (hMODEL,VELO,ND,   'VELO    ', -1, IDUMMY, IERR)
      CALL WRITMS (hMODEL,GRADI,ND,  'GRADI   ', -1, IDUMMY, IERR)

      CALL WRITMS (hMODEL,ATMEAN,1 , 'ATMEAN  ', -1, IDUMMY, IERR)
      IF (bTauUpdated .OR. bUpdateMass) THEN
        !In case of fixed g_eff => save updated Mass ang log g_grav
        IF (bUpdateMass .OR. bFixGEFF) THEN
          CALL WRITMS (hMODEL,XMSTAR,1 , 'XMSTAR  ', -1, IDUMMY, IERR)
          CALL WRITMS (hMODEL,GLOG  ,1 , 'GLOG    ', -1, IDUMMY, IERR)
        ELSEIF (.NOT. bFixGEFF) THEN
          !GEDD fixed, but not GEFF => save new g_eff, but with flexible mark (negative)
          GFLSAV = -1. * GEFFLOG
          CALL WRITMS (hMODEL,GFLSAV,1 , 'GEFFLOG ', -1, IDUMMY, IERR)
        ENDIF
      ENDIF
      
C      IF (bTauUpdated .AND. bFULLHYDROSTAT) THEN
C        NCOLIP = -1     !force COLI++ (i.e. new EDDIS) due to changes and a_rad-dependency
C        CALL WRITMS(hMODEL, NCOLIP, 1,'NCOLIP  ', -1, IDUMMY, IERR)
C      ENDIF

      IF (NEXTHYDRO > 0) THEN
        !If used, update hydro iteration counter
        CALL WRITMS (hMODEL,NEXTHYDRO,1,'NXTHYDRO',-1,IDUMMY,IERR)
      ENDIF

      IF (bForceColiPP) THEN
        NCOLIP = -1     !force COLI++ (i.e. new EDDIS) due to new RADIUS GRID
        CALL WRITMS(hMODEL, NCOLIP, 1,'NCOLIP  ', -1, IDUMMY, IERR)
      ENDIF
      
      IF (bHydroDone) THEN
        !This is only written after a successful hydro iteration
        CALL WRITMS (hMODEL,LASTHYDRO,1,'LSTHYDRO',-1,IDUMMY,IERR)
        CALL WRITMS (hMODEL,XMDOT,1 , 'XMDOT   ', -1, IDUMMY, IERR)        
        !update continous radiation field (has been interpolated)
        DO K=1, NF
          KV=1+ND*(K-1)
          WRITE (UNIT=NAME, FMT='(A3,I5)') 'XJC',K
          CALL WRITMS(hMODEL,XJC(KV,1),ND,NAME,-1, IDUMMY, IERR)
        ENDDO
        !force recalculation of inner boundary opacity
        OPARND = 0.
        CALL WRITMS(hMODEL, OPARND, 1,'OPARND  ', -1, IDUMMY, IERR)        
        CALL WRITMS(hMODEL, OPARND, 1,'OPARNDM ', -1, IDUMMY, IERR)
C          //..
      ELSEIF (LASTHYDRO > -1) THEN
        !If used, update post hydro iteration counter (+1 is in REMOST)
        CALL WRITMS (hMODEL,LASTHYDRO,1,'LSTHYDRO',-1,IDUMMY,IERR)
      ENDIF
      IF (LASTTAU >= 0) THEN
        CALL WRITMS (hMODEL, LASTTAU, 1, 'LASTTAU ', -1, IDUMMY, IERR)
      ENDIF

      IF (LASTBACKUP >= NBACKUP) THEN
        bBackupNow = .TRUE.
        LASTBACKUP = 1
      ELSE
        LASTBACKUP = LASTBACKUP + 1
      ENDIF
      CALL WRITMS (hMODEL,LASTBACKUP,1, 'LASTBAK ', -1, IDUMMY, IERR)
      
C***  UPDATING THE MODEL HISTORY
      CALL       STHIST (MODHIST, LAST, MAXHIST, GAMMAC, GAMMAL, GAMMAR, 
     >                   GAMMAD, MODHEAD, JOBNUM, CORMAX, RTCMAX,
     $                   REDUCE, NSCHAR, BUNLU, DUNLU_LOC, DUNLU_INT, 
     >                   DUNLU_RMAX, DUNLU_TB,
     >                   BXJLAPPNEW, BXJCAPPNEW, bBLOCKINVERSION,
     >                   XLAM_FINE_START, XLAM_FINE_END, bHydroDone,
     >                   bTauUpdated)

      IF (.NOT. bHLPHIST) THEN
        !write model history entry into explicit history file
        CALL GETHISTENTRY(HISTENTRY,JOBNUM,MODHIST,MAXHIST)
        OPEN (hHIST, FILE='MODHIST', STATUS='UNKNOWN',
     >              ACTION='READWRITE', POSITION='APPEND')
        WRITE (hHIST,FMT='(A)') TRIM(ADJUSTL(HISTENTRY))
        CLOSE(hHIST)      

        CALL WRITMS (hMODEL,MODHIST,MAXHIST,'MODHIST ',-1, IDUMMY, IERR)
      ENDIF

C***  Store GAMMA HISTORY
      CALL WRITMS(hMODEL,GAHIST,26*MAXGAHIST,'GAHIST  ',-1,IDUMMY,IERR)

C***  Set GF to zero
      GF = 0.
      CALL WRITMS (hMODEL, GF, 1, 'GF      ',-1, IDUMMY, IERR)


C***  DECIDE UPON THE NEXT JOB
      CALL NEXTJOB (JOBNUM, JOBMAX, MODHIST, LAST, CORMAX, EPSILON,
     >              NEWWRC, MOREJOBS, CONVERG, NOEXTRAP, NOCON, 
     >              BPRICORR, COREX, BCOREX, NCOLIP, BAG, BGFIN, 
     >              BAUTO_ABORT, FLUXEPS, ND, FLUXERR, IHSSTATUS)

C***  SECOND ENTRY FOR OUTPUT-ONLY (STEAL-HELP)
  26  CONTINUE

      IF (bHLPHIST) THEN
        !write history entry also for short steal (i.e. at the first jobs)
        WRITE(UNIT=BUFFER16, FMT='(A1,I7,A8)') '/', JOBNUM, '. STEAL '
        CALL ADDHISTENTRY(MODHIST,-1,MAXHIST,16,BUFFER16)
        WRITE(UNIT=BUFFER8,FMT='(A8)') '<BSHORT>'
        CALL ADDHISTENTRY(MODHIST,-1,MAXHIST,8,BUFFER8)

        !write model history entry into explicit history file
        CALL GETHISTENTRY(HISTENTRY,JOBNUM,MODHIST,MAXHIST)
        OPEN (hHIST, FILE='MODHIST', STATUS='UNKNOWN',
     >              ACTION='READWRITE', POSITION='APPEND')
        WRITE (hHIST,FMT='(A)') TRIM(ADJUSTL(HISTENTRY))
        CLOSE(hHIST)      

        CALL WRITMS (hMODEL,MODHIST,MAXHIST,'MODHIST ',-1, IDUMMY, IERR)
      ENDIF
   
      IF (JOBNUM == 1) GOTO 27

C********* final PLOTs and printouts *******************
C      (these plots are considering all changes from TAUMAX, HYDRO and/or NEXTJOB)

C***  THE FOLLOWING PRINTOUT IS ENFORCED IF A MODEL IS FINALLY CONVERGED:
      IF (CONVERG) THEN
         IHIST=1
         IFLUX=1
         IF (.NOT. NODATOM) IDAT=1
         IF (NOPOP) LSPOP=1
         IF (LSPOP.NE.1) THEN
            LSPOP=1
            CALL PRIPOP (LSPOP,WEIGHT,NCHARG,NOM,TNEW,BUNLU,
     >           ND,N,RNE,ITNE,LEVEL,POPNUM,DEPART,JOBNUM,MODHEAD,
     >           ABXYZ,NFIRST,NLAST,NATOM,NATOUT, SMALLPOP)
         ENDIF
         IF (LSTAU <= 0) THEN
           !Tauscal might have never been calculated => failsafe
           CALL TAUSCAL (RSTAR,ND,RADIUS,RNE,
     >                   ENTOT,
     >                   TNEW,POPNUM,NDIM,N,EN,LEVEL,NCHARG,WEIGHT,
     >                   ELEVEL,EION,EINST,ALPHA,SEXPO,
     >                   ADDCON1, ADDCON2, ADDCON3, 
     >                   IGAUNT,NOM,NF,
     >                   XLAMBDA,FWEIGHT,TAUTHOM,TAUROSScont,
     >                   MAXATOM,SIGMATHK,SEXPOK,EDGEK,KODAT,
     >                   KONTNUP,KONTLOW,LASTKON, DENSCON, FILLFAC)
         ENDIF
         LSTAU=1 !print all depth points after model is converged
         CALL PRITAU (LSTAU,MODHEAD,JOBNUM,ND,RADIUS,RNE,
     >                   ENTOT,TNEW,TAUTHOM,TAUROSS,TAUROSScont,
     >                   DENSCON,FILLFAC)
         XLOGL = -99.
         !Print model parameter summary
         IF (.NOT. bNoARAD) THEN         
           CALL PRINTMODELSUMMARY (MODHEAD, ND, TEFF, TNEW, TAUROSS, 
     >                           RSTAR, XLOGL, XMDOT, RTRANS, VELO, 
     >                           VDOP, DENSCON, FILLFAC, GLOG, GEFFLOG,
     >                           GEDD, RADIUS, RCON, INCRIT, VCRIT,
     >                           VTURB, bThinWind, bHYDROSOLVE, 
     >                           WORKRATIO, RNE, ENTOT, XMU, XMSTAR,
     >                           GRADI, NATOM, ABXYZ, ATMASS, ELEMENT,
     >                           SYMBOL, bFULLHYDROSTAT, 
     >                           GEDDRAD, ARAD, APRESS, AGRAV, 
     >                           GAMMARADMEAN, QIONMEAN, bFixGEFF)
         ELSEIF (bNoARAD) THEN
           WRITE (hCPR,'(A)') 'WARNING: *** MODEL SUMMARY NOT POSSIBLE:'
     >       // ' RADIATIVE ACCELERATIONS NOT YET ON MODEL FILE ***'
         ENDIF

         !save final XMU and TAUROSScont values
         CALL WRITMS (hMODEL,XMU  ,ND,  'XMU     ', -1, IDUMMY, IERR)
         IF (TAUROSScont(ND) > 0.) THEN
           CALL WRITMS(hMODEL,TAUROSScont,ND,'TAURCONT',-1,IDUMMY,IERR)
         ENDIF
         IF (bFixGEFF) THEN
           !ensure saving of correct mass in the MODEL file for converged models
           CALL WRITMS (hMODEL,XMSTAR,1 , 'XMSTAR  ', -1, IDUMMY, IERR)
           CALL WRITMS (hMODEL,GLOG  ,1 , 'GLOG    ', -1, IDUMMY, IERR)
         ENDIF
      ENDIF 

      IF (LSEXPO > 0)
     >      CALL PRIEXPO (POPNUM,POP1,POP2,LEVEL,N,ND,MODHEAD,JOBNUM,
     >      LSEXPO )
 
      IF ((.NOT. CONVERG) .AND. (LSTAU .GT. 0)) THEN
          CALL TAUSCAL (RSTAR,ND,RADIUS,RNE,ENTOT,
     >                  TNEW,POPNUM,NDIM,N,EN,LEVEL,NCHARG,WEIGHT,
     >                  ELEVEL,EION,EINST,ALPHA,SEXPO,
     >                  ADDCON1, ADDCON2, ADDCON3, 
     >                  IGAUNT,NOM,NF,
     >                  XLAMBDA,FWEIGHT,TAUTHOM,TAUROSScont,
     >                  MAXATOM,SIGMATHK,SEXPOK,EDGEK,KODAT,
     >                  KONTNUP,KONTLOW,LASTKON, DENSCON, FILLFAC)
          CALL PRITAU (LSTAU,MODHEAD,JOBNUM,ND,RADIUS,RNE,
     >                 ENTOT,TNEW,TAUTHOM,TAUROSS,TAUROSScont,
     >                 DENSCON,FILLFAC)
      ENDIF


C***  THE FOLLOWING PRINT/PLOT-OPTIONS ARE INACTIVE, IF THIS JOB
C***    IS AUTOMATICALLY FOLLOWED BY A SUBSEQUENT ONE!
C      IF (MOREJOBS .AND. .NOT. CONVERG) THEN
C        GOTO 30
C      ENDIF
      IF ((.NOT. MOREJOBS) .OR. (CONVERG)) THEN

C***    PRINTOUT OF MODEL HISTORY (IF REQUESTED)
        IF (IHIST == 1) CALL PRIHIST(MODHIST,LAST,MODHEAD,JOBNUM)

C***    Printout of GAMMA HISTORY
        IF (BPGAHIST) THEN
C          WRITE (hCPR,*) 'Hier Output GAMMA HISTORY'
          CALL PRIGAHIST(GAHIST, MAXGAHIST, LEVEL, N, AG, BAG, 
     >                   BPGAHISTE)
        ENDIF

C***  PRINTOUT OF EMERGENT CONT. FLUX (IF REQUESTED)
        IF (IFLUX == 1) THEN
          IF (TOTOUT == 6HUNDEF. ) THEN
            PRINT 7
    7       FORMAT (//' INVALID PRINT OR PLOT OPTION - ',
     >      'EMERGENT CONT. FLUX NOT YET CALCULATED ',//)
          ELSE
            CALL PRIFLUX (NF,XLAMBDA,EMFLUX,TOTIN,TOTOUT,RSTAR,
     >                      JOBNUM,FWEIGHT,MODHEAD,KEY,EMCOLI)
          ENDIF
        ENDIF

C***  PRINTOUT OF ATOMIC DATA (IF REQUESTED)
        IF (IDAT == 1)
     >    CALL PRIDAT(NDIM,N,LEVEL,NCHARG , WEIGHT,ELEVEL,EION,EINST,
     >                  KEY,NF,ALPHA,SEXPO,
     >                  ADDCON1, ADDCON2, ADDCON3, 
     >                  IGAUNT,SIGMATHK,SEXPOK,EDGEK, MAXATOM,
     >                  COCO,KEYCBB,ALTESUM,
     >                  NATOUT,NATOM,ELEMENT,NOM,ABXYZ,ATMASS,
     >                  NAUTO,MAXAUTO,LOWAUTO,WAUTO,EAUTO,AAUTO,IONAUTO,
     >                  INDNUP,INDLOW,LASTIND,IONGRND,
     >                  KEYCBF,KONTNUP,KONTLOW,LASTKON,
     >                  MAXNDR,ELEVDR,NCHARGDR,INDNUPDR, 
     >                  LEVELDR,NOMDR,IONDR,WEIGHTDR,KRUDAUT) 

C***    PRINTOUT OF CHEMICAL COMPOSITION (IF REQUESTED)
        IF (COMPO)
     >    CALL  PRICOMP (NDIM,EINST,N,NCHARG,NOM,NATOM,ABXYZ,ATMASS,
     >                   STAGE,NFIRST,NLAST,ELEMENT,SYMBOL,LASTIND,
     >                   INDLOW,INDNUP,NAUTO,LOWAUTO,EAUTO,
     >                   KONTNUP,KONTLOW,LASTKON,XMU(1),KRUDAUT)

C***    PRINTOUT of UNLU T-Corrections
        IF (BPRIUNLU .AND. BUNLU) THEN
          CALL PRIUNLU (KANAL, ND, DTLOCAL, DTINT, DTRMAX)
        ENDIF

      ENDIF

C********* End of final PLOTs and printouts *******************


C***  THIRD ENTRY FOR OUTPUT-ONLY (STEAL-HELP)
   27 CONTINUE

      CALL CLOSMS (hMODEL, IERR)

      !Backup Model if enforced (CARDS option)
      IF (NBACKUP > 0 .AND. bBackupNow) THEN
        WRITE(UNIT=JOBSTR, FMT='(I7,X)') JOBNUM
        WRITE(UNIT=COMMAND, FMT='(A,I1,A,A)') 
     >    "cp fort.", hMODEL, ' backup_MODEL', ADJUSTL(JOBSTR)
        CALL SYSTEM(COMMAND)        
      ENDIF
      
                                 
C***  PROGRAM STOP
      CALL JSYMSET ('G0','0')
      CALL STAMP (OPSYS, 'STEAL', TIM1)

      STOP 'O.K.'
      END
