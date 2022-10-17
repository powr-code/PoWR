C***  MAIN PROGRAM WRSTART  ****************************************************
      SUBROUTINE WRSTART
C***********************************************************************
C***  THIS PROGRAM IS TO INITIALIZE THE MODEL FILE FOR SUBSEQUENT
C***  CALCULATION OF THE NON-LTE MULTI-LEVEL LINE FORMATION.
C***  IT MAKES USE OF THE ATOMIC DATA (FILE DATOM)
C***    AND (IF NOT TAKEN FROM OLD MODEL) THE FREQUENCY GRID (FILE FGRID)
C***********************************************************************

      IMPLICIT NONE
      INCLUDE 'interfacebib.inc'
 
C***  DEFINE ARRAY DIMENSIONS
      INTEGER, PARAMETER :: MAXATOM =          26 
      INTEGER, PARAMETER :: NDIM    =         1560 
      INTEGER, PARAMETER :: NFDIM   =  2*NDIM + 400 
      INTEGER, PARAMETER :: MAXIND  =        20000 
      INTEGER, PARAMETER :: MAXFEIND  =       1500 
      INTEGER, PARAMETER :: MAXKONT =      NFDIM/2 
      INTEGER, PARAMETER :: MAXKODR =         NDIM 
      INTEGER, PARAMETER :: NDDIM   =           89 
      INTEGER, PARAMETER :: NPDIM   =           94 
      INTEGER, PARAMETER :: MAXHIST =         4000 
      INTEGER, PARAMETER :: MAXXDAT =           10 
 
C***  MAXIMUM ION CHARGE WHICH MAY OCCUR (SEE ALSO SUBR. GAUNTFF)
      INTEGER, PARAMETER :: MAXION  =   27 

C***  COMMON /VELPAR/ TRANSFERS VELOCITY-FIELD PARAMETERS (HERE USED: VMIN)
      COMMON /VELPAR/ VFINAL, VMIN, BETA, VPAR1, VPAR2, RCON, HSCALE,
     >     BETA2, BETA2FRACTION, VPAR1_2, VPAR2_2
C***  OLD RADIUS AND VELOCITY ARE TRANSFERRED TO FUNCTION WRVEL 
C***     BY SPECIAL COMMON BLOCKS:
      COMMON /COMRADI/ OLDRADI(NDDIM)
      COMMON /COMVELO/ NEWVELO, NDv, OLDVELO(NDDIM)

      CHARACTER(8), DIMENSION(NDDIM) :: INCRIT, VCRIT, VCRITold      
      !INCRIT in common block has been replaced by IDUMP as it is never used there

C***  HANDLING OF DIELECTRONIC RECOMBINATION / AUTOIONIZATION (SUBR. DATOM)
      INTEGER, PARAMETER :: MAXAUTO = 2850
      INTEGER, DIMENSION(MAXAUTO) :: LOWAUTO, IONAUTO, KRUDAUT
      REAL, DIMENSION(MAXAUTO) :: WAUTO, EAUTO, AAUTO

      INTEGER, DIMENSION(NDIM) :: NCHARG, MAINQN, NOM, IONGRND
      INTEGER, DIMENSION(MAXKONT) :: KONTNUP, KONTLOW, NFEDGE
      CHARACTER*8 IGAUNT(MAXKONT), KEYCBF(MAXKONT)
      INTEGER, DIMENSION(MAXATOM) :: KODAT, NFIRST, NLAST
      INTEGER, DIMENSION(NDDIM) :: IDUMP
      INTEGER, DIMENSION(NFDIM) :: KEY
      INTEGER, DIMENSION(MAXIND) :: INDNUP, INDLOW
      INTEGER, DIMENSION(MAXKODR) :: KODRNUP, KODRLOW

      REAL, DIMENSION(NDIM) :: WEIGHT, ELEVEL, EION, ENLTE
      REAL, DIMENSION(MAXKONT) :: ALPHA, SEXPO,
     >                            ADDCON1, ADDCON2, ADDCON3
      REAL, DIMENSION(MAXATOM) :: ABXYZ, ATMASS, STAGE
      REAL, DIMENSION(NDDIM) :: RADIUS, ENTOT, T, TOLD, VELO, GRADI,
     >                          ROLD, RNE, XJC, XJL, TAUROSSOLD, 
     >                          TAURCONT, TAURCONTOLD,
     >                          RHO, XMU, VELOold, RADIUSold, GEFFL,
     >                          ARAD, APRESS, VTEMP, DR, RI, OLDGRADI,
     >                          GAMMARAD, EN, TAUTHOM, TAUGREY, 
     >                          GRSTATIC,
C***  PROVIDE VARIABLES FOR THE READ OF THE SECOND MODEL FOR TEMPERATURE 
C***    INTERPOLATION
     >                          T2, RADIUS2, TOLD2, ROLD2, TEFFOLD2
      REAL, DIMENSION(NPDIM) :: P
      REAL, DIMENSION(NDDIM,NPDIM) :: Z
      REAL, DIMENSION(NFDIM) :: XLAMBDA, XLAMBDA2, FWEIGHT,
     >                          EXPFAC, OPAC, ETAC
      REAL, DIMENSION(NFDIM,MAXKONT) :: SIGMAKI
      REAL, DIMENSION(NFDIM,0:MAXION) :: SIGMAFF
      REAL, DIMENSION(4,NDIM) :: ALTESUM
      REAL, DIMENSION(4,MAXIND) :: COCO
      REAL, DIMENSION(NDIM,NDIM) :: EINST
      REAL, DIMENSION(NDDIM,NDIM) :: POPNUM

      REAL, DIMENSION(MAXXDAT) :: XDATA
      REAL, DIMENSION(MAXATOM,MAXATOM) :: SIGMATHK, SEXPOK, EDGEK

C***  IRON: COMMON BLOCK FOR IRON-SPECIFIC DATA
C***  include "dimblock"
      INTEGER, PARAMETER :: INDEXMAX = 1E7, NFEREADMAX = 3E5    !std
C      INTEGER, PARAMETER :: INDEXMAX = 4E7, NFEREADMAX = 5E5    !vd20
C      INTEGER, PARAMETER :: INDEXMAX = 1E8, NFEREADMAX = 6E5     !xxl

      COMMON /IRON/ FEDUMMY(NFEREADMAX),
     >              INDRB(MAXFEIND),INDRF(MAXFEIND), SIGMAFE(INDEXMAX),
     >              IFRBSTA(MAXFEIND), IFRBEND(MAXFEIND),
     >              IFENUP(MAXFEIND), IFELOW(MAXFEIND),
     >              SIGMAINT(MAXFEIND)

      INTEGER :: N, ND, NDold, JOBNUM, NATOM, LASTIND, IDX, IERR,
     >           NAUTO, IDUMMY, MAXITTAU, ITTAU,
     >           LASTKON, LASTKDR, LASTFE, LAST, IFLAG, L,
     >           NA, NF, NF2, MASSORIGIN, NP, JOBNOLD, JOBNOLD2,
     >           INDRB, INDRF, IFEUP, IFELOW, IFRBSTA, IFRBEND,
     >           NEXTK, IFENUP, LRTinput, NC, NDv, MDOTINPUT

      REAL :: GLOG, GEFFLOG, GEDD, TROLD, VMINOLD, VMINOLD2, FM,
     >        qpLINK_USER, BLACKEDGE, VA, RCSAVE, Vfac, GFLSAV,
     >        SIGMAFE, FEDUMMY, VDOPFE, DXFE, XLAM0FE, SIGMAINT,
     >        TAUMAX, TAUACC, VFINAL, XMASS, ATMEAN, STAPEL, AMIN,
     >        RSTAR, RMAX, VDOP, XMSTAR, TFAC, XMDOT, XLOGL, RTRANS,
     >        DENSCON_FIX, VMIN, TEFF, TEFFOLD, TMIN, TMIN2, TS, RL,
     >        R23, OLDRADI, OLDVELO, TROLD2, RCON, TOTOUT, BETA,
     >        VPAR1, VPAR2, BETA2, BETA2FRACTION, HSCALE, VMINhydro,
     >        VPAR1_2, VPAR2_2, TMODIFY, SPHERIC, XMDOTold, 
     >        VTURB, RCONold, fHYDROSTART, XMG, XMGold,
     >        dummy, GRADLAST, GRADIL, STEPDAMP, RINT, PL, PLP,
     >        RHOINT, VINMAX, GAMMAL, RSTARold,
     >        RADGAMMASTART, GEDDRAD, GEDDPRINT, XLAMBLUE,
     >        VFINALSCALE, DTDRIN_OLD

      CHARACTER(MAXHIST*8) :: MODHIST
      CHARACTER(80), DIMENSION(3) :: RadiusGridParameters

      CHARACTER(10), DIMENSION(NDIM) :: LEVEL
      CHARACTER(100) :: MODHEAD, MODOLD, MODOLD2
      CHARACTER(10), DIMENSION(MAXATOM) :: ELEMENT
      CHARACTER(8), DIMENSION(NDDIM) :: VELOCRITERION
      CHARACTER(4), DIMENSION(MAXIND) :: KEYCBB
      CHARACTER(2) :: WRTYPE
      CHARACTER(2), DIMENSION(MAXATOM) :: SYMBOL
      CHARACTER(80) :: DENSCON_LINE, ThinCard
      CHARACTER(6) :: BUFFER6
      CHARACTER(8) :: BUFFER8, GEFFKEY
      CHARACTER(9) :: MLRELATION
      CHARACTER(144) :: BUFFER144

      LOGICAL :: TTABLE, OLDTEMP, TEXIST, VPLOT, THIN, NEWVELO, BTWOT,
     >           OLDFGRID, LTESTART
      LOGICAL :: BFEMODEL
      LOGICAL :: BTAUR
      LOGICAL :: bTauFix,             !true if fix option has been set in the TAUMAX CARDS line
     >           bHydroStat,          !true if hydrostatic domain encountered in velocity law
     >           bNoRGrid,            !Parameter for GEOMESH, defines if radius grid should be (re)made or not
     >           bNoDetails,          !decides whether PRIMOD should print all values per depth point or not
     >           bSaveGEFF,            !if true geff is used for hydrostatic scale height instead of g
     >           bOLDSTART,           !true if OLDSTART CARDS option has been set
     >           bThinImprove,        !true if THIN artistic has already been run in this iteration
     >           bOldStratification,  !true if OLD STRATIFICATION option has been set in the CARDS file
     >           bHYDROSOLVE,         !true if hydrodynamic consistent model should be obtained
     >           bHScaleOnly,         !determines if INITVEL calulates only HSCALE or full static law
     >           bFULLHYDROSTAT,      !if true, VELTHIN uses full GAMMA instead of EDDINGTON GAMMA
     >           bGAMMARADMEAN,       !if true, a mean value for GAMMARAD is used instead of individuals
     >           bOLDRAD,             !use old radiation field for HYDROSTATIC INTEGRATION if possible
     >           bGREYSTART,
     >           bNDfirst
      INTEGER :: GEddFix              ! > 0 if GEDD should kept fixed in all calculations

      REAL, EXTERNAL :: WRVEL

C***  Tiefenabhaengiges Clumping nach Goetz Graefener
      REAL, DIMENSION(NDDIM) :: DENSCON, FILLFAC, 
     >                          DENSCON_OLD, FILLFAC_OLD

C***  COMMON /COMTEFF/  TRANSFERS THE EFF. TEMPERATURE FROM SUBR. DECSTAR
      COMMON /COMTEFF/ TEFF,TMIN,TMODIFY,SPHERIC

C***  CONTROL PARAMETER FOR TABULATED INPUT OF T(R) AND V(R)
      INTEGER :: ITAB = 0
       
      !Konstanten
      REAL, PARAMETER :: AMU = 1.66E-24     !atomic mass unit (constant)
      REAL, PARAMETER :: RSUN = 6.96E10     !SOLAR RADIUS ( CM )
      REAL, PARAMETER :: BOLTZ = 1.38E-16   !BOLTZMANN CONSTANT (ERG/DEG)
      REAL, PARAMETER :: RGAS = 8.3145E7    !Gas Constant (CGS)
      REAL, PARAMETER :: GCONST = 6.670E-8  !GRAVITATION CONSTANT (CGS UNITS)
      REAL, PARAMETER :: XMSUN = 1.989E33   !XMSUN = Solar Mass (g)

      !File and channel handles (=KANAL)
      INTEGER, PARAMETER :: hOUT = 6        !write to wruniqX.out (stdout)
      INTEGER, PARAMETER :: hCPR = 0        !write to wruniqX.cpr (stderr)
      INTEGER, PARAMETER :: hMODEL = 3      !(new) MODEL file
      INTEGER, PARAMETER :: hOldMODEL = 9   !old MODEL file (used in OLDSTART)
      INTEGER, PARAMETER :: hOldMODEL2 = 8  !second old "model" file (used if TWOTEMP cards option is set)

C***  Operating system:
      COMMON / COMOS / OPSYS
      CHARACTER(8) :: OPSYS

      CHARACTER(10) :: TIM1, TIM2

C***  Link data to identify program version
      CHARACTER(30) :: LINK_DATE
      CHARACTER(10) :: LINK_USER
      CHARACTER(60) :: LINK_HOST
      COMMON / COM_LINKINFO / LINK_DATE, LINK_USER, LINK_HOST

      bHydroStat = .TRUE.      !default value: hydrostatic domain is achieved


C***  Write Link Data (Program Version) tp CPR file
      WRITE(hCPR,'(2A)') '>>> WRSTART started: Program Version from ', 
     >                LINK_DATE
      WRITE(hCPR,'(4A)') '>>> created by ', LINK_USER(:IDX(LINK_USER)),
     >      ' at host ', LINK_HOST(:IDX(LINK_HOST))

      CALL INSTALL

      IF (OPSYS .EQ. 'CRAY') THEN
        CALL CLOCK(TIM1)
      ELSE
        CALL TIME(TIM1)
      ENDIF

      !Put dummy stuff into COMMON block data to ensure storage reservation
      NDv = NDDIM
      NEWVELO = .FALSE.
      DO L=1, NDDIM
        OLDRADI(L) = 1337.
        OLDVELO(L) = 1338.
      ENDDO
      
C***  JOB NUMBER OF THIS JOB
      JOBNUM=0

C***  INITIALIZE SOME VARIABLES (TODT 04.05.2010)
      TROLD = 0.

C***  READ ATOMIC DATA FROM FILE "DATOM"
      CALL       DATOM (NDIM,N,LEVEL,NCHARG , WEIGHT,ELEVEL,EION,MAINQN,
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
     >            'WRSTART', INDEXMAX, NFEREADMAX, MAXFEIND,
     >             LASTFE, SIGMAFE, INDRB, INDRF,
     >             IFENUP, IFELOW, IFRBSTA, IFRBEND, FEDUMMY,
     >             VDOPFE, DXFE, XLAM0FE, SIGMAINT, BFEMODEL)
 
C***  DECODING INPUT DATA
C***  VDENS1,VDENS2,DENSCON_FIX und CLUMP_CRIT nach goetz...
      CALL DECSTAR (MODHEAD, FM, RSTAR, VDOP, RMAX, TTABLE, NDDIM,
     >              OLDTEMP, MAXATOM, NATOM, ABXYZ, KODAT, VPLOT, 
     >              ATMASS, XDATA, MAXXDAT, OLDFGRID, THIN, ThinCard,
     >              GLOG, GEFFLOG, GEDD, bSaveGEFF, XMSTAR, WRTYPE, 
     >              TAUMAX, TAUACC, bTauFix, BTWOT, TFAC, DENSCON_LINE, 
     >              BLACKEDGE, bOLDSTART, RadiusGridParameters, XMDOT,
     >              XLOGL, RTRANS, BTAUR, DENSCON_FIX, MASSORIGIN,
     >              LRTinput, ELEMENT, bOldStratification, bHYDROSOLVE,
     >              GEddFix, MLRELATION,
     >              NC, VTURB, bOLDRAD, RADGAMMASTART, fHYDROSTART, 
     >              bFULLHYDROSTAT, bGAMMARADMEAN, GEFFKEY, bGREYSTART,
     >              XLAMBLUE, MDOTINPUT, LTESTART)


C***  GENERATION OF THE CONTINUOUS FREQUENCY GRID
      CALL       FGRID (NFDIM, NF, XLAMBDA, FWEIGHT, KEY, NOM, SYMBOL, 
     $                  N, NCHARG, ELEVEL, EION, EINST, NDIM,
     $                  EDGEK,KODAT,MAXATOM,
     $                  INDNUP, INDLOW, LASTIND, KONTNUP, KONTLOW, 
     >                  LASTKON, OLDFGRID, NF2, XLAMBDA2, VDOP, XLAMBLUE)
 
C***  ATMEAN = MEAN ATOMIC WEIGHT ( IN AMU )
      !What does ABXYZ contain?
      ATMEAN=0.
      DO NA=1,NATOM
          ATMEAN = ATMEAN + ABXYZ(NA) * ATMASS(NA)
      ENDDO

C***  STAPEL: NUMBER OF FREE ELECTRONS PER ATOM
C***   S T A R T   A P P R O X I M A T I O N  
      !STAPEL = RNE start approx
      STAPEL=0.
      DO NA=1,NATOM
          STAPEL = STAPEL + ABXYZ(NA) * (STAGE(NA)-1.)
      ENDDO

C***  MEAN MASS PER PARTICLE (ATOMS AND ELECTRONS) IN AMU
C***  -  HAND INPUT, AS IN OLDER PROGRAM VERSIONS, IS IGNORED
      XMASS = ATMEAN / (1. + STAPEL)
 
C***  LOOP FOR AUTOMATICAL ADJUSTMENT OF MAX. ROSSELAND DEPTH (OPTINAL)
      ITTAU = 0

      VMINOLD  = 0.
      VMINOLD2 = 0.

C***  REQUIRED ACCURACY (IN TAU): 
      !TAUACC = 1.E-4  !(This is now a CARDS option with this default value)
C***  MAXIMUM NUMBER OF ITERATIONS
      MAXITTAU = 30
      
      IF (bOLDSTART .OR. bOLDRAD .OR. bOldStratification) THEN
        !use old vmin from old MODEL instead of CARDS value if OLDSTART has been set
        CALL OPENMS (hOldMODEL, IDUMMY, IDUMMY, 1, IERR)
        CALL READMS (hOldMODEL, NDold  ,     1, 'ND      ' , IERR)
        CALL READMS (hOldMODEL, ARAD, NDold-1,   'ARAD    ', IERR)
        CALL READMS (hOldMODEL,RADIUSold,NDold,'R       ',IERR)
        IF (IERR == -10 .AND. bOldStratification) THEN
          WRITE (hCPR, '(A)') 'Error: Cannot find old radius grid!'
          STOP '*** FATAL ERROR IN WRSTART'
        ENDIF
        CALL READMS (hOldMODEL,RSTARold,1,'RSTAR   ',IERR)
        CALL READMS (hOldMODEL, XMGold, 1, 'XMSTAR  ',IERR)
        CALL READMS (hOldMODEL,RCONold,NDold,'RCON    ',IERR)
        CALL READMS (hOldMODEL,TAUROSSOLD,NDold,'TAUROSS ',IERR)
        CALL READMS (hOldMODEL,DENSCON_OLD,NDold,'DENSCON ',IERR)
        DO L=1,NDold
          IF (DENSCON_OLD(L) <= 0. ) THEN
              IF (DENSCON_OLD(1) <= 0.) THEN
                CALL REMARK (
     >            'Zero or negative depth-dep. clumping in OLD model!')
                STOP 'Error in reading depth-dep. clumping '
     >            // 'from OLD model during WRSTART'
              ENDIF
              DENSCON_OLD(L) = DENSCON_OLD(1)
              FILLFAC_OLD(L) = 1. / DENSCON_OLD(L)
          ELSE
              FILLFAC_OLD(L) = 1. / DENSCON_OLD(L)
          ENDIF
        ENDDO
        
        XMGold = XMGold * GCONST * XMSUN
        IF (IERR /= -10) THEN
          CALL READMS (hOldMODEL, VELOold, NDold, 'VELO    ' , IERR)
          IF (IERR /= -10) THEN
            IF (bOldStratification) THEN
              VMIN = VELOold(NDold)
              CALL READMS (hOldMODEL,VCRITold,NDold,'VCRIT   ' ,IERR)
              IF (IERR == -10) THEN
                VCRITold = '        '
              ENDIF
              WRITE (hCPR,'(A)') 'Using old velocity field...'
              DO L=1, NDold
                OLDRADI(L) = RADIUSold(L)
                OLDVELO(L) = VELOold(L)
              ENDDO              
              NDv = NDold
              IF (ABS(VFINAL-VELOold(1)) .GT. 0.01) THEN
                WRITE (hCPR,'(A,/,2(A,F8.2),/,A)') 
     >           '*** WARNING: VFINAL differs from the OLDSTART model:',
     >           '*** VFINAL (old): ', VELOold(1), '  (new): ', VFINAL,
     >           '*** this enforces a scaling!'
                VFINALSCALE = VFINAL / VELOold(1)                              
                DO L=1, NDold
                  OLDVELO(L) = VFINALSCALE * VELOold(L)
                ENDDO
                VFINAL = OLDVELO(1)
                VMIN = OLDVELO(NDold)
              ENDIF

              IF (RADIUSold(1) < RMAX) THEN
                RMAX = RADIUSold(1)
                WRITE (hCPR,'(A,F8.2)') 
     >            '**WARNING: OLD V forces lower RMAX = ', RMAX
              ENDIF
C***          We Use old RCON if inside the old radius grid
C***          However, this can only be done after the GEOMESH call,
C***          otherwise WRVEL would fail as we do not know the old 
C***          BETA law parameters. Therefore ensure interpolation 
C***          in WRVEL by setting RCON > RMAX
              RCON = 1.5 * RADIUSold(1)

              IF (fHYDROSTART > 0.) THEN
                CALL READMS (hOldMODEL, T, NDold,     'T       ', IERR)
                CALL READMS (hOldMODEL, RNE, NDold,   'RNE     ', IERR)
                DO L=1, NDold
                  XMU(L) = ATMEAN / (1. + RNE(L))
                ENDDO
              ENDIF
            ENDIF
          ELSEIF (bOldStratification) THEN
            !No stratification data found => create new one
            WRITE (hCPR,*) '**WARNING: OLD STRATIFICATION NOT FOUND **'
            bOldStratification = .FALSE.
          ENDIF
        ENDIF
        CALL CLOSMS (hOldMODEL, IERR)
      ELSE 
c        WRITE (hCPR,*) 'Fresh start: No old model is used!'
        NDold = 1  !must be set for DIMENSION-statements in PREP_GAMMARAD
      ENDIF
      
      
C***  Preparation of GAMMARAD for the hydrostatic equation
      CALL PREP_GAMMARAD (bOLDRAD, bFULLHYDROSTAT, GEDD,
     >      GEddFix, GAMMARAD, NDold, ARAD, GLOG, RSTAR,
     >      RADIUSold, RSTARold, XMGold, TAUROSSOLD, RCONold, GEFFLOG,
     >      STAPEL, ATMEAN, XLOGL, XMSTAR, RADGAMMASTART, 
     >      GEDDRAD, bOldStratification, bGAMMARADMEAN, bSaveGEFF )      
      
      bNDfirst = .TRUE.
      XMG = GCONST * XMSTAR * XMSUN
      VMINhydro = -99.      !init negative => not yet calculated
      
C    ------- Main TAU iteration loop starts here ----------------------
      tauit: DO
          ITTAU = ITTAU + 1
C          WRITE (hCPR,*) " ITTAU=", ITTAU

C***      CHECK WETHER NEW VMIN IS INBETWEEN 0. AND VFINAL
C          IF (VMIN < 0.) VMIN = 1.E-4
          IF (VMIN > VFINAL) VMIN = 0.9 * VFINAL
C         IF (VMIN == VMINOLD) VMIN = 0.9 * VMIN        
          IF (VMINhydro > 0. .AND. VMIN > VMINhydro) THEN
            VMIN = VMINhydro
            WRITE (hCPR,*) '*** WARNING: No hydrodynamic solution ***'
          ENDIF

C***      INITIALISATION OF THE VELOCITY-FIELD PARAMETERS
C***      Turbulence pressure added, 13-Mar-2014
          IF (.NOT. bOldStratification) THEN
            bHScaleOnly = (THIN .AND. (ITTAU > 1))
            CALL INITVEL (RMAX,TEFF,GEFFLOG,RSTAR,XMASS,
     >                    VTURB, bHScaleOnly, bHydroStat)
            IF (ITTAU == 1 .OR. (.NOT. THIN)) NEWVELO=.TRUE.
          ELSE
C***        Old stratification is used          
            ND = NDold
            DO L=1, ND
              RADIUS(L) = OLDRADI(L)
              VELO(L) = OLDVELO(L)                
            ENDDO
            NEWVELO=.FALSE.
          ENDIF

C***      LOOP FOR IMPROVED HYDROSTATIC EQUATION ("THIN WIND" OPTION)
C***       (MUST BE DONE TWICE, FIRST TO ESTABLISH A RADIUS MESH AND THE 
C***        TEMPERATURE STRUCTURE, SECOND FOR THE EXACT HYDROSTATIC EQ.)

          bThinImprove = .FALSE.

C         ------- THIN WIND loop ----------------------
          thinwind: DO

C***  GENERATION OF THE RADIUS GRID, P-GRID, Z-GRID, AND ND
            bNoRGrid = .FALSE.
            CALL GEOMESH (RADIUS,INCRIT,P,Z,ND,NDDIM,NP,NPDIM,RMAX, 
     >                    RadiusGridParameters, bNoRGrid, NC)
            IF (bNDfirst) THEN
              WRITE (hCPR,*) ' ND, NDold: ', ND, NDold
              IF (bGREYSTART) THEN
                WRITE (hCPR,*) 'TAUROSS scale used: GREY'
              ELSE
                WRITE (hCPR,*) 'TAUROSS scale used: CONT'
              ENDIF
              bNDfirst = .FALSE.
            ENDIF

C***  START APPROXIMATION FOR EL. DENSITY PUT INTO ARRAY (NEEDS ND)
            DO L=1, ND
                RNE(L) = STAPEL
            ENDDO
   
C***  READ TEMPERATURE STRUCTURE FROM OLD MODEL, IF REQUESTED
C     (note: REQUIRES GEOMESH first for interpolation)      
            IF (OLDTEMP) THEN
              CALL READOLDT (hOldMODEL, ND, NDDIM,T,RADIUS,
     >                      TOLD,ROLD,MODOLD,JOBNOLD,TEFF,TEFFOLD,
     >                      TAURCONT, TAURCONTOLD, BTAUR, DTDRIN_OLD)
              IF (BTWOT) THEN
                CALL READOLDT (hOldMODEL2, ND,NDDIM,T2,RADIUS,
     >                        TOLD2,ROLD2,MODOLD2,JOBNOLD2,TEFF,
     >             TEFFOLD2, TAURCONT, TAURCONTOLD, BTAUR, DTDRIN_OLD)
                IF (TMIN > 6000.) THEN
                  TMIN2 = TMIN
                ELSE
                  TMIN2 = 6000.
                ENDIF
                DO L=1, ND
                  TS = T(L)
                  T(L) = T(L) - TFAC*(T(L) - T2(L))
                  IF (T(L) > 1.2 * TS) THEN
                    T(L) = 1.2 * TS
                  ELSE IF (T(L) < 0.8 * TS) THEN
                    T(L) = 0.8 * TS
                  ENDIF
                ENDDO
              ENDIF
            ENDIF

C***  READ TEMPERATURE OR VELOCITY (OR BOTH) FROM FILE 'TABLE'
C***  GRADI IS OVERWRITTEN AGAIN BY GRADIFF
            IF (TTABLE) CALL TABREAD (ND,RADIUS,VELO,GRADI,T,ITAB)
 
C***  ENTOT = TOTAL NUMBER DENSITY OF ALL ATOMS
            DO L=1,ND
              RL = RADIUS(L)
              IF (.NOT.TTABLE .OR. ITAB < 2) THEN
                IF (NEWVELO) THEN
                  VELO(L) = WRVEL(RL)
                ELSE
                  CALL SPLINPOX(VELO(L), RL, OLDVELO, OLDRADI, NDv)
                ENDIF
              ENDIF
              RHO(L) = FM / RL / RL / VELO(L) / 1.E5
              ENTOT(L) = RHO(L) / AMU / ATMEAN
            ENDDO            

            CALL       GRADIFF  (ND,VELO,GRADI,RADIUS)
            DO L=1,ND
              XMU(L) = ATMEAN / (1. + RNE(L))
              IF (bFULLHYDROSTAT) THEN
                IF (bOLDRAD) THEN
                  IF (RADIUS(L) > RADIUSold(1)) THEN
                    GAMMAL = GAMMARAD(1)
                  ELSEIF (RADIUS(L) < RADIUSold(NDold)) THEN
                    GAMMAL = GAMMARAD(NDold)
                  ELSE
                    CALL SPLINPOX(GAMMAL, RADIUS(L), 
     >                           GAMMARAD, RADIUSold, NDold)
                  ENDIF
                ELSEIF (RADGAMMASTART >= 0.) THEN
                  GAMMAL = RADGAMMASTART
                ELSE 
                  GAMMAL = GEDD
                ENDIF
                GEFFL(L) = (10.**GLOG) * (1. - GAMMAL)
              ELSE
                !either Thompson only or fixed => no depth-dependent value
                GEFFL(L) = (10.**GLOG) * (1. - GEDD)     
              ENDIF
            ENDDO
C***  Definition of the Depth-dependent Clumping Factor
            IF (bOldStratification) THEN
              DO L=1, ND
                IF (RADIUS(L) > RADIUSold(1)) THEN
                  DENSCON(L) = DENSCON_OLD(1)
                  FILLFAC(L) = FILLFAC_OLD(1)
                ELSEIF (RADIUS(L) < RADIUSold(ND)) THEN
                  DENSCON(L) = DENSCON_OLD(ND)
                  FILLFAC(L) = FILLFAC_OLD(ND)
                ELSE                
                  CALL SPLINPOX(DENSCON(L), RADIUS(L), 
     >                          DENSCON_OLD, RADIUSold, NDold)
                  CALL SPLINPOX(FILLFAC(L), RADIUS(L), 
     >                          FILLFAC_OLD, RADIUSold, NDold)
                ENDIF
              ENDDO
            ELSE
              CALL CLUMP_STRUCT (DENSCON, FILLFAC, ND, DENSCON_FIX, 
     >                           VELO, TAURCONT, DENSCON_LINE, 
     >                           RADIUS, T, XMU)
            ENDIF

C***  ITAB = 2: ONLY TABULATED INPUT OF V(R), I.E. T(R) MUST BE CALCULATED
            IF (ITAB .EQ. 2) TTABLE=.FALSE.
            TEXIST=TTABLE .OR. OLDTEMP

C***  TEMPERATURE STRATIFICATION AND INITIAL POPNUMBERS (LTE)            
            CALL GREY(ND,T,RADIUS,XLAMBDA,FWEIGHT,NF,ENTOT,RNE,RSTAR,
     >                ALPHA,SEXPO,
     >                ADDCON1, ADDCON2, ADDCON3, 
     >                IGAUNT,POPNUM,TAUGREY,R23,TEXIST,NDIM,N,
     >                LEVEL,NCHARG,WEIGHT,ELEVEL,EION,EINST,ENLTE,
     >                KODAT,ABXYZ,NOM,NFIRST,NLAST,NATOM,
     >                EXPFAC,SIGMAKI,NFEDGE,OPAC,ETAC,SIGMAFF,
     >                MAXION,MAXATOM,SIGMATHK,EDGEK,
     >                SEXPOK,KONTNUP,KONTLOW,LASTKON,XDATA, DENSCON,
     >                FILLFAC)
            IF (.NOT. bGREYSTART) THEN
C***          improve TAUROSS scale by using OPAROSS instead of using OPAGREY
C               This has an influence on the velocity and the radiation field 
C               (JSTART uses TAUROSS scale, note that R23 is always from GREY)
              CALL TAUSCAL(RSTAR, ND, RADIUS, RNE, ENTOT,
     >                     T, POPNUM, NDIM, N, EN, LEVEL, NCHARG, 
     >                     WEIGHT, ELEVEL, EION, EINST, ALPHA, SEXPO,
     >                     ADDCON1, ADDCON2, ADDCON3, 
     >                     IGAUNT, NOM, NF,
     >                     XLAMBDA, FWEIGHT,
     >                     TAUTHOM, TAURCONT,
     >                     MAXATOM, SIGMATHK, SEXPOK, EDGEK, KODAT,
     >                     KONTNUP, KONTLOW, LASTKON, 
     >                     DENSCON, FILLFAC
     >        ) 
            ELSE
              DO L=1, ND
                TAURCONT(L) = TAUGREY(L)
              ENDDO
            ENDIF

            !Added XMU in WRSTART for use in VELTHIN and to avoid problems 
            !  with first steal run if hydro or TAUFIX is enabled
            !  XMU is recalculated here because RNE has been changed by GREY
            DO L=1,ND
                XMU(L) = ATMEAN / (1. + RNE(L))
            ENDDO

            IF (bHYDROSOLVE) THEN
              AMIN = SQRT(RGAS * T(ND) / XMU(ND))  !sound speed at inner boundary
              IF (VMIN*1.E5 >= AMIN) THEN
                !inner boundary values for v must be less than sound speed
                ! otherwise G~ would never cross zero
                VMINhydro = AMIN / SQRT(2.) / 1.E5
                CYCLE tauit
              ENDIF
            ENDIF
            
C***  NEW VELOCITY FIELD FROM INTEGRATION OF HYDROSTATIC EQUATION
            IF (((.NOT.TTABLE) .OR. (ITAB .LT. 2)) .AND. THIN) THEN
              CALL VELTHIN(T, RADIUS, VELO, ND, RSTAR, RMAX,
     >                    GEFFL, XMU, VTURB, ThinCard, VELOCRITERION)

C***          STORE RADIUS AND VELO FOR WRVEL, BEFORE A NEW RADIUS 
C***            GRID WILL BE DEFINED
              DO L=1, ND
                OLDRADI(L) = RADIUS(L)
                OLDVELO(L) = VELO(L)
              ENDDO
              NDv = ND
              NEWVELO = .FALSE.

              DO L=1, ND
                RL = RADIUS(L)
                RHO(L) = FM / RL / RL / VELO(L) / 1.E5
                ENTOT(L) = RHO(L) / AMU / ATMEAN
              ENDDO

              IF (.NOT. bThinImprove) THEN
                bThinImprove = .TRUE.
                CYCLE thinwind
              ELSE
                EXIT thinwind
              ENDIF
              
              EXIT
            ELSE
              EXIT
            ENDIF
            
          ENDDO thinwind
C         ------- THIN WIND loop ends here ----------------------


C***  AUTOMATICAL ADJUSTMENT OF VMIN (OPTIONAL: IF TAUMAX SPECIFIED)
          IF (bOldStratification) THEN
C***        Use old RCON in case of OLD STRATIFICATION 
C***        if inside of both, old and new grid
            IF (RCONold < RADIUSold(1) .AND. RCONold < RADIUS(1)) THEN
              RCON = RCONold
            ENDIF          
            WRITE (0,FMT='(A,F10.4,A,F11.6)') 
     >          ' *** OLD STRATIFICATION USED: TAURCONT=', TAURCONT(ND), 
     >          '   VMIN=', VMIN
            EXIT tauit
          ELSEIF (TAUMAX > .0) THEN

            IF (VMIN == VMINhydro) THEN
              TAUMAX = REAL(CEILING(TAURCONT(ND)))   !taumax = next larger integer 
              WRITE (hCPR,*) '*** TAUMAX has been adjusted for hydro'
              WRITE (hCPR,FMT='(A,F8.2)') '*** New TAUMAX = ', TAUMAX
            ENDIF
          
            !Enforce monotonic TAUROSS scale
            DO L=2, ND
              IF (TAURCONT(L) < TAURCONT(L-1)) THEN
                TAURCONT(L) = TAURCONT(L-1) + 1.E-10
              ENDIF
            ENDDO

            WRITE (hCPR,FMT=77) ITTAU, TAURCONT(ND), VMIN
   77       FORMAT (' *** TAUMAX ITERATION: ITTAU=', I3,
     >          '   TAURCONT=', F10.4, '  VMIN=', G12.4)

C***        make an extrapolation at every ten steps
            VMINOLD2 = VMINOLD
            VMINOLD = VMIN
            TROLD2 = TROLD
            TROLD = TAURCONT(ND)
            IF (ABS(TAURCONT(ND)-TAUMAX) > TAUACC  .AND.  
     >        ITTAU/10*10 == ITTAU .AND. ITTAU > 2) THEN
C***            VMIN = 0.5 * (VMINOLD + VMINOLD2)
                VMIN = (VMINOLD2*(TROLD-TAUMAX) - VMINOLD * 
     >          (TROLD2-TAUMAX)) /  (TROLD - TROLD2)
                IF (VMIN > .0) THEN
                   IF (ITTAU < MAXITTAU) CYCLE
                ELSE
                   VMIN = VMINOLD
                ENDIF
            ENDIF

C***        TAUROSS TOO LARGE:
            IF (TAURCONT(ND) > TAUMAX+TAUACC) THEN
              CALL SPLINPOX(VMIN, TAUMAX, VELO, TAURCONT, ND)
              IF (ITTAU < MAXITTAU) CYCLE
            ENDIF

C***        TAUROSS TOO SMALL:
            IF (TAURCONT(ND) < TAUMAX-TAUACC) THEN
              VMIN = VMIN * (TAURCONT(ND) / TAUMAX)**.5
              IF (ITTAU < MAXITTAU) CYCLE
            ENDIF

          ENDIF

          EXIT !Exit the TAU iteration loop if no cycling criteria was met
        
      ENDDO tauit
C     ------- Main TAU iteration loop ends here ----------------------


C***  VELOCITY GRADIENT BY DIFFERENTIATION (final calculation)
      CALL       GRADIFF  (ND,VELO,GRADI,RADIUS)

      !Store VELOCRITERION as two-character vector in model file
      VCRIT = '        '
      IF (bOldStratification) THEN
        !if old v is used, take also old VCRIT as starting approach
        ! (note: this might not be exact if taumax changes signifcantly!)
        DO L=1, ND
          VCRIT(L) = VCRITold(L)
        ENDDO
      ELSE      
        DO L=1, ND
          !Note: static identifier "ST" is already set above if used
          SELECTCASE (VELOCRITERION(L))
            CASE ('HYSTINT ')
              VCRIT(L) = 'HS      '
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
      ENDIF

C***  MODEL-FILE: MASS STORAGE FILE IN NAME-INDEX MODE ****************
      CALL OPENMS (3,IDUMMY,IDUMMY,1, IERR)
      IFLAG = -1  !'default' flagging by default ;)
      CALL WRITMS (3,ND,1,         'ND      ', IFLAG, IDUMMY, IERR)
      CALL WRITMS (3,RADIUS,ND,    'R       ', IFLAG, IDUMMY, IERR)
      CALL WRITMS (3,NP,1,         'NP      ', IFLAG, IDUMMY, IERR)
      CALL WRITMS (3,P,NP,         'P       ', IFLAG, IDUMMY, IERR)
      CALL WRITMS (3,Z,ND*NP,      'Z       ', IFLAG, IDUMMY, IERR)
      CALL WRITMS (3,ENTOT,ND,     'ENTOT   ', IFLAG, IDUMMY, IERR)
      CALL WRITMS (3,DENSCON,ND,   'DENSCON ', IFLAG, IDUMMY, IERR)
      CALL WRITMS (3,TEFF,1,       'TEFF    ', IFLAG, IDUMMY, IERR)
      CALL WRITMS (3,GLOG,1,       'GLOG    ', IFLAG, IDUMMY, IERR)
      IF (bSaveGEFF) THEN
        !Save geff only if this value should be fixed in the following
        !calculations (meaning geff or EddGamma was given in CARDS)
        CALL WRITMS (3,GEFFLOG,1,    'GEFFLOG ', IFLAG, IDUMMY, IERR)
      ELSE
        GFLSAV = -1. * GEFFLOG  !negative value indicating flexible value
        CALL WRITMS (3, GFLSAV,1,    'GEFFLOG ', IFLAG, IDUMMY, IERR)
      ENDIF
      CALL WRITMS (3,GEFFKEY,1,    'GEFFKEY ', IFLAG, IDUMMY, IERR)
      IF (GEddFix > 0) THEN
        CALL WRITMS (3,GEDD,1,       'GEDD    ', IFLAG, IDUMMY, IERR)
      ELSE
        CALL WRITMS (3,-1.0   ,1,    'GEDD    ', IFLAG, IDUMMY, IERR)
      ENDIF
      IF (GEDDRAD > 0. .AND. .NOT. bOldStratification) THEN
        CALL WRITMS (3,GEDDRAD,1,    'GEDDRAD ', IFLAG, IDUMMY, IERR)
      ENDIF
      CALL WRITMS (3,T,ND,         'T       ', IFLAG, IDUMMY, IERR)
      CALL WRITMS (3,VELO,ND,      'VELO    ', IFLAG, IDUMMY, IERR)
      CALL WRITMS (3,VMIN,1,       'VMIN    ', IFLAG, IDUMMY, IERR)
      CALL WRITMS (3,GRADI,ND,     'GRADI   ', IFLAG, IDUMMY, IERR)
      CALL WRITMS (3,TAURCONT,ND,  'TAURCONT', IFLAG, IDUMMY, IERR)
      CALL WRITMS (3,RSTAR,1,      'RSTAR   ', IFLAG, IDUMMY, IERR)
      CALL WRITMS (3,RCON,1,       'RCON    ', IFLAG, IDUMMY, IERR)
      CALL WRITMS (3,INCRIT,ND,    'INCRIT  ', IFLAG, IDUMMY, IERR)
      CALL WRITMS (3,VCRIT, ND,    'VCRIT   ', IFLAG, IDUMMY, IERR)
      CALL WRITMS (3,XMSTAR,1,     'XMSTAR  ', IFLAG, IDUMMY, IERR)
      CALL WRITMS (3,XMDOT,1,      'XMDOT   ', IFLAG, IDUMMY, IERR)
      CALL WRITMS (3,RHO,ND,       'RHO     ', IFLAG, IDUMMY, IERR)
      CALL WRITMS (3,VDOP,1,       'VDOP    ', IFLAG, IDUMMY, IERR)
      CALL WRITMS (3,NF,1,         'NF      ', IFLAG, IDUMMY, IERR)
      CALL WRITMS (3,XLAMBDA,NF,   'XLAMBDA ', IFLAG, IDUMMY, IERR)
      CALL WRITMS (3,NF2,1,        'NF2     ', IFLAG, IDUMMY, IERR)
      CALL WRITMS (3,XLAMBDA2,NF2, 'XLAMBD2 ', IFLAG, IDUMMY, IERR)
      CALL WRITMS (3,FWEIGHT,NF,   'FWEIGHT ', IFLAG, IDUMMY, IERR)
      CALL WRITMS (3,KEY,NF,       'KEY     ', IFLAG, IDUMMY, IERR)
      CALL WRITMS (3,POPNUM,ND*N,  'POPNUM  ', IFLAG, IDUMMY, IERR)
      CALL WRITMS (3,POPNUM,ND*N,  'POPLTE  ', IFLAG, IDUMMY, IERR)
      CALL WRITMS (3,RNE,ND,       'RNE     ', IFLAG, IDUMMY, IERR)
      CALL WRITMS (3,XMU,ND,       'XMU     ', IFLAG, IDUMMY, IERR)
      CALL WRITMS (3,ABXYZ,NATOM,  'ABXYZ   ', IFLAG, IDUMMY, IERR)
      CALL WRITMS (3,MODHEAD,13,   'MODHEAD ', IFLAG, IDUMMY, IERR)
      NEXTK=1
      CALL WRITMS (3,NEXTK,1,      'NEXTK   ', IFLAG, IDUMMY, IERR)
      CALL WRITMS (3,JOBNUM,1,     'JOBNUM  ', IFLAG, IDUMMY, IERR)
      BUFFER6 = 'UNDEF.'
      READ(UNIT=BUFFER6, FMT='(A6)') TOTOUT
      CALL WRITMS (3,TOTOUT,1,     'TOTOUT  ', IFLAG, IDUMMY, IERR)
      CALL WRITMS (3,XDATA,MAXXDAT,'XDATA   ', IFLAG, IDUMMY, IERR)
      CALL WRITMS (3,VTURB,1,      'VTURB   ', IFLAG, IDUMMY, IERR)
      IF (THIN) THEN
        DO L=1, ND
          GRSTATIC(L) = 1. - GEFFL(L) / (10.**GLOG)
        ENDDO
        CALL WRITMS (3,GRSTATIC, ND, 'GRSTATIC', IFLAG, IDUMMY, IERR)
      ENDIF
      IF (OLDTEMP .AND. bOldStratification) THEN
         CALL WRITMS (3,DTDRIN_OLD, 1, 'DTDRIN  ', IFLAG, IDUMMY, IERR)
      ENDIF
      

C***  WRITE THE BEGINNING OF THE MODEL HISTORY      
      MODHIST(9:32) = '/      0. WRSTART       '
      LAST=4    !3 * 8 +  offset (8)
      WRITE(UNIT=BUFFER8, FMT='(A8)') LAST 
      MODHIST(1:8)=BUFFER8 !First entry contains the number of used integers (=chars / 8)


      IF (OLDTEMP) THEN
         WRITE(UNIT=BUFFER144, FMT=2) MODOLD, JOBNOLD
    2    FORMAT ('T(R) FROM OLD MODEL: ',A,'  AFTER JOB NO.',I4,4X)
         CALL ADDHISTENTRY(MODHIST, LAST, MAXHIST, 144, BUFFER144)
         
         IF (BTWOT) THEN
           WRITE (UNIT=BUFFER144, FMT=22) TFAC, MODOLD2, JOBNOLD2
   22      FORMAT ('SECOND (TFAC=',F4.1',) : ',A,
     >             '  AFTER JOB NO.',I4,4X)
           CALL ADDHISTENTRY(MODHIST, LAST, MAXHIST, 144, BUFFER144)
         ENDIF
      ENDIF

      CALL WRITMS (3,MODHIST,MAXHIST,'MODHIST ', IFLAG, IDUMMY, IERR)
 
C***  START APPROXIMATION FOR THE RADIATION FIELD
      CALL JSTART (NF,XLAMBDA,KEY,ND,RADIUS,T,XJC,XJL,ELEVEL,
     >             N,EINST,NDIM,INDNUP,INDLOW,LASTIND,R23,TAUGREY, 
     >             LTESTART,BLACKEDGE)
 
      CALL CLOSMS (3, IERR)
C**********************************************************************

C***  PRINTOUT OF THE FUNDAMENTAL STELLAR PARAMETERS
      IF (bFULLHYDROSTAT) THEN
        GEDDPRINT = GEDDRAD
      ELSE
        GEDDPRINT = GEDD
      ENDIF
      CALL PRIPARAM (MODHEAD, TEFF, RSTAR, XMDOT, XLOGL, RTRANS, 
     >        VFINAL, VDOP, DENSCON, FILLFAC, GLOG, GEFFLOG, GEDDPRINT, 
     >        GEddFix, RMAX, XMSTAR, WRTYPE, MASSORIGIN, LRTinput, ND, 
     >        MLRELATION, VTURB, MDOTINPUT)



C***  PRINTOUT OF THE CHEMICAL COMPOSITION
      CALL PRICOMP (NDIM, EINST, N, NCHARG, NOM, NATOM, ABXYZ, ATMASS,
     $              STAGE, NFIRST, NLAST, ELEMENT, SYMBOL, LASTIND,
     $              INDLOW, INDNUP, NAUTO, LOWAUTO,
     $              EAUTO, KONTNUP, KONTLOW, LASTKON, XMASS, KRUDAUT)

C***  PRINTOUT OF THE X-RAY DATA
      CALL PRIXDAT(XDATA,MAXXDAT)

C***  PRINTOUT OF VARIOUS MODEL SPECIFICATIONS
      bNoDetails = .FALSE.  !print all details at the start
      CALL PRIMOD (ND,RADIUS,INCRIT,ENTOT,T,VELO,GRADI,NP,OLDTEMP,
     $             MODHEAD,JOBNUM,MODOLD,JOBNOLD,TTABLE,TAURCONT,R23,
     $             TEFFOLD,THIN, ITTAU, MAXITTAU, RCON, BTWOT, MODOLD2, 
     >             JOBNOLD2, TFAC, BETA, VPAR1, VPAR2, DENSCON, BETA2, 
     >             BETA2FRACTION, HSCALE, bNoDetails,.TRUE., VTURB)
 
C***  PLOT OF THE VELOCITY LAW
      IF (VPLOT) CALL PLOTV (ND,RADIUS,VELO,MODHEAD,JOBNUM)
      

      CALL JSYMSET ('G0','0')

      CALL STAMP (OPSYS, 'WRSTART', TIM1)

      STOP 'O.K.'
      END
