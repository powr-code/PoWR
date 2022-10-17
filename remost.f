      SUBROUTINE REMOST (ND,NDDIM,NP,NPDIM,P,Z,
     >                 ENTOT,T,NF,NFDIM,XLAMBDA,FWEIGHT,
     >                 RADIUS,RSTAR,RCON,VELO,GRADI,VDOP,INCRIT,VCRIT,
     >                 TOTOUT,TOTIN,EMFLUX, MAXAUTO,LASTIND,INDNUP,
     >                 INDLOW,KEY,RNE,MODHEAD,JOBNUM,XJC,XJL,MODHIST,
     >                 MAXHIST,LAST,EINST,NDIM,ABXYZ,NATOM,TEFF,
     >                 GLOG, GEFFLOG, bFixGEFF, GEDD, bGEddFix, 
     >                 ENTOTDENS, XMDOT, Rcritical, GEDDRAD,
     >                 XMSTAR, NAUTO, DRXJL, KRUDAUT, MAXXDAT,XDATA, 
     >                 HTOTL, BHTOTERR, HTOTND, DTDRIN,
     >                 DENSCON, FILLFAC, 
     >                 LASTFE, ARAD, bNoARAD, WFELOW, WFENUP, NCOLIP,
     >                 FTCOLI, EMCOLI, GAHIST, 
     >                 MAXGAHIST, MAXIND, MAXINDE, MAXKONT, LASTKON,
     >                 XJTOTL, XKTOTL, XNTOTL, WJC, GF, OPARND,
     >                 OPASMEAN, OPASMEANTC, QFJMEAN, OPAJMEAN, SMEAN,
     >                 OPAJMEANTC, OPAPMEAN, QOPAHMEAN, EDDIHOUTJMEAN, 
     >                 XLAMBDA2, NF2, HTOTOUTMINUS, 
     >                 FF_INFO, IFF_DK, IFF_MAX_MS, IFF_N_MS, 
     >                 OPALAMBDAMEAN, TAUROSS, TAUROSScont, VMIN, N, 
     >                 POPNUM, SIGMAKI, ATMASS, ACONT, ATHOM, OPAROSS,
     >                 RHO, XMU, ALPHAF, NEXTHYDRO, LASTHYDRO, VTURB,
     >              GRSTATIC, LASTBACKUP, LASTTAU, GEFFKEY, ZERO_RATES)  
C*******************************************************************************
C***  READING THE MODEL FILE FOR MAIN PROGRAM "STEAL"
C*******************************************************************************
 
ccc      IMPLICIT NONE

      INTEGER, PARAMETER :: TINYINT = SELECTED_INT_KIND(2)

      INTEGER :: ND, NDDIM, NP, NPDIM, NF, NFDIM, N, NDDENSCON,
     >           JOBNUM, NF2, NCOLIP, MAXHIST, LAST, MAXINDE,
     >           IND, NUP, INDNUP, LOW, INDLOW, INDFE, LASTFE,
     >           LASTIND, MAXXDAT, NATOM, NDIM, MAXAUTO, NAUTO,
     >           INDDR, IFF_N_MS, IFF_MAX_MS, MAXGAHIST, MAXIND,
     >           MAXKONT, LASTKON, KON, LASTBACKUP, LASTTAU
 
      REAL, DIMENSION(NDDIM), INTENT(OUT) :: RNE
      CHARACTER(8), DIMENSION(NDDIM) :: INCRIT, VCRIT
      CHARACTER(8), DIMENSION(NFDIM) :: KEY
      REAL, DIMENSION(NDIM,NDIM) :: EINST
      REAL, DIMENSION(NFDIM,MAXKONT) :: SIGMAKI
      REAL, DIMENSION(NATOM) :: ABXYZ, ATMASS
      INTEGER, DIMENSION(MAXAUTO) :: KRUDAUT
      REAL, DIMENSION(MAXXDAT) :: XDATA
      DIMENSION INDNUP(LASTIND),INDLOW(LASTIND)
      REAL, DIMENSION(ND) :: HTOTL, ENTOT, ENTOTDENS, ARAD,
     >                       XJTOTL, XKTOTL, XNTOTL, SMEAN,
     >                       RADIUS, VELO, GRADI, T, FTCOLI,
     >                       OPASMEAN, OPASMEANTC, QFJMEAN,
     >                       OPAJMEAN, OPAJMEANTC, OPAPMEAN, 
     >                       QOPAHMEAN, OPALAMBDAMEAN, OPAROSS, 
     >                       GRSTATIC
      REAL, DIMENSION(NF) :: EMCOLI, EMFLUX, FWEIGHT
      REAL, DIMENSION(NDDIM * MAXIND) :: WFELOW, WFENUP
      REAL, DIMENSION(NDDIM * MAXINDE) :: XJL
      REAL, DIMENSION(NDDIM * NFDIM) :: XJC, WJC
      REAL, DIMENSION(26,MAXGAHIST) :: GAHIST
      LOGICAL, DIMENSION(MAXAUTO) :: DRXJL
      LOGICAL :: STHLP, BHTOTERR, bGEddFix, bFixGEFF, bNoARAD
      REAL, DIMENSION(10) :: FF_INFO
      INTEGER (KIND=TINYINT), DIMENSION(IFF_MAX_MS) :: IFF_DK
      CHARACTER(8*MAXHIST) :: MODHIST

      INTEGER :: I, J, K, L, IERR, IDUMMY, NA
      REAL, INTENT(OUT) :: GF, RSTAR, VDOP, XLAMBDA, XLAMBDA2, GEDD,
     >                     TEFF, XMDOT, Rcritical
      REAL, INTENT(INOUT) :: GLOG, GEFFLOG, VMIN, RCON, XMSTAR
      REAL :: ATMEAN, OPARND, EDDIHOUTJMEAN, HTOTOUTMINUS, 
     >        TOTIN, TOTOUT, FM, RSTARSU, GEDDRAD,
     >        VFINAL, VMINorg, BETA, VPAR1, VPAR2, RCONorg, HSCALE, 
     >        BETA2, BETA2FRACTION, VPAR1_2, VPAR2_2, VTURB,
     >        HTOTND, DTDRIN

      REAL, DIMENSION(NDDIM), INTENT(INOUT) :: TAUROSS, TAUROSScont
      REAL, DIMENSION(NPDIM) :: P
      REAL, DIMENSION(NDDIM,NPDIM) :: Z
      REAL, DIMENSION(NDDIM, NDIM) :: POPNUM
      LOGICAL ZERO_RATES (NDIM*NDDIM)

      !Transfer RCON from MODEL file to routines needing the velocity field parameters (e.g. PLOTV)
      COMMON /VELPAR/ VFINAL,VMINorg,BETA,VPAR1,VPAR2,RCONorg,HSCALE,
     >                BETA2, BETA2FRACTION, VPAR1_2, VPAR2_2

C***  hydro stuff
      REAL, DIMENSION(NDDIM-1) :: ATHOM, ACONT
      REAL, DIMENSION(NDDIM) :: RHO, ALPHAF, XMU
      INTEGER, INTENT(OUT) :: NEXTHYDRO, LASTHYDRO

      CHARACTER(8) :: NAME, GEFFKEY
      CHARACTER(100) :: MODHEAD

c***  tiefenabh. clumping nach goetz
      REAL, DIMENSION(NDDIM) :: DENSCON, FILLFAC

C***  Physical constants
      REAL, PARAMETER :: AMU = 1.66E-24        !Atomic mass unit (gramm)     
      REAL, PARAMETER :: GCONST = 6.670E-8     !GRAVITATION CONSTANT (CGS UNITS)
      REAL, PARAMETER :: TEFFSUN = 5780.       !EFFECTIVE TEMPERATURE OF THE SUN
      REAL, PARAMETER :: XMSUN = 1.989E33      !XMSUN = Solar Mass (g)
      REAL, PARAMETER :: RSUN = 6.96E10        !SOLAR RADIUS ( CM )

      !File and channel handles (=KANAL)
      INTEGER, PARAMETER :: hOUT = 6        !write to wruniqX.out (stdout)
      INTEGER, PARAMETER :: hCPR = 0        !write to wruniqX.cpr (stderr)

C***  Preset VMIN (if none read)
      VMIN = -1.0

C***  Preset NCOLIP
      NCOLIP = 0

C***  Preset GAMMA HISTORY
      DO I=1, MAXGAHIST
        GAHIST(1,I) = -1.
      ENDDO

C***  Preset GF
      GF = 0.

      CALL OPENMS(3, IDUMMY, IDUMMY, 1, IERR)
      CALL READMS(3,ND,1,       'ND      ', IERR)
      IF(ND > NDDIM) THEN
            CALL REMARK ('TOO MANY DEPTH POINTS')
            STOP 'ERROR'
            ENDIF
      CALL READMS(3,NP,1,       'NP      ', IERR)   !fuer GEOMESH-Aufrufe noetig
      CALL READMS(3,NF ,1,'NF      ', IERR)
      CALL READMS(3,NF2,1,'NF2     ', IERR)
      
C***  CHECK WHETHER DIMENSIONING IS SUFFICIENT
      IF (ND > N DDIM) THEN
            CALL REMARK ('TOO MANY DEPTH POINTS')
            STOP 'ERROR'
      ENDIF
      IF (NP > NPDIM) THEN
            CALL REMARK ('TOO MANY IMPACT-PARAMETER POINTS')
            STOP 'ERROR'
      ENDIF      
      IF (NF > NFDIM) THEN
            CALL REMARK ('TOO MANY FREQUENCY POINTS')
            STOP 'ERROR'
      ENDIF
      
      CALL READMS(3,P, NP,      'P       ', IERR)
      CALL READMS(3,Z, ND*NP,   'Z       ', IERR)
      
      CALL READMS(3,ENTOT,ND,   'ENTOT   ', IERR)
      DENSCON(1:ND) = 0.0
      CALL LENGTHMS(3, NDDENSCON, 'DENSCON ', IERR)
      IF (NDDENSCON == 1) THEN
        !old MODEL file => only one DENSCON line
        CALL READMS(3,DENSCON(1),1, 'DENSCON ', IERR)
        CALL CHANGE('DENSCON ', 'DCOLD   ', 3)
        DO L=2, ND
          DENSCON(L) = DENSCON(1)
        ENDDO
        CALL WRITMS(3,DENSCON,ND,'DENSCON ',-1, IDUMMY, IERR)
      ENDIF
      CALL READMS(3,DENSCON,ND, 'DENSCON ', IERR)
      IF (IERR .EQ. -10) DENSCON(1:ND) = 1.
      DO L=1,ND
         IF (DENSCON(L) .LE. 0. ) THEN
            IF (DENSCON(1) .LE. 0. ) THEN
               CALL REMARK ('Zero or negative depth-dep. clumping!')
               STOP 'Error in depth-dep. clumping during STEAL'
            ENDIF
            DENSCON(L) = DENSCON(1)
            FILLFAC(L) = 1. / DENSCON(L)
         ELSE
            FILLFAC(L) = 1. / DENSCON(L)
         ENDIF
      ENDDO
c      DENSCON(ND+1:NDDIM) = 1.
c     CALL READMS(3,DENSCON,1, 'DENSCON ', IERR)
c      IF (IERR .EQ. -10) DENSCON = 1.
c      FILLFAC = 1. / DENSCON
      CALL READMS(3,T,ND,       'T       ', IERR)
      CALL READMS (3,RADIUS,ND, 'R       ', IERR)
      CALL READMS (3,VELO ,ND,  'VELO    ', IERR)
      CALL READMS (3,GRADI,ND,  'GRADI   ', IERR)
      CALL READMS (3,RNE,ND,    'RNE     ', IERR)
      CALL READMS (3,ARAD,ND-1, 'ARAD    ', IERR)
      IF (IERR == -10) THEN
        bNoARAD = .TRUE.
      ELSE 
        bNoARAD = .FALSE.
      ENDIF
      CALL READMS (3,FTCOLI,ND, 'FTCOLI  ', IERR)
      CALL READMS (3,OPAROSS,ND,'OPAROSS ', IERR)
      CALL READMS (3,TAUROSS,ND,'TAUROSS ', IERR)
      CALL READMS (3,TAUROSScont,ND,'TAURCONT', IERR)
      IF (IERR == -10) THEN
        !not saved yet => initialize with negative value
        TAUROSScont = -1.0
      ENDIF
      
      CALL READMS (3,VMIN   ,1 ,'VMIN    ', IERR)
      IF ((IERR /= -10) .AND. (VMIN > 0.)) THEN
        VMINorg = VMIN      !transfer to VELPAR common block if saved in the MODEL
      ENDIF
      CALL READMS (3,RCON   ,1 ,'RCON    ', IERR)
      IF (IERR == -10) THEN
        !not saved yet => initialize with negative value
        RCON = -1.0
      ELSE
        RCONorg = RCON      !transfer to VELPAR common block if saved in the MODEL
      ENDIF
      CALL READMS (3,INCRIT,ND ,'INCRIT  ', IERR)
      IF (IERR == -10) THEN
        DO L=1,ND
          INCRIT(L) = ''
        ENDDO
      ENDIF
      CALL READMS (3,VCRIT,ND ,'VCRIT   ', IERR)
      IF (IERR == -10) THEN
        DO L=1,ND
          VCRIT(L) = ''
        ENDDO
      ENDIF
      
      CALL READMS (3,POPNUM,ND*N,'POPNUM  ', IERR)
      IERR=1
      CALL READMS (3,SIGMAKI,NF*LASTKON, 'SIGMAKI ', IERR)
      IF (IERR == -10) THEN
        DO K=1, NF
          DO KON=1, LASTKON
            SIGMAKI(K,KON) = 0.
          ENDDO
        ENDDO
      ENDIF

C***  Read opacity at inner boundary from MODEL file
C**   (After the first COLI, this is the full Rosseland opacity incl. lines)
      OPARND = 0.
      CALL READMS (3,OPARND,   1, 'OPARND  ', IERR)
      IF (IERR == -10) THEN
        OPARND = -99.
      ENDIF
C***  Read Eddington flux at inner boundary (HTOTND), calculated in COLI->FREQUINT
      HTOTND = 0.
      CALL READMS (3,HTOTND,   1, 'HTOTND  ', IERR)
C***  Read temperature gradient at inner boundary;
C***   this might not exist in the first iteration or in old models
      CALL READMS (3,DTDRIN,1, 'DTDRIN  ', IERR)
      IF (IERR .EQ. -10) DTDRIN = -99.

C***  READ 'XDATA' AND CHECK WHETHER THE RECORD EXISTS
      IERR=1
      CALL READMS (3,XDATA,MAXXDAT,'XDATA   ',IERR)
      IF (IERR .LT. 0) THEN
         CALL REMARK ('ERROR WHEN READING XDATA FROM MODEL FILE')
C***     XFILL EQ 0.
         XDATA(1) = 0.  
      ENDIF

C***  READ 'ABXYZ' (REL. ABUNDANCES OF ALL ELEMENTS) AND CHECK WHETHER
C***  THE RECORD EXISTS
      IERR=1
      CALL READMS (3,ABXYZ,NATOM,'ABXYZ   ',IERR)
      IF ((IERR .LT. 0) .AND. (IERR .NE. -10)) THEN
         CALL REMARK ('ERROR WHEN READING ABXYZ FROM MODEL FILE')
         STOP 'ERROR'
      ENDIF

C***  NOT EXISTING RECORD 'ABXYZ': DEFAULT IS AN ATOMIC DATA FILE "DATOM"
C***  CONTAINING "HELIUM" AS THE ONLY ELEMENT
      IF (IERR .EQ. -10) THEN
         IF (NATOM .EQ. 1) THEN
            ABXYZ(1)=1.
         ELSE
            CALL REMARK ('NOT EXISTING RECORD ABXYZ')
            STOP 'ERROR'
         ENDIF
      ENDIF
      ATMEAN=0.
      DO NA=1, NATOM
          ATMEAN = ATMEAN + ABXYZ(NA) * ATMASS(NA)
      ENDDO

      CALL READMS(3,XLAMBDA ,NF ,'XLAMBDA ', IERR)
      CALL READMS(3,XLAMBDA2,NF2,'XLAMBD2 ', IERR)
      CALL READMS(3,FWEIGHT,NF,'FWEIGHT ', IERR)
      CALL READMS (3,KEY,NF,   'KEY     ', IERR)
      CALL READMS (3,TOTOUT, 1,'TOTOUT  ', IERR)
      IF (TOTOUT .NE. 6HUNDEF. ) THEN
            CALL READMS (3,TOTIN , 1,'TOTIN   ', IERR)
            CALL READMS (3,EMFLUX,NF,'EMFLUX  ', IERR)
            CALL READMS (3,EMCOLI,NF,'EMCOLI  ', IERR)
      ENDIF
      CALL READMS (3,RSTAR , 1, 'RSTAR   ', IERR)
      RSTARSU = RSTAR / RSUN
      CALL READMS (3,XMSTAR, 1, 'XMSTAR  ', IERR)
      IF (IERR .EQ. -10) XMSTAR = .0
      CALL READMS (3,VDOP,1,   'VDOP    ', IERR)
      CALL READMS (3,VTURB,1,  'VTURB   ', IERR)
      IF (IERR .EQ. -10) VTURB = .0

C***  Variables for the hydrodynamics branch (Andreas Sander)

      CALL READMS (3,NEXTHYDRO,1,'NXTHYDRO', IERR)
      !If NXTHYDRO does not exist in the MODEL (=old model)
      ! set value to -1 to deactivate any hydro-related stuff
      IF (IERR == -10) THEN
        NEXTHYDRO = -1
      ENDIF
      CALL READMS (3,LASTHYDRO,1,'LSTHYDRO', IERR)
      IF (IERR == -10) THEN
        LASTHYDRO = -1
      ELSE
        LASTHYDRO = LASTHYDRO + 1
      ENDIF
      CALL READMS (3,ATHOM  ,ND-1, 'ATHOM   ', IERR)
      CALL READMS (3,ACONT  ,ND-1, 'ACONT   ', IERR)
      CALL READMS (3,ALPHAF ,ND  , 'ALPHAF  ', IERR)
      IF ((IERR == -10) .AND. (NEXTHYDRO == 1)) THEN
        !Failsafe if HYDRO calculation would be this turn
        ! but no force multipliers have been calculated so far
        NEXTHYDRO = 2
      ENDIF
C      CALL READMS (3,XMU    ,ND  , 'XMU     ', IERR)      
      !IF XMU does not exist in the MODEL (=old model)
      ! recalculate it here (uses ATMASS read in DATOM and ABXYZ read here)
C      IF (IERR == -10) THEN
        DO L=1, ND
          XMU(L) = ATMEAN / (1. + RNE(L))
        ENDDO
C      ENDIF
C      CALL READMS (3,RHO    ,ND  , 'RHO     ', IERR)
C      IF (IERR == -10) THEN
      !IF RHO does not exist in the MODEL (=old model)
        DO L=1, ND
          RHO(L) = AMU * ATMEAN * ENTOT(L)
        ENDDO
C      ENDIF
      CALL READMS (3,XMDOT, 1, 'XMDOT   ', IERR)
      IF (IERR == -10) THEN
        !IF MDOT does not exist in the MODEL (=old model)
        FM=ENTOT(1)*RADIUS(1)*RADIUS(1)*VELO(1)*1.E5*AMU*ATMEAN        
        XMDOT=LOG10(FM*RSTARSU*RSTARSU)-3.02
      ENDIF
      CALL READMS (3,Rcritical, 1, 'RCSAVE  ', IERR)
      IF (IERR == -10) THEN
        Rcritical = -1.
      ENDIF

      CALL READMS (3,MODHEAD,13,'MODHEAD ', IERR)
      CALL READMS (3,JOBNUM,1, 'JOBNUM  ', IERR)
      JOBNUM=JOBNUM+1
C      IF (JOBNUM .GE. 1000) JOBNUM=JOBNUM-1000
      CALL WRITMS (3,JOBNUM,1,'JOBNUM  ',-1, IDUMMY, IERR)
 
      DO K=1,NF
        WRITE (NAME, '(A3, I4, A1)') 'XJC', K, ' '
        CALL READMS(3,XJC(1+ND*(K-1)),ND,NAME, IERR)
        WRITE (NAME, '(A3, I4, A1)') 'WJC', K, ' '
        CALL READMS(3,WJC(1+ND*(K-1)),ND,NAME, IERR)
      ENDDO

      DO IND=1, LASTIND
        NUP=INDNUP(IND)
        LOW=INDLOW(IND)
        IF (EINST(LOW,NUP) == -2.) CYCLE
        IF (IND <= 9999) THEN
          WRITE (NAME, '(A3, I4, A1)') 'XJL', IND, ' '
        ELSE
          WRITE (NAME, '(A3, I5)') 'XJL', IND
        ENDIF
        CALL READMS (3,XJL(1+ND*(IND-1)),ND,NAME, IERR)
      ENDDO
 
C***  IRON: READIN XJF'S 
      DO INDFE=1,LASTFE
        IND = LASTIND-LASTFE+INDFE
        WRITE (NAME, '(A3, I4, A1)') 'WFL', INDFE, ' '
        CALL READMS (3,WFELOW(1+ND*(IND-1)),ND,NAME, IERR)
        WRITE (NAME, '(A3, I4, A1)') 'WFU', INDFE, ' '
        CALL READMS (3,WFENUP(1+ND*(IND-1)),ND,NAME, IERR)
      ENDDO
 
      CALL READMS (3,LAST,1,'MODHIST ', IERR)
      IF (LAST.GT.MAXHIST) THEN
        WRITE (hCPR,*) 'LAST=',LAST, ' MAXHIST=',MAXHIST
        CALL REMARK ('MODEL HISTORY DIMENSION INSUFFICIENT')
        STOP 'ERROR'
      ENDIF
      CALL READMS (3, MODHIST, LAST, 'MODHIST ', IERR)

C***  Read GAMMA HISTORY
      CALL READMS (3, GAHIST, 26*MAXGAHIST, 'GAHIST  ', IERR)

C***  Read GF
      CALL READMS (3, GF, 1, 'GF      ', IERR)

      CALL READMS (3, TEFF   ,  1  , 'TEFF    ', IERR)
      CALL READMS (3, GLOG   ,  1  , 'GLOG    ', IERR)
      IF (IERR == -10) THEN
        IF (XMSTAR <= 0.) THEN
          CALL REMARK('NEITHER MSTAR NOR LOG G IS STORED -' // 
     >                'TRYING TO READ FROM CARDS')
          XMSTAR = -99.
          GLOG = -99.
        ELSE
          !GEFFLOG not in model file => calculate from XMSTAR
          GLOG = ALOG10(GCONST * XMSTAR * XMSUN / RSTAR / RSTAR)
        ENDIF
      ELSEIF (XMSTAR <= 0.) THEN
        XMSTAR = 10**(GLOG) * RSTAR * RSTAR / GCONST / XMSUN
      ENDIF
      CALL READMS (3, GEFFLOG,  1  , 'GEFFLOG ', IERR)
      IF (IERR == -10) THEN
        !GEFFLOG not in model file => mark as not found (negative value)
        GEFFLOG = -1.0  
        bFixGEFF = .FALSE.
      ELSEIF (GEFFLOG < 0.) THEN
        IF (GEFFLOG /= -1.0) THEN
          GEFFLOG = -1. * GEFFLOG
        ENDIF
        bFixGEFF = .FALSE.
      ELSE
        bFixGEFF = .TRUE.
      ENDIF
      CALL READMS (3, GEFFKEY,  1  , 'GEFFKEY ', IERR)
      IF (IERR == -10) THEN
        GEFFKEY = '        '
      ENDIF
      
      CALL READMS (3, GEDD,  1  , 'GEDD    ', IERR)
      IF (IERR == -10) THEN
        !GEDD not in model file => mark as not found (negative value)
        GEDD = -1.0  
      ENDIF
      IF (GEDD < 0) THEN
        bGEddFix = .FALSE.
      ELSE
        bGEddFix = .TRUE.
      ENDIF

      CALL READMS (3, GEDDRAD,  1  , 'GEDDRAD ', IERR)
      IF (IERR == -10) THEN
        !GEDDRAD not in model file => mark as not found (negative value)
        GEDDRAD = -1.0  
      ENDIF
      !Depth-dependend GEDDRAD is called GRSTATIC
      ! (because it is GEDDRAD used for the hydrostatic part)
      CALL READMS (3, GRSTATIC, ND  , 'GRSTATIC', IERR)
      IF (IERR == -10) THEN
        GRSTATIC = -1.0
      ENDIF
      

C***   TRYING TO READ LINE RADIATION FIELD FOR DR-TRANSITIONS 
C***   IF ONE RECORD DOES NOT EXIST (IERR=-10), THE LOGICAL DRXJL(INDDR) 
C***   IS SET TO .FALSE.
C***   RUDIMENTAL DR-TRANSITIONS ARE EXCLUDED
C***
      DO 17 INDDR=1,NAUTO
      DRXJL(INDDR) = .TRUE.
      IF (KRUDAUT(INDDR) .EQ. 1) GOTO 17
      IERR=1
      IF (LASTIND+INDDR <= 9999) THEN
        WRITE (NAME, '(A3, I4, A1)') 'XJL', LASTIND+INDDR, ' '
      ELSE
        WRITE (NAME, '(A3, I5)') 'XJL', LASTIND+INDDR
      ENDIF
      CALL READMS (3,XJL(1+ND*(INDDR+LASTIND-1)),ND,NAME,IERR)
      IF ((IERR .LT. 0) .AND. (IERR .NE. -10)) THEN
         CALL REMARK ('ERROR WHEN READING DR-XJL FROM MODEL FILE')
         STOP 'ERROR'
      ELSE IF (IERR .EQ. -10) THEN
         DRXJL(INDDR)=.FALSE.
      ENDIF
   17 CONTINUE


C***  TRY TO READ HTOTL ARRAY gives errormessage if there is no HTOTL
      IERR = 0
      CALL READMS(3, HTOTL, ND, 'HTOTL   ', IERR)
      BHTOTERR = (IERR .LT. 0)
      CALL READMS(3,XJTOTL, ND, 'JTOTL   ', IERR)
      CALL READMS(3,XKTOTL, ND, 'KTOTL   ', IERR)
      CALL READMS(3,XNTOTL, ND, 'NTOTL   ', IERR)

C***  Read NCOLIP
      CALL READMS(3, NCOLIP,        1,        'NCOLIP  ', IERR)

C***  Unsoeld-Lucy-Terms for TEMPEQ
      CALL READMS (3,OPASMEAN,      ND, 'OPASMEAN', IERR)
      CALL READMS (3,OPASMEANTC,    ND, 'OPASMTC ', IERR)
      CALL READMS (3,QFJMEAN,       ND, 'QFJMEAN ', IERR)
      CALL READMS (3,OPAJMEAN,      ND, 'OPAJMEAN', IERR)
      CALL READMS (3,OPAJMEANTC,    ND, 'OPAJMTC ', IERR)
      CALL READMS (3,SMEAN,         ND, 'SMEAN   ', IERR)

C***  Newly introduces Planck mean opacity, see PoWR-Memo 20160708
      CALL READMS (3,OPAPMEAN,      ND, 'OPAPMEAN', IERR)
      IF (IERR .EQ. -10) THEN
         DO L=1, ND
            OPAPMEAN(L) = OPASMEAN(L)
         ENDDO
      ENDIF

      CALL READMS (3,QOPAHMEAN,   ND-1, 'QOPAHMEA', IERR)
      CALL READMS (3,EDDIHOUTJMEAN,  1, 'EDDIHOJM', IERR)
      CALL READMS (3,HTOTOUTMINUS ,  1, 'HTOTOUTM', IERR)
      CALL READMS (3,OPALAMBDAMEAN, ND, 'OPALMEAN', IERR)

C***  Flags for the POPMIN levels
      CALL READMS (3,ZERO_RATES,  N*ND, 'ZERO_RAT', IERR)
C*    Default if variable does not exist yet 
      IF (IERR .EQ. -10) THEN
        DO I=1, N*ND
         ZERO_RATES(I) = .FALSE.
        ENDDO
      ENDIF

C***  Read further counters
      CALL READMS(3, LASTBACKUP, 1, 'LASTBAK ', IERR)
      IF (IERR == -10) THEN
        LASTBACKUP = 0
      ENDIF
      CALL READMS(3, LASTTAU   , 1, 'LASTTAU ', IERR)
      IF (IERR == -10) THEN
        LASTTAU = -1
      ELSE
        LASTTAU = LASTTAU + 1      
      ENDIF      
      
C***  Prepare vector with clump density stratification
      DO L=1, ND
        ENTOTDENS(L) = ENTOT(L) * DENSCON(L)
      ENDDO

      IFF_N_MS = IFF_MAX_MS     !overwrite IFF_N_MS with IFF_MAX_MS to avoid any problems in other routines

      RETURN
      END
