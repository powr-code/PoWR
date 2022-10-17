      SUBROUTINE LINPOP (T,RNE,ENTOT,ITNE,POPNUM,DEPART,POPLTE,POP1,
     >   N,ENLTE,WEIGHT,NCHARG,EION,ELEVEL,EN,EINST,LEVEL,
     >   FWEIGHT,XJC,NF,XJL,WCHARM,SCOLD,
     >   XJLAPP,
     >   DM, V1,V2,GAMMAC,GAMMAL,EPSILON,
     >   NOTEMP,DELTAC,GAMMAR,IPRICC,IPRILC,MODHEAD,JOBNUM,IFRRA,ITORA,
     >   RADIUS,RSTAR,OPA,ETA,THOMSON,IWARN,MAINPRO,MAINLEV,
     >   VELO,GRADI, VDOP,NFLDIM,PHI,PWEIGHT,
     >   LASTIND,LAINDHI,
     >   SIGMAKI,
     >   ND,LSRAT,CRATE,RRATE,RATCO,ALPHA,SEXPO, 
     >   ADDCON1, ADDCON2, ADDCON3, 
     >   IGAUNT, KEYCBB, NRB_CONT, ZERO_RATES, IPRINTZERORATES, POPMIN, 
     >   LINE,ALTESUM,NFEDGE,NOM,NATOM,ABXYZ,KODAT,NFIRST,
     >   SCRATCH,ATEST,BTEST,DTEST,BRS,BRSP,BRY,XMAX,IBLENDS,MAXLAP,
     >   XLAMZERO,ITBR,NBRCANC,NLAST,SIGMAFF,MAXION,NFL,IONGRND,IONAUTO,
     >   NAUTO,MAXAUTO,LOWAUTO,WAUTO,EAUTO,AAUTO,DRRATEN,IADR16,RDIEL,
     >   RAUTO,DRJLW,DRJLWE,DRLJW,NSCHAR,MAXATOM,BETA,PHIL,NBLENDS,
     >   SIGMATHK,SEXPOK,EDGEK,XDATA,
     >   KONTNUP,KONTLOW,LASTKON,KODRNUP,KODRLOW,LASTKDR,KEYCBF,NATOUT,
     >   BRRESET,NOCON,ENOLD,TEFF,
     >   KRUDAUT, ELEMENT, 
     >   BAUTO, BAUTO_ABORT, CKONVER, SMALLPOP, BUNLU, 
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
     >   NF2, NDDIM, 
     >   XKMIN, XKMAX, XKMID, XKRED, XKRED_CORE, XKBLUE_CORE, 
     >   BXJLAPPNEW, BXJCAPPNEW, BNEWOPER, BXJLAPPCORE,
     >   XJL_PLOTDATA, XJC_PLOTDATA_I, 
     >   IPLOT_XJLAPP, IPLOT_XJCAPP, LPLOT_XJCAPP, NITER_PLOT_JAPP, 
     >   BPLOTAPP, PWEIGHTCL, WS, BRDIVFAILED, BNEWTONRESET, 
     >   IWARN_NEG_XJCAPP, IWARN_NEG_XJLAPP,
     >   XJLAPPNEW, 
     >   XLAM_FINE_START, XLAM_FINE_END, RUDLINE, IMAXPOP,
C***  New Fine-spaced WCHARM handling
     >   IFF_MAX, IFF_MAX_MS, FF_INFO, IFF_DK, IFF_WCHARM, WCHARM_FINE, 
     >   IFF_N_MS,
C***  Split artistic
     >   bBLOCKINVERSION, AOFF, CORRLAST, bRELSCHARMERACC,
C***  Settings for ZERO_RATES (not yet implemented), BROYDEN TWOPOINT, and AUTO_MODIFY
     >   iZRType, bUseTWOPNT, iAMT,
C***  Fe rate temperature correction
     >   CORRS, DEXFAC, bFeTCORR,
C***  Contants consitence check
     >   NDIM_steal, MAXIND_steal, NFDIM_steal, MAXFEIND_steal)

C*******************************************************************************
C***  CALCULATION OF NEW NLTE POPULATION NUMBERS (ARRAY POPNUM)
C***  RADIATIVE RATES ARE CALCULATED WITH THE SCHARMER RADIATION FIELD
C***  THE (HENCE NON-LINEAR) RATE EQUATIONS ARE SOLVED BY LINEARIZATION
C***  POP1 = OLD POPULATION NUMBERS
C***  POPNUM = NEW POPULATION NUMBERS
C***  RNE = REL. ELECTRON DENSITY - UPDATED IN THIS SUBROUTINE
C***  EN(J) = NEW NLTE POP. NUMBERS  AT CURRENT DEPTH POINT
C***  ENLTE(J) = LTE POP. NUMBERS AT CURRENT DEPTH POINT
C***  --- BROYDEN METHOD USED IF REQUESTED !  ----------
C***  The option for solving energy equation as N+2nd equation became  
C***  obsolete and was removed (28-Oct-2002, wrh). 
C***  Throughout LINPOP, T(L) is the actual (TNEW) temperature after it 
C***  was updated by TEMPCORR. This might be considered to be changed, 
C***  i.e. SOLD might be better estimated with TOLD. (???)
C*******************************************************************************

      IMPLICIT NONE

      INTEGER, PARAMETER :: MAXIND  =       20000 
      INTEGER, PARAMETER :: NDIM    =        1560 
      INTEGER, PARAMETER :: NFDIM   = 2*NDIM + 400

      INTEGER, PARAMETER :: MAXFEIND  =       2500 

      INTEGER, INTENT(IN) :: ND, NDDIM, N, NF, JOBNUM, MAXION, NFLDIM,
     >                       MAXLAP, NATOM, LASTKDR, MAXFINE, INDEXMAX,
     >                       IFF_MAX, IFF_MAX_MS, MAXATOM, LASTKON,
     >                       LASTIND, MAXAUTO, LAINDHI, NKONV_THRESH

      INTEGER, DIMENSION(MAXFEIND) :: INDRB, INDRF, IFRBSTA, IFRBEND,
     >                                IFENUP, IFELOW, INDFEACT
      REAL, DIMENSION(MAXFEIND) :: SIGMAACT, SIGMAINT, FERATLU,
     >                             FERATUL, FERATLU0, FERATUL0

      COMMON /IRON/ INDRB, INDRF, IFRBSTA, IFRBEND, IFENUP, IFELOW,
     >              INDFEACT, SIGMAACT, SIGMAINT, FERATLU, FERATUL, 
     >              FERATLU0, FERATUL0

      INTEGER, DIMENSION(MAXIND) :: INDNUP, INDLOW
      REAL, DIMENSION(MAXIND) :: XRED, XBLUE, SLOLD, DETAL, OPAL, 
     >                           SLNEW, DOPAL, SCOLIND, SCNEIND,
     >                           OPACIND, XRED0, XBLUE0, OPALOLD
      REAL, DIMENSION(4, MAXIND) :: COCO

      COMMON /COMIND/  INDNUP, INDLOW, XRED,XBLUE, SLOLD, DETAL, OPAL,
     > SLNEW, DOPAL,SCOLIND,SCNEIND,OPACIND,XRED0,XBLUE0,COCO,OPALOLD

      REAL, DIMENSION(NFDIM) :: XLAMBDA, XLAMBDA2, XJCAPP, OPAC, SCNEW,
     >                          DOPA, DETA, ETAC, EXPFAC, XJCLP1, 
     >                          OPAC1, XKC, XKC2, XJCAPPNEW, FWTEST
      REAL, DIMENSION(NFDIM,4) :: XJC_PLOTDATA_L

      COMMON /COMNF/  XLAMBDA, XLAMBDA2, XJCAPP, OPAC, SCNEW, DOPA, 
     >                DETA, ETAC, EXPFAC, XJCLP1, OPAC1, XKC, XKC2, 
     >                XJCAPPNEW, XJC_PLOTDATA_L, FWTEST

      REAL, DIMENSION(ND) :: RADIUS, ENTOTDENS
      REAL, DIMENSION(N) :: ENOLD, WEIGHT
      INTEGER, DIMENSION(1) :: IGAUNT   !dies ist eigentlich ein CHARACTER-Array

      REAL, DIMENSION(NDDIM) :: OPAFE, ETAFE
      REAL, DIMENSION(ND, N) :: POPNUM, POP1, POPLTE, DEPART
      REAL, DIMENSION(ND) :: ENLTE, RNE, T, ENTOT, OPA, ETA, VELO,
     >                       GRADI, THOMSON, CORRS
      INTEGER, DIMENSION(ND) :: ITNE, IWARN
      REAL, DIMENSION(NDIM) :: ELEVEL, EION, DRRATEN, RDIEL, RAUTO,
     >                         DRJLW, DRJLWE, DRLJW
      REAL, DIMENSION(NDIM+2) :: V1, V2, EN, BRS, BRSP, BRY, DTEST
      REAL, DIMENSION(4, NDIM) :: ALTESUM
      INTEGER, DIMENSION(NDIM) :: NCHARG, IONGRND
      REAL, DIMENSION(NDIM, NDIM) :: EINST, CRATE, RRATE
      REAL, DIMENSION(NDIM+2,NDIM+2) :: RATCO, DM, AOFF, 
     >                                  ATEST, BTEST, SCRATCH

      REAL, DIMENSION(ND,LAINDHI) :: XJL
      REAL, DIMENSION(ND, 5) :: XJL_PLOTDATA
      REAL, DIMENSION(ND, 4) :: XJC_PLOTDATA_I

      REAL, DIMENSION(MAXAUTO) :: KRUDAUT, WAUTO, EAUTO, AAUTO
      INTEGER, DIMENSION(MAXAUTO) :: IONAUTO, LOWAUTO
      INTEGER, DIMENSION(LASTKON) :: KONTLOW, KONTNUP, KONTHLP, 
     >                               NFEDGE, KEYCBF, SIGMA1I

      REAL, DIMENSION(NFDIM) :: XJLAPPNEW
      REAL, DIMENSION(NF,LASTKON) :: SIGMAKI
      REAL, DIMENSION(NF,ND) :: SCOLD
      REAL, DIMENSION(NF,0:MAXION) :: SIGMAFF
      REAL, DIMENSION(NF) :: FWEIGHT
      REAL, DIMENSION(MAXATOM,MAXATOM) :: EDGEK, SIGMATHK, SEXPOK
      REAL, DIMENSION(MAXATOM) :: ABXYZ
      INTEGER, DIMENSION(MAXATOM) :: IMAXPOP
      INTEGER, DIMENSION(NATOM) :: NFIRST, NLAST, KODAT

      REAL, DIMENSION(N+2) :: EN1
      INTEGER, DIMENSION(N) :: NOM
      REAL, DIMENSION(ND,LASTIND) :: WFELOW, WFENUP
      LOGICAL :: BRDIVFAILED, BMONOTONIC, bFeTCORR
      LOGICAL, DIMENSION(LASTIND) :: BDIAG, RUDLINE
      LOGICAL, DIMENSION(LASTKON) :: NRB_CONT
      LOGICAL, DIMENSION(MAXIND) :: LINE, BLASERL
      INTEGER, DIMENSION(MAXIND) :: LINEINDEX, NUPACT, LOWACT
      REAL, DIMENSION(MAXIND) :: XLAMZERO, XLAMSOR, XLAMMIN, XLAMMAX,
     >                           ETAL, XKMIN, XKMAX, XKMID, XKRED, 
     >                           DEXFAC
      REAL, DIMENSION(LASTIND) :: XJLAPP
      INTEGER, DIMENSION(LASTIND) :: NBLENDS
      INTEGER, DIMENSION(LASTKDR) :: KODRNUP, KODRLOW
      REAL, DIMENSION(ND,NF) :: WJC, WCHARM, XJC

      REAL, DIMENSION(NFLDIM) :: PHI, PWEIGHT
      REAL, DIMENSION(NFLDIM, MAXLAP) :: PHIL
      REAL, DIMENSION(MAXLAP) :: PWEIGHTCL, WS, BETA,         !Note: BETA /= velocity beta
     >                           XKRED_CORE, XKBLUE_CORE   
      INTEGER, DIMENSION(MAXLAP) :: LIND, LINDS
      INTEGER, DIMENSION(MAXLAP,MAXIND) :: IBLENDS

      REAL, DIMENSION(MAXFINE) :: SFINE_OLD, SFINE_NEW
      REAL, DIMENSION(INDEXMAX) :: SIGMAFE

      REAL, DIMENSION(1) :: ALPHA, SEXPO, ADDCON1, ADDCON2, ADDCON3, 
     >                      XDATA

      CHARACTER(4) :: VERSION
      CHARACTER(4), DIMENSION(LASTIND) :: KEYCBB
C      CHARACTER*(*), DIMENSION(N) :: LEVEL
      CHARACTER(10), DIMENSION(N) :: LEVEL
      CHARACTER(10), DIMENSION(ND) :: MAINPRO, MAINLEV
      CHARACTER*10 ELEMENT(NATOM)

      LOGICAL, DIMENSION(N,ND) :: ZERO_RATES
      INTEGER, DIMENSION(NDIM) :: NZERORATES, LMINZERORATES, LMAXZERORATES

      REAL, DIMENSION(10) :: FF_INFO
      INTEGER, DIMENSION(IFF_MAX) :: IFF_DK, IFF_WCHARM
      REAL, DIMENSION(IFF_MAX) :: WCHARM_FINE
      
      INTEGER :: L, I, J, K, ITMAX, ITWARN, NRANK, IND, IPRICC, IPRILC,
     >           NPLUS1, NPLUS2, NOUT, ITBR, NBRCANC, NKONVER, IADR16,
     >           INEGMIN, INEGMAX, LNEGMIN, LNEGMAX, NEGINTL, NEGINTC,
     >           INEGMINC, INEGMAXC, LNEGMINC, LNEGMAXC, NFL, KANAL, LL,
     >           KON, NUP, LOW, LPLOCC, KPLOCC, ICOUNT, NA, NFIRNA, NF2,
     >           NLANA, NSCHAR, ION, NAUTO, LASTFE, NLINE, MAXFEACT,
     >           IPLOT_XJCAPP, IPLOT_XJLAPP, LPLOT_XJCAPP, ITABS, IERR,
     >           NITER_PLOT_JAPP, IWARN_NEG_XJCAPP, IWARN_NEG_XJLAPP,
     >           IDUMMY, ldepth, LSRAT, IFRRA, ITORA, NETTO, LM1, 
     >           MINEX, MAXEX, IPRINTZERORATES, NATOUT,
     >           NDIM_steal, MAXFEIND_steal, NFDIM_steal, MAXIND_steal,
     >           iPopMinViolations, iZRType, iAMT
      REAL :: RSTAR, GAMMAC, GAMMAD, GAMMAL, GAMMAR, DELTAC, EPSILON, 
     >        EPSDN, TL, TOLD, EDGE, EDGELAM, ERXMIN, WJCMIN, DELTAX,
     >        VDOP, VDOPFE, WAVENUM, XMAX, OPARL, FN, FNOLD, 
     >        FTN, FTNOLD, TSAVE, TMID, DTDR, TEFF, TLAST,  
     >        ENE, BRRESET, BRRES2, TLOG, ROOTTL, XLAM, XLAMLOG, W, W3,
     >        PRESIG, GIII, POPMIN, OPATHOM, DXFE, XLAM0FE, DFEINDR,
     >        SUM, ENMIN, ENMAX, ENNEW, XLAM_FINE_START, XLAM_FINE_END,
     >        IFF_N_MS, SMALLPOP, CORRLAST      

C***  Iteration Sequence output ("Scharmer-Tapete")
      CHARACTER(1), DIMENSION(100) :: ITCODE        !@todo: warum nicht auf ND dimensionieren?

C***  tiefenabh. clumping nach goetz
      REAL, DIMENSION(NDDIM) :: DENSCON, FILLFAC

C***  AUTO_MODIFY FACILITY
      CHARACTER(3) :: CHR
      CHARACTER(1), DIMENSION(ND) :: CKONVER

      CHARACTER(3) :: STRING3
      CHARACTER(8) :: NAME, LASTIT, LASTIT2, LASTIT3
      CHARACTER(15) :: STRING15
      CHARACTER(60) :: TEXT
      CHARACTER(80) :: CHR1, CHR2
      CHARACTER(100) :: MODHEAD
      LOGICAL :: KONVER, NOTEMP, BROYDEN, TWOPNT, NOCON, BNUEFE,
     >           BRIMP, BRIMP2, BAUTO, BAUTO_ABORT, BNEWTONRESET,
     >           BSUM, BUNLU, BPLOCC, BFECHECK, BFEWING, BFEMODEL,
     >           BXJLAPPNEW, BXJCAPPNEW, BXJLAPPCORE, BNEWOPER,
     >           BPLOTAPP, bFFASSET, bRELSCHARMERACC, bUseTWOPNT 

C***  Split switch
      LOGICAL :: bBLOCKINVERSION, bUseDMFILE

C***  Warning counter for negative line opacities --> PRICORR
      COMMON / COMNEGI / NEGINTL,INEGMIN,INEGMAX,LNEGMIN,LNEGMAX
C***  Warning counter for negative continuum opacities --> PRICORR
      COMMON / COMNEGC / NEGINTC,INEGMINC,INEGMAXC,LNEGMINC,LNEGMAXC
      COMMON / COMITWA / ITWARN, ITMAX

      CHARACTER(8) :: OPC

      INTEGER, EXTERNAL :: IDX, ISMAX, ISRCHFGT
      REAL, EXTERNAL :: ERF

      !File and channel handles (=KANAL)
      INTEGER, PARAMETER :: hCPR = 0        !write to wruniqX.cpr (stderr)
      INTEGER, PARAMETER :: hOUT = 6        !write to wruniqX.out (stdout)
      INTEGER, PARAMETER :: hDMFILE = 16    !DMFILE file (fort.17)
      INTEGER, PARAMETER :: hFF = 27        !file handle for FFASSET
      INTEGER, PARAMETER :: hMODIFY = 52    !MODIFY_INPUT file

C***  Physical Constants
      REAL, PARAMETER :: C1 = 1.4388        !C1 = H * C / K    ( CM * KELVIN )
      REAL, PARAMETER :: CFF = 1.370E-23    !CFF = COEFFICIENT FOR FREE-FREE CROSS SECTION ( ALLEN P.100 )

      bUseDMFILE = .TRUE.                   !Default: DMFILE is used
      iPopMinViolations = 0

C***  Check for consitency of the constants
      IF (NDIM /= NDIM_steal) THEN
        WRITE (hCPR,*) 'LINPOP> Constant NDIM not consistent!'
        WRITE (hCPR,*) ' STEAL vs. LINPOP: ', NDIM_steal, NDIM
        STOP 'FATAL ERROR IN STEAL'
      ENDIF
      IF (MAXIND /= MAXIND_steal) THEN
        WRITE (hCPR,*) 'LINPOP> Constant MAXIND not consistent!'
        WRITE (hCPR,*) ' STEAL vs. LINPOP: ', MAXIND_steal, MAXIND
        STOP 'FATAL ERROR IN STEAL'
      ENDIF
      IF (NFDIM /= NFDIM_steal) THEN
        WRITE (hCPR,*) 'LINPOP> Constant NFDIM not consistent!'
        WRITE (hCPR,*) ' STEAL vs. LINPOP: ', NFDIM_steal, NFDIM
        STOP 'FATAL ERROR IN STEAL'
      ENDIF
      IF (MAXFEIND /= MAXFEIND_steal) THEN
        WRITE (hCPR,*) 'LINPOP> Constant MAXFEIND not consistent!'
        WRITE (hCPR,*) ' STEAL vs. LINPOP: ', MAXFEIND_steal, MAXFEIND
        STOP 'FATAL ERROR IN STEAL'
      ENDIF
      
C***  Operator Versions:
      WRITE (hCPR,*)
      WRITE (hCPR,'(5A12)')
     >  'Operator:', 'Continuum', 'Lines', 'R-Lines', 'Fe-Lines'
      WRITE (hCPR,'(2A12)')
     >  ' ', OPC
      WRITE (hCPR,*)


C***  DEMANDED ACCURACY OF THE NEWTON-RAPHSON ITERATION
      IF (bRELSCHARMERACC) THEN
        EPSDN = MAX(EPSILON, CORRLAST) * 0.03
      ELSE
        EPSDN=EPSILON * 0.01
      ENDIF
C     goetz uses this method to determine the required ACCURACY
C      EPSDN = MAX(EPSILON, CORRLAST) * 0.03
C     but CORRLAST does not exist in the wrh branch

C--------------- additional stuff from goetz ---------------------
C***  Increase SMPOP for outer Depth Points
      NOUT = 20
      IF (NOUT .GT. ND/2) NOUT = ND/2
      
CC***  Set T and TNEW:  TEMPEQ Branch is dead --> NOTEMP = .TRUE.
CC***  ---------------
C      NOTEMP = .TRUE.
CC***  No T-Corrections:
C      IF (.NOT. BUNLU) TNEW = T

CC***  For TNEW:
CC***  RECOVER THE DIFFUSION APPROX. AT THE INNER BOUNDARY
C      IF (BTND) THEN
C         TNDM  = TNEW(ND-1)
C         TND   = TNEW(ND)
C         TSAVE = TND
C
C         IF (OPARND .GT. 0.) THEN
C            OPARL = OPARND
C            WRITE (0,*) 'LINPOP: OPARND=', OPARND
C         ELSE
C            CALL OPAROSS(OPARL,EN,TND,RNE(ND),ENTOTDENS(ND),RSTAR,
C     >                  NDIM,N,LEVEL,NCHARG,WEIGHT,ELEVEL,EION,EINST,
C     >                  ALPHA,SEXPO,ADDCON1, ADDCON2, ADDCON3, 
C     >                  IGAUNT,NF,XLAMBDA,FWEIGHT,NOM,
C     >                  MAXATOM,SIGMATHK,SEXPOK,EDGEK,KODAT,RADIUS(ND),
C     >                  KONTNUP,KONTLOW,LASTKON,BXSFINE,NXAPP,XSAPP,X0)
CC***        Scale opacity down with filling factor
C            OPARL = OPARL * FILLFAC(ND)
C         ENDIF
C
C         TBAR = TNDM / TEFF
C         DR   = RADIUS(ND-1)-RADIUS(ND)
C            
C         PP   = TBAR / DR / 3.
C         QQ   = OPARL / DR / TBAR / TBAR / 8.
C         DTDR = TEFF * (-PP +SQRT(PP*PP + QQ))
C            
C         TND  = TNDM + DR * DTDR
C            
C         TNEW(ND)= TND
C         TNDCORR = (TND-TSAVE) / TND
C      ELSE
C
C         TNDCORR = 0.
C      ENDIF
C--------------- end of additional stuff from goetz --------------


      NPLUS1=N+1
      NPLUS2=N+2
C***  RANK OF THE SYSTEM IS NOWALWAYS N+1 (no energy equation)
      NRANK=NPLUS1

C***  Note:wrh SPLIT option does not reducde the rank (unlike goetz!)
C***  implemented since May 2013 
C***  bBLOCKINVERSION is set via CARDS line SPLITINVERSION
 
C***  ITERATION NUMBER WHERE BROYDEN UPDATING STARTS (MUST BE .GE. 1)
      IF (ITBR .EQ. 0) ITBR=ITMAX+1
      IF (ITBR .LT. 1) ITBR = 1
C***  COUNTER OF DEPTH POINTS WHERE BROYDEN TREATMENT IS CANCELLED
      NBRCANC = 0

C***  INITIALIZE THE VARIABLES FOR AUTO_MODIFY
      NKONVER = 0
      DO I=1, ND
        CKONVER(I) = 'C'
      ENDDO

C***  CHECK AND OPEN "DMFILE". IF NOT USEFULL, SET ITBR=2
C      CALL DMOPEN (IADR16, ND, MODHEAD, JOBNUM, NRANK, ITBR)
 
C***  goetz version:
      IF (ITBR > 1) bUseDMFILE = .FALSE.
      IF (bUseDMFILE)
     >   CALL DMOPEN (IADR16, ND, MODHEAD, JOBNUM, NRANK, ITBR)

 
C***  INITIALIZE COUNTER AND LOGICAL FOR NOT CONVERGED 
C***      NEWTON-RAPHSON ITERATIONS
      ITWARN=0 
      NOCON = .FALSE.

C***  REMOVE NEGATIVE LINE INTENSITIES
      NEGINTL=0
      INEGMIN=LAINDHI
      INEGMAX=0
      LNEGMIN=ND
      LNEGMAX=0
      DO L=1,ND
        DO IND=1,LASTIND
          IF (EINST(INDLOW(IND),INDNUP(IND)) .EQ. -2.) EXIT
          IF (XJL(L,IND) .LT. .0) THEN
            XJL(L,IND) = .0
            NEGINTL=NEGINTL+1
            IF (IND .LT. INEGMIN) INEGMIN = IND
            IF (IND .GT. INEGMAX) INEGMAX = IND
            IF (L   .LT. LNEGMIN) LNEGMIN = L
            IF (L   .GT. LNEGMAX) LNEGMAX = L
          ENDIF
        ENDDO
        DO IND=LASTIND+1,LAINDHI
          IF (KRUDAUT(IND) .EQ. 1) EXIT
          IF (XJL(L,IND) .LT. .0) THEN
            XJL(L,IND) = .0
            NEGINTL=NEGINTL+1
            IF (IND .LT. INEGMIN) INEGMIN = IND
            IF (IND .GT. INEGMAX) INEGMAX = IND
            IF (L   .LT. LNEGMIN) LNEGMIN = L
            IF (L   .GT. LNEGMAX) LNEGMAX = L
          ENDIF
        ENDDO
      ENDDO
 
C***  REMOVE NEGATIVE CONTINUUM INTENSITIES
C     (nicht im goetz branch vorhanden)
      NEGINTC=0
      INEGMINC=NF+1
      INEGMAXC=0
      LNEGMINC=ND+1
      LNEGMAXC=0
      DO L=1,ND
         DO K=1, NF
            IF (XJC(L,K) .LT. .0) THEN
               XJC(L,K) = .0
               NEGINTC=NEGINTC+1
               IF (K .LT. INEGMINC) INEGMINC = K
               IF (K .GT. INEGMAXC) INEGMAXC = K
               IF (L .LT. LNEGMINC) LNEGMINC = L
               IF (L .GT. LNEGMAXC) LNEGMAXC = L
            ENDIF
         ENDDO
      ENDDO
 
C***  GENERATE ONCE FOR ALL PHOTOCROSSSECTIONS AT ALL FREQUENCIES
C***  SIGMAKI(K,LOW) IN CM**2
      CALL       BFCROSS (SIGMAKI,NF,N,ELEVEL,EION,EINST,NDIM,
     >                    XLAMBDA,ALPHA,SEXPO,
     >                    ADDCON1, ADDCON2, ADDCON3, 
     >                    IGAUNT,
     >                    KONTNUP,KONTLOW,LASTKON)
 
C***  DETERMINE THE WEIGHT FUNCTION WCHARM FOR THE CONTINUUM, AND SCOLD
C***  AT ALL DEPTH POINTS
      CALL CCORE (WCHARM,NF,GAMMAC,DELTAC,IPRICC,MODHEAD,JOBNUM,
     >            SCOLD,RADIUS,XLAMBDA,ND,T,RNE,POP1,ENTOTDENS,RSTAR,
     >            OPA,ETA,THOMSON,IWARN,MAINPRO,MAINLEV,NOM,
     >            NDIM,N,LEVEL,NCHARG,WEIGHT,ELEVEL,EION,EINST,SIGMAKI,
     >            MAXATOM,SIGMATHK,SEXPOK,EDGEK,KODAT,XDATA,
     >            KONTNUP,KONTLOW,LASTKON,FILLFAC, WJC, OPC, WJCMIN, 
     >            BPLOCC, LPLOCC, KPLOCC, KANAL)
 
C***  GENERATE LINE FREQUENCY GRID AND INTEGRATION WEIGHTS
      CALL FLGRID (NFLDIM,NFL,PHI,PWEIGHT,DELTAX,XMAX)
      ERXMIN=ERF(-XMAX)
 
C***  PRE-CALCULATE FREQUENCY INDICES OF IONIZATION EDGES
      DO KON=1, LASTKON
        NUP=KONTNUP(KON)
        LOW=KONTLOW(KON)
        EDGE=ELEVEL(NUP)+EION(LOW)-ELEVEL(LOW)
        EDGELAM=1.E8/EDGE
        NFEDGE(KON)=ISRCHFGT (NF,XLAMBDA,1,EDGELAM) - 1
      ENDDO
 
C***  GENERATE LOGICAL VARIABLE AS SWITCH FOR RUDIMENTAL LINES
C***         --- OPTIMIZATION OF SUBR. DERIV ---
      DO IND=1, LASTIND
        NUP=INDNUP(IND)
        LOW=INDLOW(IND)
        RUDLINE(IND) = EINST(LOW,NUP) .EQ. -2.
      ENDDO

C***  Initialize counters for ZERO_RATES
      DO J=1, N
        NZERORATES(J) = 0
        LMINZERORATES(J) = ND+1
        LMAXZERORATES(J) = 0
      ENDDO

      IF (bBLOCKINVERSION) THEN
        WRITE (hCPR,'(A)') ' DM Inversion: SPLIT'
      ELSE
        WRITE (hCPR,'(A)') ' DM Inversion: TOTAL'
      ENDIF
      IF (bFeTCORR) THEN
        WRITE (hCPR,'(A)') ' FERAT T Correction: ON'
      ELSE
        WRITE (hCPR,'(A)') ' FERAT T Correction: OFF'
      ENDIF

      CALL LOAD_FF(ND, NDDIM, MODHEAD, JOBNUM, bFFASSET,
     >             FF_INFO, IFF_DK, IFF_N_MS, IFF_MAX_MS)
      
      
C***  LOOP OVER ALL DEPTH POINTS  ---------------------------------------
C***    NOTE: OUTWARD DIRECTION NECESSARY FOR ALTERNATIVE, FLUX-CONSTANT 
C***          FORMULATION OF RADIATIVE EQUILIBRIUM IN TEMPEQ
      depthpoints: DO L=1, ND
       
      WRITE (NAME, '(A2,I3,A3)' ) 'DM', L, ' '

C***  NEWTON - RAPHSON - ITERATION  ********************************************
C***  THE ITERATION STARTS FROM THE OLD POP. NUMBERS
      DO J=1,N
        EN(J)=POP1(L,J)
        EN1(J)=EN(J)
      ENDDO
      EN(NPLUS1)=RNE(L)
      EN1(NPLUS1)=EN(NPLUS1)
      TL = T(L)
 
C***  ERROR STOP, IF NEGATIVE TEMPERATURES ON THE MODEL FILE
      IF (TL .LE. .0) THEN
        PRINT 2,L,TL
    2   FORMAT (' ERROR DISCOVERED BY LINPOP: NEGATIVE TEMPERATURE: '
     >    'T(',I2,') = ',E10.0)
        CALL REMARK ('TEMPERATURE NEGATIVE')
        STOP 'ERROR'
      ENDIF
ccc   Here one might consider to apply TOLD -- now it is identical to TNEW
ccc   TOLD is parameter of CALL COMA, goes to SETXJC and SETXJL
      TOLD=TL 

C***  IRON: IRON-LINES ARE CALCULATED WITH DIAGONAL OPERATOR
      DO IND=1, LASTIND
         BDIAG(IND) = NOM(INDLOW(IND)) .EQ. KODAT(26)
         IF (GAMMAD .GT. 0.) THEN
ccc   Der Sinn BEIDER nachstehenden Daempfungs-Varienten ist dunkel ... 
ccc     wrh, 18-Sep-2002 17:34:07
cc            WFELOW(L,IND) = WFELOW(L,IND)*EXP(1.-GAMMAD)
cc            WFENUP(L,IND) = WFENUP(L,IND)*EXP(1.-GAMMAD)
            WFELOW(L,IND) = WFELOW(L,IND) / (1. + ALOG10(GAMMAD))
            WFENUP(L,IND) = WFENUP(L,IND) / (1. + ALOG10(GAMMAD))
         ELSE
            WFELOW(L,IND) = 0.
            WFENUP(L,IND) = 0.
         ENDIF
      ENDDO

C***  DETERMINE THE CORE-CONFINING FREQUENCIES XRED, XBLUE OF THE LINES
C***  AT THE CURRENT DEPTH POINT
C***  Note for clumping: Line opacities calculated with average opacity. 
C***  This only holds if DENSCON = 1. / FILLFAC !!!
      IF (GAMMAL .GT. .0 .OR. GAMMAR .GT .0 .OR. L .EQ. 1)
     >      CALL LCORE (XRED0,XBLUE0,GAMMAL,LASTIND,INDLOW,INDNUP,
     >      GAMMAR,IPRILC,MODHEAD,JOBNUM,OPALOLD,
     >      L,ND,VELO,GRADI,RADIUS,SLOLD,XMAX,ERXMIN,
     >      VDOP,RSTAR,ENTOT,EN,NF,XLAMBDA,ELEVEL,NOM,
     >      NDIM,N,NCHARG,WEIGHT,EINST,LINE,NSCHAR,BDIAG,GAMMAD)
 
C***  IN CASE OF NON-ZERO LINE CORES: 
C***  CALCULATE BACKGROUND CONT. SOURCE FUNCTION BY INTERPOLATION
      DO IND=1,LASTIND
        IF (XRED0(IND) >= XBLUE0(IND)) EXIT
        WAVENUM=ELEVEL(INDNUP(IND))-ELEVEL(INDLOW(IND))
        CALL XRUDI (SCOLIND(IND),WAVENUM,SCOLD(1,L),XLAMBDA,1,NF,1)
      ENDDO

C***  Location of old TDIFFUS block (now before LINPOP call)
 
C***  BEGINNING OF THE ITERATION LOOP  *********************************
      ITNE(L)=0
      LASTIT = 'NONE'
      LASTIT2 = 'NONE'
      LASTIT3 = 'NONE'
      BRIMP = .TRUE.
      BRIMP = .FALSE. !auskommentiert bei goetz
      FN = 0.
      FTN = 0.
      FNOLD = 0.
      FTNOLD = 0.
C***  Define Matrix columns which will be replaced by the number conservation 
      DO NA=1,NATOM
         NFIRNA=NFIRST(NA)
         NLANA=NLAST(NA)
         IMAXPOP(NA) = ISMAX(NLANA-NFIRNA+1,EN(NFIRNA),1)+NFIRNA-1
      ENDDO
      
      !- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
   10 ITNE(L)=ITNE(L)+1     !Schleifen-Sprungmarke

C***  PREPARE BROYDEN SWITCH
C!!!      IF (ITNE(L) .EQ. ITBR) BROYDEN = .TRUE.
C!!!  THIS VERSION: CANCELLING HOLDS ONLY FOR THE PRESENT ITERATION
      IF (ITNE(L) .GE. ITBR) BROYDEN = .TRUE.
      IF (ITNE(L) .LT. ITBR) BROYDEN = .FALSE.

C***  PREPARE BROYDEN RESET
      IF (ITNE(L) .LE. 25) THEN
        BRRES2 = BRRESET
      ELSE
        BRRES2 = 1.
      ENDIF

C***  ENTRY, IF BROYDEN STEP HAS BEEN CANCELLED
107   CONTINUE 

C***  RESTORE THE ORIGINAL CORE-CONFINING FREQUENCIES FOR THE PRESENT
C***  NEWTON-RAPHSON ITERATION
C     goetz hat diese Schleife auskommentiert
      DO IND=1,LASTIND
        XRED(IND)=XRED0(IND)
        XBLUE(IND)=XBLUE0(IND)
      ENDDO

C***  REPLACE ZERO OR NEGATIVE POP. NUMBERS BY 
C***  MAX OF 1.E-100 AND 0.5 TIMES THEIR ABSOLUTE VALUE
      DO J=1, N
        IF ( EN(J) .LE. 0.0 ) EN(J) = MAX(1.E-100, -0.5*EN(J))
      ENDDO

C***  PRE-CALCULATE EXPONENTIAL FACTORS FOR THE TEMPERATURE OF THE
C***  CURRENT DEPTH POINT
C***  THIS MUST BE REPEATED, WHEN THE TEMPERATURE HAS BEEN UPDATED
      IF (ITNE(L) .EQ. 1) THEN
         DO K=1,NF
            WAVENUM=1.E8/XLAMBDA(K)
            EXPFAC(K)=EXP(-C1*WAVENUM/TL)
         ENDDO
      ENDIF
 
C***  CALCULATE LTE POP. NUMBERS
      ENE=EN(NPLUS1)*ENTOTDENS(L)

C***  Note: ENTOTL is automatically enhanced in COMA via ENE

      CALL       LTEPOP (N,ENLTE,TL,ENE,WEIGHT,NCHARG,EION,ELEVEL,NOM,
     >                  ABXYZ,NFIRST,NLAST,NATOM)

C***  PRECALCULATE FREE-FREE CROSS SECTIONS (FOR ALL FREQUENCY POINTS
C***  AND ALL ION CHARGES) AT CURRENT DEPTH POINT
C***     ----- OPTIMIZATION OF SUBR. DCOOP (CALLED FROM COMA) -----
      TLOG=ALOG10(TL)
      ROOTTL=SQRT(TL)
      DO K=1,NF
        XLAM=XLAMBDA(K)
        W=1.E8/XLAM
        W3=W*W*W
        PRESIG=CFF/W3/ROOTTL
        XLAMLOG=ALOG10(XLAM)
        SIGMAFF(K,0)=0.0
        DO ION=1,MAXION
          CALL GFFLOG (GIII,ION,XLAMLOG,TLOG)
          SIGMAFF(K,ION)=PRESIG*FLOAT(ION*ION)*GIII
        ENDDO
      ENDDO


C***  SETUP COEFFICIENT MATRICES
      CALL       COMA (CRATE,RRATE,RATCO,DM,N,NRANK ,NDIM,V1,ABXYZ,
     > ENLTE,TL,ENE,NCHARG,ELEVEL,EINST,EION,WEIGHT,ALTESUM, XLAMBDA,
     > FWEIGHT,XJC,NF,L,XJL,ND,XJLAPP,SLOLD,LASTIND,INDLOW,
     >      INDNUP,NOM,NATOM,KODAT,NFIRST,NLAST,PHI,PWEIGHT,DELTAX,XMAX,
     >      NFL,OPAC,SCNEW,DOPA,DETA,OPAL,SLNEW,DOPAL,DETAL,SIGMAKI,
     >      ETAC,NFEDGE,EXPFAC,SCOLIND,SCNEIND,OPACIND,SIGMAFF,MAXION,
     >      NOTEMP,TOLD,KONTLOW,KONTNUP,LASTKON,RUDLINE,IONGRND,
     >      XRED,XBLUE,WCHARM,EN,RSTAR,SCOLD,XJCAPP,VDOP,COCO, KEYCBB,
     >      NRB_CONT, ZERO_RATES, POPMIN,
     >      IONAUTO,NAUTO,MAXAUTO,LOWAUTO,WAUTO,EAUTO,AAUTO,DRRATEN,
     >      RDIEL,RAUTO,DRJLW,DRJLWE,DRLJW,IBLENDS,MAXLAP,XLAMZERO,
     >      KODRNUP,KODRLOW,LASTKDR,KEYCBF, OPALOLD,
     >      BETA,PHIL,NBLENDS,BROYDEN,ATEST,BTEST,MAXATOM,
     >      SIGMATHK,SEXPOK,EDGEK,XDATA,RADIUS(L),
     >      XJCLP1, OPAC1, RADIUS, ITNE(L), TEFF, OPATHOM, 
     >      LAINDHI, KRUDAUT, LEVEL, 
     >   WFELOW, WFENUP, EN1, BDIAG, 
     >   FERATLU, FERATUL, LASTFE, FERATLU0, FERATUL0, 
C*** Quantities for fine-frequency grid 
     >   SFINE_OLD, SFINE_NEW, MAXFINE, KONTHLP, MAXIND, XKC, XKC2,
     >   SIGMA1I, NLINE, LINEINDEX, XLAMSOR, XLAMMIN, XLAMMAX,
     >   NUPACT, LOWACT, BLASERL, ETAL, LIND, LINDS,
C***  IRON
     >   INDRB, IFRBSTA, IFRBEND, 
     >   VDOPFE, DXFE, XLAM0FE,
     >   INDFEACT, MAXFEACT, BFECHECK, BFEWING,
     >   DFEINDR, SIGMAFE, OPAFE, ETAFE, IFENUP, IFELOW,
     >   INDEXMAX, BFEMODEL, BNUEFE, 
C***
     >   XLAMBDA2, NF2, NDDIM, VELO(1), XKMIN, XKMAX, XKMID, XKRED, 
     >   ALPHA, SEXPO, ADDCON1, ADDCON2, ADDCON3, IGAUNT, 
     >   XKRED_CORE, XKBLUE_CORE, 
     >   BXJLAPPNEW, BXJCAPPNEW, BNEWOPER, BXJLAPPCORE,
     >   XJL_PLOTDATA, XJC_PLOTDATA_I, XJC_PLOTDATA_L, 
     >   IPLOT_XJLAPP, IPLOT_XJCAPP, LPLOT_XJCAPP, NITER_PLOT_JAPP, 
     >   BPLOTAPP, PWEIGHTCL, WS, FWTEST, 
     >   IWARN_NEG_XJCAPP, IWARN_NEG_XJLAPP, 
     >   XJCAPPNEW, XJLAPPNEW, 
     >   GAMMAC, GAMMAL, GAMMAR, 
     >   XLAM_FINE_START, XLAM_FINE_END, IMAXPOP, bBLOCKINVERSION,
C***  New Fine-spaced WCHARM handling
     >   IFF_MAX, IFF_MAX_MS, FF_INFO, IFF_DK, IFF_WCHARM, WCHARM_FINE, 
     >   IFF_N_MS, bFFASSET,
C***  FERAT correction     
     >   CORRS, DEXFAC, bFeTCORR)


C***  Update counters for ZERO_RATES
      IF (ITNE(L) .EQ. 1) THEN
         DO J=1, N
            IF (ZERO_RATES(J,L)) THEN
               NZERORATES(J) = NZERORATES(J) + 1
               LMINZERORATES(J) = MIN0 (LMINZERORATES(J),L) 
               LMAXZERORATES(J) = MAX0 (LMAXZERORATES(J),L)
            ENDIF 
         ENDDO
      ENDIF
 
C***  ALGEBRA OF ONE NEWTON ITERATION STEP:
C***  ------------------------------------

C***  V1 := RIGHT-HAND SIDE VECTOR
      CALL VMF (V2,EN,RATCO,NRANK,NRANK)
      CALL VSUB (V1,V2,NRANK)

C***  FN = THE NORM OF THE RIGHT-HAND SIDE VECTOR MEASURES HOW ACCURATE
C***       THE RESULT OF THE LAST (!) ITERATION FULLFILLS THE EQUATIONS 
      IF (ITNE(L) .EQ. 2 .OR.
     >    ITNE(L) .GT. 2 .AND. FN .LT. FNOLD) FNOLD = FN
      CALL BRNORM2 (FN, V1, NRANK)

      FTNOLD = .0
      FTN    = .0
      BSUM = .FALSE.

C***  ENTRY POINT IF SUM OF EN IS OUT OF VALID RANGE AND LASTIT = BROYDEN
  410 CONTINUE

C***    TEST WHETHER BROYDEN STEP (PRECEEDING ITERATION!) HAS REALLY 
C***    IMPROVED THE SOLUTION 
        IF (BNEWTONRESET) THEN
          IF (LASTIT .EQ. 'BROYDEN'  .AND.  
     >       (FN > BRRES2 * FNOLD  .OR. FTN > BRRES2 * FTNOLD) .OR.
     >        BSUM) THEN
C***    ... BROYDEN STEP WAS BAD: RESORE OLD POPNUMBERS, REPEAT AS NEWTON STEP
            DO J=1, NPLUS1
              EN(J) = ENOLD(J)
            ENDDO
            BROYDEN = .FALSE.
            NBRCANC = NBRCANC + 1
            LASTIT = 'CANCEL'
            LASTIT2 = 'NONE'
            LASTIT3 = 'NONE'
            ITNE(L) = ITNE(L) - 1
C!!!        WRITE(*,*) 'STEP FAILED  LASTIT=',LASTIT
            GOTO 107
          ENDIF
        ENDIF

        LASTIT3 = LASTIT2
        LASTIT2 = LASTIT

C***  ALTERNATIVE BRANCHES: BROYDEN OR NEWTON!  *****************************


C*************************
        IF (BROYDEN) THEN
C*************************

          LASTIT = 'BROYDEN'
C***      AT THE FIRST ITERATION, DM IS READ FROM DMFILE
          IF (ITNE(L) .EQ. 1) THEN
            IF (bUseDMFILE) THEN
              CALL READMS (hDMFILE, DM, NRANK*NRANK, NAME, IERR)
            ENDIF
            IF (BRIMP) THEN
              itcode(itne(l)) = 'B'
            ELSE
              itcode(itne(l)) = 'b'
            ENDIF

C***    SUBSEQUENT ITERATIONS: DM IS IMPROVED BY THE UPDATE FORMULA  
          ELSE
            CALL VSUB (BRY, V1, NRANK)
            IF ((LASTIT3 .EQ. 'BROYDEN' .OR. LASTIT3 .EQ. 'NEWTON') 
     >               .AND. bUseTWOPNT) THEN
              TWOPNT = .TRUE.
              CALL BRTPUP (DM, BRS, BRSP, BRY, BTEST, DTEST, 
     >                    NRANK, NRANK, TWOPNT, BRDIVFAILED)
            ELSE 
              TWOPNT = .FALSE.
              CALL BRTPUP (DM, BRS, BRS, BRY, BTEST, DTEST, 
     >                    NRANK, NRANK, TWOPNT, BRDIVFAILED)
            ENDIF
            IF (TWOPNT) THEN
              IF (BRIMP) THEN
                itcode(itne(l)) = 'T'
              ELSE
                itcode(itne(l)) = 't'
              ENDIF
            ELSE
              IF (BRIMP) THEN
                itcode(itne(l)) = 'B'
              ELSE
                itcode(itne(l)) = 'b'
              ENDIF
            ENDIF
          ENDIF

C***  THE CORRECTION VECTOR IS CALCULATED:
          CALL VMF (V2, V1, DM, NRANK, NRANK)

C***************
        ELSE
C***************

C***      NEWTON-RAPHSON BRANCH
          LASTIT = 'NEWTON'
          itcode(itne(l)) = 'N'
C***      SOLVE THE LINEAR SYSTEM: >>  V2 * M = V1  << FOR V2

          IF (bBLOCKINVERSION) THEN
C***        matrix M (alias DM) is inverted block by block
C           electron line will be ignored
            CALL LINSOL_SPLIT (V2, DM, V1, NRANK, NRANK,
     >                   NFIRST, NLAST, NATOM, AOFF,
     >                   SCRATCH, ATEST, BTEST, DTEST, VERSION, ELEMENT)
          ELSE
C           matrix M (alias DM) is inverted as a whole
            CALL LINSOL (V2, DM, V1, NRANK, NRANK, 
     >                   SCRATCH, ATEST, BTEST, DTEST, VERSION)
          ENDIF

          IF (VERSION .EQ. 'SING') THEN
            ITCODE(ITNE(L)) = 'I'
            KONVER = .FALSE.
            GOTO 420
          ENDIF        

C*****************
        ENDIF
C*****************
c      write (0,*) 'itne,itcode=',itne(l),itcode(itne(l))

C!!!      WRITE(*,*) 'NORMAL STEP   LASTIT=',LASTIT

C**     STORE CORRECTIONS IN EN --> V2, AND CORRECTIONS IN F(EN) --> V1,
C**       FOR THE (POSSIBLY) FOLLOWING BROYDEN ITERATION
        DO J=1,NRANK 
          BRSP(J) = BRS(J)
          BRS(J) = V2(J)
          BRY(J) = V1(J)
        ENDDO

        BRIMP2 = BRIMP
        IF ((LASTIT2 == 'NEWTON') .OR. (ITNE(L) == 1) .OR.
     >      ((FN <= FNOLD) .AND. (FTN <= FTNOLD)) ) THEN
          BRIMP = .TRUE.
          DO J=1,NPLUS1
            ENOLD(J)=EN(J)
          ENDDO
        ELSE
          BRIMP=.FALSE.
C***      goetz setzt hier BRIMP = .TRUE.
        ENDIF

C***    ADD THE CORRECTION VECTOR V2:
        DO J=1,NRANK
          EN(J) = EN(J) + V2(J)
          !DEBUG TEST:
          IF (EN(J) < POPMIN) THEN
C            EN(J) = POPMIN
            iPopMinViolations = iPopMinViolations + 1
          ENDIF
        ENDDO
 
C***    RE-ADJUST THE ELECTRON DENSITY (new suggestion from wrh: 11. Apr 2011)
        !moved here at 10.10.2011
        IF (bBLOCKINVERSION .AND. LASTIT == 'NEWTON') THEN
          ENNEW=0.0
          DO j=1,N
            ENNEW = ENNEW + EN(j) * NCHARG(j)
          ENDDO
          EN(NPLUS1)=ENNEW
        ENDIF


        SUM = 0.
        ENMAX = 0.
        ENMIN = 1000.
        DO I= 1, N
          SUM = SUM + EN(I)
C***    die folgenden zwei Zeilen sind bedeutungslos, weil nicht genutzt
C           IF (EN(I) .GT. ENMAX) ENMAX = EN(I)
C           IF (EN(I) .LT. ENMIN) ENMIN = EN(I)
        ENDDO
C***  disable this stuff for testing
        if (.true.) then
c      write (*,*) 'SUM: ', SUM
        IF (SUM .GT. 10. .OR. SUM .LT. 0.1 .OR. BRDIVFAILED) THEN
          KONVER = .FALSE.
          IF (LASTIT .EQ. 'NEWTON' .OR. .NOT. BNEWTONRESET) THEN
            IF (.NOT. BRDIVFAILED) THEN
              ITCODE(ITNE(L)) = 'S'
C***  TEST PRINTOUT: Correction Vector
C            write (*,*) 'V2: '
C            DO nTest=1, NRANK
C               write (*,*) ' ', V2(nTest)
C            ENDDO
            ELSE
              ITCODE(ITNE(L)) = 'D'
            ENDIF
            GOTO 420
          ELSE IF (LASTIT .EQ. 'BROYDEN') THEN
            BSUM = .TRUE.
            ITNE(L) = ITNE(L) + 1
            GOTO 410
          ENDIF
        ENDIF
        endif

C***    CONVERGENCE CHECK OF NEWTON ITERATION
C***    (TEMPERATURE CORRECTIONS ARE NOT TAKEN INTO ACCOUNT)
        IF (BRIMP2) THEN
          KONVER = .TRUE.
          DO J=1,NPLUS1
C          DO J=1,NRANK
            KONVER=KONVER .AND.
     >       ( ABS(V2(J)) .LT. EPSDN*ABS(EN(J)) .OR.
     >         ABS(EN(J)) .LT. SMALLPOP)
          ENDDO
        ELSE
          KONVER = .FALSE.
        ENDIF
 
        IF ( (.NOT. KONVER) .AND. (ITNE(L) < ITMAX) ) GOTO 10
      !- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 

C***  ENTRY POINT IF SUM OF EN IS OUT OF VALID RANGE AND LASTIT = NEWTON
C***  or if Matrix-Inversion (OWNINV) has failed (Singularity discovered)
  420 CONTINUE

      IF (.NOT. KONVER) THEN
        DO LL=ITNE(L)+1, ITMAX
          ITCODE(LL:LL) = ' '
        ENDDO
c        ITNE(L)=-ITMAX
        ITWARN=ITWARN+1
        NOCON = .TRUE.
      ENDIF

C***  CHANGE AUTO_MODIFY VARIABLES
      IF (BAUTO .AND. .NOT. KONVER) THEN
        NKONVER = NKONVER + 1
        CKONVER(L) = 'D'
      ENDIF

C***  END OF NEWTON-ITERATION LOOP   *******************************************

C!!!  test printout
      IF (.true.) THEN
        itabs = iabs(itne(l))
        WRITE (0,77) l, itabs, (itcode(i), i=1, itabs)
   77   FORMAT (5x, 'L=', i3, 3x, 'Niter=',i4, 5x, 1000A1)
        !write (0,*) 'LINPOP> PRIMAT of DM after inversion'
        !CALL PRIMAT(DM, NRANK, NRANK, 'Output of Matrix DM (inv)')
        !write (0,*) 'LINPOP> end of PRIMAT'
      ENDIF

C***  UPDATING ARRAY "POPNUM" AND CALCULATING DEPARTURE COEFFICIENTS
      DO J=1,N
        IF (ZERO_RATES(J,L)) THEN
           POPNUM(L,J) = POPMIN
        ELSE
           POPNUM(L,J)=EN(J)
        ENDIF
      ENDDO
C***  UPDATING THE ELECTRON DENSITY
      RNE(L)=EN(NPLUS1)
C***  Departure Coefficients are not calculated, when depths point was not 
C***    converged. DEPART is for Output only!
      DO J=1,N
        IF (KONVER) THEN
          DEPART(L,J)=EN(J)/ENLTE(J)
        ELSE
          DEPART(L,J)=1.
        ENDIF
        POPLTE(L,J)=ENLTE(J)
      ENDDO

C***  UPDATING THE DERIVATIVE MATRIX ON DMFILE (FOR BROYDEN METHOD)
C***  (NPLUS2 = MAX. RANK TO AVOID MULTIPLE INCREASE OF DMFILE)
C***  08.06.2010> NPLUS2 -> NRANK due to split
      IF (bUseDMFILE) THEN
        CALL WRITMS (hDMFILE, DM, NRANK*NRANK, NAME, -1, IDUMMY, IERR)
C      CALL WRITMS (16, DM, NPLUS2*NPLUS2, NAME, -1, IDUMMY, IERR)
      ENDIF

C***********************************************************************
C***  test output: "image" of specified matrix (NOTE:  ndim <= 126 )
      ldepth = 0
      if (l .eq. ldepth) then
        call primat (rrate,n,ndim,'RRATE')
        call primat (crate,n,ndim,'CRATE')
        call primat (ratco,nrank,nrank,'RATCO')
      endif
C***********************************************************************

C***  PRINTOUT OF RATE COEFFICIENTS ETC.  ------------------------------
      IF (LSRAT.NE.-1) THEN
        IF ((L.GE.IFRRA.AND.L.LE.ITORA).OR.ITORA.EQ.0) THEN
          NETTO=1
c!!!      NETTO=2
          LM1=L-1
          IF (IFRRA.GT.0) LM1=L-IFRRA
          IF  (((LM1)/LSRAT)*LSRAT.EQ.(LM1).OR.L.EQ.ND) THEN
            CALL PRIRAT (ITNE(L),N,LEVEL,NDIM,L,CRATE,RRATE,RATCO,EN,
     $           IFRRA,MODHEAD,JOBNUM,NETTO,NFIRST,NLAST,NATOM,NATOUT,
     $           NAUTO,RDIEL,RAUTO,IONGRND,KODRLOW,LASTKDR)
          ENDIF
        ENDIF
      ENDIF
 
C***  STORE APPROXIMATE CONTINUUM RADIATION FIELD AND OPACITIES 
C***  OF LAST-TREATED DEPTH POINT  FOR USE IN THE 
C***  ALTERNATIVE, FLUX-CONSERVING TEMPERATURE EQUATION
      DO K=1, NF 
        XJCLP1(K) = XJCAPP(K)
        OPAC1(K) = OPAC(K) + OPATHOM
      ENDDO

      ENDDO depthpoints
C***  ENDLOOP (over all depth points) ----------------------------------------------

C***  TEST OUTPUT OF NEW ELECTRON DENSITY
      !WRITE (*,*) "STEAL->LINPOP: RNE (rel. electron density)"
      !DO j=1,L
      !  WRITE (*,*) j, RNE(j)
      !ENDDO
 
c      NKONV_THRESH = 30   !NKONV_THRESH is now a CARDS option

C***  Abort execution if 
C***  Number of diverged points is greater threshold

      IF (NKONVER .GT. NKONV_THRESH) THEN
        WRITE (0,*) 'Execution of STEALCL is Stopped'
        WRITE (0,*) 'NKONVER, NKONV_THRESH=', NKONVER, NKONV_THRESH
        STOP 'ABORT in Subr. STEALCL'
      ENDIF
       

C***  CREATE THE AUTO_MODIFY FILE
      IF (BAUTO .AND. 
     >    NKONVER .GT. 0 .AND. 
     >    NKONVER .LE. NKONV_THRESH) THEN
        IF ((BUNLU .AND. (iAMT /= -1)) .OR. iAMT == 1) THEN
          CHR2 = 'TEMPERATURE CORRECTIONS'
        ELSE
          CHR2 = 'NO TEMPERATURE CORRECTIONS'
        ENDIF
        OPEN (hMODIFY, FILE='MODIFY_INPUT', STATUS='UNKNOWN')
        WRITE (hMODIFY, '(A)') 'AUTO MODIFY'
        WRITE (hMODIFY, '(A)') CHR2(:IDX(CHR2))
        WRITE (*, '(1X,A)') 'LINPOP> Next Job: MODIFY'
        WRITE (*, '(1X,A)') 'LINPOP> ' // CHR2(:IDX(CHR2))
C***  LOOP OVER ALL DEPTH POINTS
        L = 0
        dploop: DO !- - - - - - - - - - - - - - - - - - - - - - -
          L = L + 1
          IF (CKONVER(L) .EQ. 'D') THEN
            MINEX = L
  122       L = L + 1
            IF (L .LE. ND) THEN
              IF (CKONVER(L) .EQ. 'D') GOTO 122
            ENDIF
            MAXEX = L - 1
C***  DO THE INTERPOLATION
            IF (MAXEX .NE. ND .AND. MINEX .NE. 1) THEN
              WRITE (CHR,'(I3)') MINEX - 1
              CHR1 = 'INTERPOLATE FROM POINT ' // CHR
              WRITE (CHR,'(I3)') MAXEX + 1
              CHR1 = CHR1(:IDX(CHR1)) // '  TO POINT  ' // CHR
            ELSE 
              IF (MAXEX .EQ. ND) THEN
                CHR1 = 'INNER'
                WRITE (CHR,'(I3)') MINEX - 1
              ELSE IF (MINEX .EQ. 1) THEN
                CHR1 = 'OUTER'
                WRITE (CHR,'(I3)') MAXEX + 1
              ENDIF
              CHR1 = CHR1(:IDX(CHR1)) //' EXTRAPOLATION  FROM POINT '// 
     >               CHR(:IDX(CHR)) // '  SECOND POINT  ' // CHR
            ENDIF
            WRITE (hMODIFY, '(A)') CHR1(:IDX(CHR1))
            WRITE (*, '(1X,A)') 'LINPOP> ' // CHR1(:IDX(CHR1))
          ENDIF          
        
          IF (L < ND) THEN
            CYCLE dploop
          ELSE
            EXIT dploop
          ENDIF
        ENDDO dploop !- - - - - - - - - - - - - - - - - - - - - - -

        CLOSE (hMODIFY)
      ENDIF

C***  CLOSE DMFILE and FFASSET
      IF (bUseDMFILE) CALL CLOSMS (hDMFILE, IERR)
      IF (bFFASSET) CALL CLOSMS(hFF, IERR)
      
      WRITE (hCPR,*) 'Number of POPMIN VIOLATIONS: ', iPopMinViolations
      
C***  Printout of ZERO_RATES (Option: PRINT ZERORATES) 
      IF (IPRINTZERORATES .GT. 0) THEN
        WRITE (hCPR,'(/,A)') 'Levels for which the equations have been removed' 
        WRITE (hCPR,'(  A)') '------------------------------------------------' 
        DO J=1,N
          IF (NZERORATES(J) .EQ. 0) CYCLE
          IF (NZERORATES(J) .EQ. ND) THEN
            STRING15 = ' = all depths!'
          ELSE 
            STRING15 = ''
          ENDIF
          IF (NZERORATES(J) .EQ. 
     >      LMAXZERORATES(J) + 1 - LMINZERORATES(J)) THEN
            STRING3 = 'all'
          ELSE
            WRITE (STRING3, '(I3)') NZERORATES(J)
          ENDIF
          WRITE (hCPR, '(A,I4,A,I3,A,I3,A )') 
     >      'Level', J,': ' // LEVEL(J) // ' at ' // STRING3 //
     >      ' depth points between L=', LMINZERORATES(J), 
     >      ' and', LMAXZERORATES(J), STRING15
        ENDDO
      ENDIF

      RETURN
      END
