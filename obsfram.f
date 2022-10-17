      SUBROUTINE OBSFRAM (LTOT,CORE,XMAX,XMAXLIN,EMINT,CEMINT,
     >                   BCOREL,DBDRL,TAUMAX,PJPJ, 
     >                   ZFINE, OPAFINE,OPAFC,OPALFIN,
     >                   ETAFINE,ETAFC,ETALFIN,SFINE,CSFINE,RRAY,
     >                   OPARAY,OPALRAY,
     >                   ETARAY,ETACRAY,ETALRAY, ZRAY,XCMF,
     >                   MAXXN,NDADDIM,NDDIM,DELXLAP,NBLINE,IVERSION, 
     >                   LINPRO,AVOIGT,TAU,TAUC,DTAU,DTAUC,WTAU,WTAUC,
     >                   XCMFFINE, POROLENGTHFINE, MAXLAP, DXMAX,
     >                   BIRONLINES, OPAFE, ETAFE, NFLDIM, 
     >                   XCMFBLUE, XCMFRED, DXCMF, POROLENGTHRAY, 
     >                PHITAB, NFDIMPHITAB, NLDIMPHITAB, IPOINTERPHITAB,
     >                RADIUS, ND, DD_VDOPDU_RAY, NATOM, NBFIRST, NBLAST,
     >                IND_ELLINE, DD_VDOPDU_FINE_NORMFAC,
     >                DD_VDOPDU_FINE_SQRD, DD_VDOPDU_FINE,
     >                GRIEMPAR, KODAT, VDOP, INDCUT, 
     >                ZINTER, NMOD, MAXMOD, KINDEX)

C***********************************************************************
C***  INTEGRATION OF THE EMERGENT INTENSITY IN THE OBSERVER'S FRAME
C***********************************************************************
  
      INTEGER, INTENT(IN) :: NATOM
      REAL, DIMENSION(NDADDIM, NATOM), INTENT(IN) :: DD_VDOPDU_RAY
      REAL, DIMENSION(MAXXN, NATOM) :: DD_VDOPDU_FINE, DD_VDOPDU_FINE_SQRD, 
     >                                 DD_VDOPDU_FINE_NORMFAC
      INTEGER, DIMENSION (MAXLAP) :: IND_ELLINE
      DIMENSION XCMF(NDADDIM)
      DIMENSION RRAY(NDADDIM), OPARAY(NDADDIM)
      DIMENSION ETARAY(NDADDIM), ETACRAY(NDADDIM)
      DIMENSION ZRAY(NDADDIM)
      REAL, DIMENSION(NDADDIM,NBLINE) :: OPALRAY, ETALRAY
      DIMENSION DELXLAP(NBLINE)
      
      DIMENSION ZFINE(MAXXN), OPAFINE(MAXXN)
      DIMENSION OPAFC(MAXXN), OPALFIN(MAXXN)
      DIMENSION ETAFINE(MAXXN), ETAFC(MAXXN)
      DIMENSION ETALFIN(MAXXN), SFINE(MAXXN)
      DIMENSION CSFINE(MAXXN)

      LOGICAL CORE, BIRONLINES
      CHARACTER LINPRO(MAXLAP)*8
      REAL, DIMENSION(MAXLAP) :: XMAXLIN

C***  WPI = SQRT(PI)
      DATA WPI /1.772454/

C***  INITIALIZATION
      EMINT = .0
      CEMINT = .0
      TAUSUM = .0
      TAUSUMC = .0

      LRED = 1
      LBLUE = LTOT

      CALL ZONEINT (LRED, LBLUE, EMINT, CEMINT, TAUSUM, TAUSUMC,
     >                PJPJ, ZFINE, KINDEX,
     >                OPAFINE, OPAFC, OPALFIN, ETAFINE, ETAFC,
     >                ETALFIN, SFINE, CSFINE,
     >                RRAY, ZRAY, 
     >                XCMF, OPARAY,
     >                OPALRAY, ETARAY,
     >                ETACRAY, ETALRAY,
     >                MAXXN, LTOT, DELXLAP, NBLINE, NBFIRST, NBLAST, 
     >                NDADDIM, IVERSION, LINPRO, AVOIGT,
     >                TAU, TAUC, DTAU, DTAUC, WTAU, WTAUC,
     >                XCMFFINE, POROLENGTHFINE, TAUMAX, XMAX, XMAXLIN,
     >                DXMAX, BIRONLINES, OPAFE, ETAFE, NDDIM, NFLDIM, 
     >                XCMFBLUE, XCMFRED, DXCMF, CORE, POROLENGTHRAY, 
     >                PHITAB, NFDIMPHITAB, NLDIMPHITAB, IPOINTERPHITAB,
     >                RADIUS, ND, MAXLAP, DD_VDOPDU_RAY, NATOM, 
     >                IND_ELLINE, DD_VDOPDU_FINE_NORMFAC,
     >                DD_VDOPDU_FINE_SQRD, DD_VDOPDU_FINE, 
     >                GRIEMPAR, KODAT, VDOP, ZINTER, NMOD, MAXMOD)
 
C***  FOR CORE RAYS, ADD INCIDENT RADIATION
      IF (CORE) THEN
         X = OPARAY(LTOT)  ! X means kappa

C***     Without lines: PLUSIC
         PLUSIC = BCOREL + DBDRL * ZRAY(LTOT) / X

C***     With lines: add line opacities
         DO 45 NLOC = NBFIRST, NBLAST
            XLBLTOT=XCMF(LTOT)-DELXLAP(NLOC)
            IF (XLBLTOT .GT. 300) THEN
               PHI = .0
            ELSE
               PHI=EXP(-XLBLTOT*XLBLTOT)/WPI
            ENDIF            
            X=X+PHI*OPALRAY(LTOT,NLOC)
   45    CONTINUE
         PLUSI  = BCOREL + DBDRL * ZRAY(LTOT) / X
         EMINT  = EMINT  + PLUSI  * EXP(-TAUSUM)
         CEMINT = CEMINT + PLUSIC * EXP(-TAUSUMC)
      ENDIF
 
      RETURN
      END
