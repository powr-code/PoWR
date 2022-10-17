      SUBROUTINE PREPRAY (Z, COSPHI, P, ND, NDDIM, NP, JP, XO, LTOT, 
     >                   PWEIGHT, CORE, VDU, R,
     $                   OPAL, ETAL, RRAY, OPARAY, ETARAY,
     $                   ETACRAY, OPALRAY, ETALRAY, NFLDIM, ZRAY,
     $                   XCMF, NDADDIM, NBLINE, MAXLAP, 
     $                   REDIS, ETACK, ETACCK, OPACK, 
     $                   XCMFRED, XCMFBLUE, DXCMF,
     $                   BCORE,BCOREL,DBDR,DBDRL, K,
     >                   POROLENGTH, POROLENGTHRAY,
     >                   RCOROT, VSIDU, 
     >                   DD_VDOPDU_RAY, DD_VDOPDU, NATOM, NBFIRST,
     >                   NBLAST, MAXATOM, 
     >                   BVSINI_AT_RCOROT, NMOD, MAXMOD, 
     >                   ZINTER, XMAX, DELXLAP)
C***********************************************************************
C***  DEFINING WHOLE RAYS INCLUDING THE BACKWARD HALF-SPHERE
C***  All names of the rat vectors end on ...RAY
C***   Note: XCMF should also be called XCMFRAY! 
C***  THE CONTINUUM EMISSIVITY IS CALCULATED TWICE: 
C***  ETARAY:  EMISSIVITY WITH REDISTRIBUTION FROM ETACK
C***  ETACRAY: EMISSIVITY OF THE CONTINUUM FROM ETACCK
C***  If a second model is active (i.e. NMOD=2), the quantities that are
C***  copied to the RAY vectors are taken from IMOD=2 if Z falls in the
C***  interval specified by ZINTER   
C***********************************************************************

      INTEGER, INTENT(IN) :: NATOM      
      REAL, DIMENSION (NDDIM, MAXATOM, MAXMOD) :: DD_VDOPDU
      REAL, DIMENSION (NDADDIM, NATOM) :: DD_VDOPDU_RAY

      REAL, DIMENSION (NDADDIM) :: ZRAY, XCMF, RRAY, 
     >                             OPARAY, ETARAY, ETACRAY
      REAL, DIMENSION (NDDIM,MAXLAP,MAXMOD) :: OPAL, ETAL
      REAL, DIMENSION (NDADDIM,MAXLAP) :: OPALRAY, ETALRAY
      DIMENSION DELXLAP(NBLINE)

      DIMENSION VDU(NDDIM,MAXMOD), R(ND)
      DIMENSION Z(ND),P(NP), ZINTER(2)
      REAL, DIMENSION (NDDIM,NFLDIM,MAXMOD) :: ETACK, ETACCK, OPACK
      REAL, DIMENSION(NFLDIM,MAXMOD) :: BCORE, DBDR
      DIMENSION POROLENGTH(ND), POROLENGTHRAY(NDADDIM)
      LOGICAL CORE, REDIS
C***  BCOROT: local boolean, .true. if in co-rotating domain
      LOGICAL :: BCOROT, BVSINI_AT_RCOROT

C***  in case VSINI refers to RSTAR: compensate for the scaling with 1/RCOROT
      IF (BVSINI_AT_RCOROT) THEN
        VROTPROJ1 = - P(JP) * COSPHI * VSIDU   
      ELSE 
        VROTPROJ1 = - P(JP) * COSPHI * VSIDU * RCOROT
      ENDIF

C***  Check if ray intersects core 
      LMAX = MIN0(NP+1-JP,ND)
      CORE=LMAX .EQ. ND    

C***  LTOT = number of grid points along ray including back hemisphere
      IF (CORE) THEN
        LTOT = LMAX
      ELSE
        LTOT = 2*LMAX - 1
      ENDIF

C***  LL runs over both hemispheres, 
C***     while L is only within the range for the front hemisphere
      DO 1 LL = 1, LTOT
        IF (LL .LE. LMAX) THEN
           L = LL
           ZRAY(LL) = Z(L)
        ELSE
           L = 2 * LMAX - LL
           ZRAY(LL) = -Z(L)
        ENDIF

C***    Check if current point lies in second-model domain
        IMOD=1
        IF (NMOD .EQ. 2) THEN
           IF ( (ZRAY(LL)-ZINTER(1))*(ZRAY(LL)-ZINTER(2)) .LT. .0) THEN
              IMOD=2
           ENDIF
        ENDIF

        RRECIP = 1. / R(L)
        VRAD = VDU(L,IMOD)
        VRADPROJ = VRAD*ZRAY(LL)*RRECIP
        BCOROT = (R(L) .LE. RCOROT)
        IF (BCOROT) THEN
          VROTPROJ = VROTPROJ1 / RCOROT
        ELSE
          VROTPROJ = VROTPROJ1 * RCOROT * RRECIP * RRECIP
        ENDIF

        XCMF(LL) = XO - (VRADPROJ + VROTPROJ)

        RRAY(LL) = R(L)
        POROLENGTHRAY(LL) = POROLENGTH(L) 
        IF (XCMF(LL) .GT. XCMFBLUE) THEN
          XC = XCMFBLUE
        ELSE IF (XCMF(L) .LE. XCMFRED) THEN
          XC = XCMFRED
        ELSE
          XC = XCMF(LL)
        ENDIF
        KCMF = INT ( (XCMFBLUE - XC) / DXCMF) + 1
        Q = (XCMFBLUE - XCMF(LL))/DXCMF - (KCMF-1)
        PW= 1.-Q
        ETARAY(LL)  = PW*ETACK (L,KCMF,IMOD) + Q*ETACK (L,KCMF+1,IMOD)
        ETACRAY(LL) = PW*ETACCK(L,KCMF,IMOD) + Q*ETACCK(L,KCMF+1,IMOD)
        OPARAY(LL)  = PW*OPACK (L,KCMF,IMOD) + Q*OPACK (L,KCMF+1,IMOD)

        DO NA=1, NATOM
           DD_VDOPDU_RAY(LL,NA) = DD_VDOPDU(L, NA, IMOD)
        ENDDO
    1 CONTINUE    

C***  FIND THE FIRST AND LAST INVOLVED LINE FOR THAT RAY
C***  due to rotation and/or second-model, XCMF might not be monotonic
      XCMFMIN = MINVAL (XCMF(1:LTOT))
      XCMFMAX = MAXVAL (XCMF(1:LTOT))
      NBFIRST = ISRCHFGT (NBLINE, DELXLAP, 1, XCMFMIN-XMAX)
      NBLAST  = ISRCHFGE (NBLINE, DELXLAP, 1, XCMFMAX+XMAX) - 1

C***  Note: this loop over LL cannot be unified with the above one, 
C***        bacause it needs NBFIRST, NBLAST to be defined before      
      DO LL = 1, LTOT
        IF (LL .LE. LMAX) THEN
           L = LL
        ELSE
           L = 2 * LMAX - LL
        ENDIF
        DO NBL=NBFIRST, NBLAST
            OPALRAY(LL,NBL)=OPAL(L,NBL,IMOD)
            ETALRAY(LL,NBL)=ETAL(L,NBL,IMOD)
        ENDDO
      ENDDO

C***  STORE THE INTERPOLATED BCORE AND DBDR IN LOCAL VARIABLES 
C***  WHEN RAY IS A CORE RAY
C***  Since last depth index was LMAX, IMOD has still the correct value
      IF (CORE) THEN
        BCOREL = PW * BCORE(KCMF,IMOD) + Q * BCORE(KCMF+1,IMOD)
        DBDRL  = PW * DBDR (KCMF,IMOD) + Q * DBDR (KCMF+1,IMOD)
      ENDIF
      
C***  PWEIGHT = WEIGHT p*dp at current P(JP) FOR THE FLUX INTEGRAL
      IF (JP .EQ. 1) THEN
         PWEIGHT = P(2)*P(2)/3.
      ELSE 
         PWEIGHT = (P(JP-1)+P(JP)+P(JP+1))*(P(JP+1)-P(JP-1))/3.
      ENDIF

      RETURN
      END
