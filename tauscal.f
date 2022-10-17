      SUBROUTINE TAUSCAL (RSTAR,ND,RADIUS,RNE,
     $                   ENTOT,T,POPNUM,NDIM,N,EN,LEVEL,NCHARG,WEIGHT,
     $                   ELEVEL,EION,EINST,ALPHA,SEXPO,
     $                   ADDCON1, ADDCON2, ADDCON3, 
     $                   IGAUNT,NOM,NF,
     $                   XLAMBDA,FWEIGHT,TAUTHOM,TAUROSS,
     $                   MAXATOM,SIGMATHK,SEXPOK,EDGEK,KODAT,
     $                   KONTNUP,KONTLOW,LASTKON, DENSCON, FILLFAC)
C***********************************************************************
C***  CALCULATION OF THE NLTE OPTICAL DEPTH SCALES (ROSSELAND, THOMSON)
C***   for continuum opacities
C***
C***  note: this routine updates the EN vector
C***********************************************************************
 
      IMPLICIT NONE

      INTEGER, INTENT(IN) :: ND, NDIM, N, MAXATOM, LASTKON, NF
      REAL, INTENT(IN) :: RSTAR

      REAL, DIMENSION(ND), INTENT(IN) :: RADIUS, RNE, ENTOT, T
      REAL, DIMENSION(ND), INTENT(INOUT) :: TAUTHOM, TAUROSS
      REAL, DIMENSION(ND, N), INTENT(IN) :: POPNUM
      REAL, DIMENSION(N), INTENT(INOUT) :: EN
      
      REAL, DIMENSION(N), INTENT(IN) ::  WEIGHT, ELEVEL, EION
      REAL, DIMENSION(N, N), INTENT(IN) :: EINST
      REAL, DIMENSION(LASTKON), INTENT(IN) :: 
     >                ALPHA, SEXPO, ADDCON1, ADDCON2, ADDCON3
      REAL, DIMENSION(MAXATOM, MAXATOM), INTENT(IN) :: 
     >                SIGMATHK, SEXPOK, EDGEK
      REAL, DIMENSION(NF), INTENT(IN) :: FWEIGHT, XLAMBDA
      CHARACTER(10), DIMENSION(N), INTENT(IN) :: LEVEL(NDIM)

      INTEGER, DIMENSION(N), INTENT(IN) :: NCHARG, NOM
      INTEGER, DIMENSION(LASTKON), INTENT(IN) :: 
     >                IGAUNT, KONTNUP, KONTLOW
      INTEGER, DIMENSION(LASTKON), INTENT(IN) :: KODAT

C*** tiefenabh. clumping nach goetz
      REAL, DIMENSION(ND), INTENT(IN) :: DENSCON, FILLFAC
 
      INTEGER :: L, I
      REAL :: TL, RL, ENTOTDENSL, RNEL, OPAMEAN, DR, RM1, 
     >        ENELM1, ENEL, OPARM1, OPARL, ENEMEAN, aOPARL
      LOGICAL :: bUseOPAROSSpre

C***  SIGMAE = ELCTRON SCATTERING CROSS SECTION  ( CM**2 )
      REAL, PARAMETER :: SIGMAE = 0.6652E-24
 
C***  INITIALIZATION
      TAUTHOM(1) = 0.0
      TAUROSS(1) = 0.0

C***  LOOP OVER ALL DEPTH POINTS  --------------------------------------
      DO L=1,ND
        RL=RADIUS(L)
        TL=T(L)
C***    Calculate opacities for clump density
        ENTOTDENSL=ENTOT(L) * DENSCON(L)
        RNEL=RNE(L)
        ENEL=RNEL*ENTOTDENSL
        DO I=1,N
          EN(I)=POPNUM(L,I)
        ENDDO
        !note that OPAROSS only calculates continuum opacities
        !(this is used in steal for taumax fixing)
        CALL OPAROSS (OPARL,EN,TL,RNEL,ENTOTDENSL,RSTAR,NDIM,N,
     >                  LEVEL,NCHARG,WEIGHT,ELEVEL,EION,EINST,
     >                  ALPHA,SEXPO,
     >                  ADDCON1, ADDCON2, ADDCON3, 
     >                  IGAUNT,NF,XLAMBDA,FWEIGHT,NOM,
     >                  MAXATOM,SIGMATHK,SEXPOK,EDGEK,KODAT,RL,
     >                  KONTNUP,KONTLOW,LASTKON)
        IF (OPARL <= 0.) THEN
          WRITE(0,*) " WARNING: negative Opacity ", OPARL
          WRITE(0,*) '   taking absolute value... '
          OPARL = ABS(OPARL)
        ENDIF
        OPARL = OPARL * FILLFAC(L) !downscale opacity with filling factor
        IF (L > 1) THEN
          DR = RM1 - RL
          ENEMEAN    = 0.5 * (ENEL+ENELM1)
          TAUTHOM(L) = TAUTHOM(L-1) + 
     >                  FILLFAC(L) * SIGMAE * ENEMEAN * DR * RSTAR
          OPAMEAN    = 0.5*(OPARL+OPARM1) 
          TAUROSS(L) = TAUROSS(L-1)+OPAMEAN*DR
        ENDIF
        RM1=RL
        ENELM1=ENEL
        OPARM1=OPARL
      ENDDO
C***  ENDLOOP  ---------------------------------------------------------
 
      RETURN
      END
