      SUBROUTINE DIFDTDR (DTDR,TEFF,XJC,HEDDI,TRND,RND,ND,EN,POPNUM,
     $                   RNEND,ENTOTND,RSTAR,NDIM,N,LEVEL,NCHARG,WEIGHT,
     $                   ELEVEL,EION,EINST,ALPHA,SEXPO,
     $                   ADDCON1, ADDCON2, ADDCON3, 
     $                   IGAUNT,NOM,NF,
     $                   XLAMBDA,FWEIGHT,
     $                   MAXATOM,SIGMATHK,SEXPOK,EDGEK,KODAT,
     $                   KONTNUP,KONTLOW,LASTKON, DENSCON, FILLFAC, 
     >                   BKUDRITZKI,OPARND)
C***********************************************************************
C***  CALCULATES "DTDR", THE TEMPERATURE GRADIENT (-DT/DR !!) 
C***    AT THE INNER BOUNDARY 
C***  - ONLY CALLED, IF TEMPERAURE CORRECTIONS ARE TO BE APPLIED
C***  - CALLED FROM: WRCONT, COMO, ETL
C***  - "DTDR" IS USED BY: DIFFUS, MOMO
C***  FORMULATION: KUDRITZKI 1973, THESIS TU BERLIN
C***  CONCEPT: 
C***  INCIDENT RADIATION I+ = B-NUE + MUE * D(BNUE)/DTAU-ROSS.
C***                          (DIFFUSION APPROXIMATION)
C***  EMERGENT RADIATION (FLUX H- LEAVING ACROSS THE INNER BOUNDARY):
C***     CALCULATED FROM THE ACTUAL SOLUTION J, AND HEDDI (EDDINGTON FACT.)
C***********************************************************************
 
      DIMENSION FWEIGHT(NF),XLAMBDA(NF),XJC(ND,NF),HEDDI(NF)
      DIMENSION EN(N),POPNUM(ND,N)
      LOGICAL BKUDRITZKI
 
C***  tiefenabh. clumping nach goetz
      DIMENSION DENSCON(ND),FILLFAC(ND)

C***  DATA: STEBOL = STEFAN-BOLTZMANN CONSTANT (CGS-UNITS) / PI
      DATA STEBOL / 1.8046E-5 /
 
C***  CORRECT TOTAL FLUX AS DETERMINED BY THE GIVEN EFFECTIVE TEMPERATURE
      HTOT=0.25*STEBOL*TEFF*TEFF*TEFF*TEFF

C      WRITE (0,*) 'OPAROSS-COLI:', OPARND
 
      IF (OPARND .LE. 0.) THEN
C***  CALCULATION OF THE ROSSELAND MEAN OPACITY AT R=R*
        DO 20 I=1,N
   20   EN(I)=POPNUM(ND,I)
C***  Clump density
        ENTOTNDDENS = ENTOTND * DENSCON(ND)
        CALL OPAROSS (OPARND,EN,TRND,RNEND,ENTOTNDDENS,RSTAR,NDIM,N,
     $              LEVEL,NCHARG,WEIGHT,ELEVEL,EION,EINST,
     $              ALPHA,SEXPO,
     $              ADDCON1, ADDCON2, ADDCON3, 
     $              IGAUNT,NF,XLAMBDA,FWEIGHT,NOM,
     $              MAXATOM,SIGMATHK,SEXPOK,EDGEK,KODAT,RND,
     $              KONTNUP,KONTLOW,LASTKON)
C***  Scale opacity down with filling factor
        OPARND = OPARND * FILLFAC(ND)
C        WRITE (0,*) 'OPAROSS-COMO:', OPARND
      ENDIF

C***  INTEGRAL "HINT" ACCOUNTS FOR THE POSSIBLE DEVIATION
C***  OF I- FROM THE DIFFUSION APPROXIMATION:
      HINT=0.0
      IF (BKUDRITZKI) THEN
        DO K=1,NF
          HINT = HINT + 
     >      (.5*BNUE(XLAMBDA(K),TRND)-HEDDI(K)*XJC(ND,K))*FWEIGHT(K)
        ENDDO
      ENDIF

      DTDR=0.75*OPARND*(HTOT-HINT)/STEBOL/TRND/TRND/TRND
 
      RETURN
      END
