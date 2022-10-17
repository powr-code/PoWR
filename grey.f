      SUBROUTINE GREY (ND,T,RADIUS,XLAMBDA,FWEIGHT,NF,ENTOT,RNE,RSTAR,
     $            ALPHA,SEXPO,
     $            ADDCON1, ADDCON2, ADDCON3, 
     $            IGAUNT,POPNUM,TAUROSS,R23,TEXIST,NDIM,N,
     $            LEVEL,NCHARG,WEIGHT,ELEVEL,EION,EINST,ENLTE,KODAT,
     $            ABXYZ,NOM,NFIRST,NLAST,NATOM,EXPFAC,SIGMAKI,NFEDGE,
     $            OPAC,ETAC,SIGMAFF,MAXION,MAXATOM,SIGMATHK,EDGEK,
     $            SEXPOK,KONTNUP,KONTLOW,LASTKON,XDATA, DENSCON, 
     >            FILLFAC)

C***********************************************************************
C***  COMPUTATION OF THE TEMPERATURE STRUCTURE
C***********************************************************************
C***  IF (TEXIST=.TRUE.), THE TEMPERATURE STRUCTURE IS NOT CALCULATED BUT
C***  IS ASSUMED TO BE EXISTING.
C***  ELSE IF (SPHERIC=.TRUE.) THE TEMPERATURE STRUCTURE IS CALCULATED
C***  AS FOR A SPHERICAL, GREY ATMOSPHERE. THE EDDINGTON FACTOR F IS
C***  APPROXIMATED BY ASSUMING A GEOMETRICALLY DILUTED RADIATION FIELD
C***  EMERGING FROM A SPHERE OF RADIUS R23 (I.E. THE RADIUS WHERE
C***  TAU-ROSSELAND= 2/3).
C***  THE SOLUTION OF THE 1. MOMENT EQUATION IS OBTAINED BY USING THE
C***  IMSL-ROUTINE DVERK. DTDR IS A USER-PROVIDED COEFFICIENT FUNCTION.
C***  ELSE
C***  THE TEMPERATURE STRUCTURE IS CALCULATED AS FOR A PLANE-PARALLEL,
C***  GREY ATMOSPHERE.
C***  ENDIF
C***
C***  IF (TMODIFY .NE. .0), THE CALCULATED TEMPERATURE STRUCTURE IS FINALLY
C***  MODIFIED BY A FACTOR R**TMODIFY .
C***  IF (TMIN .GT. 0), THE TEMPERATURE STRUCTURE IS MODIFIED NOT TO FALL
C***  BELOW TMIN.
C***  
C***  Clumping: LTE popnumbers is calculated with clump (electron) density, 
C***  Rosseland opacity is calculated with clump density, and then
C***  scaled down with the filling factor
C***********************************************************************
 
      COMMON /COMDTDR/ OPAMEAN,R23COM,R1COM,R13COM
      COMMON /COMTEFF/ TEFF,TMIN,TMODIFY,SPHERIC

C***  Operating system:
      COMMON / COMOS / OPSYS
      CHARACTER*8 OPSYS

      DIMENSION XLAMBDA(NF),FWEIGHT(NF),EXPFAC(NF)
      DIMENSION SIGMAKI(NF,LASTKON)
      DIMENSION T(ND),RADIUS(ND),RNE(ND),ENTOT(ND),TAUROSS(ND)
      DIMENSION NCHARG(NDIM),ENLTE(NDIM),EION(NDIM),ELEVEL(NDIM)
      DIMENSION POPNUM(ND,N)
      DIMENSION KONTNUP(LASTKON),KONTLOW(LASTKON),NFEDGE(LASTKON)
      DIMENSION NOM(N)
      DIMENSION KODAT(NATOM),ABXYZ(NATOM),NFIRST(NATOM),NLAST(NATOM)
      DIMENSION COMVEC(24),WORKSP(1,9)
      LOGICAL TEXIST, SPHERIC
      CHARACTER*10 LEVEL(N)
      EXTERNAL DTDR

c*** tiefenabh. clumping nach goetz...
      DIMENSION DENSCON(ND), FILLFAC(ND)

C***  GENERATE ONCE FOR ALL PHOTOCROSSSECTIONS AT ALL FREQUENCIES
C***  SIGMAKI(K,KON) IN CM**2
      CALL BFCROSS (SIGMAKI,NF,N,ELEVEL,EION,EINST,NDIM,
     $              XLAMBDA,ALPHA,SEXPO,
     $              ADDCON1, ADDCON2, ADDCON3, 
     $              IGAUNT,
     $              KONTNUP,KONTLOW,LASTKON)

C***  PRE-CALCULATE FREQUENCY INDICES OF IONIZATION EDGES
      DO 30 KON=1,LASTKON
      NUP=KONTNUP(KON)
      LOW=KONTLOW(KON)
      EDGELAM=1.E8/(ELEVEL(NUP)+EION(LOW)-ELEVEL(LOW))
      NFEDGE(KON)=ISRCHFGT(NF,XLAMBDA,1,EDGELAM) - 1
   30 CONTINUE

C***  ABUNDANCES OF HELIUM AND HYDROGEN FOR THE MODIFIED START APPROXIMATION
            ABHE=0.
            ABH=0.
            IF (KODAT(1) .GT. 0) ABH =ABXYZ(KODAT(1))
            IF (KODAT(2) .GT. 0) ABHE=ABXYZ(KODAT(2))
C***  changed according to WR-Memo.dir/070327.txt --  wrh 26-Sep-2008 

C***  FIRST RGRID-POINT:  R=RMAX (L=1)
      TAUROSS(1)=.0
      IF (.NOT. TEXIST) THEN
         IF (SPHERIC) THEN
            T(1)=TEFF*0.7071068/SQRT(RADIUS(1))
            ELSE
            T(1)=TEFF*0.8112
            ENDIF
         ENDIF
 
C***  ITERATION LOOP  **************************************************
C***  MAIN PURPOSE IS THE ITERATION OF R23 IN THE SPHERICAL CASE
C***  START VALUE OF THE RADIUS R23  (ONLY USED IN THE SPHERICAL CASE)
      R23COM=1.
      R1COM=1.
      R13COM=1.
      ITER=0
  100 ITER=ITER+1
      DTMAX=.0
C***  INITIALIZATION OF VARIABLES USED IN THE IMSL-SUBROUTINE DVERK
      TOL=0.0001
      IND=1
      NW=1
 
C***  LOOP OVER ALL DEPTH POINTS  *******************************************
c      write (*,*) 'iter=',iter,'!!!!!!!!!!!!!!'
      DO 10 L=1,ND
      TL=T(L)
C***  Clump density!
      ENTOTL = ENTOT(L) * DENSCON(L) 
      IF (ITER .GT. 1) TOLD=T(L+1)
C***  USE MODIFIED TEMPERATURE TMIN FOR THE CALCULATION OF THE EL.DENSITY
      TP=TL
      IF (TP.LT.TMIN) TP=TMIN
C***  COMPUTATION OF THE ROSSELAND MEAN OPACITY AT POINT L
C***  FIRST: LTE POPNUMBERS, ITERATION OF ELECTRON DENSITY
      RNEL=RNE(L)
C***  MODIFIED START APPROXIMATION FOR TEMPERATURE .LT. 10000K: 
C***  SAHA EQUATION FOR HEI/HEII - IONIZED HYDROGEN
      IF (ITER .EQ. 1 .AND. TP .LT. 10000.) THEN
            T32=TP*SQRT(TP)
            RNEL=ABH/2.+SQRT(ABH*ABH/4.+2.*EXP(-285645./TP)*T32/2.07E-16
     $            *ABHE/ENTOTL)
            ENDIF

C***  Iteration of electron density
    3 ENE=RNEL*ENTOTL
      CALL       LTEPOP (N,ENLTE,TP,ENE,WEIGHT,NCHARG,EION,ELEVEL,NOM,
     $                  ABXYZ,NFIRST,NLAST,NATOM)
      RNEOLD=RNEL
      RNEL=.0
      DO 2 J=1,N
    2 RNEL=RNEL+NCHARG(J)*ENLTE(J)
      RNEDIF=RNEL-RNEOLD
      IF (ABS(RNEDIF/RNEL).GT.0.01 .AND. ABS(RNEDIF).GT.0.001) GOTO 3

C***  STORE LTE POPNUMBERS TO BE WRITTEN AT THE MODEL FILE (START APPROXIMAT9ON@
      DO 5 J=1,N
    5 POPNUM(L,J)=ENLTE(J)
C*** ELECTRON DENSITY ALSO STORED
      RNE(L)=RNEL
 
      IF (L .EQ. ND) GOTO 10
 
      CALL OPAGREY (OPARL,ENLTE,TL,RNEL,ENTOTL,RSTAR,N,
     $              NCHARG,WEIGHT,ELEVEL,EION,NF,XLAMBDA,FWEIGHT,NOM,
     $              EXPFAC,SIGMAKI,NFEDGE,OPAC,ETAC,SIGMAFF,MAXION,
     $              SIGMATHK,SEXPOK,EDGEK,KODAT,MAXATOM,
     $              KONTNUP,KONTLOW,LASTKON,RADIUS(L),XDATA)
      OPARL = OPARL * FILLFAC(L)
 
C***  COMPUTATION OF THE ROSSELAND MEAN OPACITY  AT POINT L+1
C***  IN THE FIRST ITERATION USING T(L)
C***  IN THE FOLLOWING ITERATIONS USING TOLD
      IF (ITER .EQ. 1) THEN
         TL1=TL
         ELSE
         TL1=TOLD
         ENDIF
C***  USE MODIFIED TEMPERATURE TMIN FOR THE CALCULATION OF THE EL.DENSITY
         TP1=TL1
         IF (TP1.LT.TMIN) TP1=TMIN
C***  FIRST: LTE POPNUMBERS, ITERATION OF ELECTRON DENSITY
      RNEL=RNE(L+1)
C***  Clump density!
      ENTOTL1 = ENTOT(L+1) * DENSCON(L+1)
C***  MODIFIED START APPROXIMATION FOR TEMPERATURE .LT. 10000K: 
C***  SAHA EQUATION FOR HEI/HEII - IONIZED HYDROGEN
      IF (ITER .EQ. 1 .AND. TP1 .LT. 1.E4) THEN
            T32=TP1*SQRT(TP1)
           RNEL=ABH/2.+SQRT(ABH*ABH/4.+2.*EXP(-285645./TP1)*T32/2.07E-16
     $            *ABHE/ENTOTL1)
            ENDIF
   13 ENE=RNEL*ENTOTL1
      CALL       LTEPOP (N,ENLTE,TP1,ENE,WEIGHT,NCHARG,EION,ELEVEL,NOM,
     $                  ABXYZ,NFIRST,NLAST,NATOM)
      RNEOLD=RNEL
      RNEL=.0
      DO 12 J=1,N
   12 RNEL=RNEL+NCHARG(J)*ENLTE(J)
      RNEDIF=RNEL-RNEOLD
      IF (ABS(RNEDIF/RNEL).GT.0.01 .AND. ABS(RNEDIF).GT.0.001) GOTO 13
 
      CALL OPAGREY (OPARL1,ENLTE,TL1,RNEL,ENTOTL1,RSTAR,N,
     $              NCHARG,WEIGHT,ELEVEL,EION,NF,XLAMBDA,FWEIGHT,NOM,
     $              EXPFAC,SIGMAKI,NFEDGE,OPAC,ETAC,SIGMAFF,MAXION,
     $              SIGMATHK,SEXPOK,EDGEK,KODAT,MAXATOM,
     $              KONTNUP,KONTLOW,LASTKON,RADIUS(L),XDATA)
      OPARL1 = OPARL1 * FILLFAC(L+1)
 
C***  ARITHMETIC MEAN OF OPARL AND OPARL1
      OPAMEAN=0.5*(OPARL+OPARL1)
      TAUROSS(L+1)=OPAMEAN*(RADIUS(L)-RADIUS(L+1))+TAUROSS(L)
 
      IF (.NOT. TEXIST) THEN
         IF (SPHERIC) THEN
            RL=RADIUS(L)
            RL1=RADIUS(L+1)
c            IF (OPSYS .EQ. 'CRAY') THEN
c              CALL DVERK(1,DTDR,RL,TL,RL1,TOL,IND,COMVEC,NW,WORKSP,IER)
c            ELSE
C!!!              CALL RUKU (RL, RL1, TL, 10000)
              CALL RUKU (RL, RL1, TL, 100)
c            ENDIF
c            stop 'stop in grey nach first step!!!!!!!!!'
c            write (*,*) 'grey: rl1, tl1=',rl1,tl1
            T(L+1)=TL
            ELSE
C***        HOPF FUNCTION, C.F. UNSOELD P. 138
            Q=0.6940-0.1167*EXP(-1.9720*TAUROSS(L+1))
            T(L+1)=TEFF*(0.75*(TAUROSS(L+1)+Q))**0.25
            ENDIF
         ENDIF
 
C***  MAXIMUM TEMPERATURE CORRECTION
      IF (ITER .GT. 1) DTMAX=AMAX1(DTMAX,ABS(TOLD-T(L+1)))
 
   10 CONTINUE
c      stop 'stop in grey'
C*****************************************************************************
 
C***  CALCULATE RADIUS R23 WHERE TAUROSS=2/3
      TAU23=0.666666666666
      IF (TAUROSS(ND) .LT. TAU23) THEN
         R23=1.
         ELSE
         CALL LIPO (R23,TAU23,RADIUS,TAUROSS,ND)
         ENDIF
      R23COM=R23
      TAU1=1.
      IF (TAUROSS(ND) .LT. TAU1  ) THEN
      R1 COM=1.
         ELSE
         CALL LIPO (R1 COM,TAU1 ,RADIUS,TAUROSS,ND)
         ENDIF
      TAU13=0.333333333333
      IF (TAUROSS(ND) .LT. TAU13 ) THEN
      R13COM=1.
         ELSE
         CALL LIPO (R13COM,TAU13,RADIUS,TAUROSS,ND)
         ENDIF
 
      IF (ITER .LE. 1 .OR. 
     >    (DTMAX .GT. 10. .AND. ITER .LE. 20)) GOTO 100
C      IF (DTMAX .GT. 10.) GOTO 100
      IF (ITER .GT. 20) THEN
        WRITE (0,*) 'Max Number of Iterations exceeded in Subr. GREY'
      ENDIF
 
      IF (.NOT. TEXIST) THEN
         IF (T(ND) .LE. TEFF ) THEN
            T(ND)=TEFF
            T(ND-1)=TEFF
            ENDIF
         ENDIF
 
C***  MODIFY THE TEMPERATURE STRATIFICATION BY A FACTOR RADIUS**TMODIFY
C***  AND/OR BY REQUIRING A MINIMUM TEMPERATURE TMIN
      DO 4 L=1,ND
      IF (TMODIFY .EQ. .0 .AND. T(L) .GT. TMIN ) GOTO 4
      IF (TMODIFY .NE. .0) T(L)=T(L)*RADIUS(L)**TMODIFY
      IF (T(L) .LT. TMIN) T(L)=TMIN
 
C***     CALCULATE LTE POP.NUMBERS WITH MODIFIED TEMPERATURE]
            TL=T(L)
            T32=TL*SQRT(TL)
            ENTOTL = ENTOT(L) * DENSCON(L)
            RNEL=RNE(L)
C***     MODIFIED START APPROXIMATION FOR TEMPERATURE .LT. 10000K: 
C***  SAHA EQUATION FOR HEI/HEII - IONIZED HYDROGEN
            IF (TL.LT.1.E4) THEN
            RNEL=ABH/2.+SQRT(ABH*ABH/4.+2.*EXP(-285645./TL)*T32/2.07E-16
     $            *ABHE/ENTOTL)
            ENDIF
   23       ENE=RNEL*ENTOTL
      CALL       LTEPOP (N,ENLTE,TL,ENE,WEIGHT,NCHARG,EION,ELEVEL,NOM,
     $                  ABXYZ,NFIRST,NLAST,NATOM)
            RNEOLD=RNEL
            RNEL=.0
            DO 24 J=1,N
   24       RNEL=RNEL+NCHARG(J)*ENLTE(J)
            RNEDIF=RNEL-RNEOLD
            IF (ABS(RNEDIF/RNEL).GT.0.01 .AND. ABS(RNEDIF).GT.0.00001)
     $         GOTO 23
            DO 25 J=1,N
   25       POPNUM(L,J)=ENLTE(J)
            RNE(L)=RNEL
    4 CONTINUE
 
      RETURN
      END