      SUBROUTINE TEMPEQ (I,L,DM,V1,NRANK,N,OPAC,ETAC,XLAMBDA,TL,XJLAPP,
     $   LASTIND,INDLOW,INDNUP,EINST,NDIM,OPAL,DELTAX,XMAX,
     $   XRED,XBLUE,DETAL,DOPAL,SLNEW,NFL,PHI,PWEIGHT,OPACIND,SCNEIND,
     $   VDOP,ELEVEL,RRATE,EN,ENTOTL,RSTAR,NF,ND,SCNEW,XJCAPP,
     $   FWEIGHT,
     $   WCHARM,DETA,DOPA,XJC,LASTKON,DRJLW,DRJLWE,DRLJW,IONGRND,
     $   KODRLOW,LASTKDR,IBLENDS,MAXLAP,XLAMZERO,BETA,PHIL,NBLENDS,
     $   ATEST,BTEST,BROYDEN, 
     $   OPAC1,RADIUS,ITNEL,TEFF,FEDDI,QEDDI,OPATHOM,
     >   WFLUX, WFLUXR, BCOLIUNLU, 
     >   HEDDI, HTOT, HTOTL, XJTOTL, 
     >   DTLOCAL, DTLOCC, DTLOCL, DTLOCD, DTINT, DTRMAX, NDCORR,
     >   WFELOW, WFENUP, BDIAG, HTOTCMF0, FTCOLI, 
     >   OPASMEAN, QFJMEAN, OPAJMEAN, QOPAHMEAN, EDDIHOUTJMEAN, 
     >   HTOTOUTMINUS, FTFE, WCORREC,WEIGHT)
C***********************************************************************
C***  THIS SUBROUTINE FILLS THE ADDITIONAL COLUMN (N+2)
C***  IN THE DERIVATIVE MATRIX DM AND THE RIGHT-HAND SIDE VECTOR V1
C***  NOTE: THE ADDITIONAL COLUMN IN RATCO REMAINS ZERO
C***  THIS DESCRIBES THE COMPLETELY-LINEARIZED ENERGY EQUATION
C***  NOTE: THIS VERSION ACCOUNTS FOR CONTINUA + DIEL.REC. +  LINES (!)
C***  NOTE CONCERNING THE LINE CONTRIBUTION:
C***  INSTEAD OF WRITING: ETAL - OPAL * XJL, WE MAKE USE OF THE
C***  RADIATIVE RATE COEFFICIENTS RRATE.
C***  THIS HAS THE NUMERICAL ADVANTAGE, THAT IT IMPLIES THE
C***  NET-RADIATIVE-BRACKET FORMULATION (IN CASE OF NON-ZERO LINE CORES)
C***  OPTIONALLY, A DIFFERENTIAL EQUATION ENFORCING FLUX CONSERVATION 
C***  IS ADDED WITH A PRE-SPECIDIED WEIGHT WFLUX >0. 
C***  THE EQUATIONS ARE IN CGS UNITS.
C***  CALLED FROM: SUBROUTINE COMA
C***********************************************************************

      DIMENSION SCNEW(NF),XJCAPP(NF),FWEIGHT(NF),DETA(NF),DOPA(NF)
      DIMENSION RADIUS(ND), OPAC1(NF)
      DIMENSION OPAC(NF),ETAC(NF),XLAMBDA(NF)
      DIMENSION WCHARM(ND,NF),XJC(ND,NF)
      DIMENSION DM(NRANK,NRANK),V1(NRANK)
      DIMENSION EINST(NDIM,NDIM),ELEVEL(NDIM)
      DIMENSION INDLOW(LASTIND),INDNUP(LASTIND)
      DIMENSION SLNEW(LASTIND),OPAL(LASTIND),XJLAPP(LASTIND)
      DIMENSION DETAL(LASTIND),DOPAL(LASTIND)
      DIMENSION XRED(LASTIND), XBLUE(LASTIND) 
      DIMENSION RRATE(NDIM,NDIM),EN(NRANK)
      DIMENSION KODRLOW(LASTKDR),IONGRND(NDIM)
      DIMENSION DRJLW(NDIM),DRJLWE(NDIM),DRLJW(NDIM)
      DIMENSION FEDDI(ND,NF), QEDDI(ND,NF), HEDDI(NF)
      DIMENSION HTOT(ND), HTOTL(ND), XJTOTL(ND)
      DIMENSION DTLOCAL(ND), DTLOCC(ND), DTLOCL(ND)
      DIMENSION DTLOCD(ND), DTINT(ND), DTRMAX(ND)
      DIMENSION FTCOLI(ND)
      DIMENSION WFELOW(ND,LASTIND), WFENUP(ND,LASTIND)
      DIMENSION HTOTCMF0(ND)
      DIMENSION OPASMEAN(ND), QFJMEAN(ND), OPAJMEAN(ND)
      DIMENSION QOPAHMEAN(ND)
      DIMENSION WEIGHT(NDIM)

      LOGICAL BROYDEN, BDIAG(LASTIND), BCOLIUNLU, BFIRSTITER
      DATA LASTL / 0 /

      SAVE HNULL, DHINTSUM, DHRMAX1, LASTL


C***  C1 = H * C / K        (CM * KELVIN)
      DATA C1 / 1.4388 /
C***  C2 = 2 * H * C     ( CGS UNITS )
      DATA C2 / 3.9724E-16 /
C***  DATA: STEBOL = STEFAN-BOLTZMANN CONSTANT (CGS-UNITS) / PI
      DATA STEBOL / 1.8046E-5 /
C***  PI8 = 8*PI
      DATA PI8 /25.1327412288 /

C***  BCOLIUNLU IST NOW TAKEN FROM CARDSFILE!!!!!!!!!! 
C***  if BCOLIUNLU = .true. OPASMEAN, QFJMEAN, OPAJMEAN, QOPAHMEAN
C***  are taken from ROUTINE COLI  


C******************************************************
C***  BRANCH FOR INNER BOUNDARY POINT ONLY
C******************************************************
      IF (L .NE. ND) GOTO 22

C***  OMIT CALCULATION OF DERIVATIV MATRIX
      IF (BROYDEN) GOTO 51

C***  MATRIX ELEMENT RATCO(I,N+2): DERIVATIVE OF THE ENERGY EQ. WITH RESPECT
C***  TO EN(I)
      DFDNI=.0
      DO 11 K=1,NF
      BRACKET = DOPA(K) * (BNUE(XLAMBDA(K),TL)-XJCAPP(K)) +
     +      WCHARM(L,K) * (DOPA(K)*SCNEW(K)-DETA(K)) 
      IF (I .EQ. NRANK) 
     $  BRACKET = BRACKET + OPAC(K)*DBNUEDT(XLAMBDA(K),TL)
      DFDNI=DFDNI + FWEIGHT(K) * BRACKET
   11 CONTINUE
C***  CONVERT INTO CGS UNITS 
      DFDNI = DFDNI / RSTAR

      DM(I,NRANK) = DFDNI

51    CONTINUE

      IF (I .EQ. NRANK) THEN
C***    FILL THE RIGHT-HAND SIDE VECTOR ELEMENT, V1(NRANK)
        FT =.0
        FT2=.0
        DO K=1,NF
          FT = FT + 
     >       FWEIGHT(K) * OPAC(K) * (BNUE(XLAMBDA(K),TL) - XJCAPP(K))
        ENDDO
C***    CONVERT INTO CGS UNITS
        FT = FT / RSTAR
        V1(NRANK)=-FT
C***  END OF LOOP "IF (I .EQ. NRANK) ..."
      ENDIF

      RETURN
C***  END OF SPECIAL PART FOR INNERMOST DEPTH POINT (L=ND) ONLY
C*********************************************************************

   22 CONTINUE

C***********************************************************
C***  BRANCH FOR ALL DEPTH POINTS EXCEPT THE INNER BOUNDARY:
C***********************************************************

C***  OMIT CALCULATION OF DERIVATIV MATRIX IF NOT REQUIRED
      IF (BROYDEN) GOTO 50

C***  CLASSICAL EQUATION: LOCAL ENERGY BALANCE

C***  CONTINUA: INTEGRAL OVER FREQUENCY  ============================
      DFCDNI=.0
      DO K = 1,NF
         DFCDNI = DFCDNI + ( DETA(K)*(1-WCHARM(L,K)) - 
     >         DOPA(K)*(XJCAPP(K)-WCHARM(L,K)*SCNEW(K)) )*FWEIGHT(K)
      ENDDO
C***  CONVERT TO CGS UNITS
      DFCDNI = DFCDNI / RSTAR

C***  DIEL. RECOMBINATION/AUTOIONIZATION: LOOP OVER ALL DR CONTINUA  
C***  NOTE: ALL DR-RATES IN COMBINATION WITH WAVENUMBER OF STABILIZING 
C***        TRANSITIONS OR ENERGY OF AUTOIONIZATION LEVELS ARE
C***        PRECALCULATED IN SUBR. RADNET
      TLINV=1./TL
      ENIINV = 1./ EN(I)
      DFDRDNI=.0
      DO 14 KDR=1, LASTKDR
       LOW = KODRLOW(KDR)
       NUP = IONGRND(LOW)
       IF ( I .EQ. NUP ) THEN
          DFDRDNI = DFDRDNI + DRJLW(LOW)
       ELSE IF ( I .EQ. LOW ) THEN
          DFDRDNI = DFDRDNI -DRLJW(LOW)
       ELSE IF ( I .EQ. N+1 ) THEN
          DFDRDNI = DFDRDNI + EN(NUP) * DRJLW(LOW) * ENIINV
       ELSE IF ( I .EQ. NRANK ) THEN
          DFDRDNI = DFDRDNI + EN(NUP) * 
     *                     (-1.5*DRJLW(LOW)+C1*TLINV*DRJLWE(LOW))*TLINV
       ENDIF
   14 CONTINUE
C***  END-OF-LOOP OVER ALL DR CONTINUA  ==============================
C***  CONVERT INTO CGS UNITS: RATES * h nue / 4 pi, ABSOLUTE POPNUMS      
      DFDRDNI = DFDRDNI * ENTOTL * C2 / PI8 

C***  LINES: LOOP OVER ALL BOUND-BOUND TRANSITIONS  ====================
      DFLDNI=.0
      DO 4 IND=1,LASTIND
       LOW=INDLOW(IND)
       NUP=INDNUP(IND)

C***   THE DERIVATIVES OF THE CONT. BACKGROUND-OPACITY AND EMISSIVITY
C***   ARE NEGLECTED. HENCE:
       IF ((I .NE. NUP) .AND. (I .NE. LOW)) GOTO 4
       W = ELEVEL(NUP) - ELEVEL(LOW)

C***   ZERO LINE CORE (ALSO VALID FOR RUDIMENTAL LINES):
C***   ZERO DERIVATIVE OF JL, NO NET RADIATIVE BRACKET
       IF (XRED(IND) .GE. XBLUE(IND)) THEN
         SUM = 0.
         IF (I .EQ. NUP) SUM= RRATE(NUP,LOW)
         IF (I .EQ. LOW) SUM=-RRATE(LOW,NUP)      

C***   Iron line: no Net Radiative Bracket
         IF (BDIAG(IND)) THEN
C***        DERIVATIVES FOR IRON-LINES ARE TAKEN FROM DIAGONAL WEIGHTS
           IF (I .EQ. LOW .OR. I .EQ. NUP) THEN
              IF(I .EQ. LOW) DJLDNI = WFELOW(L, IND)
              IF(I .EQ. NUP) DJLDNI = WFENUP(L, IND)
              WAVENUM=ELEVEL(NUP)-ELEVEL(LOW)
              W3=WAVENUM*WAVENUM*WAVENUM
              DLOWUP = (EN(NUP) - EN(LOW) * WEIGHT(NUP)/WEIGHT(LOW) ) *
     >     	       EINST(NUP,LOW) * DJLDNI / (C2*W3)
              SUM = SUM + DLOWUP
           ENDIF
         ENDIF

       ELSE
C***   Non-iron lines with core: Net Radiative Bracket
C***     KETTENREGEL: 1. DERIVATIVE OF APPROXIMATE LINE RADIATION FIELD:
         CALL JLDERIV (DJLDNI,IND,DELTAX,XMAX,XRED,XBLUE,
     $     DETAL,DOPAL,SLNEW,OPAL,NFL,PHI,PWEIGHT,OPACIND,SCNEIND,
     $     IBLENDS,MAXLAP,XLAMZERO,BETA,PHIL,NBLENDS,
     $     VDOP,I,INDNUP,INDLOW,ATEST,BTEST)
C***     (NOTE THAT DJLDNI IS THE NEGATIVE (!!) DERIVATIVE OF XJLAPP!
         DZDNI =  DJLDNI / SLNEW(IND)
C***     ...  2. DERIVATIVE OF 1/SLNEW:
         ETAL=OPAL(IND)*SLNEW(IND)
         DSLINV=(DOPAL(IND)*ETAL-OPAL(IND)*DETAL(IND))/(ETAL*ETAL)
         DZDNI=DZDNI-XJLAPP(IND)*DSLINV
         SUM = EN(NUP) * DZDNI * EINST(NUP,LOW)
C***     DERIVATIVE WITH RESPECT TO NUP REPRODUCES THE NET RATE (UP-LOW):
         IF (I .EQ. NUP) SUM = SUM + RRATE(NUP,LOW)
       ENDIF

       DFLDNI = DFLDNI + W*SUM
    4 CONTINUE
C***  END-OF-LOOP OVER ALL BOUND-BOUND TRANSITIONS  ====================

C***  CONVERT INTO CGS UNITS: RATES * h nue / 4 pi, ABSOLUTE POPNUMS      
      DFLDNI = DFLDNI * ENTOTL * C2 / PI8

C***  ADD ALL CONTRIBUTIONS: CONTINUUM + DIEL.RECOMBINATION + LINE
      DFDNI = DFCDNI + DFDRDNI + DFLDNI

C***  Note: The right-hand side with the flux conservation terms is assumed to 
C***        have no derivatives!

C***  MATRIX ELEMENT DM(I,N+2): DERIVATIVE OF THE ENERGY EQ. WITH RESPECT
C***  TO EN(I)
      DM(I,NRANK) = DFDNI


50    CONTINUE
C***  END OF CALCULATION OF THE DERIVATIVES (LEFT_HAND SIDE MATRIX)
C***  THIS PART WAS JUMPED OVER IF "BROYDEN"



C**********************************************************************
C***  FILL THE RIGHT-HAND SIDE VECTOR ELEMENT, V1(NRANK)
C***  Note: as long as I .LT. NRANK, only derivative terms are
C***        preperaed for the left-hand side matrix 
      IF (I .NE. NRANK) RETURN
C**********************************************************************



C***  CHECK IF CURRENT DEPTHPOINT IS ENCOUNTERED FOR THE FIRST TIME
      BFIRSTITER = LASTL .NE. L 
      LASTL = L

C***  PREPARATIONS WHICH MUST BE DONE ONLY IN THE FIRST ITERATION AT 
C***  EACH DEPTH POINT
      RL2 = RADIUS(L) * RADIUS(L)
      HNULL = 0.25 * STEBOL * TEFF*TEFF*TEFF*TEFF

C**************************************************************
C***  FLUX CONSERVATION TERMS following the Unsoeld-Lucy method
C**************************************************************

C *** calculate OPASCMEAN also for OUTPUT ONLY in case of BCOLIUNLU Branch 
      IF (BCOLIUNLU) THEN
        SUM   = 0.
        SMEAN = 0.
        DO K=1, NF
           SUM   = SUM + OPAC(K) * SCNEW(K) * FWEIGHT(K)
           SMEAN = SMEAN + SCNEW(K) * FWEIGHT(K)
        ENDDO
        OPASCMEAN = SUM / (SMEAN * RSTAR)
      ENDIF

C *** IF BCOLIUNLU =.TRUE. take values for OPASMEAN,OPAJMEAN,
C *** QFJMEAN and QOPAHMEAN from coli
 
      IF (.NOT. BCOLIUNLU) THEN

C***  OPASMEAN = Source-function weighted mean of opacity (CGS UNITS)
C***             - only used for output: conversion into temperatures
        SUM   = 0.
        SMEAN = 0.
        DO K=1, NF
           SUM   = SUM + OPAC(K) * SCNEW(K) * FWEIGHT(K)
           SMEAN = SMEAN + SCNEW(K) * FWEIGHT(K)
        ENDDO
C***  Conversion into cgs-units
        OPASMEAN(L) = SUM / (SMEAN * RSTAR)


C***  QFJMEAN = J-weighted mean of (qf) at RL
        SUM   = 0.
        XJMEAN = 0.
        DO K=1, NF
           SUM = SUM + QEDDI(L,K) * FEDDI(L,K) * XJCAPP(K) * FWEIGHT(K)
           XJMEAN = XJMEAN + XJCAPP(K) * FWEIGHT(K)
        ENDDO
          QFJMEAN(L) = SUM / XJMEAN


C***  OPAJMEAN = J-weighted mean of OPACITY at RL
C***  NOTE: OPAJMEAN without Thomson-Scatter
        SUM   = 0.
        DO K=1, NF
           SUM   = SUM + OPAC(K) * XJCAPP(K) * FWEIGHT(K)
        ENDDO
C***  Conversion into cgs-units
        OPAJMEAN(L) = SUM / (XJMEAN * RSTAR)

      ENDIF


C***  CLASSICAL FORMULATION: LOCAL ENERGY BALANCE

C***  CONTINUA:  INTEGRAL OVER FREQUENCY  ==============================
      FTCONT=.0
      DO 2 K=1,NF
        FTCONT = FTCONT + (ETAC(K) - OPAC(K)*XJCAPP(K)) * FWEIGHT(K)
    2 CONTINUE
C***  CONVERT INTO CGS UNITS
      FTCONT = FTCONT / RSTAR


C***  DIEL. RECOMBINATION/AUTOIONIZATION: LOOP OVER ALL DR CONTINUA
      FTDR = .0
      DO 33 KDR=1, LASTKDR
        LOW = KODRLOW(KDR)
        NUP = IONGRND(LOW)
        FTDR = FTDR + EN(NUP) * DRJLW(LOW) - EN(LOW) * DRLJW(LOW)
   33 CONTINUE
C***  CONVERT INTO CGS UNITS: RATES * h nue / 4 pi, ABSOLUTE POPNUMS      
      FTDR= FTDR * ENTOTL * C2 / PI8 

C***  LINES: LOOP OVER ALL BOUND-BOUND TRANSITIONS
      FTLINES = .0
      DO 3 IND=1, LASTIND
        LOW = INDLOW(IND)
        NUP = INDNUP(IND)
        W = ELEVEL(NUP) - ELEVEL(LOW)
        FTLINES = FTLINES + W * 
     *    (EN(NUP) * RRATE(NUP,LOW) - EN(LOW) * RRATE(LOW,NUP))
    3 CONTINUE
C***  CONVERT INTO CGS UNITS: RATES * h nue / 4 pi, ABSOLUTE POPNUMS      
      FTLINES = FTLINES * ENTOTL * C2 / PI8 

C***  ADD ALL CONTRIBUTIONS: CONTINUUM + DIEL.RECOMBINATION + LINE 
      FT = FTCONT + FTDR + FTLINES


C**********************************************************************
C***  if option "FLUXCORR = ndcorr" is set, 
C***  FT is corrected such that for the first inner iteration
C***  (ITNEL = 1, i.e. still with the formal-solution radiation field)
C**   FTCONT + FTLINES becomes = FTCOLI , which was calculated in COLI.
C***  (FTDR is left aside, because dielectronic lines are not accounted 
C***  for in COLI as well.) 
C***********************************************************************

C***  Store FTCONT + FTLINES from first inner Iteration 
      IF (BFIRSTITER) THEN
         FTNULL = FTCONT + FTLINES
      ENDIF


C***  Add Difference of FTCOLI - FTNULL to FT
C***  The correction is smoothly switched off for depth indices below NDCORR
      IF (NDCORR .LT. ND .AND. L .LT. ND) THEN
         IF (L .LT. NDCORR) THEN
            FCORRFAC = EXP(-(FLOAT(NDCORR)/FLOAT(L))**2)
         ELSE
            FCORRFAC = 1. 
         ENDIF
         FT = FT + (FTCOLI(L) - FTNULL) * FCORRFAC
      ENDIF

C***  Jump calculation of further UNLU quantities if taken from COLI
      IF (BCOLIUNLU) GOTO 71

      IF (L .EQ. 1) THEN
C***********************************************************************
C***    Prepare third Term of Unsoeld-Lucy-Procedure: outer Boundary
C***********************************************************************
C***    EDDIHOUTJMEAN = J-weighted mean of (H) at R1 (Outer Boundary)
        SUM   = 0.
        DO K=1, NF
           SUM   = SUM + HEDDI(K) * XJCAPP(K) * FWEIGHT(K)
        ENDDO
        EDDIHOUTJMEAN = SUM / XJMEAN

        DHRMAX1 = HTOTCMF0(1) / EDDIHOUTJMEAN - XJMEAN * RL2

        DHINTSUM = .0
      ELSE
C***********************************************************************
C***    Second Term of Unsoeld-Lucy-Procedure: contribution to depth integral
C***    H-total at radius-interstice (L,L-1)
C***    Note the different indexing of HTOT... from 1...ND-1
C***    Note: The continuum flux is taken from COMO, because the 
C***          predicted approximate XJCAPP is not suitable for taking a 
C***          depth derivative. 
C***          For the same reason, QOPAHMEAN is weighted with HNUE obtained 
C***          from differencind the formal-solution XJC 
C***          HTOTCLM1  guarantees exact normalization of the weighting
C***********************************************************************
        HTOTCLM1 = 0.
        SUM = .0
        DO 6 K=1, NF
          DR = RADIUS(L) - RADIUS(L-1)
          QL   = QEDDI(L  ,K)
          QLP1 = QEDDI(L-1,K)
          FL   = FEDDI(L  ,K)
          FLP1 = FEDDI(L-1,K)
          QFRL = QL * FL * RL2
          QFRLP1 = QLP1 * FLP1 * RADIUS(L-1) * RADIUS(L-1)
          X    = 0.5 * (OPAC(K) + OPATHOM + OPAC1(K))
C***      NOTE: OPAC1 IS THE OPACITY FROM THE PREVIOUS DEPTH POINT L-1
C***          AND CONTAINS ALREADY THE THOMSON CONTRIBUTION
C***          (CF. SUBR. LINPOP)
          Q = 0.5 * (QL + QLP1)
          DTQ = X * Q * DR
          HNUE = - (QFRL * XJC(L,K) - QFRLP1 * XJC(L-1,K)) / DTQ
          HTOTCLM1 = HTOTCLM1 + HNUE * FWEIGHT(K)

C***      Occasional computation of QOPAHMEAN
C***        = FLUX-weighted mean of q*opacity (CGS UNITS) at interstice (L,L-1)
          SUM = SUM + Q * X * HNUE * FWEIGHT(K)

    6   CONTINUE

C***    Normalization to H

        QOPAHMEAN(L-1) = SUM / HTOTCLM1

        DH = HTOTL(L-1) - HTOTCMF0(L-1)

        WRADIUS = RADIUS(L) - RADIUS(L-1)

C***    Save integral contribution of current depth pount from last inner
C***      iteration in order to replace it by the current iteration value  
        IF (BFIRSTITER) THEN
           DHLOLD = .0
        ELSE
           DHLOLD = DHL
        ENDIF

C***    Integral summand for actual L,L-1 interval
        DHL = DH * WRADIUS * QOPAHMEAN(L-1)

        DHINTSUM = DHINTSUM + DHL - DHLOLD

C***    Factors outside integral
        DHINT = OPAJMEAN(L) * DHINTSUM / (RL2 * QFJMEAN(L))
 
      ENDIF

C***********************************************************************
C***  Third Term of Unsoeld-Lucy-Procedure: Flux at outer Boundary
C***********************************************************************

      DHRMAX = OPAJMEAN(L) * DHRMAX1 * QFJMEAN(1) / (QFJMEAN(L) * RL2) 


C******************************************************************
C***  END OF BRANCH FOR BCOLIUNLU= .FALSE.
C******************************************************************
      GOTO 72



C*******************************************************************
C *** NOW THE NEW BRANCH WITH VALUES TAKEN FROM MAINPROGRAM: COLI
C *** i.e. if BCOLIUNLU = .TRUE
C*******************************************************************
 71   CONTINUE 

C***  In the BCOLIUNLU version, the flux conservation terms are not 
C***  updated in the inner iteration. Therefore:
      IF (.NOT. BFIRSTITER) GOTO 72

C *** DEFINE FACTOR FOR Unsoeld-Lucy-Terms
      ODQF= OPAJMEAN(L) / QFJMEAN(L)
 

      IF (L .EQ. 1) THEN
C***********************************************************************
C***    Prepare third Term of Unsoeld-Lucy-Procedure: outer Boundary
C***********************************************************************
C *** New treatment with h+ and H-;   EDDIHOUTJMEAN = <h+>_JMEAN


        DHRMAX1 = (HTOTCMF0(1) - HTOTOUTMINUS) / EDDIHOUTJMEAN - 
     >             XJTOTL(1) 
        DHINTSUM = .0
      ELSE
C***********************************************************************
C***    Second Term of Unsoeld-Lucy-Procedure: integrated flux deviation
C***********************************************************************

        DH = HTOTL(L-1) - HTOTCMF0(L-1)
        WRADIUS = RADIUS(L) - RADIUS(L-1)


        DHINTSUM = DHINTSUM + DH * WRADIUS * QOPAHMEAN(L-1)

C***    Factors outside integral
        DHINT = ODQF * (DHINTSUM /RL2)
      ENDIF

C***********************************************************************
C***  Third Term of Unsoeld-Lucy-Procedure: Flux at outer Boundary
C***********************************************************************
      DHRMAX = ODQF * (DHRMAX1 * QFJMEAN(1) / RL2) 


C**********************************************************************
C***  Final composition of the right-hand side
C**********************************************************************
   72 CONTINUE

C***  Add the contributions into the right-hand side vector
C***  NOTE : Not active for inner boundary
      IF (L .EQ. ND) THEN
        V1(NRANK) = -FT
      ELSE
        V1(NRANK) = -FT + WFLUX * DHINT + WFLUXR * DHRMAX
      ENDIF
 

C***  Temperature correction (output only!)
C***  Under the assumption S = Planck, the error term FT is converted
C***  into a temperature correction (in kK)

      DTB     = 1. / (TL * TL * TL * 4. * STEBOL * OPASMEAN(L) * 1000.)
      IF (BCOLIUNLU) THEN
       DTBC    = 1. / (TL * TL * TL * 4. * STEBOL * OPASCMEAN * 1000.)
      ELSE
       DTBC = DTB
      ENDIF
      DTLOCAL(L) = FT      * DTB
c      DTLOCC(L)  = FTCONT  * DTBC
      DTLOCC(L)  = FTCONT  * DTB
      DTLOCL(L)  = FTLINES * DTB
      DTLOCD(L)  = FTDR    * DTB
      DTINT(L)   = DHINT   * DTB
      DTRMAX(L)  = DHRMAX  * DTB

      RETURN

      END
