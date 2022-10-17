      SUBROUTINE FREQUINT (K, FWEIGHTL, XLAMK, XLAMKOLD,    
     >                   XLAMBDA, 
     >                   XJL, XJLOLD, XHL, XHLOLD, 
     >                   OPAKOLD, NF, ND, NDDIM,
     >                   KONTACT, KONTAUP, DFKONT, LIND, PWEIGHT, 
     >                   MAXLIN, INDFEACT, SIGMAACT, MAXFEACT, LASTFE,
     >                   OPAK, ETAK, OPAFEI, ETAFEI, IFENUP, IFELOW,
     >                   DJDS, DJDO, DJDS_OLD, DJDO_OLD,
     >                   WEIGHT, N, POPNUM, XKLOLD,XNLOLD,
     >                   OPAKNOTH, ETAKNOTH, OPAKNOTHO, ETAKNOTHO,
     >                   EDDIFO, EDDIGO, RADIUS, BCOLIRAY,
     >                   ETANOTH, OPA, THOMSON, T, ETAKOLD,
     >                   QLFOLD, QLHOLD, OPAKHOLD,
     >                   EDDIHOUTOLD, EDDIHINOLD, XHID,
     >                   FWEIGHT, OPAO, THOMSONO,
C***  Temporary vectors
     >                   SUMJ, SUMJW, SUMDJDSC, SUMDJDSCW, SC, SCO, 
C***             OUTPUT:
     >                   XJCINT, FWTEST, XJLMEAN, XJFEMEAN,
     >                   HTOTL, HTOTND, HTOTNDCOR,
     >                   ARAD, ACONT, ATHOM,
     >                   XJTOTL, XKTOTL, XNTOTL, WFELOW, WFENUP, FTCOLI,
     >                   DSDFELOW, DSDFENUP, WJC, WJC_MIN, 
     >                   DSDSC, DJDSC, DJDSCOLD,
     >                   DBDTINT, DBDTOPAINT, DBDTINT_M, DBDTOPAINT_M,
C***       Unsoeld-Lucy:
     >                   OPASMEAN, OPASMEANTC, SMEAN, QFJMEAN, OPAJMEAN,
     >                   OPAJMEANTC, OPAPMEAN, QOPAHMEAN, HMEAN, 
     >                   EDDIHOUTJMEAN, XHOM, HTOTOUTMINUS,
     >                   RADIUS2, OPC, 
     >                   FERATUL, FERATLU, ELEVEL, 
     >                   EMCOLI, FTFE, IVERS_FE_EXPFAC, LPLOT_WCHARM, 
     >                   GAMMACOLI, OPAROSS, OPALAMBDAMEAN, 
     >                   GAMMAT, UNLU_TAUMAX, UNLU_TAUMAX2, TEFF,
C*** die folgenden SKALAREN Parameter werden ausgereicht, um den
C*** Compiler-Bug seit Compaq Tru64 UNIX V5.0A (Ende 2000) zu umgehen
     >                   XNUEK, XNUEKOLD, XNUEKOLDOLD)
C***  CLIGHT = SPEED OF LIGHT IN ANGSTROEM/SECOND
      DATA CLIGHT / 2.99792458E18 /
      DATA PI8 /25.1327412288 /
C***  C1 = H * C / K        (CM * KELVIN)
      DATA C1 / 1.4388 /
C***  C2 = 2 * H * C     ( CGS UNITS )
      DATA C2 / 3.9724E-16 /

      DATA WJMAX / 0.9999999999 /

      DIMENSION XLAMBDA(NF),XJL(ND),XJLOLD(ND),XHL(ND),XHLOLD(ND)
      DIMENSION OPAKOLD(ND), ETAKOLD(ND)
      DIMENSION XKLOLD(ND),XNLOLD(ND)
      DIMENSION LIND(MAXLIN), PWEIGHT(MAXLIN)
      DIMENSION INDFEACT(MAXFEACT), SIGMAACT(MAXFEACT)
      DIMENSION RADIUS(ND), OPAKNOTHO(ND), ETAKNOTHO(ND)

      DIMENSION XJCINT(ND,NF), FWTEST(NF)
      DIMENSION HTOTL(ND), XJTOTL(ND), XKTOTL(ND), XNTOTL(ND)
      DIMENSION ARAD(ND), ACONT(ND), ATHOM(ND)
      DIMENSION XJFEMEAN(ND, LASTFE), XJLMEAN(NDDIM,MAXLIN)

      DIMENSION OPAK(ND), ETAK(ND)
      DIMENSION OPAKNOTH(ND), ETAKNOTH(ND)
      DIMENSION OPAFEI(ND,LASTFE), ETAFEI(ND,LASTFE)
      DIMENSION WFELOW(ND,LASTFE), WFENUP(ND,LASTFE)
      DIMENSION DSDFELOW(ND,LASTFE), DSDFENUP(ND,LASTFE)
      DIMENSION DSDETA(ND), DSDOPA(ND)
      DIMENSION IFENUP(LASTFE), IFELOW(LASTFE)
      DIMENSION WEIGHT(N), POPNUM(ND,N), DJDS(ND), DJDO(ND), FTCOLI(ND)
      DIMENSION SC(ND), SCO(ND)

      DIMENSION DJDS_OLD(ND), DJDO_OLD(ND)
      DIMENSION WJC(ND,NF), DSDSC(ND), DJDSC(ND), DJDSCOLD(ND)
      DIMENSION WJC_MIN(ND,NF)
      DIMENSION ETANOTH(ND), OPA(ND), THOMSON(ND), T(ND)
      DIMENSION OPAO(ND), THOMSONO(ND)      
      DIMENSION RADIUS2(ND)
      DIMENSION FWEIGHT(NF), EMCOLI(NF),  FTFE(ND,LASTFE)

C***  Unsoeld-Lucy
      DIMENSION OPASMEAN(ND), SMEAN(ND), QFJMEAN(ND)
      DIMENSION OPAJMEAN(ND), OPAPMEAN(ND)
      DIMENSION QOPAHMEAN(ND), HMEAN(ND)
      DIMENSION OPASMEANTC(ND), OPAJMEANTC(ND) 

      DIMENSION QLFOLD(ND), QLHOLD(ND), OPAKHOLD(ND)
      DIMENSION EDDIFO(ND), EDDIGO(ND)
      LOGICAL BCOLIRAY, BNEWINTERVAL, BNEWINTERVALOLD, BEMPTY
      DIMENSION SUMJ(ND), SUMJW(ND), SUMDJDSC(ND), SUMDJDSCW(ND)      
      CHARACTER(LEN=8) OPC
      DIMENSION FERATLU(LASTFE,ND), FERATUL(LASTFE,ND)
      DIMENSION ELEVEL(LASTFE)
C***  ALI Acceleration terms
      DIMENSION OPAROSS(ND), OPALAMBDAMEAN(ND)

      LOGICAL :: BWCELAB, BCUTOFF, BCUTOFF2
      REAL :: HTOTND, XHID, EDDIHINOLD, HTOTNDCOR

      SAVE 

cccccc!!!!!!!!!! LTE test option: set radiation field to Blackbody!
ccc   Consistent change AT TWO OTHER PLACES in this subroutine !
ccc      dimension xjlsave(100), xjloldsave(100)
ccc      DO L=1, ND
ccc      xjlsave(l) = xjl(l)      
ccc      xjloldsave(l) = xjlold(l)      
ccc      XJL(L)    = amax1 ( 1.e-100, BNUE(XLAMK   , T(L)) * RADIUS2(L))
ccc      if (k .gt. 0) XJLOLD(L) = 
ccc     >      amax1 ( 1.e-100, BNUE(XLAMKOLD, T(L)) * RADIUS2(L))
ccc      ENDDO
cccccccccccccccccccccccccccccccccccc

C***  FIXED Switch for the 'elaborated' Calculation of WJC
      BWCELAB = .TRUE.

C***  Application of GAMMACOLI to DJDS
C***  --------------------------------
      DJDS = MIN(DJDS, WJMAX)
      DJDS = MAX(DJDS, 0.)

C***  DJDS -> TAU
      DJDS = -LOG(1.-DJDS)
C***  DJDS -> TAU/GAMMA
      DJDS = DJDS / GAMMACOLI
C***  Back Transformation
      DJDS = 1.-EXP(-DJDS)

C********************************************************************
C***  Lines: Integration of the scattering integrals XJLMEAN
C********************************************************************

C***  ADDING THE NEW XJL TO THE MEAN INTENSITY XJLMEAN (NORMAL LINES)
      DO NL=1, MAXLIN
        IF (LIND(NL) .EQ. 0) CYCLE
        DO L=1, ND
           XJLMEAN(L,NL) = XJLMEAN(L,NL) + XJL(L)*PWEIGHT(NL)
        ENDDO
      ENDDO

C***  ADDING THE NEW XJL TO THE MEAN INTENSITY XJFEMEAN (IRON LINES)
C***     AND CALCULATION OF DIAGONAL WEIGHTS

C***  DERIV. OF S WITHOUT THOMSON TERMS
      DSDETA = 0.
      DSDOPA = 0.
      WHERE (OPAKNOTH .GT. 0.)
         DSDETA = 1. / OPAKNOTH
         DSDOPA = -ETAKNOTH / (OPAKNOTH * OPAKNOTH)
      END WHERE

      XLAMKCM = XLAMK * 1.E-8
      WAVK = 1. / XLAMKCM
      WAVK3 = WAVK * WAVK * WAVK  
      DO INDACT=1, MAXFEACT
         IND = INDFEACT(INDACT)
         SIGMA = SIGMAACT(INDACT)

         NUP = IFENUP(IND)
         LOW = IFELOW(IND)
         IF (LOW .EQ. NUP) CYCLE
         WLU = WEIGHT(LOW)/WEIGHT(NUP)
         WAV0 = ELEVEL(NUP) - ELEVEL(LOW)
         XLAM0CM = 1. / WAV0

         DO L=1, ND
C***       INTEGRATION OF XJFEMEAN
            XJFEMEAN(L,IND) =XJFEMEAN(L,IND) + XJL(L) * FWEIGHTL * SIGMA

C***       CORRECTION FOR ENERGY EQUATION
            FTFE(L,IND) = FTFE(L,IND)
     >            + (ETAFEI(L,IND) - OPAFEI(L,IND)*XJL(L)/RADIUS2(L))
     >              * FWEIGHTL

C***        Radiative rates of iron superlines
            XJLPLAIN = XJL(L) / RADIUS2(L)

C***        optional de-activation of exp term - see comment in CMFFEOP
            IF (IVERS_FE_EXPFAC .EQ. 0) THEN
               WLUEXP = WLU 
            ELSEIF (IVERS_FE_EXPFAC .EQ. 1) THEN
               WLUEXP = WLU * EXP(C1*(WAV0-WAVK)/MAX(T(L),TEFF))
            ELSEIF (IVERS_FE_EXPFAC .EQ. 2) THEN
               WLUEXP = WLU * EXP(C1*(WAV0-WAVK)/T(L))
            ELSE
              STOP '*** FATAL INTERNAL ERROR 1 in subr. FREQUINT'
            ENDIF

C*** This version XLAMKCM --> XLAM0CM in the next 2 eqns.
            DRATLU =  
     >          PI8/C2 * XJLPLAIN * XLAM0CM * SIGMA * FWEIGHTL 
            DRATUL =  
     >          WLUEXP * PI8/C2 * 
     >          (C2 * WAVK3 + XJLPLAIN) * XLAM0CM * SIGMA * FWEIGHTL

            FERATLU(IND,L) = FERATLU(IND,L) + DRATLU
            FERATUL(IND,L) = FERATUL(IND,L) + DRATUL

C***       CALCULATION OF ACTUAL DERIVATIVES OF S

            POPNUP = POPNUM(L,NUP)
            POPLOW = POPNUM(L,LOW)

C***       DERIV. OF ETAI
            IF (OPAFEI(L,IND) .GT. 0. .AND. POPNUP .NE. .0) THEN
               DETAIDNUP = ETAFEI(L,IND) / POPNUP
            ELSE
               DETAIDNUP = 0.
            ENDIF

C***       DERIV. OF OPAI
            IF (OPAFEI(L,IND) .GT. 0. .AND. POPLOW .NE. .0) THEN
               OPANULL  = OPAFEI(L,IND) / (POPLOW - WLUEXP*POPNUP)
            ELSE
               OPANULL = 0.
            ENDIF
            DOPAIDLOW = OPANULL
            DOPAIDNUP = -OPANULL*WLUEXP

C***       DERIV. OF S WITH RESPECT TO LOW AND NUP
            DSDFELOW(L,IND) = DSDOPA(L)*DOPAIDLOW
            DSDFENUP(L,IND) = DSDOPA(L)*DOPAIDNUP + DSDETA(L)*DETAIDNUP

C***       ADD DERIV. AT THE NEW FREQUENCY TO DERIV. OF J
            DJDLOW = DJDS(L) * DSDFELOW(L,IND)
            DJDNUP = DJDS(L) * DSDFENUP(L,IND)

C***       INTEGRATION OF DIAGONAL WEIGHTS
            WFELOW(L,IND) = WFELOW(L, IND) + DJDLOW * FWEIGHTL * SIGMA
            WFENUP(L,IND) = WFENUP(L, IND) + DJDNUP * FWEIGHTL * SIGMA

         ENDDO
      ENDDO


C********************************************************************
C***  Averaging the intensity from fine to coarse frequency grid
C********************************************************************

C***  Note: The current (fine) frequency point bears the index K; the 
C***        integration step trails and hits the frequency index K-1 = OLD
C***  KOLD defines the current coarse interval (KCONT, KCONT+1)

C***  Check whether the current frequency point is the first which 
C***  has entered a new coarse interval
      IF (K .EQ. 0) THEN 
         BNEWINTERVALOLD = .FALSE.
         KCONT = 1
         BEMPTY = .FALSE.
      ELSE
         BNEWINTERVALOLD =  BNEWINTERVAL
      ENDIF
      BNEWINTERVAL = XLAMK .GE. XLAMBDA(KCONT+1) 

C***  Frequency Integration Weight for old finegrid point
C***  Integration weight nue*dnue for test function 
C***   Note: First and last interval get additional terms - see below!
      IF (K .EQ. 0) THEN

         IF (BWCELAB) THEN
            FWTEST = 1.
         ELSE
            FWTEST = 0.
         ENDIF
         DO L=1, ND
            SUMJ     (L) = .0
            SUMJW    (L) = .0
            SUMDJDSC (L) = .0
            SUMDJDSCW(L) = .0
         ENDDO
         SUMEMF = 0.
         SUMEMFW= 0.
         SUMFW  =0.
         SUMFWW =0.
         XNUEK = CLIGHT / XLAMK

         IF (OPC .EQ. 'DIAGMIN') THEN
           DO KK=1, NF
             DO L=1, ND
               WJC_MIN(L,KK) = 1.
             ENDDO
           ENDDO
         ENDIF

         DLEFT  = 0.
         DRIGHT = 0.

cccccccccccccccccccccc LTE-Test Facility cccccccccccccccc
ccc      DO L=1, ND
ccc         xjl(l) = xjlsave(l) 
ccc         xjlold(l) = xjloldsave(l) 
ccc      ENDDO
cccccccccccccccccccccccccccccccccccccccccccccccccccccccc

         RETURN

      ELSE IF (K .EQ. 1) THEN
         XNUEKOLD  = XNUEK
         XNUEK     = CLIGHT / XLAMK
         XNUE2 = CLIGHT / XLAMBDA(1)
         XNUE3 = CLIGHT / XLAMBDA(2)

         DLEFT = (XNUEKOLD - XNUE2) / (XNUE3 - XNUE2)
         DRIGHT = 1. - DLEFT
         FWKOLD      = 0.5 * (XNUEKOLD - XNUEK)
         FWNUEKOLD   = (2. * XNUEKOLD + XNUEK) *
     >                 (XNUEKOLD - XNUEK) / 6.        

      ELSE
         XNUEKOLDOLD = XNUEKOLD
         XNUEKOLD  = XNUEK
         XNUEK     = CLIGHT / XLAMK
         IF (BNEWINTERVALOLD) THEN
            XNUELEFT = CLIGHT / XLAMBDA(KCONT)
         ELSE
            XNUELEFT = XNUEKOLDOLD
         ENDIF
         IF (BNEWINTERVAL) THEN
            XNUERIGHT = CLIGHT / XLAMBDA(KCONT+1)
         ELSE
            XNUERIGHT = XNUEK
         ENDIF

         DLEFT = (XNUEKOLD - XNUE2) / (XNUE3 - XNUE2)
         DRIGHT = 1. - DLEFT
         FWKOLD      = 0.5 * (XNUELEFT-XNUERIGHT)
         FWNUEKOLD   = (XNUELEFT + XNUEKOLD + XNUERIGHT) *
     >                 (XNUELEFT-XNUERIGHT) / 6.        

C***       Test Output for the current Problems with Compaq Fortran
C           IF ((K .GT. 60900) .AND. (K .LT. 61100)) THEN
C              WRITE(0,*) FWNUEKOLD, XNUELEFT, XNUERIGHT, XNUEKOLD
C           ENDIF

      ENDIF

C***  Determine minimum of 1.-EXP(-TAU), if requested
      IF (OPC .EQ. 'DIAGMIN') THEN
        DLEFT_MIN  = (XNUEKOLD - XNUE2) / (XNUE3 - XNUE2)
        DRIGHT_MIN = 1. - DLEFT_MIN
C        IF (DLEFT_MIN .GT. DRIGHT_MIN) THEN
        IF (DLEFT_MIN .LT. DRIGHT_MIN) THEN
          DO L=1, ND
            WJC_MIN(L,KCONT) = AMIN1(WJC_MIN(L,KCONT), DJDS_OLD(L))
          ENDDO
        ELSE
          DO L=1, ND
            WJC_MIN(L,KCONT+1) = AMIN1(WJC_MIN(L,KCONT+1), DJDS_OLD(L))
          ENDDO
        ENDIF
c        write (0,'(a,4(1x,e12.5), 2(1x,f7.3), 3(1x,e12.5))')
c     >    '===', XLAMKOLD, xnuekold, xnue2, xnue3, dleft, dright, 
c     >    WJC_MIN(50,KCONT), WJC_MIN(50,KCONT+1), DJDS_OLD(50)
      ENDIF

C***  Here is the entry for an empty interval:
   10 CONTINUE

      IF (KCONT .LT. 1) THEN
         WRITE (0,*) '== WARNING from FREQUINT: KCONT .LT. 1'
         KCONT = 1
      ENDIF
      IF (KCONT .GT. NF-1) THEN
         WRITE (0,*) '== WARNING from FREQUINT: KCONT .GT. NF-1'
         KCONT = NF-1
      ENDIF

      IF (.NOT. BEMPTY) THEN
         DO L=1, ND
C***        Average the Jnue's from the fine frequency step for the coarse
C***        continuum frequency mesh used by STEAL
            SUMJ (L)  = SUMJ (L) + XJLOLD(L) * FWKOLD
C***        The same is done after multiplication with the testfunction NUE
            SUMJW(L)  = SUMJW(L) + XJLOLD(L) * FWNUEKOLD

C***        Old Dioagonal Weight for Continuum is saved
            DJDSCOLD(L) = DJDSC(L)

C***        Dioagonal Weight for Continuum is calculated:
            DJDSC(L) = 0.

C***        Deriv. of J with resp. to Sc at the current Frequency
            IF (OPAKNOTH(L) .GT. 0.) THEN
              OPAC  = OPA(L)*(1. - THOMSON(L))
              DSDSC(L) = OPAC/OPAKNOTH(L)
            ELSE
              DSDSC(L) = 0.
            ENDIF

            IF (OPC .EQ. 'DIAGTAU') THEN
C***          In this version DIAGTAU, the Scharmer-weight factor DJDSC is first 
C***            convertet into a kinf of TAU before averaged over the fine grid
C***            The back conversion is made in SUBR. FREQUNORM
              DJFO = DJDS(L)*DSDSC(L)
C!!!              Operator without Opacity Correction
C!!!              DJFO = DJDS(L)

C***          Special WCHARM plot
              IF (L .EQ. LPLOT_WCHARM) THEN
                WRITE (105, '(2(E20.10,1X))') XLAMK, DJFO
              ENDIF
              IF (DJFO .GE. 1.) THEN
                DJDSC(L) = 0.
              ELSE
                DJDSC(L) = -ALOG(1. - DJFO)
              ENDIF

C***        In this version the averaging is made directly with DJDSC
C***        WARNING: DERIVATIVE TERM UNCLEAR !!!
            ELSE IF (OPC .EQ. 'DIAG') THEN
C***          Scale DJDSC with ratio of Source funktions
              DJDSC(L) = DJDS(L)*DSDSC(L)
C!              IF (DJDSC(L) .LT. 0.) DJDSC(L)=0.

            ELSE IF (OPC .EQ. 'DIAGNC') THEN
C***          B) Add Deriv. at the actual Frequency
              DJDSC(L) = DJDSC(L) + DJDS(L)*DSDSC(L)
              IF (DJDSC(L) .LT. 0.) DJDSC(L)=0.
            ELSE IF (OPC .EQ. 'DIAGMIN') THEN
              DJFO = DJDS(L)
C***          Special WCHARM plot
              IF (L .EQ. LPLOT_WCHARM) THEN
                WRITE (105, '(2(E20.10,1X))') XLAMK, DJFO
              ENDIF
            ELSE IF (OPC .EQ. 'NONE') THEN
              DJDSC(L) = 0.
            ELSE
              WRITE (0,*) 'Unknown Version OPC = ', OPC
            ENDIF

C***        DJDSCOLD is the diagonal element for the continuum 
C***        for optional use in STEAL
            IF (BWCELAB) THEN
C***        Average DJDSCOLD from the fine frequency step for the coarse
C***        continuum frequency mesh used by STEAL
               SUMDJDSC (L)  = SUMDJDSC (L) + DJDSCOLD(L) * FWKOLD
C***        Same is done after multiplication with the testfunction NUE
               SUMDJDSCW(L)  = SUMDJDSCW(L) + DJDSCOLD(L) * FWNUEKOLD
            ELSE
C***        Integration of WJC ('non elaborated', trapezregel)
               WJC(L, KCONT)  = WJC(L, KCONT)
     >                         + DRIGHT * DJDSCOLD(L) * FWKOLD
               WJC(L, KCONT+1)= WJC(L, KCONT+1)
     >                         + DLEFT  * DJDSCOLD(L) * FWKOLD
            ENDIF

         ENDDO
C***     Integration of Emergent Flux
         SUMEMF  = SUMEMF + XHLOLD(1) * FWKOLD
         SUMEMFW = SUMEMFW+ XHLOLD(1) * FWNUEKOLD
         SUMFW   = SUMFW  + FWKOLD
         SUMFWW  = SUMFWW + FWNUEKOLD
C***     Frequency Weight for Normalization
         IF (.NOT. BWCELAB) THEN
            FWTEST(KCONT)  = FWTEST(KCONT)  + DRIGHT * FWKOLD
            FWTEST(KCONT+1)= FWTEST(KCONT+1)+ DLEFT  * FWKOLD
         ENDIF

      ENDIF

      IF (BNEWINTERVAL) THEN 
C***     Now the (half-)interval (KCONT, KCONT+1) is completed:
C***     XLAMKOLD to XLAMBDA(KCONT+1), i.e. XNUEKOLD ... XNUE2
C***     and subtract (!) the excess part of the last integration step
         XNUE1 = CLIGHT / XLAMBDA(KCONT  )
         XNUE2 = CLIGHT / XLAMBDA(KCONT+1)
         XNUE3 = CLIGHT / XLAMBDA(KCONT+2)

         QSPECIAL = (XNUEKOLD - XNUE2) / (XNUEKOLD - XNUEK)
         QSPECIAL2 = 1. - QSPECIAL

         IF (BEMPTY) THEN
            XNUELEFT = CLIGHT / XLAMBDA(KCONT)
         ELSE
            XNUELEFT = XNUEKOLD
         ENDIF

C***     Complete the integration with the part of the step from
C***     MAX(XLAMKOLD,XLAMBDA(KCONT) ... XLAMBDA(KCONT+1),
C***     i.e. XNUELEFT ... XNUE2
         FWSPECIAL    = 0.5 * (XNUELEFT - XNUE2)
         FWNUESPECIAL = (XNUELEFT - XNUE2) * (2.*XNUE2 + XNUELEFT) / 6.

         DO L=1, ND
            XJSPECIAL = QSPECIAL2 * XJLOLD(L) + QSPECIAL * XJL(L)
            SUMJ (L)  = SUMJ (L) + XJSPECIAL * FWSPECIAL
            SUMJW(L)  = SUMJW(L) + XJSPECIAL * FWNUESPECIAL
            DJDSCSPECIAL = QSPECIAL2 * DJDSCOLD(L) 
     >                    + QSPECIAL * DJDSC(L)

            IF (BWCELAB) THEN
               SUMDJDSC (L)  = SUMDJDSC (L) + DJDSCSPECIAL * FWSPECIAL
               SUMDJDSCW(L)  = SUMDJDSCW(L)
     >                         + DJDSCSPECIAL * FWNUESPECIAL
            ELSE
C***        Integration with DLEFT=0., DRIGHT=1.
               WJC(L, KCONT+1)= WJC(L, KCONT+1)
     >                          + DJDSCSPECIAL * FWSPECIAL
            ENDIF

         ENDDO

         XHSPECIAL = QSPECIAL2 * XHLOLD(1) + QSPECIAL * XHL(1)
         SUMEMF    = SUMEMF + XHSPECIAL * FWSPECIAL
         SUMEMFW   = SUMEMFW+ XHSPECIAL * FWNUESPECIAL
         SUMFW     = SUMFW  + FWSPECIAL
         SUMFWW    = SUMFWW + FWNUESPECIAL

         IF (.NOT. BWCELAB) THEN
            FWTEST(KCONT+1) = FWTEST(KCONT+1) + FWSPECIAL
         ENDIF

C***     Generate intergration weights at coarse frequency points
C***     which incorporate the test function, i.e. nue*dnue
         FWNUE1 = (2.*XNUE1 + XNUE2) * (XNUE1-XNUE2) / 6.
         FWNUE2 = (2.*XNUE2 + XNUE1) * (XNUE1-XNUE2) / 6.

         FW1 = 0.5 * (XNUE1 - XNUE2)         
         FW2 = FW1
         IF (KCONT .EQ. 1) THEN 
            FP = 1.
            DO L=1, ND
               XJCINT(L,1) = 0.
            ENDDO 
            EMCOLI(1) = 0.
         ELSE
            XNUE0 = CLIGHT / XLAMBDA(KCONT-1) 
            FP = (XNUE1 - XNUE2) / (XNUE0 - XNUE2) 
         ENDIF
         FQ = 1. - FP

         DO L=1, ND
C***       Note: At the left point, J has already the contribution 
C***                from the right end of the last interval
C***       Note: The contributions are now weighted 
            XJCINT(L,KCONT+1) =
     >               ( FWNUE1 * SUMJ(L) - FW1 * SUMJW(L) ) / 
     >               ( FW2 * FWNUE1 - FW1 * FWNUE2)
            XJCINT(L,KCONT) = FQ * XJCINT(L,KCONT) + FP * 
     >               ( FWNUE2 * SUMJ(L) - FW2 * SUMJW(L))/ 
     >               (FW1 * FWNUE2 - FW2 * FWNUE1)

C***       **** ELABORATED Integration of WJC ******************
            IF (BWCELAB) THEN
               WJC(L,KCONT+1) =
     >               ( FWNUE1 * SUMDJDSC(L) - FW1 * SUMDJDSCW(L) ) / 
     >               ( FW2 * FWNUE1 - FW1 * FWNUE2)
               WJC(L,KCONT) = FQ * WJC(L,KCONT) + FP * 
     >               ( FWNUE2 * SUMDJDSC(L) - FW2 * SUMDJDSCW(L))/ 
     >               (FW1 * FWNUE2 - FW2 * FWNUE1)
C***       ******************************************************
            ENDIF

         ENDDO
         EMCOLI(KCONT+1) = 4.*
     >               ( FWNUE1 * SUMEMF - FW1 * SUMEMFW ) /
     >               ( FW2 * FWNUE1 - FW1 * FWNUE2)
         EMCOLI(KCONT)   = FQ * EMCOLI(KCONT) + 4.*FP *
     >               ( FWNUE2 * SUMEMF - FW2 * SUMEMFW)/
     >               (FW1 * FWNUE2 - FW2 * FWNUE1)

C***     Now increment KCONT and check for empty coarse interval 
C***     (i.e. which contains no fine frequency point).
C***     In that case the fine Interval ends at XLAMBDA(KCONT+1)
         KCONT = KCONT + 1
         IF (XLAMK .GT. XLAMBDA(KCONT+1)) THEN
            BEMPTY = .TRUE.
            XNUERIGHT = CLIGHT / XLAMBDA(KCONT+1)
         ELSE
            BEMPTY = .FALSE.
            XNUERIGHT = XNUEK
         ENDIF

C***     Start the next integration (KCONT, KCONT+1) with the old
C***     (half-)interval:
C***     XLAMBDA(KCONT) to MIN(XLAMK,XLAMBDA(KCONT+1),
C***     i.e. XNUE2 ... XNUERIGHT
         FWSPECIAL    = 0.5 * (XNUE2 - XNUERIGHT)
         FWNUESPECIAL = (XNUE2 - XNUERIGHT)
     >                  * (2.*XNUE2 + XNUERIGHT) / 6.

         DO L=1, ND
            XJSPECIAL = QSPECIAL2 * XJLOLD(L) + QSPECIAL * XJL(L)
            SUMJ (L)  = XJSPECIAL * FWSPECIAL
            SUMJW(L)  = XJSPECIAL * FWNUESPECIAL
            DJDSCSPECIAL = QSPECIAL2 * DJDSCOLD(L) 
     >                    + QSPECIAL * DJDSC(L)

            IF (BWCELAB) THEN
               SUMDJDSC (L)  = DJDSCSPECIAL * FWSPECIAL
               SUMDJDSCW(L)  = DJDSCSPECIAL * FWNUESPECIAL
            ELSE
C***           Integration with DLEFT=1., DRIGHT=0.
               WJC(L, KCONT) = WJC(L, KCONT) + DJDSCSPECIAL * FWSPECIAL
            ENDIF

         ENDDO
         XHSPECIAL = QSPECIAL2 * XHLOLD(1) + QSPECIAL * XHL(1)
         SUMEMF  = XHSPECIAL * FWSPECIAL
         SUMEMFW = XHSPECIAL * FWNUESPECIAL
         SUMFW   = FWSPECIAL
         SUMFWW  = FWNUESPECIAL
C***     Normalization
         IF (.NOT. BWCELAB) THEN
            FWTEST(KCONT) = FWTEST(KCONT) + FWSPECIAL
         ENDIF

C***     In case of empty Interval go again to 'New Interval':
         IF (BEMPTY) GOTO 10

      ENDIF


C**********************************************************************
C***  INTEGRATE TOTAL FLUX AND RADIATIVE FORCE
C**********************************************************************
      IF (K .EQ. 0) THEN
         FWMID = 0.5 * (XNUEKOLD - XNUEK)
      ELSE
         FWMID = 0.5 * (XNUEKOLDOLD - XNUEK)
      ENDIF

      DO L=1, ND
         IF (XJLOLD(L) .GT. 1.E-15)
     >   XJTOTL(L) = XJTOTL(L) + XJLOLD(L)*FWMID
         IF (BCOLIRAY) THEN
              XKTOTL(L) = XKTOTL(L) + XKLOLD(L)*FWMID
         ELSE
              XKTOTL(L) = XKTOTL(L) + XJLOLD(L)*EDDIFO(L)*FWMID
         ENDIF

         IF (L == ND) THEN
           IF (K > 0) THEN
C***         Calculate Eddington Flux at inner boundary 
C***          via frequency integration of
C***         HTOTND: H_nu,ND = H_nu,diff + B_nu/2 - h_nu,in * J_nu,ND
             HTOTND = HTOTND + ( XHID  + 
     >          BNUE(XLAMKOLD,T(ND))/2.- EDDIHINOLD*XJLOLD(L)) * FWMID
C***         HTOTNDCOR: frequency integration without the diffusion term H_nu,diff
C***              (i.e. of the terms accounting for LTE/diffusion deviations)
             HTOTNDCOR = HTOTNDCOR + 
     >         (BNUE(XLAMKOLD,T(ND))/2.- EDDIHINOLD*XJLOLD(L)) * FWMID
           ENDIF
           
           CYCLE
         ENDIF

C***    Radiation Pressure
         ARAD(L) = ARAD(L)
     >               + 0.5*(OPAKOLD(L)+OPAKOLD(L+1))*XHLOLD(L)*FWMID
C+C***    Contributions of the Continuum
         OPTH = 0.5 * (OPAO(L)*THOMSONO(L) + OPAO(L+1)*THOMSONO(L+1))
         OPAC = 0.5 * (OPAO(L) + OPAO(L+1))
         ACONT(L) = ACONT(L) + OPAC*XHLOLD(L)*FWMID
         ATHOM(L) = ATHOM(L) + OPTH*XHLOLD(L)*FWMID

         HTOTL(L) = HTOTL(L) + XHLOLD(L)*FWMID
         IF (BCOLIRAY) THEN
              XNTOTL(L) = XNTOTL(L) + XNLOLD(L)*FWMID
         ELSE
              XNTOTL(L) = XNTOTL(L) + XHLOLD(L)*EDDIGO(L)*FWMID
         ENDIF

      END DO

C***  Calculation of OPAROSS at L=ND
      DBDT   = DBNUEDT(XLAMKOLD, T(ND))
      DBDTINT    = DBDTINT    + DBDT*FWMID
      DBDTOPAINT = DBDTOPAINT + DBDT*FWMID/OPAKOLD(ND)
C- ***********************


C***  Unsoeld-Lucy Procedure for COLI
C***  Note: Integrands are only incremented here for the OLD frequency point!
C***        Final calculation is performed in SUBR. FREQUNORM

C***  NOTE: cut off long wavelength -- de-activated 22-Apr-2002 17:00:37 wrh
ccc      IF (XLAMKOLD .GT. 20000.) GOTO 70

C***  Cut-Off of lauge optical depths governed by UNLU_TAUMAX und UNLU_TAUMAX2
C***     is switched off if ZERO ( = infinity)
      BCUTOFF  = UNLU_TAUMAX  .NE. 0. 
      BCUTOFF2 = UNLU_TAUMAX2 .NE. 0. 

C***  Special Eddi for outer boundary
      EDDIHOUTJMEAN  = EDDIHOUTJMEAN + EDDIHOUTOLD * XJLOLD(1) * FWMID 


C***  Begin of Loop ******************************************
      DO L=1, ND

C***   Rosseland Opacity
       DBDT       = DBNUEDT(XLAMKOLD, T(L))
       OPAROSS(L) = OPAROSS(L) + DBDT * FWMID / OPAKOLD(L)

       IF (BCUTOFF2 .AND.
     >    OPAKNOTHO(L) * (RADIUS(L)-1.) .GT. UNLU_TAUMAX2) GOTO 45
          IF (XJLOLD(L) .LT. 1.E-50) GOTO 45
       FTCOLI(L) = FTCOLI(L)
     >  + (ETAKNOTHO(L) - OPAKNOTHO(L)*XJLOLD(L)/RADIUS2(L)) * FWMID
 45    continue

C***  CUT OFF VERY LARGE OPTICAL DEPTHS ("detailed balance") 
       IF (BCUTOFF .AND.
     >    (OPAKNOTHO(L) * (RADIUS(L)-1.) .GT. UNLU_TAUMAX)) GOTO 44
C***  Special treatment for negative 'true' Opacity
C***   (Ist dieser Fallback schlau oder eher problematisch?)
c            => greift meist nur aussen!
       IF (OPAKNOTHO(L) .LE. 1.E-30) THEN
           OPAJMEAN(L) = OPAJMEAN(L) + OPAKOLD(L)*XJLOLD(L)      * FWMID
           OPASMEAN(L) = OPASMEAN(L) + ETAKOLD(L)                * FWMID
           OPAPMEAN(L) = OPAPMEAN(L) 
     >                          + OPAKOLD(L)*BNUE(XLAMKOLD,T(L)) * FWMID
           SMEAN(L)    = SMEAN(L)    + ETAKOLD(L)/OPAKOLD(L)     * FWMID
       ELSE
          OPAJMEAN(L) = OPAJMEAN(L) + OPAKNOTHO(L)*XJLOLD(L)    * FWMID
          OPASMEAN(L) = OPASMEAN(L) + ETAKNOTHO(L)              * FWMID
          OPAPMEAN(L) = OPAPMEAN(L) 
     >                       + OPAKNOTHO(L)*BNUE(XLAMKOLD,T(L)) * FWMID
          SMEAN(L)    = SMEAN(L)    + ETAKNOTHO(L)/OPAKNOTHO(L) * FWMID
       ENDIF

C***    Option for ALI-like amplification of the temperature correction 
C***      of the "first term" in the Unsoeld-Lucy procedure
C***      De-activated if GAMMAT=0.
        IF (GAMMAT .NE. .0) THEN
c            IF (L .EQ. 1) THEN
c               DR = 0.
c            ELSE IF (L .EQ. ND) THEN
c               DR = 0.
c            ELSE
c               DR = AMIN1( RADIUS(L) - RADIUS(ND), 
c     >                     RADIUS(1) - RADIUS(L))
c            ENDIF

C***  A more sophisticated estimate of TAU = Delta Tau to nearest border:
C***  Assume that the opacity dilutes with r**(-2). Integration yields:
            DROUT = RADIUS(L) * ( 1. - RADIUS(L) / RADIUS(1))
            DRIN  = RADIUS(L) * ( RADIUS(L) - 1.)
            DR = AMIN1 ( DROUT, DRIN )

            TAU = OPAKOLD(L) * DR / GAMMAT
            IF (TAU .GT. 1.E-6) THEN
               ALONUE = 1. - (1. - EXP(-TAU)) / TAU
            ELSE
               ALONUE = 0.
            ENDIF

            OPALAMBDAMEAN(L) = OPALAMBDAMEAN(L) + 
     >                        ETAKNOTHO(L) * ALONUE * FWMID
        ENDIF

   44  CONTINUE

C***    The "TC" quantities are for TEMPCORR (Unsoeld-Lucy)
        OPASMEANTC(L) = OPASMEAN(L)
        OPAJMEANTC(L) = OPAJMEAN(L)

        QFJMEAN(L)  = QFJMEAN(L)  + 
     >       QLFOLD(L)*EDDIFO(L)*XJLOLD(L) * FWMID


C***    Flux quantities defined at interstices - skip for L=ND !
        IF (L .EQ. ND) CYCLE

CCC     "EFFECTIVE OPACITY" ACCOUNTING FOR THE INFLUENCE OF  
CCC        THOMSON SCATTERING ON THE THEMALIZATION DEPTH 
C***    ADDITIONALLY, CUT OFF VERY LARGE OPTICAL DEPTHS ("detailed balance") 
        OPATRUE    =  0.5 * (OPAKNOTHO(L) + OPAKNOTHO(L+1)) 
        IF (BCUTOFF) THEN
           IF (OPATRUE * (RADIUS(L)-1.) .GT. UNLU_TAUMAX2) GOTO 46
        ELSE
           IF (XHLOLD(L) .LT. 1.E-50) GOTO 46
        ENDIF

ccc  merkwuerdige Artistik, war wohl auch anders gemeint. 
ccc testweise abgeschaltet wrh 28-Nov-2003 14:35:58
c        OPASCATTER = OPAKHOLD(L) - OPATRUE 
c        OPAOPA     = OPATRUE * OPASCATTER
c        IF (OPAOPA .GT. .0) THEN
c           OPAEFFECTIVE = SQRT (OPAOPA)
c           QOPAHMEAN(L) = QOPAHMEAN(L) + 
c     >          QLHOLD(L) * OPAEFFECTIVE * XHLOLD(L) * FWMID
c        ENDIF
ccc dafuer angeschaltet:

ccc    noch eine Variante: nach Ansicht von Goetz kommt hier: 
ccc    in der 1. Momentengleichung, aus der in TEMPCORR der "2. Term" 
ccc    ausgerechnet wird, die *volle* Opazitaet (incl. Thomson) hin:
ccc        wrh  7-Dec-2005 11:56:43
ccc    Also seine Empfehlung: Ersetze das folgende Statement
        QOPAHMEAN(L) = QOPAHMEAN(L) +
     >       QLHOLD(L) * OPATRUE * XHLOLD(L) * FWMID
ccc    durch:
cc        OPAFULL    =  0.5 * (OPAKOLD(L) + OPAKOLD(L+1)) 
cc        QOPAHMEAN(L) = QOPAHMEAN(L) +
cc     >       QLHOLD(L) * OPAFULL * XHLOLD(L) * FWMID

        HMEAN(L)= HMEAN(L) + XHLOLD(L) * FWMID

 46     CONTINUE

      END DO
C***  End of L-Loop

C***  Jump label: integrations skipped for very long wavelengths 
   70 CONTINUE



cccccccccccccccccccccc LTE-Test Facility cccccccccccccccc
ccc      DO L=1, ND
ccc         xjl(l) = xjlsave(l) 
ccc         xjlold(l) = xjloldsave(l) 
ccc      ENDDO
cccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      RETURN
      END
