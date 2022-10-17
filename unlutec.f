      SUBROUTINE UNLUTEC (L, NRANK, N, OPAC, ETAC, XLAMBDA, TL, TEFF,  
     $   LASTIND, INDLOW, INDNUP, NDIM, OPAL,
     $   NF, ND, XJCAPP, SCNEW, FWEIGHT,
     $   LASTKON, DRJLW, DRLJW, IONGRND,
     $   KODRLOW, LASTKDR, FEDDI, QEDDI, HEDDI, 
     >   LEVEL, RNEL, KONTLOW, SIGMAKI,
     $   SIGMAFF, MAXION, NCHRG, RSTAR, EN, ENTOTL, LSTEQ, 
     >   ELEVEL, RRATE, OPAC1, XJCLP1, RADIUS, OPATHOM, DUNLU, 
     >   HTOT, HTOTC, HTOTL, SLNEW, VDOP,
     >   DTLOCAL, DTLOCC, DTLOCL, DTLOCD, DTINT, DTRMAX, TMIN, NEGT, 
     >   DUNLUR)

C***********************************************************************
C***  
C***  Unsoeld-Lucy Temperature Correction Procedure 
C***  
C***  The local equation version accounts for continua + diel.rec. + lines
C***     Note concerning the line contribution:
C***     Instead of writing: ETAL - OPAL * XJL, we make use of the
C***     radiative rate coefficients RRATE.
C***     This has the numerical advantage, that it implies the
C***     net-radiative-bracket formulation (in case of non-zero line cores)
C***     the equations are in cgs units.
C***  CALLED FROM: SUBROUTINE LINPOP
C***
C***  List of terms under consideration:
C***  Term I (local Term):
C***    Continuum with approximate XJC (XJCAPP)   ==> FTCONT
C***    Dielectronic Recombination from the rates ==> FTDR
C***    Lines form the rates                      ==> FTLINES
C***    Note that Kappa_S is 
C***      corrected for the lines                 ==> OPASMEAN2
C***  Term II (Integral Term):
C***    HNUE from differentiation of XJCAPP
C***    Flux of the lines added                   ==> DHINT
C***  Term III (Boundary Term):
C***    HNUE from XJCAPP and HEDDI
C***    Flux of the lines added:                  ==> DHRMAX
C***    
C***    
C***    
C***********************************************************************

      DIMENSION XJCAPP(NF), SCNEW(NF), FWEIGHT(NF), RADIUS(ND) 
      DIMENSION OPAC(NF), ETAC(NF), XLAMBDA(NF)
      DIMENSION ELEVEL(NDIM), INDLOW(LASTIND), INDNUP(LASTIND)
      DIMENSION RRATE(NDIM,NDIM), EN(NRANK)
      DIMENSION KODRLOW(LASTKDR), IONGRND(NDIM)
      DIMENSION DRJLW(NDIM), DRJLWE(NDIM), DRLJW(NDIM)
      DIMENSION FEDDI(ND,NF), QEDDI(ND,NF), HEDDI(NF)
      DIMENSION OPAC1(NF), XJCLP1(NF)
      DIMENSION HTOT(ND), HTOTC(ND), HTOTL(ND)
      DIMENSION SLNEW(LASTIND), OPAL(LASTIND)

      SAVE HNULL, DHINT, HJMEAN, DHRMAX

C***  C1 = H * C / K        (CM * KELVIN)
      DATA C1 / 1.4388 /
C***  C2 = 2 * H * C     ( CGS UNITS )
      DATA C2 / 3.9724E-16 /
C***  DATA: STEBOL = STEFAN-BOLTZMANN CONSTANT (CGS-UNITS) / PI
      DATA STEBOL / 1.8046E-5 /

C***  PI8 = 8*PI
      DATA PI8 /25.1327412288 /

C***  PREPARATIONS AT EACH DEPTH POINT
      RL2 = RADIUS(L) * RADIUS(L)

C***  OPASMEAN = Source-function weighted mean of opacity (CGS UNITS)
      SUM   = 0. 
      SMEAN = 0.
      DO K=1, NF
         SUM   = SUM + OPAC(K) * SCNEW(K) * FWEIGHT(K)
         SMEAN = SMEAN + SCNEW(K) * FWEIGHT(K)
      ENDDO
C***  Conversion into cgs-units
      OPASMEAN = SUM / (SMEAN * RSTAR)

      VDOPCM = VDOP * 1.E5
      ETAL = 0.
      DO IND=1, LASTIND
         LOW = INDLOW(IND)
         NUP = INDNUP(IND)
         W = ELEVEL(NUP) - ELEVEL(LOW)
         ETAL = ETAL + (OPAL(IND) * SLNEW(IND) * W)
      ENDDO
C***  Conversion into cgs-units
      ETAL = ETAL * VDOPCM / (SMEAN * RSTAR)
      OPASMEAN2 = OPASMEAN + ETAL

C***  QFJMEAN = J-weighted mean of (qf) at RL
      SUM   = 0.
      XJMEAN = 0.
      DO K=1, NF
         SUM   = SUM + QEDDI(L,K) * FEDDI(L,K) * XJCAPP(K) * FWEIGHT(K)
         XJMEAN = XJMEAN + XJCAPP(K) * FWEIGHT(K)
      ENDDO
      QFJMEAN = SUM / XJMEAN 

C***  OPAJMEAN = J-weighted mean of OPACITY at RL
C***  NOTE: OPAJMEAN without Thomson-Scatter
      SUM   = 0.
      DO K=1, NF
         SUM   = SUM + OPAC(K)              * XJCAPP(K) * FWEIGHT(K)
      ENDDO
C***  Conversion into cgs-units
      OPAJMEAN = SUM / (XJMEAN * RSTAR)


C***  CLASSICAL FORMULATION: LOCAL ENERGY BALANCE

C***  1. CONTINUA:  INTEGRAL OVER FREQUENCY  ==============================
      FTCONT=.0
      DO K=1,NF
         FTCONT = FTCONT + (ETAC(K) - OPAC(K)*XJCAPP(K)) * FWEIGHT(K)
      ENDDO
C***  CONVERT INTO CGS UNITS
      FTCONT = FTCONT / RSTAR

C***  2. DIEL. RECOMBINATION/AUTOIONIZATION: LOOP OVER ALL DR CONTINUA
      FTDR = .0
      DO KDR=1, LASTKDR
         LOW = KODRLOW(KDR)
         NUP = IONGRND(LOW)
         FTDR = FTDR + EN(NUP) * DRJLW(LOW) - EN(LOW) * DRLJW(LOW)
      ENDDO
C***  CONVERT INTO CGS UNITS: RATES * h nue / 4 pi, ABSOLUTE POPNUMS      
      FTDR= FTDR * ENTOTL * C2 / PI8 

C***  3. LINES: LOOP OVER ALL BOUND-BOUND TRANSITIONS
      FTLINES = .0
      DO IND=1, LASTIND
         LOW = INDLOW(IND)
         NUP = INDNUP(IND)
         W = ELEVEL(NUP) - ELEVEL(LOW)
         FTLINES = FTLINES + W * 
     >    (EN(NUP) * RRATE(NUP,LOW) - EN(LOW) * RRATE(LOW,NUP))
      ENDDO
C***  CONVERT INTO CGS UNITS: RATES * h nue / 4 pi, ABSOLUTE POPNUMS      
      FTLINES = FTLINES * ENTOTL * C2 / PI8 

C***  ADD ALL CONTRIBUTIONS: CONTINUUM + DIEL.RECOMBINATION + LINE 
      FT = FTCONT + FTDR + FTLINES

C***  Under the assumption S = Planck, the error term FT is converted 
C***  into a temperature correction 

      DTB     = 1. / (TL * TL * TL * 4. * STEBOL * OPASMEAN2)
      DTLOCAL = - FT      * DTB
      DTLOCC  = - FTCONT  * DTB
      DTLOCL  = - FTLINES * DTB
      DTLOCD  = - FTDR    * DTB

      IF (L .EQ. 1) THEN
        HNULL = 0.25 * STEBOL * TEFF*TEFF*TEFF*TEFF
C***    Randterm
C***    HJMEAN = J-weighted mean of (H) at R1 (Outer Boundary)
        SUM   = 0.
        DO K=1, NF
           SUM   = SUM + HEDDI(K) * XJCAPP(K) * FWEIGHT(K)
        ENDDO
        HJMEAN = SUM / XJMEAN 
C***    Prepare third Term of Unsoeld-Lucy-Procedure: 
C***      contribution of the Boundary
        DHRMAX = (HNULL - HTOTL(1)) / HJMEAN - XJMEAN * RL2
        DHINT = .0

      ELSE
C***    Second Term of Unsoeld-Lucy-Procedure: contribution to depth integral
C***    H-total at radius-interstice (L,L-1)
C***    Note the different indexing of HTOT... from 1...ND-1
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
          HNUE = - (QFRL * XJCAPP(K) - QFRLP1 * XJCLP1(K)) / DTQ
          HTOTCLM1 = HTOTCLM1 + HNUE * FWEIGHT(K) 
C***      Occasional computation of QOPAHMEAN 
C***        = FLUX-weighted mean of q*opacity (CGS UNITS) at interstice (L,L-1)
          SUM = SUM + Q * X * HNUE * FWEIGHT(K)

    6   CONTINUE

C***  Normalization to H
        QOPAHMEAN = SUM / HTOTCLM1

C***  ADD Line and Continuum fluxes
        HTOT(L-1) = HTOTC(L-1) + HTOTL(L-1)

        DH = HTOT(L-1) - HNULL
        WRADIUS = RADIUS(L) - RADIUS(L-1) 
C***  Integral summand for actual L,L-1 interval
        DHINT = DHINT + DH * WRADIUS * QOPAHMEAN 
C***  Factors outside integral
        DTINT = OPAJMEAN * DHINT / 
     >      (TL * TL * TL * 4. * STEBOL * OPASMEAN * RL2 * QFJMEAN)
      ENDIF

C***  Third Term of Unsoeld-Lucy-Procedure: 
C***    contribution of the Boundary
      DTRMAX = OPAJMEAN * DHRMAX / 
     >    (TL * TL * TL * 4. * STEBOL * OPASMEAN * RL2)

C***  TEST PRINTOUT
      IF (.false.) THEN
         IF (L .EQ. 1) PRINT 98, DUNLU, DUNLUR
   98    FORMAT (//, ' Damping of Terms II and III: ', 2(1X,F5.2), //, 
     >   '    L    FT-LOCAL  T(OLD)  DT-Local   DT-Integral   DT-Rmax')
         PRINT 99, L, FT, TL, DTLOCAL, DTINT, DTRMAX
   99    FORMAT (I6, E12.4, 4(1X,G20.6))

      ENDIF

C***  PRINT OPTION: PRESENT VERSION ONLY FOR CLASSICAL FORMULATION
      IF (LSTEQ .GT. 0) THEN
         IF (((L-1)/LSTEQ)*LSTEQ .EQ. (L-1)) 
     >       CALL PRITEQ (L,NRANK,XLAMBDA,TL,
     $                 LASTIND,INDLOW,INDNUP,NDIM,ELEVEL,RRATE,EN,N,
     $                 NF,ND,XJCAPP,FWEIGHT,LASTKON,KONTNUP,KONTLOW,
     $                 NFEDGE,ENLTE,EXPFAC,SIGMAKI,LEVEL,
     $                 DRJLW,DRLJW,IONGRND,KODRLOW,LASTKDR,FT,FT2,
     $                 RNEL,ENTOTL,SIGMAFF,MAXION,NCHARG,OPAMEAN)
         ENDIF

C***  MODIFY TL IN ORDER TO IMPROVE THE TEMPERATURE
C***  Maximum Correction 10 percent
      DTL = DTLOCAL + DUNLU * DTINT + DUNLUR * DTRMAX
      IF (ABS(DTL) .GT. 0.1 * TL) THEN
        IF (DTL .LT. 0.) THEN
          DTL = -0.1 * TL
        ELSE
          DTL = 0.1 * TL
        ENDIF
      ENDIF
      TL = TL + DTL

C***  REPLACE TEMPERATURES BELOW TMIN BY TMIN
      IF (TL .LT. TMIN) THEN
        TL = TMIN
        NEGT = NEGT + 1
      ENDIF
    
      RETURN
      END
