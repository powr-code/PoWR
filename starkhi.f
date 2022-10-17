      FUNCTION STARKHI(N,M,WAVE,WAVEH,T,XNE, LEVELLOW, LEVELNUP)
C********************************************************************
C     Stark-Broadening for hydrogen (H I) lines
C     XNE - electron density
C     M - upper , N - lower principle quantum number
C
C     The unmodified Code returns a stark-profile for lyman-series
C
C     the modified code returns stark-profiles for lyman, paschen and 
C     bracket-series
C     mario parade - !!! using a vcs-approximation-function vcse1f !
C
C     Returns the Stark broadened line profile.  The code follows Griem's 
C     theories (mostly 1960 ApJ, 132, 883) with corrections to approximate 
C     the Vidal, Cooper & Smith (1973, ApJS 25, 37) profiles.
C
C     Area normalised to unity with frequency.
C
C     by Deane Peterson & Bob Kurucz.
C     (adapted, corrected and comments added by PB)
C************************************************************************

      REAL K
      DIMENSION Y1WTM(2,2),XKNMTB(4,3)
      LOGICAL LYMANALF
      CHARACTER*10 LEVELLOW, LEVELNUP
      SAVE

C  Knm constants as defined by Griem (1960, ApJ 132, 883) for the long 
C  range Holtsmark profile (due to ions only). Lyman and Balmer series 
C  are from VCS, higher series from elsewhere.

      DATA XKNMTB/0.0001716, 0.0090190, 0.1001000, 0.5820000,
     1            0.0005235, 0.0177200, 0.1710000, 0.8660000,
     2            0.0008912, 0.0250700, 0.2230000, 1.0200000/
  
C    new data - mario parade
C
C    DATA KNMTAB/.000356,	.000523,	.00109,		.00149,
C		 .00225,	.0125,		.0177,		.028,
C    		1.0348,		.0493,		.124,		.171,
C		 .223,		.261,		.342,		.683,
C		 .866,		1.02,		1.19,		1.46/
C     DATA FSTARK/.1387,.07910,.02126,.01394,.006462,.004814,.002779,
C    1 .002216,.001443,.001201,.3921,.1193,.03766,.02209,.01139,
C    2 .008036,.005007,.003850,.002658,.002151,.6103,.1506,.04931,
C    3 .02768,.01485,.01023,.006588,.004996,.003542,.002838,.8163,.1788,
C    4 .05985,.03189,.01762,.01196,.007825,.005882,.004233,.003375/

      DATA Y1WTM/1.E18, 1.E17, 1.E16, 1.E14/
      DATA N1/0/, M1/0/

      PARAMETER (CLIGHT = 2.9979258E18)
      PARAMETER (PI = 3.14159265359, SQRTPI = 1.77245385)
      PARAMETER (H = 6.62618E-27)  !Planck in cgs
      PARAMETER (K = 1.38066E-16)  !Boltzmann in cgs

C  Variables depending on conditions

C***  Check if Principle Quantum numbers are known
      IF (N .LE. 0) GOTO 93
      IF (M .LE. 0) GOTO 94
	
cc unused! wrh 19-Feb-2010      HE1FRC=1.0
      T4 = T/10000.
      T43 = T4**0.3
cc unused! wrh 19-Feb-2010     T3NHE = T43*HE1FRC ! fraction H/He
      XNE16 = XNE**0.1666667
cc      PRINT *,(XNE16)
      PP = XNE16*0.08989/SQRT(T) ! the shielding parameter  austausch F=-->F1 mario
      F1 = (XNE16*XNE16*XNE16*XNE16)*1.25E-9      ! Holtsmark normal field strength
cc      PRINT *,(PP)
      Y1B = 2./(1.+0.012/T*SQRT(XNE/T))
      Y1S = T43/XNE16
      C1D = F1*78940./ T
      C2D = F1**2/5.96E-23/XNE
      GCON1 = 0.2+0.09*SQRT(T4)/(1.+XNE/1.E13)
      GCON2 = 0.2/(1.+XNE/1.E15)
cc	PRINT *,(N),(M),(WAVE),(XNE16),(C1D),(C2D)
C

      DELW = WAVE-WAVEH
      FREQNM = CLIGHT/WAVEH
      FREQ = CLIGHT/WAVE
      DEL = FREQ-FREQNM
C
C  Variables dependent on line - compute first time only
C

      IF((N.NE.N1).OR.(M.NE.M1)) THEN  
         N1 = N
         M1 = M
         MMN = M-N
         XN = N
         XN2 = XN*XN
         XM = M
         XM2 = XM*XM
         XMN2 = XM2*XN2
         XM2MN2 = XM2-XN2
         GNM = XM2MN2/XMN2
C
C  Knm constants not tabulated from approximate asymptotic expression 
C  (Griem 1960 eqn 33) where 1./(1.+.13/FLOAT(MMN)) appears to be a 
C  correction factor to match to the tables.
C

         IF ((MMN.LE.3).AND.(N.LE.4)) THEN
            XKNM = XKNMTB(N,MMN)
         ELSE
            XKNM = 5.5E-5/GNM*XMN2/(1.+.13/FLOAT(MMN))
         END IF
	
	
C
C  Some weighting factors which relate to y1, which is the velocity at 
C  which the minimum impact parameter (where second order perturbation 
C  theory breaks down) and the Lewis cutoff (the limit of validity of 
C  the impact approximation) are the same.
C

         IF(M.EQ.2) THEN
            Y1NUM = 550.
         ELSE IF (M.EQ.3) THEN
            Y1NUM = 380.
         ELSE
            Y1NUM = 320.
         END IF
         IF (MMN.LE.2 .AND. N.LE.2) THEN
            Y1WHT = Y1WTM(N,MMN)
         ELSE IF (MMN.LE.3) THEN
            Y1WHT = 1.E14
         ELSE
            Y1WHT = 1.E13
         END IF
C

         C1CON = XKNM/WAVEH*GNM*XM2MN2
         C2CON = (XKNM/WAVEH)**2
cc	 PRINT *,(C1CON),(C2CON)


      ENDIF

C  Compute line profile
C
C  PRQS is the quasistatic ion contribution
C  FNS  is the quasistatic electron contribution rel to PRQS
C  F    is the impact electron contribution
C
C  First compute the width of the impact electron profile roughly Griem
C  (1967, ApJ 147, 1092) eqn for w.
C

      WTY1 = 1./(1.+XNE/Y1WHT)
      Y1SCAL = Y1NUM*Y1S*WTY1+Y1B*(1.-WTY1)
      C1 = C1D*C1CON*Y1SCAL
      C2 = C2D*C2CON
      G1 = 6.77*SQRT(C1)
      BETA = ABS(DELW)/F1/XKNM
      Y1 = C1*BETA
      Y2 = C2*BETA**2

      IF ((Y2.LE.1.E-4).AND.(Y1.LE.1.E-5)) THEN
         GAM = G1*AMAX1(0.,0.2114+LOG(SQRT(C2)/C1))*(1.-GCON1-GCON2)
      ELSE
         GAM = G1*(0.5*EXP(-MIN(80.,Y1))+VCSE1F(Y1)-0.5*VCSE1F(Y2))*
     *            (1.-GCON1/(1.+(90.*Y1)**3)-GCON2/(1.+2000.*Y1))
         IF (GAM.LE.1.E-20) GAM = 0.
      END IF

C
C  Compute individual quasistatic and impact profiles.
C  PP - shielding parameter
C
	
      PRQS = SOFBET(BETA,PP,N,M) 
      IF (GAM.GT.0.) THEN
         F = GAM/PI/(GAM*GAM+BETA*BETA)
      ELSE
	 F = 0.
      ENDIF
C
C  Fraction of electrons which count as quasistatic. A fit to eqn 8 
C  (2nd term) of Griem (1967, ApJ 147, 1092).

      P1 = (0.9*Y1)**2
      FNS = (P1+0.03*SQRT(Y1))/(P1+1.)
C
C  DBETA (=dBeta/dfreq) changes the area normalisation. 
C  DSQRT(WAVE/WAVEH) corrects the long range part to dfreq**-5/2
C  asymptote, (see Stehle and Hutcheon 1999, A&AS 140, 93).

cc      DBETA = CLIGHT/FREQ/FREQ/XKNM/F1
cc    richtiger mit Zentralwellenlaenge:
      DBETA = WAVEH/FREQNM/XKNM/F1
      STARKHI = (PRQS*(1.+FNS)+F)*DBETA * SQRT(WAVE/WAVEH)

C
C  The red wing is multiplied by the Boltzmann factor to roughly account
C  for quantum effects (Stehle 1994, A&AS 104, 509 eqn 7). Assume 
C  absorption case.  If emission do for DEL.GT.0.
C

      IF (DEL.LT.0.) STARKHI = STARKHI * EXP(-ABS(H*DEL)/K/T)
C
      RETURN

C***  Error branches *******************************************

   93 WRITE (0,*) '*** ERROR: Missing MAINQN for lower level, '
     >   , ' needed for Stark broadening of line ' 
     >   , LEVELNUP, ' - ', LEVELLOW      
      WRITE (*,*) '*** ERROR: Missing MAINQN for lower level, ' 
     >   , ' needed for Stark broadening of line ' 
     >   , LEVELNUP, ' - ', LEVELLOW
      GOTO 100

   94 WRITE (0,*) '*** WARNING: Missing MAINQN for upper level, ' 
     >   , ' needed for Stark broadening of line ' 
     >   , LEVELNUP, ' - ', LEVELLOW
      WRITE (*,*) '*** WARNING: Missing MAINQN for upper level, ' 
     >   , ' needed for Stark broadening of line ' 
     >   , LEVELNUP, ' - ', LEVELLOW

      GOTO 100

  100 STOP '*** FATAL ERROR detected by Subr. STARKHI'

      END

