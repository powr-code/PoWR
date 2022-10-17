	SUBROUTINE STARKDAMP_HEI (GAMMAHE1, SHIFT, T, XNE, XLAM, 
     >                NUP, LOW, LINPRO, LEVEL)
C ***
C ******************************************************************************
C ***      
C ***	ABSORPTIONSKOEFFIZIENT ISOLIERTER HELIUMLINIEN.	
C ***   NACH GRIEM, PHYS.REV. 125,177 (1962)
C ***	STARK-PARAMETER NACH BENNETT/GRIEM, UNIVERS. MARYLAND REPORT (71)
C ***
C ***   T       > R	Temperature [K]
C ***   XNE     > R	Electron density [cm**-3]
C ***   K       > I	Line index, see below
C ***   XLAM    > R	Central wavelength of line [Å]
C ***
C ***   Modifications:
C ***
C ***   28-OCT-1997 22:58:34.11		ML / Stegaurach.  Copied 4026 & 4471
C ***                                   from BCSS.FOR.  Verify correctness...
C ***   23-May-2008 --> only gam4 berechnung / 
C***             short fast version without bcss data tables
C ***   		older version linhe1.f
C ***
C ******************************************************************************

      CHARACTER LINPRO*8, LEVEL(2)*10
      REAL TEM(4), ALF(4,16), W(4,16), SHF(4,16), EP(16)
      REAL BK, AN0, HM, PI, HC, C

      parameter( bk = 1.38046E-16, an0 = 8.85255E-13,
     >             hm = 1.6732E-24, pi = 3.14159, hc = 1.986182E-16,
     >             c = 2.997929E+10 )

      real c1, c2, c3, c4, c5, c6, c7, c8  ! DATA doesn't allow expressions
      parameter( c1 = 0.172*5.62, c2 = 0.193*5.62,
     >             c3 = 0.218*5.62, c4 = 0.249*5.62,
     >             c5 = 0.107*5.62, c6 = 0.119*5.62,
     >             c7 = 0.134*5.62, c8 = 0.154*5.62 ) 

      DATA TEM / 5000., 10000., 20000., 40000. /

******************************************************************
C**** LINES CALCULATED
C**** K = 1    LINE 4922    NOT USED, SEE HELIUM5
C**** K = 2    LINE 4713
C**** K = 3    LINE 4437
C**** K = 4    LINE 4388    NOT USED SEE HELIUM2
C**** K = 5    LINE 4121
C**** K = 6    LINE 3964    BLENDED WITH H. CALCULATED ONLY IF H.LT. .01
C**** K = 7    LINE 3888
C**** K = 8    LINE 5875    CALCULATED IF KHEM .GT. 7
C**** K = 9    LINE 6678    AS 5875
C**** K = 10   LINE 7065    AS 5875
C**** K = 11   LINE 5015   
C**** K = 12   LINE 5048
C**** K = 13   LINE 2945    CALCULTED IF KHEM .GT. 12
C**** K = 14   LINE 2829    AS 2945
C *** K = 15   LINE 4026    Griem 1974, copied from BCSS.FOR -- shift ignored
C *** K = 16   LINE 4471    BCS 1974, JQSRT 14, 1025; BCSS.FOR -- no shift
C ***                       forbidden component ignored.
C***********************************************************************

C***  In order to find the proper "line index" K, we will compare 
C***  XLAM with the approximate wavelengths:
      DIMENSION XLAMK(16)
      DATA XLAMK / 4922., 4713., 4437., 4388., 
     2             4121., 3964., 3888., 5875., 
     3             6678., 7065., 5015., 5048., 
     4             2945., 2829., 4026., 4471.  / 

      DATA ALF / 
     1           0.683, 0.773, 0.885, 1.023,
     2		 0.115, 0.103, 0.095, 0.092,
     3		 0.199, 0.184, 0.177, 0.179,
     4		 1.159, 1.321, 1.527, 1.782,
     5		 0.171, 0.155, 0.145, 0.141,
     6		 0.275, 0.290, 0.311, 0.341,
     7	         0.075, 0.070, 0.067, 0.067,
     8		 0.064, 0.061, 0.059, 0.059,
     9		 0.146, 0.157, 0.169, 0.181,
     Z		 0.067, 0.060, 0.055, 0.052,
     1           0.154, 0.160, 0.169, 0.180,
     2           0.135, 0.123, 0.117, 0.116,
     3		 0.204, 0.195, 0.195, 0.203,
     4		 0.285, 0.276, 0.279, 0.294,
     >           c1, c2, c3, c4,		! see PARAMETER above
     >           c5, c6, c7, c8/

      DATA W / 
     1         2.30,  1.96,  1.63,  1.35,
     2	       0.343, 0.394, 0.438, 0.459,
     3	       1.41,  1.57,  1.65,  1.62,
     4	       6.13,  5.15,  4.24,  3.45,
     5	       0.787, 0.900, 0.987, 1.02,
     6	       1.030, 0.966, 0.877, 0.776,
     7	       0.102, 0.112, 0.117, 0.117,
     8	       0.159, 0.170, 0.176, 0.177,
     9	       0.423, 0.386, 0.349, 0.318,
     Z	       0.180, 0.208, 0.235, 0.254,
     1         0.378, 0.359, 0.334, 0.306,
     2         0.625, 0.704, 0.755, 0.760,
     3	       0.808, 0.857, 0.857, 0.812,
     4	       1.790, 1.870, 1.840, 1.720,
     5         4.040, 3.490, 2.960, 2.470,
     6         1.460, 1.269, 1.079, 0.898/

      DATA SHF / 
     1            1.020,  0.772,  0.584,  0.439,
     2		  0.403,  0.416,  0.391,  0.335,
     3		  1.51,   1.43,    1.24,  0.995,
     4		  2.52,   1.87,    1.38,   1.01,
     5		  0.900,  0.989,  0.811,  0.673,
     6		 -0.647,  -.504,  -.374,  -.266,
     7		  .0743,  .0603,  .0463,  .0347,
     8		 -.0880, -.0552, -.0255, -.0050,
     9		  0.275,  0.233,  0.196,  0.161,
     Z		  0.215,  0.231,  0.227,  0.203,
     1           -0.250, -0.200, -0.152, -0.111,
     2            0.700,  0.685,  0.611,  0.504,
     3		  0.521,  0.411,  0.311,  0.231,
     4		  1.070,  0.830,  0.620,  0.455,
     5            0.000,  0.000,  0.000,  0.000,
     6            0.000,  0.000,  0.000,  0.000/

      DATA EP / 171135., 169087., 171135., 171135., 169087., 166278.,
     A		159856., 169087., 171135., 169087., 166278., 171135.,
     A          159856., 159856., 169087., 169087./

C***  Find K index
      DO KK = 1, 16
         IF (XLAMK(KK) .LT. 0.999*XLAM) CYCLE 
         IF (XLAMK(KK) .GT. 1.001*XLAM) CYCLE 
         K = KK
         GOTO 1
      ENDDO

      WRITE (*,11)  LEVEL(LOW), LEVEL(NUP) 
      WRITE (0,11)  LEVEL(LOW), LEVEL(NUP) 
   11 FORMAT ('*** WARNING from subr. STARKDAMP_HEI: ', 
     >        'no tabulated data for ', A, ' - ', A, 
     >        '--> using Q-STARK instead')

      GAMMAHE1 = 0.0
      LINPRO   = 'Q-STARK '
      RETURN

    1 CONTINUE

ccc test output
ccc      WRITE (0,'(A,2F8.2)') '*** STARKDAMP_HEI: ', XLAMK(K),  XLAM 

C**** Interpolation der Daten
      CALL NEWPOL2(W(1,K), TEM, WE, T, 4 )
      WE = WE * XNE * 1.E-16
      CALL NEWPOL2(ALF(1,K), TEM, ALFI, T, 4 )
      ALFI = ALFI * 1.E-4 * XNE**0.25
      CALL NEWPOL2(SHF(1,K), TEM, D, T, 4 )
      D = D * XNE * 1.E-16
      VV = SQRT(8./PI*BK*T*(1./(4.*HM) + 1./HM/16.))	
      SIG = WE * (0.75/PI/XNE)**(1./3.) / VV / 1.E-8 * C / (XLAM**2)

C**** GESAMTDAEMPFUNG UND VERSCHIEBUNG
      B = ALFI**(8./9.) * SIG**(-1./3.)
      GAMMAHE1 = WE*(1. + 1.36*B)

      IF (D .LT. 0.)  B = - B
      SHIFT = D + WE*2.36*B

C***  Beziehung zwischen gam4 und gamrad
C***  gam = gammahe1 + gamrad/4./pi  
C***  a = gam/dop      

      RETURN
      END
