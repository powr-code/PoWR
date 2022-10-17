      SUBROUTINE LIMB_INFO(JPFIRST, NP, ND, MAXXLIMB, NFOBS, K, JP,
     >                     NPBAND, LIMB_UNIT, N_INTBIN, LIMB_LINE, 
     >                     BRANCH, MODHEAD, LIMBFILENAME, XLAMREF, 
     >                     ALN, DXOBS, RANGERED, RANGEBLUE, FINCRI, 
     >                     SUMFILT, BEGLAM, JPLASTLIMB, EMINTK, 
     >                     SUMINT, P, TAUROSS, BANDLAM, WGTBAND, 
     >                     WGTBANDMTX, DXLAM)

C***************************************************************************
C*** This routine prepares and plots/writes the intensities I_nu(p) 
C*** This is available either for (1) a specific wavelength, 
C*** (2) a band-average, or (3) to write data in a file.
C*** The routine has TWO main branches called from formal via "Branch". 
C*** The first, 'DECODE', decodes the LIMB_LINE. The second branch, which 
C*** has several options, is called from the frequency loop and evaluates
C*** the needed data. Options:
C*** (1) PLOT LIMB LAM x.x [PCUT|TAUCUT x.x] [NORMALIZED] [MU]: Plots
C*** I_nu(p) at lam = x.x up to impact parameter p=pcut (or the one implied
C*** TAUCUT, if not stated then up to p(NP)). Can be normalized to I(p=0).
C*** MU plots it with MU as the X-axis (MU = p/plast). 
C*** (2) PLOT LIMB BAND U|B|V|u|b|v|y [PCUT|TAUCUT x.x] [NORMALIZED] [MU]:
C*** Plots averaged I_nu(p) over chosen band using pre-defined filter 
C*** functions.
C*** (3) WRITE LIMB [LAM1=x.x] [LAM2=x.x] [SPACING=x.x] [FILE=ssss]
C*** [PCUT|TAUCUT x.x]: writes I_nu(p) at equally spaced (spacing) wavelengths
C*** between LAM1 and LAM2. The intensities are binned over each of these
C*** intervals and written out for every impact parameter at every 
C*** mid-wavelength. If not spcified, default values are taken. 
C************************************************************************** 
!       IMPLICIT NONE
      
      INTEGER, INTENT(IN) :: JPFIRST, NP, ND, MAXXLIMB, NFOBS, K, JP
      INTEGER :: NPBAND, LIMB_UNIT, N_INTBIN, JPLASTLIMB
      CHARACTER LIMB_LINE*(*), BRANCH*(*), MODHEAD*(*)
      CHARACTER LIMBFILENAME*(*)
      REAL, INTENT(IN) :: XLAMREF, ALN, DXOBS, RANGERED, RANGEBLUE
      REAL :: FINCRI, SUMFILT, BEGLAM, EMINTK, DXLAM
      REAL, DIMENSION (NP):: SUMINT
      REAL, DIMENSION (NP), INTENT(IN) :: P
      REAL, DIMENSION (ND), INTENT(IN) :: TAUROSS
      REAL, DIMENSION (MAXXLIMB) :: BANDLAM, WGTBAND
      REAL, DIMENSION (MAXXLIMB,7):: WGTBANDMTX
      
      DIMENSION UJ(14),BJ(22),VJ(25),US(28),VS(29)
      DIMENSION BS(29),YS(29)
      CHARACTER BANDNAME(7)*2
      DIMENSION FCAL(7),NFILT(7),FBLUE(7),FINCR(7), FCENTER(7)
      
C***  LOCAL VARIABLES:
      CHARACTER ACTPAR*80
      REAL :: XLAMNEXT, XLAMLAST, FREQLAST, FREQNEXT, MU, PCUT
      REAL :: TAUCUT, DUMMY, DNUE, FREQ
      REAL :: MESHLAM, WGT
      INTEGER :: I, J, L, NPAR, PLIPO
      LOGICAL BSKIP
      
      REAL, PARAMETER :: CLIGHT2 = 2.9979E18 
      REAL, PARAMETER :: PI = 3.1415926535898 
      
C***  BAND NAMES IN THIS ROUTINE, corresponding to J=1,2,...7
      DATA BANDNAME / 'U', 'B', 'V', 'u', 'v', 'b', 'y' /


C***  CENTER WAVELENGHTS OF THE BANDS (ONLY FOR IDENTIFICATION)
      DATA FCENTER / 3500., 4350., 5550., 3500., 4100., 4700., 5500. /

C***         NUMBER OF DATA POINTS
C***       (1 - 3: JOHNSON U,B,V : 4 - 7: STROEMGREN U,V,B,Y)
      DATA NFILT / 14, 22, 25, 28, 29, 29, 29 /

C***         BLUE EDGES OF FILTERS (NUMBERS: SEE NFILT)
      DATA FBLUE / 2800., 3400., 4700., 3125., 3725., 4325., 5000. /

C***         INCREMENTS OF WAVELENGHTS (NUMBERS: SEE NFILT)
      DATA FINCR / 100., 100., 100., 25., 25., 25., 25. /

C***         CALIBRATION 
      DATA FCAL /  49.424, 48.525, 48.624, 
     -             48.324, 48.316, 48.452, 48.600 /

C***         JOHNSON FILTER FUNCTIONS

      DATA UJ/0.000,0.025,0.250,0.680,1.137,1.650,2.006,2.250,2.337,    
     -     1.925,0.650,0.197,0.070,0.000/               
      DATA BJ/0.000,0.006,0.080,0.337,1.425,2.253,2.806,2.950,3.000,    
     -     2.937,2.780,2.520,2.230,1.881,1.550,1.275,0.975,0.695,   
     -     0.430,0.210,0.055,0.000/                 
      DATA VJ/0.020,0.175,0.900,1.880,2.512,2.850,2.820,2.625,2.370,    
     -     2.050,1.720,1.413,1.068,0.795,0.567,0.387,0.250,0.160,   
     -     0.110,0.081,0.061,0.045,0.028,0.017,0.007/            

C***     STROEMGREN FILTER FUNCTIONS

      DATA US/00.00,00.60,05.92,11.80,17.80,23.80,29.12,33.32,36.08,    
     -        38.12,39.28,39.44,39.08,38.52,37.08,35.32,33.08,30.08,    
     -        26.16,21.68,16.92,11.76,07.68,04.68,02.16,00.88,00.40,    
     -        00.00/                            
      DATA VS/00.00,00.16,00.36,00.84,01.48,02.28,03.04,04.80,07.84,    
     -        13.00,20.00,29.88,40.00,47.28,49.32,48.00,43.52,37.20,    
     -        28.20,18.12,11.12,06.72,04.00,02.68,02.00,01.40,00.72,    
     -        00.36,00.00/                      
      DATA BS/00.00,00.40,00.96,01.60,02.32,03.60,05.00,08.00,12.40,    
     -        20.00,30.20,40.40,45.88,46.88,45.04,38.20,27.84,17.40,    
     -        11.00,06.88,03.92,02.52,01.52,01.24,00.92,00.72,00.40,    
     -        00.20,00.00/                      
      DATA YS/00.00,00.68,01.68,02.72,04.00,06.96,10.20,15.20,23.40,    
     -        33.00,41.00,45.40,47.92,50.12,52.20,51.32,45.44,34.80,    
     -        24.80,15.48,10.00,06.60,04.16,02.52,01.36,00.96,00.80,    
     -        00.40,00.00/    


C*****************  FIRST BRANCH *****************   
C***  Decodes the line LIMB_LINE given by the user in FORMAL_CARDS

C**** WEIGHTSMATRIXBAND(I,J) gives the weight of the J'th band at the I'th mesh wavelength
      IF (BRANCH(1:6) .EQ. 'DECODE') THEN
        FINCRI = .0
        DO J=1, NP
           SUMINT(JP) = .0
        ENDDO
        SUMFILT = .0
C**** WEIGHTSMATRIXBAND(I,J) gives the weight of the J'th band at the I'th mesh wavelength
        DO I=1, NFILT(1)
            WGTBANDMTX(I,1) = UJ(I)
        ENDDO
        DO I=1, NFILT(2)
            WGTBANDMTX(I,2) = BJ(I)
        ENDDO
        DO I=1, NFILT(3)
            WGTBANDMTX(I,3) = VJ(I)
        ENDDO
        DO I=1, NFILT(4)
            WGTBANDMTX(I,4) = US(I)
        ENDDO
        DO I=1, NFILT(5)
            WGTBANDMTX(I,5) = VS(I)
        ENDDO
        DO I=1, NFILT(6)
            WGTBANDMTX(I,6) = BS(I)
        ENDDO
        DO I=1, NFILT(7)
            WGTBANDMTX(I,7) = YS(I)
        ENDDO
C***  Decodes the line LIMB_LINE given by the user in FORMAL_CARDS
        CALL SARGC(LIMB_LINE, NPAR)
        IF ((NPAR .LT. 2) .OR. (NPAR .GT. 12)) GOTO 100
        DO I=1, NPAR
            CALL SARGV(LIMB_LINE, I, ACTPAR)
C***  Intensities taken up to impact parameter PCUT
            IF (ACTPAR .EQ. 'PCUT') THEN
              CALL SARGV(LIMB_LINE, I+1, ACTPAR)
              READ (ACTPAR, '(F12.0)') PCUT
C***  Intensities taken up to impact parameter where tauross=2/3
              JPLASTLIMB = ISRCHFGT(NP,P,1,PCUT) - 1
            ELSE IF (ACTPAR .EQ. 'TAUCUT') THEN
              CALL SARGV(LIMB_LINE, I+1, ACTPAR)
              READ (ACTPAR, '(F12.0)') TAUCUT
              JPLASTLIMB = NP - ISRCHFGT(ND,TAUROSS,1,TAUCUT) - 1
              PCUT = P(JPLASTLIMB+1)
            ENDIF
        ENDDO
        CALL SARGV(LIMB_LINE, 1, ACTPAR)
C** WRITE BRANCH:
C** True if intensities are to be written in a file:
C** User generally specifies blue+red boundaries LAM1 and LAM2
C** And desired spacing FINCRI (SPACING in FORMAL_CARDS). Otherwise
C** Default values are taken.  
C** The intensities at each impact parameters are then binned over 
C** wavelength intervals of FINCRI angstrom and printed out for each
C** mid-point of such an interval in the following format:
C** First & second line: File and model data
C** third & fourth line: Impact parameter and corresponding mu 
C** Further lines: I_nu(p) for each mesh wavelength and impact parameter
        IF (ACTPAR .EQ. 'WRITE') THEN
            BRANCH='WRITE'
C** Default file name
            LIMBFILENAME = "limb.txt"
            DO I=1, NPAR
                CALL SARGV(LIMB_LINE, I, ACTPAR)
C** blue wavelength boundary for file
                IF (ACTPAR .EQ. 'LAM1') THEN
                    IF (I+1 .GT. NPAR) GOTO 100
                    CALL SARGV(LIMB_LINE, I+1, ACTPAR)
                    READ(ACTPAR, '(F12.0)') BANDLAM(1)
C** red wavelength boundary
                ELSE IF (ACTPAR .EQ. 'LAM2') THEN
                    IF (I+1 .GT. NPAR) GOTO 100
                    CALL SARGV(LIMB_LINE, I+1, ACTPAR)
                    READ(ACTPAR, '(F12.0)') BANDLAM(2)
C** spacing between wavelengths in file
                ELSE IF (ACTPAR .EQ. 'SPACING') THEN
                    IF (I+1 .GT. NPAR) GOTO 100
                    CALL SARGV(LIMB_LINE, I+1, ACTPAR)
                    READ(ACTPAR, '(F12.0)') FINCRI
C** User may specify file name
                ELSE IF (ACTPAR .EQ. 'FILE') THEN
                    CALL SARGV(LIMB_LINE, I+1, LIMBFILENAME)
                ENDIF
            ENDDO
C**  Default values if values not specified
            IF (BANDLAM(1) .EQ. 0.) THEN 
                BANDLAM(1) = RANGEBLUE
                WRITE(0,'(A,F10.0)') "Lam1 not specified, set to: ", RANGEBLUE
            ENDIF
            IF (BANDLAM(2) .EQ. 0) THEN
                BANDLAM(2) = RANGERED
                WRITE(0,'(A,F10.0)') "Lam2 not specified, set to: ", RANGERED
            ENDIF
C*** Order of lam1 lam2 doesn't matter
            IF (BANDLAM(1) .GT. BANDLAM(2)) THEN
                DUMMY = BANDLAM(1)
                BANDLAM(1) = BANDLAM(2)
                BANDLAM(2) = DUMMY
            ENDIF   
            IF (FINCRI .EQ. 0) THEN 
                FINCRI = 100.
                WRITE(0,'(A,F10.2)') "Spacing not specified, set to: ", 100.
            ENDIF
            NPBAND = 2
C** Limb file is opened
            OPEN (UNIT=LIMB_UNIT, FILE=LIMBFILENAME, ACTION="write", 
     >            status="unknown")
C** Basic information written
            WRITE(LIMB_UNIT, '(A,F8.2,A,F8.2,A,F5.2, A, F8.2, A)') 
     >      '#Intensities I_nu from ', BANDLAM(1), ' to ', BANDLAM(2),
     >      ' [Ang] from impact parameter', P(JPFIRST), ' to ', 
     >      P(JPLASTLIMB), '[R*] in erg/s/cm^2/nu(!)'
            WRITE(LIMB_UNIT, '(A,A)') '#',MODHEAD
C** Impact parameter line
            WRITE(LIMB_UNIT, '(A12)', advance='no') 'P[R*]: '
            DO J=JPFIRST, JPLASTLIMB            
                WRITE(LIMB_UNIT, '(1P,G12.4)', advance='no')  P(J)
            ENDDO
            WRITE(LIMB_UNIT, '(A)') ''
C** Mu line
            WRITE(LIMB_UNIT, '(A12)', advance='no') 'MU: '
            DO J=JPFIRST, JPLASTLIMB            
                MU = P(J) / P(JPLASTLIMB)
                WRITE(LIMB_UNIT, '(1P,G12.4)', advance='no')  MU
            ENDDO
C** Printout continues in frequency loop.

C** PLOT BRANCH:
C** True if I_nu(p) is plotted at specific wavelength/band
        ELSE IF (ACTPAR .EQ. 'PLOT') THEN
            CALL SARGV(LIMB_LINE, 3, ACTPAR)
C** Intensities I_nu(p) are averaged over band via: 
C** I(p) = int(I_nu(p) * FILTER_WGT(nu) dnu) / int(FILTER_WGT dnu)
            IF (ACTPAR .EQ. 'BAND') THEN
                BRANCH = 'PLOT BAND'
                J = 1
                CALL SARGV(LIMB_LINE, 4, ACTPAR)
C** Band's index is searched
                DO WHILE (BANDNAME(J) .NE. ACTPAR)
                    J = J + 1
                    IF (J .EQ. 8) THEN
                        WRITE (0,'(A)') "Specified band not found"
                        STOP "Error in Subroutine LIMB_INFO"
                    ENDIF
                ENDDO
                WRITE (0,'(A,A)') "Plotting I(p) in BAND ", BANDNAME(J)
C*** Number of filter wavelengths
                NPBAND = NFILT(J)
C*** Band mesh wavelengths and corresponding weights
                DO I=1, NPBAND
                    BANDLAM(I) = FBLUE(J) + I*FINCR(J)
                    WGTBAND(I) = WGTBANDMTX(I,J)
                ENDDO
C*** Warning if range does not contain band
                IF ((RANGEBLUE .GT. BANDLAM(1)) .OR. 
     >             (RANGERED .LT. BANDLAM(NPBAND))) THEN
                    WRITE(0,'(A,F10.0,F10.0)') 
     >                  "WARNING: Range must contain band ",
     >              BANDLAM(1), BANDLAM(NPBAND)
                    STOP 'Fatal error in subroutine LIMB_INFO'
                ENDIF
C** continues in frequency loop.

C** The intensities I_nu(p) are evaluated at a single frequncy
            ELSE IF ((ACTPAR .EQ. 'LAM') .OR. (ACTPAR .EQ. 'LAMBDA')) THEN
                BRANCH = 'PLOT LAM'
                CALL SARGV (LIMB_LINE, 4, ACTPAR)
C** BANDLAM(1) holds the frequency at which I_nu(p) is to be evaluated
                READ (ACTPAR, '(F12.0)') BANDLAM(1)
                WRITE(0,'(A,F10.2)') "Plotting I(p) at wavelength ", 
     >                            BANDLAM(1)
                IF ((RANGEBLUE .GT. BANDLAM(1)) .OR. 
     >              (RANGERED .LT. BANDLAM(1))) THEN
                    WRITE(0,'(A,F10.5,A)') "Wavelength ",BANDLAM(1), 
     >              " is not in range"
                    STOP "Fatal error in subroutine LIMB_INFO"
                ENDIF
                NPBAND = 1
            ELSE 
                GOTO 100
            ENDIF
        ELSE
            GOTO 100
        ENDIF
        RETURN
C****** END OF DECODING BRANCH **********************************        
        
C****** SECOND BRANCH: CALLED FROM FREQUENCY LOOP ***************
      ELSE
C*** Nothing needs to be evaluated beyond JPLASTLIMB
        IF (JP .GT. JPLASTLIMB) RETURN
C*** WRITE BRANCH continues here!
        IF (BRANCH .EQ. 'WRITE') THEN
C*** If obs wavelength outside desired range for file - skip.
            BSKIP = (XLAMREF .LT. BANDLAM(1)) .OR. 
     >              (XLAMREF .GT. BANDLAM(2))
            IF (BSKIP) RETURN
C*** These statements are performed only for the first JP at every frequency
            IF (JP .EQ. JPFIRST) THEN
C*** N_INTBIN = number of intensities which were summed for binning
C*** N_INTBIN = 0 => New binning interval is to be defined, hence 
C*** new line in file.
               IF (N_INTBIN .EQ. 0) THEN
                WRITE(LIMB_UNIT, '(A)') ''
C*** BEGLAM = blue wavelength of current binning interval
                BEGLAM = XLAMREF
C*** MESHLAM = midpoint of the binning interval ("representative of interval")
                MESHLAM = BEGLAM + 0.5*FINCRI
                WRITE(LIMB_UNIT,'(F12.2)', advance='no') MESHLAM
               ENDIF
C*** DXLAM = interval between current wavelength and blue wavelength for 
C*** binning interval. Is always smaller than spacing!
               DXLAM = XLAMREF - BEGLAM
               N_INTBIN = N_INTBIN + 1
            ENDIF
C*** As long as interval smaller than spacing, or frequncy isn't last one,
C*** keep adding up the intensities at each JP
            IF ((DXLAM .LT. FINCRI) .AND. (K .NE. NFOBS)) THEN
               SUMINT(JP) = SUMINT(JP) + EMINTK
C*** ELSE: average sum over binning interval, and write result.
            ELSE 
               SUMINT(JP) = SUMINT(JP) / N_INTBIN
               WRITE(LIMB_UNIT,'(1PG12.4)', advance='no') SUMINT(JP)
               SUMINT(JP) = .0
               IF (JP .EQ. JPLASTLIMB) N_INTBIN=0
            ENDIF           
C*** Plotting at specific wavelength lam continues here!
        ELSE IF (BRANCH .EQ. 'PLOT LAM') THEN
C*** XLAMNEXT/LAST: Next/last observer's frame wavelengths.
            IF (JP .EQ. JPFIRST) THEN 
                XLAMNEXT=XLAMREF*EXP(-DXOBS*ALN)
                XLAMLAST=XLAMREF*EXP(DXOBS*ALN)
            ENDIF
C*** XLAMLAST >= lam => I_nu(lam,JP) already saved, skip routine.
            BSKIP = (XLAMLAST .GT. BANDLAM(1))
            IF (BSKIP) RETURN
C*** True if current wavelength is closest to lam "from the left" (bluer)
            IF ((XLAMREF .LT. BANDLAM(1)) .AND. 
     >         (XLAMNEXT .GT. BANDLAM(1))) THEN
                SUMINT(JP) = SUMINT(JP) + EMINTK
C*** True only if last statement was already true ==> for next wavelength!
            ELSE IF (SUMINT(JP) .NE. 0) THEN
C*** PLIPO = linear interpolation between lam- and lam+
                PLIPO = (BANDLAM(1) - XLAMLAST) / (XLAMREF - XLAMLAST)
                SUMINT(JP) = PLIPO * EMINTK + (1 - PLIPO) * SUMINT(JP)
            ENDIF
C*** I_nu(JP) averaged over band.
        ELSE IF (BRANCH .EQ. 'PLOT BAND') THEN
C*** Skip if current wavelength outside band
            BSKIP = (XLAMREF .LT. BANDLAM(1)) .OR. 
     >              (XLAMREF .GT. BANDLAM(NPBAND))
            IF (BSKIP) RETURN
C*** following statements only at first impact parameter of each frequency
            IF (JP .EQ. JPFIRST) THEN
C*** wavelength ---> Frequency (in hz!)
                FREQ = CLIGHT2 / XLAMREF
C*** corresponding filter function weight interpolated at current wavelength
                CALL LIPO (WGT, XLAMREF, WGTBAND, BANDLAM, NPBAND)
C*** blue edge
                IF (K .EQ. 1) THEN
                    XLAMNEXT = XLAMREF * EXP(-DXOBS*ALN)
                    FREQNEXT = CLIGHT2 / XLAMNEXT
                    DNUE = 0.5*(FREQ - FREQNEXT)
C*** red edge
                ELSE IF (K .EQ. NFOBS) THEN
                    XLAMLAST = XLAMREF * EXP(DXOBS*ALN)
                    FREQLAST = CLIGHT2 / XLAMLAST
                    DNUE = 0.5*(FREQLAST - FREQ)
                ELSE 
                    XLAMNEXT = XLAMREF * EXP(-DXOBS*ALN)
                    XLAMLAST = XLAMREF * EXP(DXOBS*ALN)
                    FREQNEXT = CLIGHT2 / XLAMNEXT
                    FREQLAST = CLIGHT2 / XLAMLAST
                    DNUE = 0.5*(FREQLAST - FREQNEXT)
                ENDIF
C*** SUMFILT is the normalization factor of the filters: int(wgt(nu) dnu)
                SUMFILT = SUMFILT + WGT*DNUE
            ENDIF
C*** SUMINT(p) = int(I_nu(p) * wgt(nu) dnu)
            SUMINT(JP) = SUMINT(JP) + EMINTK*WGT*DNUE
C*** True if band is now covered. Calculated final result:
C*** int (I_nu (p) wgt(nu) dnu) / int(wgt(nu) dnu)
            IF (XLAMNEXT .GT. BANDLAM(NPBAND)) THEN
                SUMINT(JP) = SUMINT(JP) / SUMFILT
            ENDIF
        ELSE
            WRITE (0,'(A, A)') "Unknown branch: ", BRANCH
            STOP "Fatal error in Subroutine LIMB_INFO"
        ENDIF
        RETURN
      ENDIF
            
    



C**** ERROR BRANCHES

  100 WRITE(0,'(A)') "Wrong syntax"
      WRITE(0,'(A)') "The error occured in the following line"
      WRITE(0,'(A)') LIMB_LINE
      WRITE(0,'(A)') "Possible syntax:"
      WRITE(0,'(A)') "PLOT LIMB BAND U|B|V|u|v|b|y [PCUT|TAUCUT x.x] [MU] [NORMALIZED]"
      WRITE(0,'(A)') "PLOT LIMB LAM x.x [PCUT|TAUCUT x.x] [MU] [NORMALIZED]"
      WRITE(0,'(A)') "WRITE LIMB LAM1 x.x LAM2 x.x SPACING x.x [PCUT|TAUCUT x.x.] [FILE = ...]"
      STOP 'Fatal error in subroutine LIMB_INFO'
      
      END
