      SUBROUTINE VELTHIN (T,       !Electron Temperature per depth point
     >                    RADIUS,  !Radius (in units of Rstar) per depth point
     >                    VELO,    !Velocity (in km/s) per depth point (modified here)
     >                    ND,      !Total number of depth points
     >                    RSTAR,   !Rstar in Rsun
     >                    RMAX,    !Rmax in Rstar 
     >                    GEFFL,   !g_eff per depth point
     >                    XMU,     !mean particle mass per depth point
     >                    VTURB,   !turbulence velocity (hydrostatic part only) in km/s
     >                    ThinCard,    !CARDS line with HYDROSTATIC INTEGRATION (for parameters)
     >                    CRITERION    !return velocity criterion
     >                   )
C*******************************************************************************
C***  INITIALIZATION OF THE VELOCITY-FIELD VECTOR VELO(L) 
C***  IN THE INNER PART, THE CORRECT HYDROSTATIC EQUATION IS INTEGRATED,
C***  ACCOUNTING FOR THE DEPTH-DEPENDENT SCALE HEIGHT ACCORDING TO THE
C***  GIVEN TEMPERATURE STRATIFICATION AND THE R**2 - SCALING OF GRAVITY
c***  MAJOR REVISION by Andreas Sander (Feb. 2012)
C***  detailed documentation of the first part is given in Sander et al. (2015)
C*******************************************************************************
  
      IMPLICIT NONE
      INCLUDE 'interfacebib.inc'

      INTEGER, INTENT(IN) :: ND
      REAL, INTENT(IN) :: RSTAR, RMAX, VTURB
      !Note: XMU is now an array, earlier versions used only XMASS = XMU(ND)
      REAL, DIMENSION(ND), INTENT(IN) :: T, RADIUS, XMU, GEFFL
      !Do not use INTENT(OUT) for VELO because outer part is read later on
      REAL, DIMENSION(ND), INTENT(INOUT) :: VELO     
      CHARACTER(8), DIMENSION(ND), INTENT(OUT) :: CRITERION
      CHARACTER(80), INTENT(IN) :: ThinCard

      REAL, DIMENSION(4) ::  RIP, VIP
      CHARACTER(40), DIMENSION(20) :: CURPAR

      INTEGER, PARAMETER :: NDTHINMAX = 100.            
      REAL, DIMENSION(NDTHINMAX) :: VBETA, VHELP, RHELP, 
     >                              A2SUM, VELOMINUSSONIC

      INTEGER :: I, J, IMIN, ICON, ICONOLD, K, L, MAXIT, IMAXBETA,
     >           NONMONO, NEXTRA, NSameRcon, LOszil, IE, NPAR, IERR, 
     >           LL, Lguess, LEnforcerNeeded, INMS, INME, NDHELP
      REAL :: DGP, GRAD1, GRAD2, BPL, DR, H, RM, TMID, Q, HCONST, Vcon,
     >        VA, VMAX, S, ST, GRADIN, Vdummy, VINT, VCUR, VFINAL, DXX,
     >        VMIN, BETA, VPAR1, VPAR2, RCON, HSCALE, FPLAST, FPMAX, H0,
     >        BETA2, BETA2FRACTION, VPAR1_2, VPAR2_2, RCONOLD, ROUTMAX,
     >        DVDRcon, VDIFF, RCONMIN, RCONMAX, RGRADMAX, BAD, VADD,
     >        RFA, RFB, RFX, RFFA, RFFB, RFFX, RFD, EXPO, RSTEP, RLAST,
     >        RINT, RIN, XMUCUR, Hlast, RBETAMIN, RCONSTART, P, VAM,
     >        GEFF, GEFFR, A2SUML, A2, DA2DR, BETACON, VOFF, RONSET,
     >        RCONMAXguess, RCONMAXnow, tempREAL, fsonic, DVDR, VSOUND,
     >        VMINUSA, VMINUSAM, DELTAB, VL,
     >        GRADRCONMIN, GRADRCONMAX
      LOGICAL :: bConFound, bConPointConv, bRFbi, bSmoothVelo, bDEBUG,
     >           bUseSonicConnectionCriterion, bForceMonotonic, bRUKU,
     >           bFailsafe, bPrintVelo, bFullHD, bSMOCO, bTestMono, 
     >           bNewSonic

      REAL, EXTERNAL :: DELTAGRTHIN, VELOBETA

      !Common block VELPAR needed to transfer velocity parameters
      ! This routine updates VPAR1, VPAR2, RCON
      COMMON /VELPAR/ VFINAL,VMIN,BETA,VPAR1,VPAR2,RCON,HSCALE,
     >                 BETA2, BETA2FRACTION, VPAR1_2, VPAR2_2

      !Numerical constants
      INTEGER, PARAMETER :: MAXCONIT = 99    !Maximum iterations for connection point search (was 20)
      REAL, PARAMETER :: EXPMAX = 500        !Maximum exponent for velocity scaling (prevent overflow)
      REAL, PARAMETER :: RCONACC = 1.E-6     !Accuracy for RCON
      REAL, PARAMETER :: RFACC = 1.E-12      !Accuracy for Regula Falsi Call
      REAL, PARAMETER :: EPSRCONMIN = 1.E-12 !Offset from beta law singularity

      REAL, DIMENSION(MAXCONIT) :: RCONHIST !History with calculated RCON values


      !Physical Constants
      REAL, PARAMETER :: BOLTZ = 1.38E-16   !BOLTZMANN CONSTANT (ERG/DEG)
      REAL, PARAMETER :: AMU = 1.66E-24     !ATOMIC MASS UNIT (GRAMM)
      REAL, PARAMETER :: RGAS = BOLTZ / AMU !Gas Constant (CGS) = BOLTZ / AMU

      !File and channel handles (=KANAL)
      INTEGER, PARAMETER :: hOUT = 6        !write to wruniqX.out (stdout)
      INTEGER, PARAMETER :: hCPR = 0        !write to wruniqX.cpr (stderr)

      !Defaults for Switches
      bUseSonicConnectionCriterion = .FALSE.
      bForceMonotonic = .TRUE.
      bSmoothVelo = .FALSE.
      bFailsafe = .FALSE.
      LEnforcerNeeded = 0
      bPrintVelo = .FALSE.        !Debug option: Print resulting VELO vector
      bFullHD = .FALSE.
      bRUKU = .TRUE.
      fsonic = 1.

      IF (ND > NDTHINMAX) THEN
        WRITE (hCPR,'(A)') 'VELTHIN: FATAL ERROR ******'
        WRITE (hCPR,'(A)') 'VELTHIN: NDTHINMAX INSUFFICIENT'
        WRITE (hCPR,'(2(A,I4))') 'ND = ', ND,
     >                         ', NDTHINMAX = ', NDTHINMAX
        STOP 'FATAL ERROR IN VELTHIN'      
      ENDIF
          
c      IF (BETA < .0) THEN
c        CALL REMARK ('VELTHIN: BETA OPTION INVALID')
c        STOP 'ERROR'
c      ENDIF
      
      !Read additional CARDS line parameters (if existing)
      CALL SARGC (ThinCard, NPAR)
      DO i=1, NPAR
        CALL SARGV(ThinCard,i,CURPAR(I))
      ENDDO
      IF (NPAR >= 2) THEN
        DO i=2, NPAR 
          SELECTCASE(CURPAR(i))
            CASE ('NONMONO', 'NONMONOTONIC') 
              bForceMonotonic = .FALSE.
            CASE ('SMOOTH')
              bSmoothVelo = .TRUE.
            CASE ('SONIC')
              bUseSonicConnectionCriterion = .TRUE.
              IF (NPAR >= (i+1)) THEN
                READ (CURPAR(i+1), '(F10.0)', IOSTAT=IERR) tempREAL
                IF (IERR == 0) THEN
                  fsonic = tempREAL
                ENDIF                          
              ENDIF
            CASE ('RUKU')
              bRUKU = .TRUE.
            CASE ('NORUKU')
              bRUKU = .FALSE.
            CASE ('FULLHD')
              bFullHD = .TRUE.
            CASE ('PRINTV')
              bPrintVelo = .TRUE.
          ENDSELECT
        ENDDO
      ENDIF
       
  20  IF (bFailsafe) THEN
        !failsafe branch (only reached via jump)
        bUseSonicConnectionCriterion = .TRUE. 
      ENDIF
      
     
C***  LOWER PART OF VELOCITY LAW: 
C***  INTEGRATION OF THE HYDROSTATIC EQUATION

C***  ANSATZ    p = p0 f(r) exp(-(r-1)/HSCALE)
C***            rho = AMU 
C***   
C***            f(r) is obtained numerically by integrating (f = FP in the code)
C***            df/dr = f (1/HSCALE - 1/H)

C***  FIRST STEP: INTEGRATE   dp/dr = -1/H p
C***  WITH DEPTH-DEPENDEND SCALE HEIGHT H = HCONST * T(I)
C***  1. STEP: SOLVE DIFF. EQ. FOR THE CORRECTION FACTOR TO THE
C***  EXPONENTIAL LAW WITH CONSTANT SCALE HEIGHT HSCALE
C***  - THIS FACTOR IS STORED IN VECTOR VELO

      DO L=1, ND
        !initialize VELO vector with VFINAL (just in case)
        VELO(L) = VFINAL
      ENDDO
      
      H0 = HSCALE 

      IMIN = 1
      BPL = 0.
      VELO(ND) = VMIN

C***  Calculate isothermal sound speed a (squared) 
C        and add turbulence contribution
      DO L=1, ND
        A2  = RGAS * T(L) / XMU(L)
        A2SUM(L) = A2 + VTURB*VTURB*1.E10        
      ENDDO
      
C***  Initialize the loop at l=ND
      GEFFR = GEFFL(ND) / (RADIUS(ND) * RADIUS(ND))
      H     = A2SUM(ND) / GEFFR / RSTAR

      hydrostat: DO L=ND-1,1,-1        
        DR = RADIUS(L) - RADIUS(L+1)         


        GEFFR = GEFFL(L) / (RADIUS(L) * RADIUS(L))
        
        IF (bFullHD) THEN
C***      use full hydrodynamic formula  
          IF (bRUKU) THEN
c            WRITE (hCPR,'(A,3(2X,G15.8))') 'A/V= ', 
c     >        SQRT(A2SUM(L+1))/VELO(L+1)/1.E5, 
c     >        SQRT(A2SUM(L+1))/1.E5, VELO(L+1)
            IF (SQRT(A2SUM(L+1))/VELO(L+1)/1.E5 - 1. < 1.) THEN
              RCONMAX = RADIUS(L+1)
              IMIN = L+1
              EXIT hydrostat
            ELSE
              CALL HYSTHDRUKU(RADIUS(L+1), RADIUS(L), VELO(L+1), VL, 
     >                        RADIUS, GEFFL, A2SUM, ND, RSTAR)
              VELO(L) = VL
            ENDIF
          ELSE 
            CALL SPLINPOX(A2SUML, RADIUS(L), A2SUM, RADIUS, ND, 
     >                                               DFDX=DA2DR)
            DVDR = (GEFFR - 2. * A2SUML/RADIUS(L) + DA2DR) /
     >             ( A2SUML/VELO(L+1) - VELO(L+1) )
            IF (DVDR < 0.) THEN
              RCONMAX = RADIUS(L+1)
              IMIN = L+1
              EXIT hydrostat
            ENDIF
            VELO(L) = DVDR * DR
          ENDIF
        ELSE
C***      uses hydrostatic formula with exp split        
          Hlast = H        ! from last loop, i.e. at L+1
          
          H     = A2SUM(L) / GEFFR / RSTAR

          H = 0.5 * (H + Hlast)
          IF (bRUKU) THEN
            CALL HYSTRUKU (RADIUS(L+1), RADIUS(L), DELTAB, 
     >                   RADIUS, GEFFL, A2SUM, ND, H0, RSTAR)
            BPL = BPL + DELTAB
c          WRITE (0,*) 'BLA: ', DR * (1./H0 - 1./H), DELTAB
          ELSE
            BPL = BPL + DR * (1./H0 - 1./H)
          ENDIF
          
          EXPO = (RADIUS(L)-1.)/H0 - BPL

C***    prevent overflow because of exponential growth of VELO
          IF (EXPO .LT. EXPMAX) THEN
            VELO(L) = VMIN * A2SUM(L)/A2SUM(ND) * EXP(EXPO)/RADIUS(L)**2
            IMIN = L 
            RCONMAX = RADIUS(IMIN)
C*         end the static part one point after VFINAL
            IF (VELO(L+1) >= VFINAL) THEN
              VELO(L) = MIN(10.*VFINAL, VELO(L))
              EXIT hydrostat 
            ENDIF
          ELSE
            VELO(L) = VMIN*A2SUM(L)/A2SUM(ND) * EXP(EXPMAX)/RADIUS(L)**2
            VELO(L) = MIN (VELO(L), VFINAL)
            EXIT hydrostat
          ENDIF
        ENDIF
      ENDDO hydrostat

C***  Define the index of the connection point 
C***  (only if defined by the sonic speed, other wise this will happen later)
C***  - find index interval (ICON-1, ICON) where VELO first time exceeds f*VSOUND 
      IF (bUseSonicConnectionCriterion) THEN
        IF (bFullHD) THEN
C***      For FULLHD option take last good point of integration        
C***      and check if we can still go lower from there
          ICON = IMIN
          IF (ICON < ND) THEN
            VSOUND = SQRT(RGAS * T(ICON) / XMU(ICON)) * 1.E-5   !speed of sound in km/s
            VMINUSAM = VELO(ICON) - fsonic*VSOUND
            DO L=ICON+1, ND
              VSOUND = SQRT(RGAS * T(L) / XMU(L)) * 1.E-5       !speed of sound in km/s
              VMINUSA = VELO(L) - fsonic*VSOUND
              IF (VMINUSA < 0. .OR. L == ND) THEN
                ICON = L
                EXIT
              ELSE 
                VMINUSAM = VMINUSA
              ENDIF
            ENDDO
          ENDIF
          IF (bDEBUG) THEN
            WRITE (0,*) 'DEBUG: ICON, IMIN = ', ICON, IMIN
          ENDIF          
        ELSE
          ICON = ND+1  ! No hydrostatic domain!
          DO L=ND, IMIN, -1
              IF (L .LT. ND) VMINUSA = VMINUSAM
              VSOUND = SQRT(RGAS * T(L) / XMU(L)) * 1.E-5       !speed of sound in km/s
              VMINUSAM = VELO(L) - fsonic*VSOUND
              IF (VMINUSAM > .0) EXIT
              ICON = L
          ENDDO
        ENDIF
      ENDIF

C***  Ensure monotonic increase of the velocity
C***  note: this applies only to the hydrostatic part, i.e. from index IMIN
      IF (.NOT. bForceMonotonic) GOTO 16
   15 Continue
C***  Copy all points to help vectors, omitting those which are non-monotonic      
      LL = 0
      DO L=IMIN, ND
        IF (L .EQ. ND) THEN
           bTestMono = .TRUE.
        ELSE
           bTestMono =  (VELO(L) > VELO(L+1)) 
        ENDIF
        IF (bTestMono) THEN
          LL = LL + 1
          VHELP(LL) = VELO(L)
          RHELP(LL) = RADIUS(L)
        ELSEIF (LEnforcerNeeded == 0) THEN
           LEnforcerNeeded = L
        ENDIF
      ENDDO
      NDHELP = LL

      IF (NDHELP .EQ. ND-IMIN+1) GOTO 16 ! all points are (now) monotonic
C***  If not, then replace all points by their interpolated value      
C***    (trivially, points which had been copied will not change)
      DO L=ND-1, IMIN, -1
        IF (RADIUS(L) < RHELP(1)) THEN
c          CALL SPLINPO (VELO(L), RADIUS(L), VHELP, RHELP, NDHELP)
          CALL SPLINPOX(VELO(L), RADIUS(L), VHELP, RHELP, NDHELP)
        ELSE 
C***       non-monotonic point(s) might include IMIN, thus shorteneing 
C***       the range in which interpolation is possible. 
C***       Linear extrapolation is applied for these point(s) 
           IF (NDHELP >= 2) THEN
             VELO(L) = VHELP(1) + (RADIUS(L)-RHELP(1)) * 
     >               (VHELP(2)-VHELP(1))/(RHELP(2)-RHELP(1))  
           ELSE ! extreme case that only one good point exists
             VELO(L) = VELO(L+1) * 1.001
           ENDIF
        ENDIF
      ENDDO
      GOTO 15 ! Check once more for monotony, just in case

   16 CONTINUE

      RCONMAX = RADIUS(IMIN+1)+0.9*(RADIUS(IMIN)-RADIUS(IMIN+1))

      DO I=ND, 1, -1
        RHELP(ND-I+1) = RADIUS(I)
        VHELP(ND-I+1) = VELO(I)
      ENDDO
      
      VMAX = VFINAL * (1.-BETA2FRACTION)

C***  ITERATIVE SEARCH FOR THE CONNECTION POINT INDEX ICON
C***   WHERE THE DERIVATIVES ARE CONTINUOUS
C***   - GIVEN THAT POINT, THE PARAMETERS OF THE ANALYTIC LAW ARE
C***     RE-ADJUSTED SUCH THAT THE VELOCITIES ARE CONTINUOUS.
C***  THE SEARCH PROCEEDS INWARDS.
      
      IF (bUseSonicConnectionCriterion) THEN
C***    Sonic point is between indices (ICON,ICON-1) 
        bConPointConv = .TRUE.
        MAXIT = 1
        IF (bFullHD) THEN
          VCON = VELO(ICON)
          RCON = RADIUS(ICON)
        ELSEIF (ICON .LE. ND) THEN
C*      linear interpolation for RCON
          P = VMINUSA / (VMINUSA - VMINUSAM) 
          RCON = P * RADIUS(ICON-1) + (1.-P) * RADIUS(ICON)
          VA  = VELO(ICON)   - VMINUSA
          VAM = VELO(ICON-1) - VMINUSAM
          VCON = VA + (RADIUS(ICON-1) - RCON) *
     >      (RADIUS(ICON-1)-RADIUS(ICON))/(VAM - VA) 

c          write (*,*) 'VCON, VA, ICON=', VCON, VA, ICON
cc          CALL SPLINPO (VCON, RCON, 
cc     >     VELO(ICON-1:ND), RADIUS(ICON-1:ND), ND-ICON+2) 
        ELSE
           RCON = 1.
           Vcon = VMIN
        ENDIF

      ELSE
        bConPointConv = .FALSE.
        RCON = 1.
        MAXIT = MAXCONIT
        Vcon = VMIN
      ENDIF
            
      RCONSTART = RCON
      RCONOLD = RCON
      NSameRcon = 0
c      WRITE (0,*) 'VELTHING: Searchin RCON...'

      conpointloop: DO K=1, MAXIT !------------------------------------------

        RCONOLD = RCON
        IF (bDEBUG) WRITE (hCPR,*) ' conpointloop IT, =', K

        IF (.NOT. bUseSonicConnectionCriterion) THEN
          !Obtain VCON via spline interpolation on integration result
          IF (RCON > RADIUS(IMIN)) THEN
            WRITE (hCPR,*) "WARNING RCON=", RCON, " > ", RADIUS(IMIN)
          ENDIF
          IF (RCON < RHELP(1)) THEN
            WRITE (hCPR,*) ' RCON limited to minimum value: ', 
     >                             RHELP(1)            
            WRITE (hCPR,*) ' RHELP-: Rcon=', Rcon,' Rhelp(1)=', RHELP(1)
            RCON = RHELP(1)
            Vcon = VMIN
          ELSEIF (RCON >= RHELP(ND-IMIN+1)) THEN
            WRITE (hCPR,*) ' RCON limited to maximum value: ', 
     >                             RHELP(ND-IMIN+1)
C            WRITE (hCPR,*) ' RHELP+: Rcon=',Rcon,' Rhelp(+)=', 
C     >                             RHELP(ND-IMIN)
            RCON = RHELP(ND-IMIN+1)
            Vcon = VHELP(ND-IMIN+1)
          ELSE
            CALL SPLINPOX(Vcon,RCON,VHELP,RHELP,ND-IMIN+1,DFDX=DVDRcon)
          ENDIF
          !Failsafe-Check for Vcon (in case of bad SPLINPOX result)
          IF (Vcon < VMIN) THEN
            Vcon = VMIN
          ENDIF
        ENDIF

C***    CALCULATION OF VELOCITY-FIELD PARAMETERS (ANALYTIC LAW)
        IF (bUseSonicConnectionCriterion .AND. bNewSonic 
     >                                      .AND. BETA /= 0.) THEN
C***      In this branch v(Rmax) = v_inf is not guaranteed
C***      Possibly add iteration loop to ensure v_inf again
C***      TODO: Extend this for 2-beta-laws
          CALL SPLINPOX(Vdummy, RCON, VELO, RADIUS, ND, DFDX=DVDRcon)
          WRITE (0,*) 'DEBUG:  DVDRcon = ', DVDRcon
          VPAR1 = VFINAL - Vcon
          DO
          IF (BETA > 0) THEN
            VOFF = Vcon - (DVDRcon/BETA/VPAR1)**(BETA/(BETA-1.))
            VPAR2 = 1./(1.- ((Vcon - VOFF)/(VPAR1 - VOFF))**(1./BETA) ) 
     >                                   - RCON
          ELSE
            VOFF = Vcon - (RCON*DVDRcon/BETA/VPAR1)**(BETA/(BETA-1.))
            VPAR2 = RCON*(1.- ((Vcon - VOFF)/(VPAR1 - VOFF))**(1./BETA))
          ENDIF          
C***        Iteration to ensure VMAX 
c            WRITE (0,*) 'DEBUG: VINF, WRVEL(RMAX)', VFINAL, VELOBETA(RMAX)
            IF (ABS(VFINAL - VELOBETA(RMAX)) < 1.) EXIT
            IF (VFINAL - VELOBETA(RMAX) > 0) THEN
              VPAR1 = 1.1 * VPAR1
            ELSE 
              VPAR1 = 0.9 * VPAR1
            ENDIF
          ENDDO
        ELSE
          CALL INITVELBETAPAR(RMAX, RCON, Vcon, .TRUE.)
        ENDIF
c        CALL INITVELBETAPAR(RMAX, RCON, Vcon, .FALSE.)
c        WRITE (0,*) 'VCON = ', VCON, VELOBETA(RCON)

        IF (bUseSonicConnectionCriterion) THEN
          !no connection point loop needed if sonic criterion is used
          IF (BDEBUG) 
     >       CALL PLOTVGRAD(ND ,RADIUS, VELO, 'DEBUG', 1, 1., 1.)
          EXIT conpointloop
        ENDIF

C***    Define minimum radius for the BETA law       
        IF (BETA < 0.) THEN
          RBETAMIN = VPAR2 + EPSRCONMIN
        ELSE
          RBETAMIN = 1.-VPAR2+EPSRCONMIN  !Beta law cannot be calculated inside of this point
        ENDIF
        
C***    Initialize RCON and ICON in the first run of this loop        
        IF (K <= 1) THEN
          RIN = MAX(RBETAMIN,1.)
          IF (DELTAGRTHIN(RIN,ND,RADIUS,VELO) < .0) THEN
C***        Inner gradient is already larger than outer gradient at innermost point:
C***        => no hydrostatic domain          
            RCON=RIN
            ICON = ND
            IF (BDEBUG) 
     >         CALL PLOTVGRAD(ND ,RADIUS, VELO, 'DEBUG', 1, 1., 1.)
            EXIT conpointloop
          ENDIF

          IF ((IMIN >= ND) .OR. RCONMAX > RMAX) THEN          
            RCON=RMAX
            WRITE (hCPR,'(A)') 
     >          'VELTHIN: NO BETA-LAW DOMAIN ENCOUNTERED'
            RETURN
          ENDIF
        ENDIF

C***    LOWER BOUNDARY FOR RCON         
        RCONMIN=MAX(RCONSTART, RBETAMIN) !Lower guess for RCON
        !Radius of maximum of gradient in beta laws
        RGRADMAX = 1. + (BETA-1.)/2. - VPAR2    

        CALL PLOTVGRAD(ND ,RADIUS, VELO, 'DEBUG', 1, RCONMIN, RGRADMAX)

C***    Loop: INCREASE RCONMIN STEPWISE TO ENSURE DELTAGRTHIN(RCONMIN) .GE. 0.
C***    DELTAGRTHIN returns GRAD_BETA - GRAD_HYST
C***    In the outer wind this value should be below zero (i.e. GRAD_HYST > GRAD_BETA)
C***    However, for BETA > 1, there can also be an inner range where 
C***      GRAD_HYST > GRAD_BETA, which must be avoided to find the correct solution.
C***    We therefore increase RCONMIN until we get out of this region.
        GRADRCONMIN = DELTAGRTHIN(RCONMIN,ND,RADIUS,VELO)
        DO WHILE (GRADRCONMIN < .0) 
          DXX = 0.1 * HSCALE        
          RCONMIN = RCONMIN + DXX
          GRADRCONMIN = DELTAGRTHIN(RCONMIN,ND,RADIUS,VELO)
          IF (RCONMIN > RGRADMAX) THEN
            WRITE (hCPR,'(A)') '*** VELTHIN: NO CONNECTION POINT FOUND'
            WRITE (hCPR,'(A)') '*** using sonic connection criterion...'
            bFailsafe = .TRUE.
            GOTO 20            
C            STOP               '*** ERROR IN VELTHIN'
          ENDIF
        ENDDO

C***    RCONMIN has now been increased such that we avoid the second solution for 
C***    cases with BETA > 1
C***    RCONMAX denotes the maximum radius where we have a hydrostatic solution
C***    IF RCONMAX > RCONMIN we are fine, but for RCONMIN > RCONMAX we cannot find 
C***    a solution
        
        RCONMAXguess = MAX(1., RBETAMIN)           !RCONMAXguess must be initialized
        Lguess = ND
        conmaxmin: DO L=ND-1, 1, -1
          IF (RADIUS(L) > RCONMIN) THEN
            Lguess = L
            RCONMAXguess = RADIUS(L)
            EXIT conmaxmin
          ENDIF
        ENDDO conmaxmin


C***    Try to find a good range for RCONMAX    
C       DELTAGRTHIN = BETALAWGRADIENT - HYSTGRADIENT
C***    We loop outwards through the grid until we reach the end of the 
C***     tabulated (hydrostatic) velocity field (index IMIN)
        conmaxmax: DO L=Lguess, IMIN, -1
          RCONMAXguess = RADIUS(L)
          GRADRCONMAX = DELTAGRTHIN(RCONMAXguess,ND,RADIUS,VELO)
          IF (GRADRCONMAX < .0) THEN
C***        Hydrostatic gradient is (again) larger than BETA law gradient
C***        => we have reached the outer limit for RCON
            EXIT conmaxmax
          ENDIF
          IF (L == IMIN) THEN
C***        We have never reached a part where we can find a connection point          
            RCONMAXguess = RCONMAX
          ENDIF
        ENDDO conmaxmax
        GRADRCONMAX = DELTAGRTHIN(RCONMAXguess,ND,RADIUS,VELO)

        
C***    Update vgrad (debug) plot with final boundaries
        CALL PLOTVGRAD(ND ,RADIUS, VELO, 'DEBUG', 1, 
     >                 RCONMIN, RCONMAXguess)
        
C***    Check if RCON search boundaries exclude each other
C***      or if there is no sign change between both boundaries
        IF (RCONMIN > RCONMAXguess .OR. 
     >           (GRADRCONMAX > 0. .AND. GRADRCONMIN > 0.)) THEN
C***         => no hydrostatic domain        
          RCON = 1.
          ICON = ND
          EXIT conpointloop
        ELSEIF (GRADRCONMAX < 0. .AND. GRADRCONMIN < 0.) THEN
C***         => no wind domain        
          RCON = 1.1 * RMAX
          ICON = 1
CCC       TODO: Check if ICON = 0 would be possible and better
          EXIT conpointloop
        ENDIF                        


C***    Now find new RCON where DELTAGRTHIN is zero:        
        !------------- Regula Falsi -------------------------
 
        bRFbi=.TRUE.
        RFA=RCONMIN
        RFB=RCONMAXguess
        RFFA = DELTAGRTHIN(RFA, ND, RADIUS, VELO)
        RFFB = DELTAGRTHIN(RFB, ND, RADIUS, VELO)
        IF (RFFA == 0.) THEN
          RCON = RFA
        ELSEIF (RFFB == 0.) THEN
          RCON = RFB
        ELSEIF (RFFA*RFFB > 0) THEN
          WRITE (hCPR,*) '*** INVALID ARGUMENTS FOR REGULA FALSI ***'
          WRITE (hCPR,FMT=10) RFA, RFB
   10     FORMAT ('INTERVAL:     A  =', E15.5, 5X, '  B  =', E15.5)
          WRITE (hCPR,FMT=11) RFFA, RFFB
   11     FORMAT ('FUNCTION:   F(A) =', E15.5, 5X, 'F(B) =', E15.5)
          WRITE (hCPR,'(A)') '*** using sonic connection criterion...'
          bFailsafe = .TRUE.
          GOTO 20            
C          CALL TRBK
C          STOP 'VELTHIN: ERROR IN REGULA FALSI'
        ELSE
          rfloop: DO !-  -  -  -  -  -  -  -  -  -  -
            RFD = RFA - RFB
            IF (ABS(RFD) < RFACC) THEN
              RFX = 0.5 * (RFA + RFB)
              EXIT rfloop
            ENDIF
            bRFbi = .NOT. bRFbi
            IF (bRFBi) THEN
              RFX = RFA - RFFA*RFD/(RFFA-RFFB)
            ELSE
              RFX = 0.5 * (RFA + RFB)
            ENDIF
            RFFX = DELTAGRTHIN(RFX, ND, RADIUS, VELO)
            IF (RFFX == 0.) EXIT rfloop
            IF (RFFX * RFFA > 0.) THEN
              RFA = RFX
              RFFA = RFFX
            ELSE
              RFB = RFX
              RFFB = RFFX
            ENDIF
          ENDDO rfloop !-  -  -  -  -  -  -  -  -  -  -
          RCON = RFX
        ENDIF
          
        !-------------- end of Regula Falsi -----------------

CC        WRITE (0,*) ' RCON found = ', RCON
        IF (RCON > RCONMAX) THEN
          WRITE (0,*) 'Error: RCON larger than RCONMAX '
          STOP "*** FATAL ERROR IN VELTHIN"        
        ENDIF
        RCONHIST(K) = RCON
C        WRITE (hCPR,*) ' THIN-IT:', K, ' RCON=', RCON

        IF (ABS(RCON-RCONOLD) <= RCONACC) THEN
          bConPointConv = .TRUE.
          EXIT conpointloop
        ELSEIF (K >= 2) THEN
          !Check if solution is oscillating
          NSameRcon = 0
          DO I=K-1, 1, -1
            IF (ABS(RCON-RCONHIST(I)) <= RCONACC) THEN
              NSameRcon = NSameRcon + 1
              !Cancel if more than 10 oscillations
              IF (NSameRcon > 10) THEN
                DO J=K-1, 1, -1
                  IF (ABS(RCON-RCONHIST(J)) <= RCONACC) THEN
                    LOszil = K - J
                    EXIT
                  ENDIF
                ENDDO
                !RCON = (RCONOLD + RCON) / 2.
                RCON = 0.
                DO J=K, K-LOszil+1, -1
                  RCON = RCON + RCONHIST(J)
                ENDDO
                RCON = RCON / FLOAT(LOszil)
                bConPointConv = .TRUE.
                BAD = ABS( MAXVAL(RCONHIST(K-LOszil+1:K))
     >                   - MINVAL(RCONHIST(K-LOszil+1:K)) )
     >            / ABS( 0.5 * ( MAXVAL(RCONHIST(K-LOszil+1:K))
     >                   + MINVAL(RCONHIST(K-LOszil+1:K)) ) )
C                CALL PLOTVGRAD(ND ,RADIUS, VELO, 'DEBUG', 1337)
                WRITE (hCPR,'(A,A,I2,A,F12.8,A)')
     >              'VELTHIN: WARNING - OSCILLATING SOLUTION',
     >              ' (LENGTH: ', LOszil, '  BADNESS: ', BAD, ')'
C                WRITE (hCPR,'(A,F12.8)')
C     >              '  averaged RCON: ', RCON
                EXIT conpointloop
              ENDIF
            ENDIF
          ENDDO
        ENDIF

      ENDDO conpointloop !------------------------------------------

      IMAXBETA = 0
      
      IF (ICON /= ND) THEN
        IF (.NOT. bConPointConv) THEN
          WRITE (hCPR,'(A)') 'VELTHIN: CONVERGENCE PROBLEMS'
          WRITE (hCPR,*) ABS(RCON-RCONOLD), ' ACC:',RCONACC
          STOP 'ERROR'
        ENDIF

        IF (.NOT. bUseSonicConnectionCriterion) THEN
          iconloop: DO I=ND, 1, -1  
            IF (RADIUS(I) >= RCON) THEN
              ICON = I
              EXIT iconloop
            ENDIF
          ENDDO iconloop
        ENDIF

C        WRITE (hCPR,FMT='(A,I2)') " ICON=", ICON
        IMAXBETA = ICON-1
      ELSE
        !no hydrostatic domain: use beta-law only
        IMAXBETA = ND
      ENDIF
      
c      IF (BETA > 0.) THEN
        !Calculate beta law velocities for comparison
        DO I=1, ND
          RM = RADIUS(I)
          IF (RM <= RCON) THEN
            IMAXBETA = I-1
            EXIT
          ELSE
            VBETA(I) = VELOBETA(RM)
            IF (bDEBUG) THEN
               WRITE (hCPR,FMT='(A,I2,A,F12.6)') 
     >           "DEBUG: VELOBETA(",I,")=",VBETA(I)
            ENDIF
          ENDIF
        ENDDO
  
c      ENDIF


C***  DEFINITION OF VELO(I) - OVERWRITING THE OUTER PART WITH THE
C***      BETA OR SQRT(ALOG(R)) LAW
      CRITERION = "HYSTINT "
C      WRITE (hCPR,*) " IMAXBETA=", IMAXBETA
      DO I=1, IMAXBETA
        IF (ABS(BETA) <= .0) THEN
          VELO(I) = VPAR1*SQRT(ALOG(RADIUS(I)/VPAR2))
          CRITERION(I) = "SQRT    "
        ELSE
          VELO(I) = VBETA(I)
          IF (BETA2FRACTION > .0) THEN 
            CRITERION(I) = "2BETA   "
          ELSE
            CRITERION(I) = "BETA    "
          ENDIF
        ENDIF
      ENDDO

C***  Ensure monotonic field at connection point     
C***  (This can be an issue in case of oscillating solutions)
      IF (IMAXBETA > 1 .AND. IMAXBETA < ND 
     >      .AND. VELO(IMAXBETA) < VELO(IMAXBETA+1)) THEN
        INME = IMAXBETA+1        
        nmcloop: DO J=IMAXBETA-1, 1, -1
          IF (VELO(J) > VELO(IMAXBETA+1)) THEN
            INMS = J
            EXIT nmcloop
          ENDIF
        ENDDO nmcloop
        IF (INMS > 2) THEN
          RIP(1) = RADIUS(INME)
          RIP(2) = RADIUS(INMS)
          RIP(3) = RADIUS(INMS-1)
          RIP(4) = RADIUS(INMS-2)          
          VIP(1) = VELO(INME)
          VIP(2) = VELO(INMS)
          VIP(3) = VELO(INMS-1)
          VIP(4) = VELO(INMS-2)          
        ELSE
          WRITE (hCPR,*) 'NONMONOTONIC BETA-PART IN WIND VELOCITY'
          STOP 'FATAL ERROR IN VELTHIN'
        ENDIF
        DO J=INMS+1, INME, 1
          CALL SPLINPOX(VELO(J), RADIUS(J), VIP, RIP, 4)
        ENDDO        
      ENDIF      
      
C***  Smooth velocity field if enforced (re-uses VHELP)
      IF (bSmoothVelo) THEN
        DO L=2, ND-1
          VHELP(L) = 0.25 * 
     >      (LOG10(VELO(L-1)) + 2. * LOG10(VELO(L)) + LOG10(VELO(L+1)))        
        ENDDO
        DO L=2, ND-1
          VELO(L) = 10**VHELP(L)
        ENDDO
      ENDIF

      !Debug Output
      IF (bPrintVelo) THEN
        DO I=1, ND
          WRITE (hCPR,FMT='(A,G10.3,A,I2,A,F12.6,A,A)') 
     >      "R-1=", RADIUS(I)-1.,
     >      '  VELTHIN-VELO(',I,')=',VELO(I),'  ',CRITERION(I)
        ENDDO
        WRITE (hCPR,FMT='(A,F10.5)') ' *** VELTHIN: RCON = ', RCON
        IF (BETA2FRACTION > .0) THEN 
          WRITE (hCPR,FMT='(A,4(2X,G15.8))') 'BETA VPARs ', 
     >          VPAR1, VPAR2, VPAR1_2, VPAR2_2
        ELSE 
          WRITE (hCPR,FMT='(A,2(2X,G15.8))') 'BETA VPARs ', VPAR1, VPAR2
        ENDIF
      ENDIF
C      CALL SPLINPOX(Vcon,RCON,VELO,RADIUS,ND)
      IF (IMAXBETA < 1) THEN
        WRITE (hCPR,'(A)') 'VELTHIN: NO BETA-LAW DOMAIN ENCOUNTERED'
      ENDIF
      IF (IMAXBETA >= ND-1) THEN
        WRITE (hCPR,'(A)') 'VELTHIN: NO HYDROSTATIC DOMAIN ENCOUNTERED'
      ENDIF
      IF (bForceMonotonic .AND. LEnforcerNeeded > ICON) THEN
        WRITE (hCPR,'(A)') 'VELTHIN: MONOTONIC STRATIFICATION ENFORCED'
      ENDIF

      RETURN
      END


