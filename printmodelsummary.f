      SUBROUTINE PRINTMODELSUMMARY (
     >                              MODHEAD,    !Model start header string
     >                              ND,         !number of depth points
     >                              TEFF,       !effective temperature
     >                              T,          !(Electron) Temperature per depth point
     >                              TAUROSS,    !(complete) Rosseland optical depth per depth point
     >                              RSTAR,      !Stellar radius in cm (inner boundary)
     >                              LOGL,       !logL
     >                              LOGMDOT,    !log Mdot in Msun/yr (from MODEL or calculated in REMOST)
     >                              RTRANS,     !transformed radius
     >                              VELO,       !velocity field
     >                              VDOP,       !Doppler velocity
     >                              DENSCON,    !Density contrast (D) per depth point
     >                              FILLFAC,    !Volume filling factor (1/D) per depth point
     >                              GLOG,       !log g_grav
     >                              GEFFLOG,    !log g_eff
     >                              GEDD,       !Eddington Gamma
     >                              RADIUS,     !Radius grid
     >                              RCON,       !Rcon in Rstar
     >                              INCRIT,     !Rgrid depth point criterion per depth point
     >                              VELOCRIT,   !Vgrid calculation criterion per depth point
     >                              VTURB,      !Turbulence velocity in km/s
     >                              bThinWind,  !determines if HYDROSTATIC INTEGRATION was used
     >                              bHYDROSOLVE,!determines if HYDROSOLVE was used
     >                              WORKRATIO,  !current work ratio of the model
     >                              RNE,        !relative electron density
     >                              ENTOT,      !Total number of electrons per depth points
     >                              XMU,        !mean particle mass (mu) per depth point
     >                              MSTAR,      !stellar mass
     >                              GRADI,      !velocity field gradient
     >                              NATOM,      !Number of different atoms considered
     >                              ABXYZ,      !Abundances (X, Y, Z) - fractions per number
     >                              ATMASS,     !mean atomic mass
     >                              ELEMENT,    !Element corresponding to atom integer index
     >                              SYMBOL,     !same as above but only symbol instead of full name
     >                              bFULLHYDROSTAT,
     >                              GEDDRAD,
     >                              ARAD,
     >                              APRESS,
     >                              AGRAV,
     >                              GAMMARADMEAN,
     >                              QIONMEAN,
     >                              bFixGEFF
     >           )
C*******************************************************************************
C***  PRINT A SUMMARY OF THE FINAL STELLAR PARAMETERS 
C***  CALLED FROM STEAL AFTER MODEL IS CONVERGED OR STHLP IS SET
C*******************************************************************************

      IMPLICIT NONE
      INCLUDE 'interfacebib.inc'

      INTEGER, INTENT(IN) :: ND, NATOM
      REAL, DIMENSION(ND), INTENT(IN) :: DENSCON, FILLFAC, RADIUS, VELO,
     >                                   GRADI, RNE, ENTOT, T, TAUROSS
      REAL, DIMENSION(ND-1), INTENT(IN) :: ARAD, APRESS, AGRAV
      REAL, DIMENSION(ND), INTENT(INOUT) :: XMU
      REAL, DIMENSION(ND) :: AMACH, VSCRATCH
      REAL, DIMENSION(ND-1) :: RI
      REAL, DIMENSION(NATOM), INTENT(IN) :: ABXYZ, ATMASS 
      REAL, INTENT(IN) :: TEFF, RSTAR, VDOP, RCON, LOGMDOT, VTURB, 
     >                    WORKRATIO, GEDDRAD
      REAL, INTENT(INOUT) :: LOGL, RTRANS, GEDD, GLOG, GEFFLOG, MSTAR
      
      CHARACTER(2), DIMENSION(NATOM), INTENT(IN) :: SYMBOL
      CHARACTER(8), DIMENSION(ND), INTENT(IN) :: INCRIT, VELOCRIT
      CHARACTER(10), DIMENSION(NATOM), INTENT(IN) :: ELEMENT
      CHARACTER(100), INTENT(IN) :: MODHEAD

      LOGICAL, INTENT(IN) :: bThinWind, bHYDROSOLVE, 
     >                       bFULLHYDROSTAT, bFixGEFF

      INTEGER, EXTERNAL :: IDX

      INTEGER :: NA, L, Lcand, LCON
      REAL :: RSTARSU, ATMEAN, XELEM, VFINAL, VMIN, FM, VESC, VESCFULL,
     >        Dout, T1, T13, T23, R1, R13, R23, TAU1, TAU13, TAU23,
     >        Tindex, Rtindex, RL1, RL2, LOGQ, gammaIon, MDOTNL, GEDDe,
     >        LOGMDOTTRANS, LOGDMOM, Rsonic, Vsonic, XMUsonic, Tsonic,
     >        GAMMARAD, LOGLINF, QIONMEAN, GAMMARADMEAN, fdummy, q,
     >        G23LOG, GEFF23LOG, ASTAR, ETA, Rscrit, Vscrit, VTH, VTEST,
     >        TAUsonic, TAUscrit, VMICND, LOGLINFTEST, LWINDSUB
      CHARACTER(2) :: VCRIT
      CHARACTER(100) :: CTABLINE

      REAL, PARAMETER :: PI4 = 12.5663706144    !PI4 = 4*PI
      REAL, PARAMETER :: AMU = 1.66E-24         !Atomic mass unit (gramm)     
      REAL, PARAMETER :: RGAS = 8.3145E7        !Gas Constant (CGS)
      REAL, PARAMETER :: RSUN = 6.96E10         !SOLAR RADIUS ( CM )
      REAL, PARAMETER :: TEFFSUN = 5780.        !SOLAR EFFECTIVE TEMPERATURE
      REAL, PARAMETER :: MSUN = 1.989E33        !SOLAR MASS (gramm)
      REAL, PARAMETER :: XLSUN = 3.85E33        !Solar Luminosity (CGS-Units)
      REAL, PARAMETER :: LYR = 7.4988           !LOG SECONDS in a year
      REAL, PARAMETER :: GCONST = 6.670E-8      !GRAVITATION CONSTANT (CGS UNITS)
      REAL, PARAMETER :: CLIGHT = 2.99792458E10     !Speed of Light in cm/s
      REAL, PARAMETER :: STEBOL = 5.6705E-5     !STEFAN-BOLTZMANN CONSTANT (CGS-UNITS)
      REAL, PARAMETER :: XMSUNPYR = 6.303E25    !1 SOLAR MASS / YR in CGS
      REAL, PARAMETER :: ASTARREF = 5E11        !Reference value for Astar (Chevalier & Li 1999)

      !File and channel handles (=KANAL)
      INTEGER, PARAMETER :: hOUT = 6            !write to wruniqX.out (stdout)
      INTEGER, PARAMETER :: hCPR = 0            !write to wruniqX.cpr (stderr)
      INTEGER, PARAMETER :: hSUMMARY = 77       !write to modinfo.kasdefs

      !Link data to identify program version
      CHARACTER(30) :: LINK_DATE
      CHARACTER(10) :: LINK_USER
      CHARACTER(60) :: LINK_HOST
      COMMON / COM_LINKINFO / LINK_DATE, LINK_USER, LINK_HOST
  
      OPEN (hSUMMARY, FILE='modinfo.kasdefs')
      WRITE(hSUMMARY,'(A)')  '*Model parameter summary file'
      WRITE(hSUMMARY,'(5A)') '* created by PoWR (version from ', 
     >                      LINK_USER(:IDX(LINK_USER)), ' on ', 
     >                      LINK_DATE(:IDX(LINK_DATE)) ,')'

      ATMEAN = .0
      DO NA=1, NATOM
        ATMEAN = ATMEAN + ABXYZ(NA)*ATMASS(NA)
      ENDDO
      DO L=1, ND
        XMU(L) = ATMEAN / (1. + RNE(L))
      ENDDO      
      
      RSTARSU = RSTAR / RSUN
      IF (LOGL <= 0.) THEN
        !Calculate log L from Rstar and Teff if not saved
        !(This is normal if called from STEAL)
        LOGL = ALOG10( RSTARSU**2 * (TEFF/TEFFSUN)**4)
      ENDIF
      VFINAL = VELO(1)
      VMIN = VELO(ND)

      !For fixed g_eff, ensure correct mass calculation
      IF (bFixGEFF) THEN
        DO L=1, ND-1
          RI(L) = 0.5 * (RADIUS(L) + RADIUS(L+1))
        ENDDO
        !This updates MSTAR, GLOG, GEDD, GAMMARADMEAN
        CALL CALCMASSFROMGEFF(bFixGEFF, bFULLHYDROSTAT,
     >                        MSTAR, GLOG, GEFFLOG, 
     >                        ARAD, APRESS, AGRAV, RADIUS, ND,
     >                        RSTAR, RCON, RI, TAUROSS,
     >                        GAMMARADMEAN, GEDD,
     >                        LOGL, QIONMEAN, fdummy)
      ENDIF
      
      
      !Find sonic point parameters
      Lcand = 0
      DO L=1, ND
        AMACH(L) = SQRT( RGAS * T(L) / XMU(L) ) / 1.E5      !a_mach in km/s
        VSCRATCH(L) = VELO(L) - AMACH(L)        
        IF ((VSCRATCH(L) < 0) .AND. (Lcand == 0)) THEN
          Lcand = L
        ENDIF
      ENDDO
      IF (Lcand > 0) THEN
        CALL SPLINPOX(Rsonic,0.,RADIUS,VSCRATCH,ND,.FALSE.,Lcand)
        CALL SPLINPOX(Vsonic,Rsonic,AMACH,RADIUS,ND)
        CALL SPLINPOX(Tsonic,Rsonic,T,RADIUS,ND)
        CALL SPLINPOX(XMUsonic,Rsonic,XMU,RADIUS,ND)
        CALL SPLINPOX(gammaIon,Rsonic,RNE,RADIUS,ND)
      ENDIF
            
        
      Dout = DENSCON(1)
      RTRANS = 
     >  (VFINAL/2500. *1.E-4 / 10.**LOGMDOT /SQRT(Dout))**(2./3.)
     >      * RSTARSU


      WRITE(hSUMMARY,'(A,F10.2)') '\VAR TEFF    = ', TEFF
      WRITE(hSUMMARY,'(A,1PG14.4)') '\VAR RSTAR   = ', RSTAR   !in cm
      WRITE(hSUMMARY,'(A,F8.3)') '\VAR RSTARSU = ', RSTARSU  !in Rsun
      WRITE(hSUMMARY,'(A,F5.3)') '\VAR LOGL    = ', LOGL     !log L/Lsun
      WRITE(hSUMMARY,'(A,F12.4)') '\VAR RTRANS  = ', RTRANS   !Rtrans in solar radii
      WRITE(hSUMMARY,'(A,F9.3)') '\VAR VFINAL  = ', VFINAL   !v_inf in km/s
      WRITE(hSUMMARY,'(A,F9.3)') '\VAR VMIN    = ', VMIN     !v_min in km/s
      IF (Lcand > 0) THEN
        WRITE(hSUMMARY,'(A,F9.3)') '\VAR VSONIC  = ', Vsonic   !v_sonic in km/s
        WRITE(hSUMMARY,'(A,F9.3)') '\VAR RSONIC  = ', Rsonic   !Rsonic in Rstar
        WRITE(hSUMMARY,'(A,F9.3)') '\VAR RSONICSU = ', Rsonic * RSTARSU  !Rsonic in Rsun
      ENDIF
      WRITE(hSUMMARY,'(A,F7.3)') '\VAR LOGMDOT = ', LOGMDOT  !log Mdot/(Msol/yr)
      WRITE(hSUMMARY,'(A,F6.2)') '\VAR DENSCON = ', 
     >                                 MAXVAL(DENSCON(1:ND)) !maximum D value
      WRITE(hSUMMARY,'(A,F9.3)') '\VAR VDOP    = ', VDOP
      WRITE(hSUMMARY,'(A,F8.2)') '\VAR RMAX    = ', RADIUS(1)        !in Rstar (als Kontrolle fuer RMAX)
      IF (RCON > 0.) THEN
        WRITE(hSUMMARY,'(A,F10.5)') '\VAR RCON    = ', RCON     !in Rstar
      ENDIF
      
      ETA = 10**(LOGMDOT) * XMSUNPYR * VFINAL * 1.E5 * CLIGHT 
     >        / ( 10**(LOGL) * XLSUN )
      WRITE(hSUMMARY,'(A,F7.3)') '\VAR ETA   = ', ETA        !Eta = Mdot vinf c / L
      ASTAR = 10**(LOGMDOT) * XMSUNPYR / PI4 / (VFINAL*1.E5) / ASTARREF
      WRITE(hSUMMARY,'(A,F7.3)') '\VAR ASTAR = ', ASTAR      !Astar = Mdot / (4 pi vinf) * 5.E11

      !Mass, log g, log g_eff, Eddington Gamma
      q = QIONMEAN
      GAMMARAD = GAMMARADMEAN      
      IF (GAMMARAD == 0. .AND. GEDDRAD > 0.) THEN
        !Use value from last taumax iteration as fallback if available
        GAMMARAD = GEDDRAD
      ENDIF

      IF (.NOT. bFULLHYDROSTAT) THEN
        GEDDe = GEDD
      ELSE
        GEDDe = -99.
      ENDIF

      !Eddington Gamma_e
      IF (.NOT. bFixGEFF .OR. GEFFLOG < 0.) THEN
        IF (GEDDe < 0.) THEN
          GEDDe = 10.**(-4.51) * q * 10**(LOGL) / MSTAR        
        ENDIF
        IF (bFULLHYDROSTAT) THEN
          GEFFLOG = ALOG10( (10.**GLOG) * (1. - GAMMARAD) ) 
        ELSE
          GEFFLOG = ALOG10( (10.**GLOG) * (1. - GEDDe) ) 
        ENDIF
      ELSEIF (GLOG < 0.) THEN
        IF (GEDDe < 0.) THEN
          GEDDe = 10.**(-4.51) * q * 10**(LOGL) / MSTAR
        ENDIF
        IF (bFULLHYDROSTAT) THEN
          GLOG = ALOG10( (10.**GEFFLOG) / (1. - GAMMARAD) ) 
        ELSE
          GLOG = ALOG10( (10.**GEFFLOG) / (1. - GEDDe) ) 
        ENDIF
      ELSE
        !GEFF and G seem to be fixed => Calculate GEDDe and new mass
        IF (GEDD < 0) THEN
          GEDD = 1. - 10.**( GEFFLOG - GLOG )
        ENDIF
        IF (.NOT. bFULLHYDROSTAT) THEN
          MSTAR = 10.**( GEFFLOG - 4.4371 ) * RSTARSU**(2.)
     >                   + 10.**(-4.51) * q * 10**(LOGL)
          GEDDe = GEDD
        ELSE
          MSTAR = 10.**(GLOG) * RSTAR**2 / GCONST / MSUN
          GEDDe = 10.**(-4.51) * q * 10**(LOGL) / MSTAR
        ENDIF
      ENDIF

      IF (GEDDe < 1.) THEN
        VESC = SQRT(2*GCONST*MSTAR*MSUN * (1. - GEDDe) / RSTAR) / 1.E5
      ELSE 
        VESC = -99.
      ENDIF
      IF (GAMMARAD < 1.) THEN
        VESCFULL = SQRT(2*GCONST*MSTAR*MSUN*(1.-GAMMARAD)/RSTAR) / 1.E5
      ELSE 
        VESCFULL = -99.
      ENDIF
      LWINDSUB = 10**(LOGMDOT) * XMSUNPYR/XLSUN * 1.E10 *
     >                    ( VFINAL**2 + VESC**2 ) / 2.
      IF (LWINDSUB < 10**LOGL) THEN
        LOGLINF = LOG10( 10**(LOGL) 
     >                 - LWINDSUB )
      ELSE
        WRITE(hCPR,'(A)')
     >    '*** WARNING: WIND CONSUMES MORE LUMINOSITY THAN PROVIDED'
        WRITE(hOUT,'(A)')
     >    '*** WARNING: WIND CONSUMES MORE LUMINOSITY THAN PROVIDED'
        WRITE(hCPR,'(A,3(1X,F8.3))') 'L, Lsub, Linftest ',
     >                      LOGL, LOG10(LWINDSUB)
      ENDIF

      WRITE(hSUMMARY,'(A,F8.3)') '\VAR MSTAR   = ', MSTAR
      WRITE(hSUMMARY,'(A,F8.4)') '\VAR GLOG    = ', GLOG
      IF (GEFFLOG > 0.) THEN
        WRITE(hSUMMARY,'(A,F8.4)') '\VAR GEFFLOG = ', GEFFLOG
      ENDIF
      WRITE(hSUMMARY,'(A,F6.3)') '\VAR GEDD    = ', GEDDe     !Eddington Gamma_e (Thomson only)
      WRITE(hSUMMARY,'(A,F6.3)') '\VAR GAMMARAD = ', GAMMARAD !Eddington Gamma with full a_rad
      WRITE(hSUMMARY,'(A,F9.3)') '\VAR VTURB   = ', VTURB     !Turbulence velocity
      VMICND = VTURB * SQRT(2.)
      WRITE(hSUMMARY,'(A,F9.3)') '\VAR VMIC   = ', VMICND     !literature microturbulence value
      IF (VESC >= 0.) THEN
        WRITE(hSUMMARY,'(A,F9.3)') '\VAR VESC    = ', VESC    !Escape velocity at Rstar
      ELSE 
        WRITE(hSUMMARY,'(A)')    '\VAR VESC    = invalid'     !supereddington, no escape velocity at Rstar
      ENDIF
      IF (VESCFULL >= 0.) THEN
        WRITE(hSUMMARY,'(A,F9.3)') '\VAR VESCFULL = ', VESCFULL !Escape velocity at Rstar w/ full GAMMARAD
      ELSE 
        WRITE(hSUMMARY,'(A)')    '\VAR VESCFULL = invalid'      !supereddington, no escape velocity at Rstar
      ENDIF
      IF (LWINDSUB < 10**LOGL) THEN
        WRITE(hSUMMARY,'(A,F5.3)') '\VAR LOGLINF = ', LOGLINF   !log Linf/Lsun (L as seen from Obs.)
      ENDIF
      WRITE(hSUMMARY,'(A,F6.3)') '\VAR qion    = ', q         !Ionization parameter



      !Chemical composition
      DO NA=1, NATOM
        !calculate mass fraction of current element
        !note: ELEMENT(NA) returns the Name, SYMBOL(NA) just He, C, ...
        XELEM = ABXYZ(NA)*ATMASS(NA) / ATMEAN
        WRITE(hSUMMARY,'(A8,A2,A5,F7.5)')
     >     '\VAR Xn_',ADJUSTL(SYMBOL(NA)), '   = ', ABXYZ(NA)      !by number
        WRITE(hSUMMARY,'(A8,A2,A5,F7.5)') 
     >     '\VAR Xm_',ADJUSTL(SYMBOL(NA)), '   = ', XELEM          !by mass
      ENDDO

      WRITE(hSUMMARY,'(A,F6.2)') '\VAR ATMEAN  = ', ATMEAN       !mean atomic mass
      WRITE(hSUMMARY,'(A,F6.2)') '\VAR XMU     = ', XMU(ND)      !mean particle mass (in AMU => mu) 
                                                                 ! at inner boundary
                                                                 !XMU is called XMASS in some routines

      !Calucate logQ from O-star people: Q = Mdot / ( R_star * vinf )^(3/2)                                                                 
      LOGQ = LOGMDOT + 0.5 * LOG10(Dout) 
     >                 - 3. / 2. * ( LOG10(RSTARSU) + LOG10(VFINAL) )
      WRITE(hSUMMARY,'(A,F7.3)') '\VAR LOGQ    = ', LOGQ                                                           
                                                                 
      !Calculate transformed Mass-loss rate (Bestenlehner et al. 2013)
      LOGMDOTTRANS = LOGMDOT + 0.5 * LOG10(Dout) - LOG10(VFINAL/1000)
     >                 - 0.75 * (LOGL - 6.)
      WRITE(hSUMMARY,'(A,F7.3)') '\VAR LOGMDOTTRANS = ', LOGMDOTTRANS                                                           

      !Calculate modified wind momentum
      LOGDMOM = 0.5 * LOG10(Dout) + LOGMDOT + LOG10(MSUN) - LYR + 
     >          LOG10(VFINAL) + 5.0 + 0.5 * LOG10(RSTARSU)
      WRITE(hSUMMARY,'(A,F7.3)') '\VAR LOGDMOM = ', LOGDMOM

      !Calculate MDOT from Nugis & Lamers (2002) for thick winds
      IF (Lcand > 0) THEN
        MDOTNL = 16./3. * PI4 * STEBOL / CLIGHT
     >           * SQRT( RGAS * (gammaIon + 1.) / XMUsonic )
     >           * (Rsonic*RSTAR)**3 * Tsonic**(4.5) 
     >           / (MSTAR * MSUN) / GCONST
        WRITE(hSUMMARY,'(A,F7.3)') '\VAR LOGMDOTNL = ', 
     >                                    LOG10 (MDOTNL / XMSUNPYR)
      ENDIF
      
      !Calculate R and T at TauRoss = 1/3, 2/3 and 1
      TAU13=0.333333333333
      IF (TAUROSS(ND) < TAU13) THEN
         R13=1.
         T13=T(ND)
      ELSE
         CALL LIPO (R13,TAU13,RADIUS,TAUROSS,ND)
         T13 = SQRT(1./R13) * TEFF
      ENDIF
      TAU23=0.666666666666
      IF (TAUROSS(ND) < TAU23) THEN
         R23=1.
         T23=T(ND)
      ELSE
         CALL LIPO (R23,TAU23,RADIUS,TAUROSS,ND)
         T23 = SQRT(1./R23) * TEFF
      ENDIF
      TAU1=1.
      IF (TAUROSS(ND) < TAU1 ) THEN
         R1 =1.
         T1 =T(ND)
      ELSE
         CALL LIPO (R1 ,TAU1 ,RADIUS,TAUROSS,ND)
         T1 = SQRT(1./R1) * TEFF
      ENDIF

      
      WRITE(hSUMMARY,'(A,F6.3)') '\VAR R1      = ', R1              !Radius where (Tau_Ross = 1)
      WRITE(hSUMMARY,'(A,F6.3)') '\VAR R13     = ', R13             !Radius where (Tau_Ross = 1/3)
      WRITE(hSUMMARY,'(A,F6.3)') '\VAR R23     = ', R23             !Radius where (Tau_Ross = 2/3)
      WRITE(hSUMMARY,'(A,F10.2)') '\VAR T1      = ', T1             !Eff. Temperature where (Tau_Ross = 1)
      WRITE(hSUMMARY,'(A,F10.2)') '\VAR T13     = ', T13            !Eff. Temperature where (Tau_Ross = 1/3)
      WRITE(hSUMMARY,'(A,F10.2)') '\VAR T23     = ', T23            !Eff. Temperature where (Tau_Ross = 2/3)  
      G23LOG = GLOG - 2. * LOG10(R23)
      WRITE(hSUMMARY,'(A,F8.4)') '\VAR G23LOG  = ', G23LOG
      IF (GEFFLOG > 0.) THEN
        GEFF23LOG = GEFFLOG - 2. * LOG10(R23)
        WRITE(hSUMMARY,'(A,F8.4)') '\VAR GEFF23LOG = ', GEFF23LOG
      ENDIF

      WRITE(hSUMMARY,'(A,F10.4)') '\VAR TAUMAX  = ', TAUROSS(ND)
      WRITE(hSUMMARY,'(A,F6.3)')  '\VAR WORKRATIO = ', WORKRATIO    !calculated in INITFCORR
      
      !Grid index calculation
      Tindex = (log10(TEFF) - 4.35) / 0.05
      Rtindex = (2.1 - log10(RTRANS)) / 0.1
      WRITE(hSUMMARY,'(A,F6.2)') '\VAR Tgrid   = ', Tindex
      WRITE(hSUMMARY,'(A,F6.2)') '\VAR Rtgrid  = ', Rtindex
      
      WRITE(hSUMMARY,*)
      CLOSE (hSUMMARY)

      !Write important parameters into out file:
      
      WRITE(hOUT,*) ''
      WRITE(hOUT,*) '------------------------------------'
      WRITE(hOUT,*) 'Output of converged Model Parameters'
      WRITE(hOUT,'(A, F7.0)') ' T* = ', TEFF
      WRITE(hOUT,'(A, 1PG12.3,0P,A,F7.3,A)') 
     >      ' Rt = ', RTRANS, ' = ', LOG10(RTRANS), 'dex'
      WRITE(hOUT,'(A, F7.4)') ' R* = ', RSTARSU
      WRITE(hOUT,'(A, F7.3,A)') ' Md = ', LOGMDOT, 'dex'
      WRITE(hOUT,'(A, F7.1)') ' V8 = ', VFINAL
      WRITE(hOUT,'(A, F7.2)') ' D(ND)  = ', DENSCON(ND)
      WRITE(hOUT,'(A, F7.2)') ' M/Msun = ', MSTAR
      WRITE(hOUT,'(A, F8.4)') ' log g = ', GLOG
      WRITE(hOUT,'(A, F8.4)') ' log geff = ', GEFFLOG
      WRITE(hOUT,'(A, F7.3)') ' q_ion = ', q
      IF (bFULLHYDROSTAT) THEN
        WRITE(hOUT,'(A, F7.3)') ' Ed.Gamma = ', GAMMARAD
      ELSE
        WRITE(hOUT,'(A, F7.3)') ' Ed.Gamma = ', GEDDe
      ENDIF
      WRITE(hOUT,'(A, F7.3)') ' ATMEAN = ', ATMEAN
      WRITE(hOUT,'(A, F7.3)') ' mu = ', XMU(1)
      IF (Lcand > 0) THEN
        WRITE(hOUT,'(A, F9.3)') ' Vs = ', Vsonic
      ENDIF
      WRITE(hOUT,*) '------------------------------------'

C      WRITE (hCPR,*) ' DEBUG-GRADI: '
C      DO L=1, ND
C        WRITE (hCPR,'(A,I2,A,F12.4)') 'L=', L, ' GRADI=', GRADI(L)
C      ENDDO

      !Note: TAUROSS and DENSITY are not given because the are already printed
      !      in O P T I C A L   D E P T H   S C A L E S
      WRITE (hOUT,'(A)') 'Depth-dependent quantities'
      WRITE (hOUT,*)
      WRITE (hOUT,'(2A)')
C     >   ' DEPTH     R-1     RMAX-R    CRITERION     VELOCITY        ',
     >   ' DEPTH     R-1    CRITERION      VELOCITY (TYPE)       ',
     >   'GRADIENT     CLUMPING    CLUMPING   MEAN PARTICLE'
      WRITE (hOUT,'(2A)')
C     >   ' INDEX                                      (KM/S)     (KM/',
     >   ' INDEX                            (KM/S)           (KM/',
     >   'S PER RSTAR)  DENSCON     FILLFAC    MASS (XMU)'
      WRITE (hOUT,*)
      CTABLINE = '(I3,5X,G10.3,2X,A8,3X,G11.4,1X,A2,5X,' // 
     >               'F12.4,5X,F9.4,3X,F9.4,4X,F9.4)'
      DO L=1,ND
        IF (TRIM(VELOCRIT(L)) == "") THEN
          IF (bHYDROSOLVE) THEN
            VCRIT = 'HD'          !hydrodynamic law
          ELSEIF (RADIUS(L) < RCON) THEN
            IF (bThinWind) THEN
              VCRIT = 'HS'        !correct hydrostatic integration
            ELSE
              VCRIT = 'ST'        !static (exponential) approach
            ENDIF
          ELSE
            VCRIT = 'B '          !Beta law
          ENDIF
        ELSE
          VCRIT = VELOCRIT(L)(1:2)
        ENDIF
        RL1=RADIUS(L)-1.
        RL2=RADIUS(1)-RADIUS(L)
        WRITE(hOUT,FMT=CTABLINE) L, RL1, INCRIT(L), VELO(L), 
     >     VCRIT, GRADI(L), DENSCON(L), FILLFAC(L), XMU(L)
      ENDDO
      WRITE (hOUT,*) ''
      WRITE (hOUT,*) '------------------------------------'
      WRITE (hOUT,*)

      WRITE (hCPR,*) 
      WRITE (hCPR,*) 'A model parameter summary file was written'
      WRITE (hCPR,*) 

      RETURN
      END
