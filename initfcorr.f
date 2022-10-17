      SUBROUTINE INITFCORR (TEFF, RSTAR, T, RNE, bNoARAD,
     >                      ABXYZ, ATMASS, ELEMENT, NATOM, WRTYPE,
     >                      RADIUS, VELO, GRADI, VTURB, ENTOT, ND,
     >                      HTOTM, HTOTG, HTOTOBS, HTOTL, ACONT, ATHOM,
     >                      AGRAV, AMECH, ARAD, APRESS, WORKRATIO, 
     >                      CLUMP_SEP, OPAROSS, MacroCard, DENSCON,
     >                      FILLFAC, OPALINE_SCALE,
     >                      FTCOLI, XJTOTL, XKTOTL, XNTOTL, HTOTCMF0, 
     >                      XMSTAR, MFORM, GLOG, ATMEAN, VMACH,
     >                      TAUROSS, QIONMEAN, GAMMARADMEAN, RCON)

      IMPLICIT NONE
      INCLUDE 'interfacebib.inc'
     
      INTEGER, INTENT(IN) :: NATOM, ND     
      REAL, DIMENSION(NATOM) :: ABXYZ, ATMASS
      CHARACTER(10), DIMENSION(NATOM) :: ELEMENT
      REAL, DIMENSION(ND) :: RADIUS, VELO, GRADI, ENTOT, T, RNE,
     >                       HTOTM, HTOTG, HTOTOBS, HTOTL, HTOTCMF0,
     >                       AGRAV, AMECH, ARAD, APRESS, GEFFL, RI,
     >                       ATHOM, ACONT, ALINES, DENSCON, FILLFAC, 
     >                       XJTOTL, XKTOTL, XNTOTL, FTCOLI, 
     >                       OPAROSS, ClumpScale, TAUROSS, VMACH
     
      INTEGER :: MFORM, NA, L, LCON
      REAL :: GLOG, VTURB, ATMEAN, XMSTARG, XMSTAR, TEFF, XMDOT,
     >        EWTOT, WMEAN, XHY, YHE, XC, SUM, SUMADD, HNULL, VELOINT,
     >        RINT, DELTAH, DELTAV, DELTAR, XMECH, VDR, RSTAR, RL2,
     >        HMECH, HMECH1, HGRAV, XGRAV, XKINT, XJINT, WRAD, WWIND, 
     >        RHOL, RHOLP, XMUL, XMULP, AL2, ALP2, PL, PLP, WORKRATIO,
     >        CLUMP_SEP, OPALINE_SCALE, DCmacro, fAmod, TauCL,
     >        QIONMEAN, GAMMARADMEAN, RCON
      
      CHARACTER(2) :: WRTYPE
      CHARACTER(9) :: MLRELATION
      CHARACTER(120) :: MacroCard

      LOGICAL :: BFMEC, bNoARAD

C***  Physical constants
      REAL, PARAMETER :: PI4 = 12.5663706144    !PI4 = 4*PI
      REAL, PARAMETER :: GCONST = 6.6726E-8     !GCONST = GRAVITATIONAL CONSTANT (CGS-UNITS)
      REAL, PARAMETER :: STEBOL = 1.8046E-5     !STEBOL = STEFAN-BOLTZMANN CONSTANT (CGS-UNITS) / PI
      REAL, PARAMETER :: CLIGHTKM = 2.9979E5    !CLIGHT = SPEED OF LIGHT IN KILOMETER/SECOND
      REAL, PARAMETER :: BOLTZK = 1.3807E-16    !BOLTZMANN CONSTANT on cgs units
      REAL, PARAMETER :: XMSUN = 1.989E33       !XMSUN = Solar Mass (g)
      REAL, PARAMETER :: XLSUN = 3.85E33        !Solar Luminosity (CGS-Units)
      REAL, PARAMETER :: AMU = 1.6605E-24       !Atomic mass unit (gramm) = m_H

C***  If the stellar mass XMSTAR is not yet available from the MODEL 
C***  File, it is now calculated (this should only happen for old models)

      IF (XMSTAR .EQ. .0) THEN
C***     Find type (WN or WC):
C***     Calculate Helium and Carbon Mass fraction from the number fractions
         SUM = .0
         XHY = .0
         YHE = .0
         XC = .0
         DO NA=1, NATOM
            SUMADD = ABXYZ(NA)*ATMASS(NA)
            SUM = SUM + SUMADD
            IF (ELEMENT(NA) .EQ. 'HYDROGEN  ') XHY = SUMADD
            IF (ELEMENT(NA) .EQ. 'HELIUM    ') YHE = SUMADD
            IF (ELEMENT(NA) .EQ. 'CARBON    ') XC  = SUMADD
         ENDDO
C***     Hydrogen mass fraction:
         XHY = XHY / SUM
C***     Helium mass fraction:
         YHE = YHE / SUM
C***     Carbon mass fraction; if X_C > 0.1, WRtype is assumed to be WC
         XC  = XC  / SUM
         IF (WRTYPE == '') THEN
           IF (XHY .GE. 0.45) THEN 
             WRTYPE = 'OB' 
           ELSE IF (XC .GE. 0.1) THEN 
             WRTYPE = 'WC' 
           ELSE
             WRTYPE = 'WN' 
           ENDIF
         ENDIF           

C***     M-L relations: 
C***        OB type: Goetz (mandatory) - H-burner
C***        WN type: Goetz (default) or Langer '89 (optional) - He burner
C***        WC type: Langer (mandatory)
         IF (WRTYPE == 'WC' .OR. 
     >      (WRTYPE == 'WN' .AND. MFORM == 1) .OR.
     >       WRTYPE == 'WC') THEN
            CALL MLANGER (XMSTARG, TEFF, RSTAR, YHE, WRTYPE)         
            XMSTAR = XMSTARG / XMSUN
            MLRELATION='Langer'
         ELSE
            CALL MGOETZ (XMSTAR, TEFF, RSTAR, XHY, WRTYPE)
            XMSTARG = XMSTAR * XMSUN
            MLRELATION='Graefener'
         ENDIF
         GLOG = ALOG10(GCONST * XMSTARG / RSTAR / RSTAR)
         WRITE (0,'(A, /, A, F7.3)') 
     >   'INITFCORR: WARNING: MASS (FOR ACC_PLOT) NOT ON MODEL FILE -',
     >   '   CALCULATED FROM M-L-RELATION: TYPE = ' // WRTYPE // 
     >   ' FROM ' // MLRELATION // '   M/M_sun =', XMSTAR
      ELSE
         XMSTARG = XMSTAR * XMSUN
      ENDIF

      
C***  Integrals of total work (initialization)
      WRAD  = .0
      WWIND = .0

C***   Mittlere Massenzahl
      ATMEAN = 0.
      DO NA=1, NATOM
        ATMEAN = ATMEAN + ABXYZ(NA) * ATMASS(NA)
      ENDDO
      WMEAN = AMU*ATMEAN

            
      DO L=1, ND          
         RL2 = RADIUS(L) * RADIUS(L)
         XMDOT = ENTOT(L) * VELO(L)*1.E5 * PI4*RL2*RSTAR*RSTAR * WMEAN 

C***     Kinetic Part of Wind-Flux
         XMECH = XMDOT * 0.5 * VELO(L)*VELO(L)*1.E10
         HMECH = XMECH / (RSTAR*RSTAR*PI4*PI4)

C***     Gravitational Part of Wind-Flux
         XGRAV = XMDOT*GCONST*XMSTARG/RSTAR * (1. - 1./RADIUS(L))
         HGRAV = XGRAV / (RSTAR*RSTAR*PI4*PI4)
         
         IF (L .EQ. 1) HMECH1 = HMECH

         HTOTM(L) = HMECH
         HTOTG(L) = HGRAV

C***  Gravitational and Wind acceleration at radius interstices
         IF (L .EQ. ND) EXIT

         RI(L) = 0.5 * ( RADIUS(L) + RADIUS(L+1) )
         RINT = RI(L) * RSTAR
         AGRAV(L) = XMSTARG * GCONST / RINT / RINT

         VELOINT = 0.5 * (VELO(L) + VELO(L+1)) 
         DELTAV = (VELO(L) - VELO(L+1)) 
         DELTAR = (RADIUS(L) - RADIUS(L+1)) * RSTAR
         AMECH(L) = VELOINT * DELTAV / DELTAR * 1.E10
C***  Gas Pressure term
C***     a_press = - 1/rho dp/dr
C***       including turbulence:
C***           p = rho * ( k_B * T  / ( mu * m_H ) + v_turb^2 )
C***             = rho * ( a^2 + v_turb^2 )
C***         rho = ENTOT * m_H * ATMEAN = ENTOT * WMEAN
C***          mu = ATMEAN / ( 1 + RNE )
         EWTOT = 0.5 * (ENTOT(L) + ENTOT(L+1)) * WMEAN
         RHOL  = ENTOT(L) * WMEAN
         RHOLP = ENTOT(L+1) * WMEAN
         XMUL  = ATMEAN / (1.+RNE(L))
         XMULP = ATMEAN / (1.+RNE(L+1))
         AL2   = BOLTZK / AMU * T(L) / XMUL
         ALP2  = BOLTZK / AMU * T(L+1) / XMULP
         PL  = RHOL  * (AL2  + (VTURB*1.E5)**2)
         PLP = RHOLP * (ALP2 + (VTURB*1.E5)**2)
         APRESS(L) = (PLP-PL) / (DELTAR * EWTOT)

         VMACH(L) = SQRT(AL2) / 1.E5
         IF (L == ND-1) THEN
           VMACH(ND) = SQRT(ALP2) / 1.E5
         ENDIF
C         IF (L > 50) THEN
C           WRITE (0,'(I4,6(G15.8))') L, PL, EWTOT,
C     >     DELTAR, (PLP-PL), EWTOT*DELTAR, APRESS(L)
C         ENDIF

C***     Scale line acceleration if CARDS option is set        
         IF (OPALINE_SCALE /= 1.) THEN
           ALINES(L) = ARAD(L) - ACONT(L)
           ALINES(L) = ALINES(L) * OPALINE_SCALE
           ARAD(L) = ALINES(L) + ACONT(L)
         ENDIF

C***     Adjust with macroclumping effect 
         IF (CLUMP_SEP > 0.) THEN
           DCmacro = 0.5 * (DENSCON(L) + DENSCON(L+1))
           ClumpScale(L) = CLUMP_SEP * ( 3./PI4 /DCmacro )**(1./3.)
     >         * ( (RINT/RSTAR)**2 * VELOINT/VELO(1) )**(1./3.)
           TauCL = OPAROSS(L) * DCmacro * ClumpScale(L)
           IF ( TauCL <= 1.E-4 ) THEN
             fAmod = 1.
           ELSE
             fAmod = (1 . - EXP(-1.*TauCL) ) / TauCL
           ENDIF
           ARAD(L) = ARAD(L) * fAmod
           ACONT(L) = ACONT(L) * fAmod
         ENDIF          

C***  Integrals of total work (arbitrary units)
         WRAD = WRAD  +
     >    (ARAD(L)+APRESS(L)) * VELOINT * RINT * RINT * EWTOT * DELTAR
         WWIND= WWIND +
     >    (AMECH(L)+AGRAV(L)) * VELOINT * RINT * RINT * EWTOT * DELTAR

C***  Flux in the observer's frame (only for HSUM-PLOT!)
         XJINT = 0.5*(XJTOTL(L)+XJTOTL(L+1))
         XKINT = 0.5*(XKTOTL(L)+XKTOTL(L+1))
         HTOTOBS(L) = HTOTL(L) + VELOINT/CLIGHTKM*(XJINT+XKINT)

      ENDDO

C***  Calculate mean GAMMARAD and mean ionization parameter q  
      QIONMEAN = RNE(ND) / ATMEAN
      DO L=ND-1, 1, -1
        IF (RADIUS(L) > RCON) EXIT
        QIONMEAN = 1./REAL(ND-L+1)*(REAL(ND-L)*QIONMEAN+RNE(L)/ATMEAN)
      ENDDO
      IF (.NOT. bNoARAD) THEN
        CALL CALCGAMMARADMEAN(ARAD, AGRAV, RADIUS, TAUROSS, 
     >                        ND, RCON, RI, GAMMARADMEAN)
        GAMMARADMEAN = MIN(GAMMARADMEAN, 0.9)                      
      ENDIF
      
      
C***  CMF-Flux, starting at the inner boundary and intergrating the losses (see UNLU-Paper Eq. 62)
C***  Note: Previous versions of this subroutine contained a minor bug, taking 
C***        the integration step delta-r as one-sided, instead of between 
C***        interstices L+1/2, L-1/2. 
      HNULL = 0.25 * STEBOL * TEFF*TEFF*TEFF*TEFF
      HTOTCMF0(ND-1) = HNULL
      DO L=ND-1, 2, -1
         VDR = VELO(L)/RADIUS(L)
         DELTAH = ((GRADI(L) -VDR)*XKTOTL(L) + VDR*XJTOTL(L)) *
     >            0.5 * (RADIUS(L-1)-RADIUS(L+1)) / CLIGHTKM
         HTOTCMF0(L-1) = HTOTCMF0(L) - DELTAH 
      ENDDO

C***  Total work ratio
      WORKRATIO = WRAD / WWIND


      RETURN
      END
