      SUBROUTINE INITVEL (RMAX, TEFF, GEFFLOG, RSTAR, XMASS,
     >                    VTURB, bHScaleOnly, bHydroStat)
C*******************************************************************************
C***  INITIALIZATION OF THE VELOCITY-FIELD PARAMETERS
C***  CALLED FROM: WRSTART; STEAL - ENSURETAUMAX
C***  New parametrization of beta law, wrh  5-Mar-2001 
C***  2BETA-LAW, 6-Apr-2006
C*******************************************************************************

      IMPLICIT NONE

      REAL, INTENT(IN) :: RMAX, TEFF, GEFFLOG, RSTAR, XMASS
      LOGICAL, INTENT(INOUT) :: bHydroStat,    !returns .FALSE. if no hydrostatic domain is encountered
     >                          bHScaleOnly    !if true only HSCALE is calculated in this routine

      REAL :: VFINAL, VMINCAND, BETA, BETA2, BETA2FRACTION, VMIN,
     >        VPAR1, VPAR2, HSCALE, VPAR1_2, VPAR2_2, RCON, VCON, 
     >        VTURB

      COMMON/VELPAR/ VFINAL,VMIN,BETA,VPAR1,VPAR2,RCON,HSCALE,
     >            BETA2, BETA2FRACTION, VPAR1_2, VPAR2_2
      REAL, EXTERNAL :: DELTAGR
 
      REAL :: VMAX, Q, S, T, RCONOLD, RCONMIN, RCONMAX, RGRADMAX, DXX
      INTEGER :: ITER

      REAL, PARAMETER :: BOLTZ = 1.38E-16   !BOLTZMANN CONSTANT (ERG/DEG)
      REAL, PARAMETER :: AMU   = 1.66E-24   !ATOMIC MASS UNIT (GRAMM)

C***  DETERMINATION OF THE SCALE HEIGHT (HYDROSTATIC EQ.)
C***  XMASS= MEAN MASS IN AMU (e.g.  2.0 FOR HE II, 1.33 FOR HE III )
C***  Turbulence pressure added 13-Mar-2014
      HSCALE = (BOLTZ*TEFF/(XMASS*AMU) + (VTURB*1.E5)**2)
     >                   /10.**GEFFLOG /RSTAR

      IF (bHScaleOnly) RETURN     
      
C***  In case of 2BETA-lwa, the second beta must be larger than 1
C***  in order to make sure that the contribution at the connection
C***  point vanishes 
      IF (BETA2FRACTION > 0) THEN
         IF (BETA2 .LT. 1.) THEN       
           WRITE (0, '(A)') '*** SECOND BETA MUST BE > 1' 
           STOP             '*** FATAL ERROR IN INITVEL'
         ENDIF
      ENDIF

      VMAX = VFINAL * (1.-BETA2FRACTION)
      Q = (VMIN/VMAX)**(1./BETA)
      
ccc   HSCALE=(BOLTZ*TEFF) / (XMASS*AMU*RSTAR*10.**GEFFLOG)
ccc   WRITE (0,*) " HSCALE: ", HSCALE
 

C***  ITERATIVE DETERMINATION OF THE CONNECTION POINT RCON
C***  BY REQUIRING A CONTINUOUS AND SMOOTH CONNECTION
 
C**   Initialization
      RCON=1.
      ITER = 0

C---  RCON Iteration loop ------------------------------------------
      DO
        RCONOLD=RCON
        ITER = ITER + 1
      
C        WRITE (0, '(A, I2, A, F12.6, A, F12.6)') 
C     >      'ITER=', ITER, '  VPAR2=', VPAR2, '  RCON=', RCON 
        IF (ITER >= 99) THEN
          WRITE (0, '(A)') '*** ITERATION FOR VPARs NOT CONVERGED!' 
          STOP             '*** FATAL ERROR IN INITVEL'
        ENDIF

C***    CHOSE THE BETA-LAW PARAMETERS SUCH THAT BOTH VELOCITIES AGREE
C***    AT THE CONNECTION POINT
        VCON = VMIN * EXP((RCON-1.)/HSCALE)
        CALL INITVELBETAPAR(RMAX, RCON, VCON, .TRUE.)
C        Q = (VCON/VMAX)**(1./BETA)
C        S = (RMAX - 1. + RCON) / 2.
C        T = (RMAX - Q * RCON) / (1.-Q) - RMAX * RCON
C        VPAR2 = SQRT(S*S + T) - S
C        VPAR1 = VMAX / (1.-1./(VPAR2+RMAX))**BETA
C        IF (BETA2FRACTION > 0.) THEN
C          VPAR2_2 = 1. - RCON   
C          VPAR1_2 = BETA2FRACTION * 
C     >             VFINAL/(1.-1./(VPAR2_2+RMAX))**BETA2
C        ENDIF

C***    Now it is additionally demanded that the velocity gradient
C***    is also smooth. The FUNCTION DELTAGR gives the difference 
C***    of the gradients from the outer minus the inner law 

C***    CHECK IF AT THE INNER BOUNDARY THE EXPONENTIAL LAW GRADIENT
C***    EXCEEDS ALREADY THE BETA-LAW GRADIENT
C***    Note: this check may only be performed in the first iteration
        IF (ITER <= 1) THEN
          IF (DELTAGR(1.) .LT. .0) THEN
            RCON=1.
            CALL REMARK('INITVEL: NO HYDROSTATIC DOMAIN ENCOUNTERED')
            bHydroStat = .FALSE.
            RETURN
          ENDIF
        ENDIF
 
C***    LOWER GUESS FOR RCON
        RCONMIN=AMAX1(1.,1.-VPAR2+1.E-10)
C***    Maximum of gradient in beta laws:
        RGRADMAX = 1. + (BETA-1.)/2. - VPAR2


C***    Loop: INCREASE RCONMIN STEPWISE TO ENSURE DELTAGR(RCONMIN) .GE. 0.
        !@Check: Kann dies eine Endlosschleife werden?
        DO
          IF (DELTAGR(RCONMIN) >= .0) EXIT
          DXX = 0.1 * HSCALE
          RCONMIN = RCONMIN + DXX
          IF (RCONMIN > RGRADMAX) THEN
            WRITE (0, '(A)') '*** INITVEL: NO CONNECTION POINT FOUND'
            STOP             '*** ERROR IN INITVEL'
          ENDIF
cc           WRITE (0, '(A, F12.6, A, F12.6)') 
cc     >           'RCONMIN INCREASED TO', RCONMIN, 
cc     >           '  DELTAGR=', DELTAGR(RCONMIN)
        ENDDO
 
C***    UPPER GUESS FOR RCON
        RCONMAX=RCONMIN
C**     Loop: 
        DO
           RCONMAX=RCONMAX+HSCALE
           IF (DELTAGR(RCONMAX) < .0) EXIT
           IF (RCONMAX >= RMAX) THEN
             RCON=RMAX
             CALL REMARK('INITVEL: NO BETA-LAW DOMAIN ENCOUNTERED')
             RETURN
           ENDIF
        ENDDO
 
C***    Now find new RCON where DELTAGR is zero:
        CALL REGULA (DELTAGR,RCON,.0,RCONMIN,RCONMAX,1.E-8)

        IF (ABS(RCON-RCONOLD) <= 1.E-7) EXIT
        
      ENDDO
C---  End of RCON Iteration loop --------------------------------------
 
      RETURN
      END
