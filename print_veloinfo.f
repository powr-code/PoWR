      SUBROUTINE PRINT_VELOINFO (bVELOTABLE, THIN, VCON, VTURB,
     >          FSONICTA, FSONICTA_CARDS, TEFF, XMASS, RMAX, 
     >          GEFFLOG, RSTAR)
C*******************************************************************************
C***  PRINT INFORMATION ABOUT THE VELOCITY FIELD 
C***  CALLED FROM WRSTART, STEAL
C*******************************************************************************
 
      IMPLICIT NONE
       
      REAL,    INTENT(IN) :: VCON, VTURB, FSONICTA, FSONICTA_CARDS
      REAL,    INTENT(IN) :: TEFF, XMASS, RMAX, GEFFLOG, RSTAR
      LOGICAL, INTENT(IN) :: bVELOTABLE, THIN
      REAL,    EXTERNAL   :: WRVEL

C***  COMMON /VELPAR/ TRANSFERS VELOCITY-FIELD PARAMETERS 
      COMMON /VELPAR/ VFINAL, VMIN, BETA, VPAR1, VPAR2, RCON, HSCALE,
     >     BETA2, BETA2FRACTION, VPAR1_2, VPAR2_2
      REAL VFINAL, VMIN, BETA, VPAR1, VPAR2, RCON, HSCALE,
     >     BETA2, BETA2FRACTION, VPAR1_2, VPAR2_2

C***  BOLTZMANN CONSTANT (ERG/DEG) for HSCALE
      REAL, PARAMETER :: BOLTZ = 1.38E-16   
      REAL, PARAMETER :: AMU = 1.66E-24     !ATOMIC MASS UNIT (GRAMM)

C***  INFO CONCERNING THE VELOCITY FIELD 
      PRINT *, ' '
      IF (RCON == 1.) THEN
         PRINT *, 'VELOCITY FIELD: NO HYDROSTATIC DOMAIN ENCOUNTERED'
      ELSE 
         PRINT '(2A,1PG12.3)', 
     >         'VELOCITY FIELD: HYDROSTATIC DOMAIN ',
     >         'BELOW RADIUS RCON = 1 + ', RCON-1. 
         PRINT '(17X,A,1PG12.3,A)', 
     >         'CORRESPONDING TO VCON = ', VCON, ' km/s' 
C***     HSCALE printed only from wrstart where XMASS is available 
         IF (XMASS .GT. .0) THEN
            HSCALE = (BOLTZ*TEFF/(XMASS*AMU) + (VTURB*1.E5)**2)
     >                   /10.**GEFFLOG /RSTAR
            IF (.NOT. THIN) PRINT '(A, 1PG12.3)', 
     >         'FIXED SCALE HEIGHT HSCALE = ', HSCALE
            IF (THIN) PRINT '(2A,1PG12.3)',
     >         ' HYDROSTATIC EQ. WITH DEPTH-DEPENDENT ',
     >         'SCALE HEIGHT, APPROXIMATELY HSCALE = ', HSCALE
            IF (VTURB .GT. .0) PRINT '(A, F6.1, A)',
     >         ' TURBULENCE PRESSURE TAKEN INTO ACCOUNT, VTURB=', 
     >          VTURB, ' km/s'
         ENDIF
      ENDIF

      IF (bVELOTABLE) THEN
       PRINT *, 'OUTER PART: input as VELOTABLE'
       PRINT '(A, F7.3)', 
     >     'Connection velocity VCON = VSOUND * ', FSONICTA
       IF (FSONICTA .EQ. FSONICTA_CARDS) THEN 
          PRINT *, ' as requested'
       ELSE
          PRINT '(A, F7.3)', 
     >     ' has been reduced; requested was VSOUND * ', FSONICTA_CARDS
       ENDIF          

      ELSE
C***  Re-construct the beta-law parameters from RCON, VCON
      CALL INITVELBETAPAR (RMAX, RCON, VCON, .TRUE.)
       IF (BETA .GT. .0) PRINT '(A,F4.1,A,1PG12.4,A,1PG12.4)',
     >  ' OUTER PART: BETA=', BETA, ' LAW, VPAR1=', VPAR1
     >  ,' , VPAR2=', VPAR2

       IF (BETA2FRACTION .GT. .0) 
     >    PRINT '(13X,A,F4.2,A,F5.2,A,1PG12.4,A,1PG12.4)',
     >  '2-BETA-LAW:  BETA2FRACTION=', BETA2FRACTION,
     >    '  BETA2=', BETA2, 
     >    '  VPAR1_2=', VPAR1_2, '  VPAR2_2=', VPAR2_2
      ENDIF 

      PRINT *, ' '

      RETURN
      END
