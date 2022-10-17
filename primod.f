      SUBROUTINE PRIMOD (ND,RADIUS,INCRIT,ENTOT,T,VELO,GRADI,NP,OLDTEMP,
     >          MODHEAD,JOBNUM,MODOLD,JOBNOLD,TTABLE,TAUROSS,R23,
     >          TEFFOLD, THIN, ITTAU, MAXITTAU, RCON, BTWOT, MODOLD2, 
     >          JOBNOLD2, TFAC, BETA, VPAR1, VPAR2, DENSCON, BETA2,
     >          BETA2FRACTION, HSCALE, bNoDetails,bStartCall, VTURB)
C*******************************************************************************
C***  PRINTOUT OF THE DEPTH-DEPENDENT MODEL SPECIFICATIONS
C***  CALLED FROM WRSTART, STEAL
C*******************************************************************************
 
      IMPLICIT NONE
       
      REAL, DIMENSION(ND), INTENT(IN) :: 
     >      RADIUS, ENTOT, T, VELO, GRADI, TAUROSS, DENSCON
      REAL, INTENT(IN) :: 
     >      RCON, BETA, VPAR1, VPAR2, BETA2, BETA2FRACTION, HSCALE,
     >      R23, TEFFOLD, TFAC, VTURB
      INTEGER, INTENT(IN) ::
     >      ND, NP, JOBNUM, JOBNOLD, JOBNOLD2, ITTAU, MAXITTAU
      LOGICAL, INTENT(IN) ::
     >      OLDTEMP,
     >      TTABLE,
     >      BTWOT,
     >      THIN,         !THIN-Wind option set
     >      bNoDetails,   !Do not print table with values for all depth points
     >      bStartCall    !if false, informations only known by WRSTART are skipped
      CHARACTER(8), DIMENSION(ND), INTENT(IN) :: 
     >      INCRIT
      CHARACTER(100), INTENT(IN) :: 
     >      MODHEAD, MODOLD, MODOLD2

      COMMON /COMTEFF/ TEFF,TMIN,TMODIFY,SPHERIC

      REAL :: RL1, RL2, DENS, TEFF, TMIN, TMODIFY
      INTEGER :: L
      LOGICAL :: SPHERIC
      CHARACTER(40) :: MESSAGE
     
      IF (MODHEAD /= ' ') THEN
        PRINT 1,MODHEAD,JOBNUM
    1   FORMAT (1X,  A  ,20X,'JOB NO.',I7,
     $       ///,20X,'M O D E L   S P E C I F I C A T I O N',/,20X,
     $  38('='))
      WRITE (*,*)
      ENDIF
      
      IF (.NOT. bNoDetails) THEN    !details option
        WRITE (*,*)
        WRITE (*,'(20X,A)') '   DEPTH-DEPENDENT QUANTITIES'
        WRITE (*,*)
        WRITE (*,'(3A)')
     >   ' DEPTH     R-1     RMAX-R    CRITERION     VELOCITY        ',
     >   'GRADIENT      LOG NUMB. DENS.   EL. TEMPERATURE ',
     >   'TAU-ROSS.  CLUMPING'
        WRITE (*,'(3A)')
     >   ' INDEX                                      (KM/S)     (KM/',
     >   'S PER RSTAR)  (ATOMS PER CM+3)        (KELVIN)   ',
     >   '(IN LTE)   DENSCON'
        WRITE (*,*)

C***    LOOP OVER ALL DEPTH POINTS **************************************
        DO L=1,ND
          RL1=RADIUS(L)-1.
          RL2=RADIUS(1)-RADIUS(L)
          DENS=ALOG10(ENTOT(L))
          PRINT 3, L, RL1, RL2, INCRIT(L), VELO(L), GRADI(L), DENS,
     >         T(L), TAUROSS(L), DENSCON(L)
        ENDDO
    3   FORMAT(I3, 5X, G10.3, G10.3, 3X, A8, 2X, G11.4, F14.2, F19.3,
     >       F19.0, F12.3, F10.3)
C**********************************************************************
      ENDIF !end of details option

      IF (bStartCall) PRINT 4,NP
    4 FORMAT (/,' NUMBER OF IMPACT-PARAMETER POINTS  NP =',I3)

      PRINT 9, R23
    9 FORMAT (' RADIUS WHERE TAU-ROSSELAND = 2/3: R23 =',F7.3)


C***  MESSAGES CONCERNING THE TEMPERATURE STRUCTURE
      IF (TTABLE) PRINT *,'TEMPERATURE STRUCTURE FROM TABLE INPUT'

      IF (bStartCall .AND. SPHERIC) PRINT 7,TEFF
    7 FORMAT (' TEMPERATURE STRUCTURE AS FOR A SPHERICAL, ',
     $   'GREY LTE ATMOSPHERE (APPROXIMATELY) WITH TEFF=',F7.0)

      IF (OLDTEMP) THEN
        PRINT 11, MODOLD, JOBNOLD
   11   FORMAT (' TEMPERATURE STRUCTURE TAKEN FROM THE FOLLOWING ',
     >     'OLD MODEL:',/,1X,A,5X,'AFTER JOB NO.',I3)
        IF (BTWOT) THEN
          IF (TFAC .GE. 0. .AND. TFAC .LE. 1.) THEN
            IF (TFAC .EQ. 0.5) THEN 
              MESSAGE = ' INTERPOLATION'
            ELSE IF (TFAC .LT. 0.5) THEN 
              MESSAGE = ' INTERPOLATION, MODEL 1 STRESSED'
            ELSE
              MESSAGE = ' INTERPOLATION, MODEL 2 STRESSED'
            ENDIF
          ELSE
            IF (TFAC .LT. 0.) THEN
              MESSAGE = ' EXTRAPOLATION, MODEL 1 STRESSED'
            ELSE
              MESSAGE = ' EXTRAPOLATION, MODEL 2 STRESSED'
            ENDIF
          ENDIF
          WRITE (*,'(A,/,A,A,A,I3)') 
     >           ' SECOND MODEL : ',' ', MODOLD2, 
     >           'AFTER JOB NO.', JOBNOLD2
          WRITE (*,'(A,F4.1,A)')
     >           ' COMBINE FACTOR = ', TFAC, MESSAGE
        ENDIF
      ENDIF

      IF (OLDTEMP .AND.  TEFF .NE. TEFFOLD) PRINT 12, TEFF, TEFFOLD
   12 FORMAT (' ... SCALED WITH TEFF (NEW/OLD) =', F7.0, ' / ', F7.0)

      IF (bStartCall .AND.
     >      .NOT. TTABLE .AND. .NOT. SPHERIC .AND. .NOT. OLDTEMP)
     $   PRINT 8,TEFF
    8    FORMAT (' TEMPERATURE STRUCTURE AS FOR A ',
     $           'PLANE-PARALLEL, GREY LTE ATMOSPHERE WITH TEFF=',F7.0)
          
      IF (TMODIFY .NE. .0) PRINT 5,TMODIFY
    5 FORMAT (' TEMPERATURE STRATIFICATION MODIFIED: TMODIFY=',F6.3)
 
      IF (TMIN .GT. .0) PRINT 10, TMIN
   10 FORMAT (' MINIMUM TEMPERATURE SPECIFIED: TMIN=',F7.0)



C***  MESSAGES CONCERNING THE VELOCITY FIELD 
      IF (RCON == 1.) THEN
         PRINT *, 'VELOCITY FIELD: NO HYDROSTATIC DOMAIN ENCOUNTERED'
      ELSE 
         PRINT '(2A,G10.3,A,1PG12.3,A)', 
     >         ' VELOCITY FIELD: HYDROSTATIC DOMAIN ',
     >         'BELOW RADIUS RCON = 1 + ', RCON-1. , '  (HSCALE = ',
     >         HSCALE,')'
         IF (THIN) PRINT *,'HYDROSTATIC EQ. WITH TEMPERATURE',
     >                     '-DEPENDENT SCALE HEIGHT IS USED'
         IF (VTURB .GT. .0) PRINT '(A, F6.1)',
     >      ' TURBULENCE PRESSURE TAKEN INTO ACCOUNT, VTURB=', VTURB
         ENDIF
      IF (BETA .GT. .0) PRINT '(A,F4.1,A,3PG12.6,A,1PG12.3)',
     >  ' OUTER PART: BETA=', BETA, ' LAW, VPAR1=', VPAR1
     >  ,' , VPAR2=', VPAR2
      IF (BETA2FRACTION .GT. .0) PRINT '(13X,A,F4.2,A,F5.2)',
     >  '2-BETA-LAW:  BETA2FRACTION=', BETA2FRACTION,
     >    '  BETA2=', BETA2

      IF (ITTAU .GT. 1) THEN
       IF (ITTAU .LT. MAXITTAU) THEN
        PRINT '(A,I3,A)', 
     >    ' VMIN AUTOMATICALLY ADJUSTED TO SPECIFIED TAUMAX (', 
     >    ITTAU, ' ITERATIONS)'
        ELSE                  
        PRINT '(A,I3,A)', 
     >   '*** WARNING: AUTOMATICAL ADJUSTMENT OF VMIN NOT CONVERGED IN',
     >   ITTAU, ' ITERATIONS ! ***'
        ENDIF
       ENDIF
 
      PRINT *, ' '

      RETURN
      END
