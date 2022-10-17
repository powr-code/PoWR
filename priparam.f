      SUBROUTINE PRIPARAM (MODHEAD, TEFF, RSTAR, XMDOT, XLOGL, RTRANS,
     >           VFINAL, VDOP, DENSCON, FILLFAC, GLOG, GEFFLOG, GEDD, 
     >           GEddFix, RMAX, XMSTAR, WRTYPE, MASSORIGIN, LRTinput, 
     >           ND, MLRELATION, VTURB, MDOTINPUT)
C*******************************************************************************
C***  PRINTOUT OF THE STELLAR PARAMETERS
C***  CALLED FROM WRSTART
C***   Note: Some of these values can chance significantly during the iteration
C***    process, use results of STEAL->PRINTMODELSUMMARY to obtain final values
C*******************************************************************************

      IMPLICIT NONE

      INTEGER, INTENT(IN) :: ND, MASSORIGIN, GEddFix, LRTinput
      INTEGER, INTENT(IN) :: MDOTINPUT
      REAL, DIMENSION(ND), INTENT(IN) :: DENSCON, FILLFAC
      REAL, INTENT(IN) :: TEFF, RSTAR, XMDOT, GEDD,
     >                    VFINAL, VDOP, GLOG, GEFFLOG, RMAX, XMSTAR
      REAL, INTENT(INOUT) :: XLOGL, RTRANS
      CHARACTER(100), INTENT(IN) :: MODHEAD
      CHARACTER(2), INTENT(IN) :: WRTYPE
      CHARACTER(9) :: MLRELATION
      CHARACTER(30) :: FORMSTR

      REAL :: RSTARSU, VTURB

      REAL, PARAMETER :: RSUN = 6.96E10     !SOLAR RADIUS ( CM )
      REAL, PARAMETER :: TEFFSUN = 5780.    !SOLAR EFFECTIVE TEMPERATURE
      REAL, PARAMETER :: GCONST = 6.670E-8  !GRAVITATION CONSTANT (CGS UNITS)

      !File and channel handles (=KANAL)
      INTEGER, PARAMETER :: hOUT = 6        !write to wruniqX.out (stdout)
      INTEGER, PARAMETER :: hCPR = 0        !write to wruniqX.cpr (stderr)


      RSTARSU = RSTAR/RSUN                  !STELLAR RADIUS IN SOLAR UNITS
 
      WRITE (hOUT,1) MODHEAD(:32), MODHEAD(33:)
    1 FORMAT (///,80('*'),/,'*',/,
     >        '* FUNDAMENTAL PARAMETERS',/,
     >        '* ======================',/,'*',/,
     >        '* ',A,/,'* ',A,/'*')

      SELECTCASE (LRTinput) 
        CASE (1)  !TEFF & RSTAR given => L calculated
          WRITE (hOUT,'(A,I7,A)')   '* TEFF    =', INT(TEFF), ' K'
     >          // '       (INPUT)'
          WRITE (hOUT,'(A,F7.3,A)') '* RSTAR   =', RSTARSU, ' R_SUN'
     >              // '   (INPUT)'
          WRITE (hOUT,'(A,F7.3,A)') '* LOG L   =', XLOGL,   ' L_SUN' 
     >              // '   (CALCULATED FROM RSTAR AND TEFF)'
        CASE (2)  !TEFF & L given => RSTAR calculated
          WRITE (hOUT,'(A,I7,A)')   '* TEFF    =', INT(TEFF), ' K'
     >          // '       (INPUT)'
          WRITE (hOUT,'(A,F7.3,A)') '* LOG L   =', XLOGL,   ' L_SUN' 
     >              // '   (INPUT)'
          WRITE (hOUT,'(A,F7.3,A)') '* RSTAR   =', RSTARSU, ' R_SUN'
     >              // '   (CALCULATED FROM LUMINOSITY AND TEFF)'
        CASE (3)  !L & RSTAR given => TEFF calculated
          WRITE (hOUT,'(A,F7.3,A)') '* RSTAR   =', RSTARSU, ' R_SUN'
     >              // '   (INPUT)'
          WRITE (hOUT,'(A,F7.3,A)') '* LOG L   =', XLOGL,   ' L_SUN' 
     >              // '   (INPUT)'
          WRITE (hOUT,'(A,I7,A)')   '* TEFF    =', INT(TEFF), ' K'
     >          // '       (CALCULATED FROM RSTAR AND LUMINOSITY)'
        CASE (4)  !MSTAR & GEDD & TEFF given => L & RSTAR calculated
          WRITE (hOUT,'(A,I7,A)')   '* TEFF    =', INT(TEFF), ' K'
     >          // '       (INPUT)'
          WRITE (hOUT,'(A,F7.3,A)') '* LOG L   =', XLOGL,   ' L_SUN' 
     >              // '   (CALCULATED FROM MSTAR AND EDDINGTON_GAMMA)'
          WRITE (hOUT,'(A,F7.3,A)') '* RSTAR   =', RSTARSU, ' R_SUN'
     >              // '   (CALCULATED FROM LUMINOSITY AND TEFF)'
        CASE (5)  !MSTAR & GEDD & RSTAR given => TEFF & L calculated
          WRITE (hOUT,'(A,F7.3,A)') '* RSTAR   =', RSTARSU, ' R_SUN'
     >              // '   (INPUT)'
          WRITE (hOUT,'(A,F7.3,A)') '* LOG L   =', XLOGL,   ' L_SUN' 
     >              // '   (CALCULATED FROM MSTAR AND EDDINGTON_GAMMA)'
          WRITE (hOUT,'(A,I7,A)')   '* TEFF    =', INT(TEFF), ' K'
     >          // '       (CALCULATED FROM RSTAR AND LUMINOSITY)'
        CASE DEFAULT
          WRITE(hCPR,'(A,I5)') 'INVALID VALUE FOR LRTinput: ', LRTinput
          STOP 'FATAR ERROR detected in DECSTAR'
      ENDSELECT
            

C***  calculation of Rtrans if not done in DECSTAR
      IF (RTRANS .LE. .0)  
     >    RTRANS = RSTARSU * ( VFINAL / 2500. 
     >          / 10.**(XMDOT+0.5*ALOG10(DENSCON(1))+4.0))**(2./3.)
      IF (RTRANS > 999.) THEN
         FORMSTR = '(A,F7.0,A,F6.3,A)'
      ELSE
         FORMSTR = '(A,F7.3,A,F6.3,A)'

      ENDIF

      SELECTCASE (MDOTINPUT) 
        CASE (1)
          WRITE (hOUT,'(A,F7.3,A)') '* M-DOT   =', XMDOT, 
     >           ' DEX, IN M_SUN/YR   (INPUT)'
          WRITE (hOUT,FMT=FORMSTR) '* RTRANS  =', RTRANS, 
     >           ' = ', ALOG10(RTRANS), ' DEX   (CALCULATED FROM MDOT)'

        CASE (2)
          WRITE (hOUT,'(A,F7.3,A)') '* M-DOT   =', XMDOT, 
     >           ' DEX, IN M_SUN/YR   (CALCULATED FROM RTRANS)'
          WRITE (hOUT,FMT=FORMSTR) '* RTRANS  =', RTRANS, 
     >           ' = ', ALOG10(RTRANS), ' DEX   (INPUT)'

        CASE (3)
          WRITE (hOUT,'(A,F7.3,A)') '* M-DOT   =', XMDOT, 
     >           ' DEX, IN M_SUN/YR   (CALCULATED FROM MDTRANS)'
          WRITE (hOUT,FMT=FORMSTR) '* RTRANS  =', RTRANS, 
     >           ' = ', ALOG10(RTRANS), ' DEX   (CALCULATED FROM MDTRANS)'

        CASE (4)
          WRITE (hOUT,'(A,F7.3,A)') '* M-DOT   =', XMDOT, 
     >           ' DEX, IN M_SUN/YR   (CALCULATED FROM LOG Q)'
          WRITE (hOUT,FMT=FORMSTR) '* RTRANS  =', RTRANS, 
     >           ' = ', ALOG10(RTRANS), ' DEX   (CALCULATED FROM LOG Q)'

      ENDSELECT 
     

      WRITE (hOUT,'(A,I7,A)') '* VFINAL  =', NINT(VFINAL), ' KM/S'
      WRITE (hOUT,'(A,I7,A)') '* VDOP    =', NINT(VDOP)  , ' KM/S'


c***  tiefenabh. clumping...
C***  aufpassen, evtl. ueber den L-Bereich??
      IF (DENSCON(1) .NE. 1.) THEN 
         WRITE (hOUT,6) DENSCON(1), FILLFAC(1)
    6    FORMAT ('* DENSCON(1) =', F7.2, '    FILLFAC(1) =', F6.4)  
      ELSE 
         WRITE (hOUT,'(A)') '* MODEL WITHOUT CLUMPING'
      ENDIF

      IF (MASSORIGIN == 0) THEN
         WRITE (hOUT,'(A,F7.2, A)') 
     >         '* LOG G_GRAV =', GLOG, ' [CGS]  CALCULATED FROM:'  
         WRITE (hOUT,'(A)') '*   MASS-LUMINOSITY ' //
     >                   'RELATION FOR TYPE = ' // WRTYPE
     >                   // ' from ' // MLRELATION
         WRITE (hOUT,'(A,F6.2,A,A,F5.2,A,A,F5.2)')
     >          '*   MASS =', XMSTAR, ' M_SUN',  
     >                ' -- CALCULATED LOG G_EFF = ', GEFFLOG, ' [CGS]',
     >                ' VIA EDDINGTON_GAMMA = ', GEDD
      ELSEIF (MASSORIGIN == 1) THEN
         WRITE (hOUT,'(A,F7.2, A)') 
     >         '* MSTAR   =', XMSTAR, ' M_SUN   (INPUT)'  
         IF (GEddFix == 0) THEN
           WRITE (hOUT,'(A,F5.2, A,/,A,F5.2,A,A,F5.2)')
     >         '*   IMPLIED LOG G_GRAV = ', GLOG, ' [CGS]',
     >         '*   CALCULATED LOG G_EFF = ', GEFFLOG, ' [CGS]',
     >                ' VIA EDDINGTON_GAMMA = ', GEDD 
         ELSE
           IF (GEddFix == 1) THEN
             WRITE (hOUT,'(A,F7.2, A,/,A,F5.2,A,A,F5.2,A)')
     >         '* LOG G_EFF = ', GEFFLOG, ' [CGS]   (INPUT)',
     >         '*   IMPLIED LOG G_GRAV = ', GLOG, ' [CGS]',
     >                ' AND EDDINGTON_GAMMA = ', GEDD , ' (FIXED)'
           ELSE
             WRITE (hOUT,'(A,F7.2, A,/,A,F5.2,A,A,F5.2,A)')
     >         '* EDDINGTON_GAMMA = ', GEDD, '   (INPUT)',
     >         '*   IMPLIED LOG G_GRAV = ', GLOG, ' [CGS]',
     >                ' AND LOG G_EFF = ', GEFFLOG, ' [CGS]'
           ENDIF
         ENDIF
      ELSEIF (MASSORIGIN == 2) THEN
         WRITE (hOUT,'(A,F7.2, A)') 
     >         '* LOG G_GRAV =', GLOG, ' [CGS]   (INPUT)'  
         IF (GEddFix == 0) THEN
           WRITE (hOUT,'(A,F6.2, A,/,A,F5.2,A,A,F5.2)') 
     >         '*   IMPLIED STELLAR MASS = ', XMSTAR, ' M_SUN',
     >         '*   CALCULATED LOG G_EFF = ', GEFFLOG, ' [CGS]',
     >                ' VIA EDDINGTON_GAMMA = ', GEDD
         ELSE
           IF (GEddFix == 1) THEN
             WRITE (hOUT,'(A,F7.2, A,/,A,F6.2,A,A,F5.2,A)') 
     >         '* LOG G_EFF = ', GEFFLOG, ' [CGS]   (INPUT)',
     >         '*   IMPLIED STELLAR MASS = ', XMSTAR, ' M_SUN',
     >                ' AND EDDINGTON_GAMMA = ', GEDD, ' (FIXED)'
           ELSE
             WRITE (hOUT,'(A,F7.2, A,/,A,F6.2,A,A,F5.2,A)') 
     >         '* EDDINGTON_GAMMA = ', GEDD, '   (INPUT)',
     >         '*   IMPLIED STELLAR MASS = ', XMSTAR, ' M_SUN',
     >                ' AND LOG G_EFF = ', GEFFLOG,' [CGS]'
           ENDIF
         ENDIF
      ELSEIF (MASSORIGIN == 3) THEN
         IF (GEddFix == 0) THEN
           WRITE (hOUT,'(A,F7.2, A)') 
     >         '* LOG G_EFF =', GEFFLOG, ' [CGS]   (INPUT)'  
           WRITE (hOUT,'(A,F5.2, A, A, F5.2)') 
     >         '*   CALCULATED LOG G_GRAV = ', GLOG, ' [CGS]',
     >                ' VIA EDDINGTON_GAMMA = ', GEDD
           WRITE (hOUT,'(A,F6.2, A)') 
     >         '*   IMPLIED STELLAR MASS = ',  XMSTAR, ' M_SUN'
         ELSE
           WRITE (hOUT,'(A,F7.2, A)') 
     >         '* LOG G_GRAV =', GLOG, ' [CGS]   (INPUT)'  
           IF (GEddFix == 1) THEN
             WRITE (hOUT,'(A,F7.2, A)') 
     >         '* LOG G_EFF = ', GEFFLOG, ' [CGS]   (INPUT)'           
             WRITE (hOUT,'(A,F6.2, A, A, F5.2,A)') 
     >         '*   IMPLIED STELLAR MASS = ',  XMSTAR, ' M_SUN',
     >                ' AND EDDINGTON_GAMMA = ', GEDD, ' (FIXED)'
           ELSE
             WRITE (hOUT,'(A,F7.2, A)') 
     >         '*  EDDINGTON_GAMMA = ', GEDD, '   (INPUT)'           
             WRITE (hOUT,'(A,F6.2, A, A, F5.2,A)') 
     >         '*   IMPLIED STELLAR MASS = ',  XMSTAR, ' M_SUN',
     >                ' AND LOG G_EFF = ', GEFFLOG, ' [CGS]'
           ENDIF
        ENDIF
      ELSE
         STOP 'PRIPARAM: ERROR - UNKNOWN CODE FOR MASSORIGIN'
      ENDIF

      IF (VTURB > .0) THEN
        WRITE (hOUT,'(A,F7.2,A)') '* VTURB   =', VTURB, ' KM/S'
      ENDIF

      WRITE (hOUT,'(A,F7.2,A,F8.2,A)') '* RMAX    =', RMAX, ' RSTAR =', 
     >                          RMAX*RSTARSU, ' R_SUN'


      WRITE(hOUT,80) 
   80 FORMAT ('*',/,80('*'),//)


      RETURN
      END
