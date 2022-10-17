      SUBROUTINE READ_LINECARD_PARAMETERS (KARTE_REST, XLAMBL, 
     >             BROAD, LINPROBL, AVOIGTBL)
C*************************************************************************
C***  This subroutine reads additional (optional) parameters :
C***  - Wavelength (first parameter, no keyword)
C***  - 'AIR' or 'VAC'  (second parameter if wavelenth given; defaults)
C***  - VOIGT (radiation-damping on for this line only)
C***  - VOIGT-parameter a (optional, otherwise calculated from lifetimes)
C***  CALLED FROM FORMAL - PREFORM 
C*************************************************************************

      IMPLICIT NONE

      INTEGER NPAR, I, IDX
      REAL XLAMBL, XLAM2, AVOIGTBL

      CHARACTER*80 KARTE_REST
      CHARACTER*20 ACTPAR
      CHARACTER*8  LINPROBL

      LOGICAL BAIR, BVAC, BROAD

      BAIR = .FALSE.
      BVAC = .FALSE.
      LINPROBL = ''
      AVOIGTBL = -1.
 
C***  Count supplementary parameters
      CALL SARGC (KARTE_REST,NPAR)

C***  Input empty? (Should not happen!)
      IF (NPAR .LE. 0) RETURN

C***  reference wavelength explicitely given as 1st parameter?
      CALL SARGV(KARTE_REST,1,ACTPAR)
      IF (ACTPAR(1:1) .LT. 'A') THEN 
         READ (ACTPAR,'(F10.0)', ERR=100) XLAMBL

C***     Wavelength might be followed by keyword "AIR" or "VAC" 
         IF (NPAR .GE. 2) THEN
            CALL SARGV (KARTE_REST, 2, ACTPAR)
            IF (ACTPAR .EQ. 'AIR') BAIR = .TRUE.
            IF (ACTPAR .EQ. 'VAC') BVAC = .TRUE.
         ENDIF

C***     All calculations will be done for vacuum wavelengths
C***     Therefore, explicitely given wavelenths are converted if appropriate 
         IF (BAIR .OR. (.NOT. BVAC .AND.
     >       XLAMBL .GT. 2000. .AND. XLAMBL .LT. 20000.)) THEN
             XLAM2  = XLAMBL * XLAMBL
             XLAMBL = XLAMBL * (1.0 + 2.735182E-4 + 131.4182
     >                 / XLAM2 + 2.76249E8 / (XLAM2*XLAM2))
          ENDIF
      ENDIF

C***  Voigt parameter definition by input (in fact outdated!) 
      IF (BROAD) THEN
         DO 12 I=1, NPAR
            CALL SARGV (KARTE_REST, I, ACTPAR)
            IF (ACTPAR .EQ. 'VOIGT') THEN
               IF (I .LT. NPAR) THEN
                  CALL SARGV (KARTE_REST, I+1, ACTPAR)
cc                  IF ((ACTPAR(1:1) .LT. 'A'))
                  READ (ACTPAR,'(F20.0)', ERR=110) AVOIGTBL
                  LINPROBL = 'VOIGT   '
               ELSE
                  WRITE (0,*) '*** WARNING: VOIGT keyword without ' //
     >            'value ignored' 
               ENDIF
            ENDIF
   12    CONTINUE
      ENDIF

      RETURN

C***  ERROR Branches ******************************************

  100 WRITE (0,*) '*** ERROR: Wavelength cannot be decoded as number'
      GOTO 130

  110 WRITE (0,*) '*** ERROR: AVOIGT value cannot be decoded as number'
      GOTO 130

  130 WRITE (0,*) '*** The problematic string is: ' // 
     >      KARTE_REST(:IDX(KARTE_REST))
      STOP '*** FATAL ERROR detected bu subr. READ_LINECARD_PARAMETERS'

      END
