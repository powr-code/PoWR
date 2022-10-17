      SUBROUTINE DECVELPAR(KARTE, VFINAL, VMIN, BETA, RMAX)
C***  Decodes the VELPAR line form the CARDS file      

      IMPLICIT NONE

      CHARACTER(40) :: TRYPAR
      CHARACTER(40), DIMENSION(20) :: CURPAR
      CHARACTER(100) :: KARTE
      REAL :: VFINAL, VMIN, BETA, RMAX

      INTEGER :: NPAR, i, IERR
      
      LOGICAL :: bOldDecode 
      LOGICAL, DIMENSION(4) :: bParamFound

      !File and channel handles (=KANAL)
      INTEGER, PARAMETER :: hOUT = 6        !write to wruniqX.out (stdout)
      INTEGER, PARAMETER :: hCPR = 0        !write to wruniqX.cpr (stderr)
     
      bOldDecode = .FALSE.
      CALL SARGC (KARTE, NPAR)
      IF (NPAR < 5) THEN
        WRITE (hCPR,'(A)') '*** VELPAR: NOT ENOUGH PARAMETERS'
        STOP '*** FATAL ERROR WHILE DECODING VELPAR CARDS-LINE'
      ENDIF

      !New decoding => allows flexible format and modern syntax
      bParamFound = .FALSE.
      DO i=1, NPAR
        CALL SARGV(KARTE,i,CURPAR(i))
      ENDDO
      IF (NPAR > 2) THEN
        DO i=2, NPAR 
         SELECTCASE (CURPAR(i))
          CASE ('VFINAL')          
            IF (NPAR >= (i+1)) THEN
              TRYPAR = CURPAR(i+1)
              IF (TRYPAR == '(KM/S)') THEN
                IF (NPAR >= (i+2)) THEN
                  TRYPAR = CURPAR(i+2)
                ELSE
                  GOTO 92
                ENDIF
              ENDIF
              READ (TRYPAR, '(F10.0)', IOSTAT=IERR, ERR=92) VFINAL
              IF (IERR == 0) THEN
                bParamFound(1) = .TRUE.
              ENDIF      
            ENDIF
          CASE ('VMIN')
            IF (NPAR >= (i+1)) THEN
              READ (CURPAR(i+1), '(F10.0)', IOSTAT=IERR) VMIN
              IF (IERR == 0) THEN
                bParamFound(2) = .TRUE.
              ENDIF                  
            ENDIF
          CASE ('BETA')
            IF (NPAR >= (i+1)) THEN
              READ (CURPAR(i+1), '(F10.0)', IOSTAT=IERR) BETA
              IF (IERR == 0) THEN
                bParamFound(3) = .TRUE.
              ENDIF                  
            ENDIF
          CASE ('RMAX')
            IF (NPAR >= (i+1)) THEN
              READ (CURPAR(i+1), '(F10.0)', IOSTAT=IERR) RMAX
              IF (IERR /= 0) THEN
                WRITE(hCPR,'(A)') '*** DECVELPAR: CANNOT READ RMAX'
                STOP 'ERROR'
              ELSE
                bParamFound(4) = .TRUE.
              ENDIF                  
            ENDIF
         ENDSELECT
        ENDDO
      ENDIF

      DO i=1, 4 
        !One or more parameters have not been found => switch to old decoding
        IF (.NOT. bParamFound(i)) THEN
          WRITE (hCPR,*) '*** DECVELPAR: Old VELPAR decoding used'
          bOldDecode = .TRUE.
        ENDIF
      ENDDO

      IF (bOldDecode) THEN
         READ (KARTE,19,ERR=99) VFINAL,VMIN,BETA,RMAX
   19    FORMAT(22X,F7.0,5X,F6.0,5X,F4.0,5X,F6.0)
      ENDIF
      
      RETURN
      
C***  FATAL ERROR CODES      
      
   92 WRITE(hCPR,'(A)') '*** DECVELPAR: CANNOT READ VFINAL IN:'
      WRITE (hCPR,*) KARTE
      STOP 'ERROR'
      
   99 WRITE (hCPR,*)
     >   'DECVELPAR: ERROR WHILE DECODING THE FOLLOWING CARDS-LINE:'
      WRITE (hCPR,*) KARTE
      STOP 'ERROR'
      
      END
      