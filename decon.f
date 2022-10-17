      SUBROUTINE DECON (LSOPA,LSINT,IFLUX,JOBMAX,MODHIST, 
     >                  BUNLU, IVERS)
C***********************************************************************
C***  DECODING INPUT OPTIONS, CALLED FROM WRCONT *******************************
C***********************************************************************

      IMPLICIT NONE

      INTEGER :: I, NPAR, IFLUX, JOBMAX, IVERS,
     >           LSOPA, LSINT, IDX
      INTEGER, DIMENSION(1) :: MODHIST

      CHARACTER(14), DIMENSION(5) :: ACTPAR
      CHARACTER(80) :: KARTE

      LOGICAL BUNLU
      !File and channel handles (=KANAL)
      INTEGER, PARAMETER :: hCPR = 0 !write to wruniqX.cpr (stderr)

 
C***  DEFAULT VALUES
      LSOPA=-1
      LSINT=-1
      IFLUX=-1
      JOBMAX=-1
      BUNLU = .FALSE.
      IVERS = 4
 
      OPEN (1,FILE='CARDS', STATUS='UNKNOWN')
      REWIND (1)

      DO !--- Loop over all CARDS lines ---

        READ (1, '(A)', END=1) KARTE

        CALL SARGC(KARTE,NPAR)
        IF ( NPAR .LT. 1) CYCLE
        IF ( NPAR .GT. 5) NPAR = 5

        DO I=1, NPAR
          CALL SARGV(KARTE,I,ACTPAR(I))
        ENDDO

C***    PRINT options ************************

        IF (ACTPAR(1) .EQ. 'PRINT') THEN
C                           =====
          IF (ACTPAR(2) .EQ. 'FLUX') THEN
C                             ====
             IFLUX=1

          ELSE IF (ACTPAR(2) .EQ. 'INT') THEN
C                                  ===
             READ (ACTPAR(3), *, ERR=90) LSINT
             IF (LSINT .EQ. 0) LSINT=1

          ELSE IF (ACTPAR(2) .EQ. 'OPA ') THEN
C                                  ===
             READ (ACTPAR(3), *, ERR=90) LSOPA
             IF (LSOPA .EQ. 0) LSOPA=1

          ENDIF

C***    Other (not PRINT) options *****************************

        ELSE IF (ACTPAR(1) .EQ. 'JOBMAX') THEN
C                                ======
          READ (ACTPAR(2), *, ERR=90) JOBMAX

        ELSE IF (ACTPAR(1) .EQ. 'UNLUTEC') THEN
C                                =======
          BUNLU=.TRUE.

        ELSE IF (ACTPAR(1) .EQ. 'OB-VERS') THEN
C                                =======
          IF (NPAR .EQ. 2 .OR. ACTPAR(3) .EQ. 'WRCONT') THEN
            READ (ACTPAR(2), *, ERR=90) IVERS
          ENDIF

        ENDIF  
 
      ENDDO

C***  END-OF-FILE REACHED *****************************************
    1 CONTINUE
      CLOSE (1)
      RETURN
 
C***  ERROR EXITS **************************************************

   90 WRITE (hCPR,*) '*** DECON: FATAL ERROR WHEN DECODING NUMBER'
      WRITE (hCPR,*) '*** THE ERROR OCCURED IN THE FOLLOWING LINE:'
      WRITE (hCPR,*) KARTE(:IDX(KARTE))
      STOP 'ERROR in WRCONT'
      
      END
