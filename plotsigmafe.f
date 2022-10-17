      SUBROUTINE PLOTSIGMAFE(PLOTOPT, MODHEAD, JOBNUM,
     >                       SIGMAFE, INDRB, MAXFEIND, 
     >                       LASTIND, LASTFE, IFRBSTA, IFRBEND, 
     >                       LEVEL, N, INDNUP, INDLOW,
     >                       INDEXMAX, NFEREADMAX,
     >                       VDOPFE, DXFE, XLAM0FE, bOwn)
      
      IMPLICIT NONE
      
      INTEGER, INTENT(IN) :: MAXFEIND, INDEXMAX, NFEREADMAX,
     >                       LASTIND, LASTFE, N, JOBNUM
      REAL, INTENT(IN) :: XLAM0FE, DXFE, VDOPFE
      INTEGER, DIMENSION(LASTIND), INTENT(IN) :: INDNUP, INDLOW
      REAL, DIMENSION(INDEXMAX), INTENT(IN) :: SIGMAFE
      CHARACTER(10), DIMENSION(N) :: LEVEL
      CHARACTER PLOTOPT*(*)
      
      INTEGER, DIMENSION(MAXFEIND) :: INDRB, IFRBEND, IFRBSTA
      
      REAL, DIMENSION(NFEREADMAX) :: X, Y

      CHARACTER(5) :: CIFE, CIND
      CHARACTER(10) :: LEV1, LEV2
      CHARACTER(100) :: MODHEAD
      CHARACTER(110) :: HEADLINE
      
      REAL :: XLOGSTEP
      INTEGER :: NPTS, IndexFE, IndexFELAM, I, ILAM, IND,
     >           hPLOT, ILEV1, ILEV2, IERR, NPAR, LOW, NUP
      
      LOGICAL :: bOwn
      
      INTEGER, EXTERNAL :: IDX

      !Physical constants
      REAL, PARAMETER :: CLIGHT = 2.9979E10     !Speed of Light in cm/s

      !File and channel handles (=KANAL)
      INTEGER, PARAMETER :: hOUT = 6        !write to wruniqX.out (stdout)
      INTEGER, PARAMETER :: hCPR = 0        !write to wruniqX.cpr (stderr)
      

C***  DECODE LEVELS (can be index number or name)      
      ILEV1 = 0
      ILEV2 = 0
      CALL SARGC (PLOTOPT, NPAR)
      IF (NPAR == 4) THEN
        CALL SARGV (PLOTOPT, 3, LEV1)
        READ (LEV1(:IDX(LEV1)), '(I4)', IOSTAT=IERR) ILEV1
        IF (IERR > 0) THEN
          DO I=1, N
            IF (LEVEL(I) == LEV1) THEN
              ILEV1 = I
              EXIT
            ENDIF
          ENDDO      
        ELSE
          LEV1 = LEVEL(ILEV1)          
        ENDIF
        CALL SARGV (PLOTOPT, 4, LEV2)
        READ (LEV2(:IDX(LEV2)), '(I4)', IOSTAT=IERR) ILEV2
        IF (IERR > 0) THEN
          DO I=1, N
            IF (LEVEL(I) == LEV2) THEN
              ILEV2 = I
              EXIT
            ENDIF
          ENDDO      
        ELSE
          LEV2 = LEVEL(ILEV2)          
        ENDIF
      ELSE
        WRITE (hCPR,*) 'PLOTSIGMAFE: Invalid parameters ************'
        WRITE (hCPR,*) '***** The following plot was aborted:'
        WRITE (hCPR,*) PLOTOPT(:IDX(PLOTOPT))        
        RETURN
      ENDIF            

      IF (ILEV1 == 0 .OR. ILEV2 == 0) THEN
        WRITE (hCPR,*) 'PLOTSIGMAFE: Invalid levels ************'
        WRITE (hCPR,*) '***** The following plot was aborted:'
        WRITE (hCPR,*) PLOTOPT(:IDX(PLOTOPT))        
        RETURN
      ENDIF
      
      IND = 0
      DO I=1, LASTIND
        IF (INDNUP(I) == ILEV1 .AND. INDLOW(I) == ILEV2) THEN
          LOW = ILEV2
          NUP = ILEV1
          IND = I
          EXIT
        ELSEIF (INDNUP(I) == ILEV2 .AND. INDLOW(I) == ILEV1) THEN
          LOW = ILEV1
          NUP = ILEV2
          IND = I
          EXIT
        ENDIF
      ENDDO
      
C***  Check if transition is found and really an iron bound-bound transition
      IF (IND < (LASTIND - LASTFE)) THEN
        WRITE (hCPR,*) 'PLOTSIGMAFE: Invalid transition ************'
        WRITE (hCPR,*) '***** The following plot was aborted:'
        WRITE (hCPR,*) PLOTOPT(:IDX(PLOTOPT))        
        RETURN
      ENDIF
      
      
C***  PLOT FILE INITIALIZATION
      IF (bOwn) THEN
        hPLOT = 37       !write to seperate plot file
        OPEN (hPLOT, FILE='sigmafe.plot', STATUS='UNKNOWN')
        WRITE (hCPR,*) 'SIGMAFE PLOTTED FOR IND = ', IND
      ELSE
        hPLOT = 2       !write to file PLOT (becomes steal.plot)
        OPEN (hPLOT, FILE='PLOT', STATUS='UNKNOWN')
        CALL JSYMSET ('G2','TRANSFER')
        CALL REMARK ('SIGMAFE TO BE ROUTED')
      ENDIF
      
      HEADLINE = 'SIGMAFE:'//MODHEAD(13:)
C      HEADLINE = 'SIGMAFE: to be done'
      
      WRITE (hPLOT, '(A,I8)') '*NEXTPLOT: SIGMAFE ', IND
      WRITE (hPLOT, '(A)') 'PLOT: ' // HEADLINE
      WRITE (hPLOT, '(A)') '\FONT=HELVET'      
      WRITE (hPLOT,'(A)') '\COLOR=1'
      WRITE (hPLOT,'(A)') '\PEN=1'
      
      XLOGSTEP = LOG10(1. + VDOPFE*1.E5*DXFE/CLIGHT)
      IndexFE = IND - (LASTIND - LASTFE)        !convert total IND to fe-only IND number
      IndexFELAM = INDRB(IndexFE)               !starting index inside SIGMAFE for current FeIND
      DO 
        I = IndexFELAM - INDRB(IndexFE) + 1
        ILAM = IFRBSTA(IndexFE) + I - 1
C        WRITE (hCPR, *) I, ILAM, XLAM0FE * 10.**(ILAM*XLOGSTEP)
        X(I) = XLAM0FE * 10.**(ILAM*XLOGSTEP)       !lambda from index
        Y(I) = SIGMAFE(IndexFELAM) / 1.E-15          !plot in 10^(-15) cm^2
C        Y(I) = SIGMAFE(IndexFELAM)
        IF (ILAM == IFRBEND(IndexFE)) EXIT
        IndexFELAM = IndexFELAM + 1
      ENDDO
      NPTS = I      

      WRITE (UNIT=CIFE, FMT='(I5)') IndexFE
      WRITE (UNIT=CIND, FMT='(I5)') IND
      WRITE (hPLOT,'(5A)') 
     >  '\LUN XMAX YMAX R-0.5 U-0.5 0.2 Fe-Transition ',
     >    TRIM(ADJUSTL(CIFE)), ' (', TRIM(ADJUSTL(CIND)), ')'
      WRITE (hPLOT,'(2A)') 
     >  '\LUN XMAX YMAX R-0.5 U-1.0 0.25 Low: ', LEVEL(LOW)
      WRITE (hPLOT,'(2A)') 
     >  '\LUN XMAX YMAX R-0.5 U-1.5 0.25 Up: ', LEVEL(NUP)
      
      CALL PLOTANFS (hPLOT,HEADLINE, '&E'//HEADLINE,
     >        '\CENTER\#l# / \A',
     >        '\CENTER\#s#&TLU&M / 10&H-15&M\,cm&H2&M',
     >        .0, .0, .0, .0, .0, .0,
     >        .0, .0, .0, .0, .0, .0,
     >        X, Y, NPTS, 'COLOR=4')
     
C      CALL PLOTCONS (hPLOT, X, TCORR,  ND,   'PEN=4 COLOR=2')

C      IF (bOwn) THEN
C        CLOSE (hPLOT)
C      ENDIF


      RETURN

      END
