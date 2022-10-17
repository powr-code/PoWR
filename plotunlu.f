      SUBROUTINE PLOTUNLU (KANAL, PLOTOPT, MODHEAD, JOBNUM, 
     >                    ND, T, TNEW, DUNLU_LOC, DUNLU_TB, bTDIFFUS,
     >                    DTLOCAL, DTINT, DTRMAX, DTKUBAT)
C*****************************************************************
C***  PLOT of temperature corrections from UNSOELD-LUCY procedure
C*****************************************************************

      IMPLICIT NONE

      INTEGER, PARAMETER :: NDMAX = 200
      INTEGER, INTENT(IN) :: ND, KANAL, JOBNUM

      CHARACTER PLOTOPT*(*)
      CHARACTER(100) :: MODHEAD
      CHARACTER(110) :: HEADLINE
      CHARACTER(20) :: CUROPT
      CHARACTER(8) :: CENTER

      REAL, DIMENSION(NDMAX) :: X, TCORR
      REAL, DIMENSION(ND) :: DTLOCAL, DTINT, DTRMAX, DTKUBAT, T, TNEW

      INTEGER :: L, I, NPAR, NDin
      REAL :: XLABOFF, DUNLU_LOC, DUNLU_TB
      
      LOGICAL :: bTDIFFUS, bPlotKUBAT, bPlotLOCAL

C***  STEBOL = STEFAN-BOLTZMANN CONSTANT / PI (ERG/CM**2/SEC/STERAD/KELVIN**4
      REAL, PARAMETER :: STEBOL = 1.8046E-5

      IF (ND > NDMAX) THEN
        WRITE (0, '(A)') 'DIMENSION INSUFFICIENT - UNLU-PLOT SUPPRESSED'
        WRITE (0, '(A)') 'NON-FATAL ERROR IN SUBROUTINE PLOTUNLU'
        RETURN
      ENDIF

C***  Decode possible options
      bPlotKUBAT = DUNLU_TB  .GT. .0
      bPlotLOCAL = DUNLU_LOC .GT. .0
      CALL SARGC (PLOTOPT, NPAR)
      IF (NPAR > 2) THEN
        parloop: DO I=3, NPAR
          CALL SARGV (PLOTOPT, I, CUROPT)
          SELECTCASE(CUROPT)
            CASE ('LOCAL', 'TE')
              bPlotLOCAL = .TRUE.
            CASE ('TBALANCE', 'TB', 'KUBAT')
              bPlotKUBAT = .TRUE.
            CASE ('ALLTERMS', 'ALL')
              bPlotKUBAT = .TRUE.
              bPlotLOCAL = .TRUE.
          ENDSELECT
        ENDDO parloop
      ENDIF
      
      
      CENTER = CHAR(92) // 'CENTER' // CHAR(92)
      CALL JSYMSET ('G2','TRANSFER')

      HEADLINE = 'UNLU: M'//MODHEAD(13:)
      WRITE (HEADLINE(90:), '(A8,I7)') ' JOB No.', JOBNUM

      DO L=1, ND
         X(L) = FLOAT(L)
         TCORR(L) = (TNEW(L) - T(L)) / 1000.
      ENDDO
      
      XLABOFF = 13.

      WRITE (KANAL,'(A)') 'PLOT   :' // HEADLINE
      WRITE (KANAL, '(A)') '\DEFINECOLOR=9 0. .5 0.'
      WRITE (KANAL, '(A)') '\COLOR=8'
      WRITE (KANAL, '(A)') 
     >         'KASDEF LINUN XMIN 0. XMAX 0.  0. 0.'

        IF (bPlotKUBAT .AND. bPlotLOCAL) THEN
          WRITE (KANAL,'(A)') '\COLOR=4'
          WRITE (KANAL,'(A,F8.4,A)') 
     >      '\LINUN 0. YMIN  5. YMIN  ', XLABOFF, ' 3.0'
          WRITE (KANAL,'(A)') '\COLOR=9'
          WRITE (KANAL,'(A,F8.4,A)') 
     >     '\LINUN 0. YMIN 2.5 YMIN  ', XLABOFF, ' 3.0'
          WRITE (KANAL,'(A,F8.4,A)') 
     >      '\LUN   7. YMIN ', XLABOFF, ' M3.0 0.3 &1Local (&9TB&1,&4TE&1)'
        ELSEIF (bPlotKUBAT) THEN
          WRITE (KANAL,'(A)') '\COLOR=9'
          WRITE (KANAL,'(A,F8.4,A)') 
     >      '\LINUN 0. YMIN  5. YMIN  ', XLABOFF, ' 3.0'
          WRITE (KANAL,'(A,F8.4,A)') 
     >      '\LUN   7. YMIN ', XLABOFF, ' M3.0 0.3 Local (TB)'
        ELSEIF (bPlotLOCAL) THEN
        WRITE (KANAL,'(A)') '\COLOR=4'
          WRITE (KANAL,'(A,F8.4,A)') 
     >      '\LINUN 0. YMIN  5. YMIN  ', XLABOFF, ' 3.0'
          WRITE (KANAL,'(A,F8.4,A)') 
     >      '\LUN   7. YMIN ', XLABOFF, ' M3.0 0.3 Local (TE)'
        ENDIF        
        
        WRITE (KANAL,'(A)') '\COLOR=1'
        WRITE (KANAL,'(A,F8.4,A)') 
     >    '\LINUN 0. YMIN  5. YMIN  ', XLABOFF, ' 2.4 SYMBOL=9 SIZE=0.2'
        WRITE (KANAL,'(A,F8.4,A)') 
     >    '\LUN   7. YMIN          ', XLABOFF, ' M2.4 0.3 Int'

        WRITE (KANAL,'(A)') '\COLOR=1'
        WRITE (KANAL,'(A,F8.4,A)') 
     >    '\LINUN 0. YMIN  5. YMIN  ', XLABOFF, ' 1.8 SYMBOL=10 SIZE=0.2'
        WRITE (KANAL,'(A,F8.4,A)') 
     >    '\LUN   7. YMIN          ', XLABOFF, ' M1.8 0.3 Rmax'

      WRITE (KANAL,'(A)') '\COLOR=2'
      WRITE (KANAL,'(A)') '\PEN=4'
      WRITE (KANAL,'(A,F8.4,A)') 
     >  '\LINUN 0. YMIN  5. YMIN  ', XLABOFF, ' 1.2'
      WRITE (KANAL,'(A,F8.4,A)') 
     >  '\LUN   7. YMIN          ', XLABOFF, ' M1.2 0.3 T-Corr.'
      WRITE (KANAL,'(A)') '\COLOR=1'
      WRITE (KANAL,'(A)') '\PEN=1'

      CALL PLOTANFS (KANAL,HEADLINE, '&E'//HEADLINE,
     >        CENTER//'Depth Index L',
     >        CENTER//'Temperature Correction / kK',
     >        .0, .0, .0, .0, .0, .0,
     >        .0, .0, .0, .0, .0, .0,
     >        X, TCORR, ND, 'PEN=4 COLOR=2')
      
      IF (.NOT. bTDIFFUS) THEN
        NDin = ND
      ELSE
        NDin = ND-1
      ENDIF
      
      IF (bPlotLOCAL) THEN
        CALL PLOTCONS (KANAL, X, DTLOCAL, NDin, 'SYMBOL=5 COLOR=4')
      ENDIF
      IF (bPlotKUBAT) THEN
        CALL PLOTCONS (KANAL, X, DTKUBAT, NDin, 'SYMBOL=5 COLOR=9')
      ENDIF
      CALL PLOTCONS (KANAL, X, DTINT , NDin, 'SYMBOL=9 SIZE=0.2')
      CALL PLOTCONS (KANAL, X, DTRMAX, NDin, 'SYMBOL=10 SIZE=0.2')

      RETURN
      END
