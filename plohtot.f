      SUBROUTINE PLOHTOT (HTOT, ND, JOBNUM, MODHEAD, TEFF)
C**********************************************************************
C***  PLOT OF THE CONTINUUM FLUX(EXPRESSED AS T-EFF)
C***       AS FUNCTION OF DEPTH-INDEX L
C***  THIS OUTPUT MAY BE USEFUL TO INSPECT FLUX CONSERVATION.
C***  CALLED FROM: COMO
C**********************************************************************

C***  DARR AND TARR ARE USED TO PLOT THE STRATIFICATION OF HTOT
      PARAMETER (NDMAX = 200)
      DIMENSION HTOT (ND), DARR(NDMAX), TARR(NDMAX)

C***  STEBOL = STEFAN-BOLTZMANN CONSTANT / PI (ERG/CM**2/SEC/STERAD/KELVIN**4
      DATA STEBOL /1.8046E-5/
      CHARACTER MODHEAD*100, HEADLINE*100

      IF (ND .GT. NDMAX) THEN
        WRITE (0,*) 'DIMENSION OF DARR AND TARR INSUFFICIENT!'
        WRITE (0,'(2(A,I3))') ' ND=', ND, '  NDMAX=', NDMAX
        STOP 'ERROR IN SUBR. PLOHTOT'
      ENDIF

      OPEN (UNIT=1, FILE='PLOT', STATUS='UNKNOWN')

      CALL JSYMSET ('G2', 'TRANSFER')

      DO 2 L=2,ND
        TRAD = (4.*ABS(HTOT(L))/STEBOL)**0.25
        IF (HTOT(L) .LE. .0) TRAD = -TRAD
        DARR(L-1) = FLOAT(L)
        TARR(L-1) = TRAD / 1000.
    2 CONTINUE

      HEADLINE = 'M'//MODHEAD(13:)
      WRITE (HEADLINE(90:), '(A8,I3)') ' JOB No.', JOBNUM

      XMIN = 0.
      XMAX = FLOAT(ND)
      YMIN = TEFF/ 1000. - 5.
      YMAX = TEFF/ 1000. + 5.

      CALL PLOTANF(1, 'CONTINUUM FLUX' , '&E'//HEADLINE, 
     >             'Depth Index L', 'Continuum Flux T&Trad&M / kK', 
     >             0., XMIN, XMAX, 5., 10., 0., 
     >             0., YMIN, YMAX, .2,  1., 0., 
     >             DARR, TARR, ND-1, 5)

C***  CONTINUUM
      DARR(1) = XMIN
      DARR(2) = XMAX
      TARR(1) = TEFF / 1000.
      TARR(2) = TEFF / 1000.
      CALL PLOTCON (1, DARR, TARR, 2, 5)      

      CLOSE(1)

      RETURN
      END
