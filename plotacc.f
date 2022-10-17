      SUBROUTINE PLOTACC(PLOTOPT, AGRAV, AMECH, ARAD, APRESS, ACONT, 
     >                   ATHOM, WORKRATIO, VELO, RADIUS, ND, ATMEAN,
     >                   RNE, RCON, T, TEFF, RSTAR, XMU, XMSTAR,
     >                   Rcritical, bFULLHYDROSTAT, 
     >                   MODHEAD, JOBNUM, hPLOT)
C******************************************************************************
C***  DIRECT TRANSFER OF HSUM PLOT
C***  TOTAL (FREQUENCY-INTEGRATED) FLUX versus DEPTH INDEX
C***  for both, continuum, line and sum of these
C******************************************************************************
 
      IMPLICIT NONE
      INCLUDE 'interfacebib.inc'

      INTEGER, PARAMETER :: NDMAX = 200

      INTEGER, INTENT(IN) :: ND, JOBNUM, hPLOT

      REAL, DIMENSION(ND), INTENT(IN) :: AGRAV, AMECH, ARAD, APRESS, T,
     >                                   ACONT, ATHOM, VELO, RADIUS, 
     >                                   RNE, XMU
      REAL, INTENT(IN) :: WORKRATIO, RCON, ATMEAN, XMSTAR, RSTAR, TEFF, 
     >                    Rcritical
     
      LOGICAL, INTENT(IN) :: bFULLHYDROSTAT

      CHARACTER(110) :: MODHEAD, HEADLINE
      CHARACTER(8) :: CENTER, CNORM
      CHARACTER(40) :: XTEXT
      CHARACTER(20) :: CUROPT, XAXISMODE
      CHARACTER PLOTOPT*(*)

      REAL, DIMENSION(NDMAX) :: X, Y1, Y2, Y3, Y4, Y5, 
     >                          RI, ANORM, VMACH, VSCRATCH

      INTEGER :: I, L, NPAR, NDIN, Lcand
      REAL :: XMIN, XMAX, YMIN, YMAX, YMINVAL, XOFF, VINT,
     >        GEDDL, XLSTAR, XLSTARS, RNEINT, Xcrit,
     >        XABST, XTICK, XLENGTH, XICON, Rsonic, Vsonic, Xsonic

      LOGICAL :: bNormalizeToGEFF

      INTEGER, EXTERNAL :: IDX

      !Physical constants
      REAL, PARAMETER :: PI4 = 12.5663706144    !PI4 = 4*PI
      REAL, PARAMETER :: STEBOL = 5.6705E-5    !STEFAN-BOLTZMANN CONSTANT (CGS-UNITS)
      REAL, PARAMETER :: XLSUN = 3.85E33        !Solar Luminosity (CGS-Units)
      REAL, PARAMETER :: XMSUN = 1.989E33       !XMSUN = Solar Mass (g)
      REAL, PARAMETER :: RGAS = 8.3145E7        !Gas Constant (CGS)

      !File and channel handles (=KANAL)
      INTEGER, PARAMETER :: hOUT = 6        !write to wruniqX.out (stdout)
      INTEGER, PARAMETER :: hCPR = 0        !write to wruniqX.cpr (stderr)

      IF (NDMAX < ND) THEN
         WRITE (hCPR,'(A)') 'PLOTACC: NON-FATAL ERROR ******'
         WRITE (hCPR,'(A)') 'PLOTACC: DIMENSION NDMAX INSUFFICIENT'
         WRITE (hCPR,'(2(A,I4))') '  NDMAX = ', NDMAX, ',   ND = ', ND
         WRITE (hCPR,'(A)') 'PLOTACC: PLOT SKIPPED'
         RETURN
      ENDIF

      CENTER = '\CENTER\'

C***  Decode possible options for x-axis mode
      XAXISMODE = 'DEPTHINDEX'
      bNormalizeToGEFF = .FALSE.
      CALL SARGC (PLOTOPT, NPAR)
      IF (NPAR > 2) THEN
        DO I=3, NPAR
          CALL SARGV (PLOTOPT, I, CUROPT)
          SELECTCASE (CUROPT)
            CASE ('DEPTHINDEX', 'VELOCITY', 'RADIUS',
     >            'L', 'R', 'V')
C***          simple keyword option mode            
              XAXISMODE = CUROPT
            CASE ('XAXISMODE', 'X')
              IF (NPAR >= (I+1)) THEN
                CALL SARGV (PLOTOPT, I+1, XAXISMODE)
              ENDIF
            CASE ('GEFF')
              bNormalizeToGEFF = .TRUE.
          ENDSELECT
        ENDDO
      ELSE 
C***    Fallback for unique labeling inside steal.plot / wruniq.plot
        PLOTOPT = PLOTOPT(:IDX(PLOTOPT)) // ' L'
      ENDIF

      XLSTAR = PI4 * STEBOL * RSTAR*RSTAR * TEFF*TEFF*TEFF*TEFF
      XLSTARS = XLSTAR / XLSUN

      ANORM(ND) = 0.
      DO L=1, ND-1
        RI(L) = 0.5 * ( RADIUS(L) + RADIUS(L+1) )      
        IF (bNormalizeToGEFF) THEN
          !calculate depth-dependend GEFF 
          IF (bFULLHYDROSTAT) THEN
            !calculate geff with full a_rad
            GEDDL = MIN(ARAD(L)/AGRAV(L), 0.9)
          ELSE
            !calculate geff with 
            RNEINT = 0.5 * (RNE(L) + RNE(L+1))
            GEDDL = 10.**(-4.51) * (RNEINT/ATMEAN) * XLSTARS / XMSTAR
          ENDIF
          ANORM(L) = AGRAV(L) * ( 1. - GEDDL )
          CNORM = 'g&Teff&M'
        ELSE
          ANORM(L) = AGRAV(L)
          CNORM = 'g'
        ENDIF          
      ENDDO

      !Find sonic point parameters
      Lcand = 0
      DO L=1, ND
        VMACH(L) = SQRT( RGAS * T(L) / XMU(L) ) / 1.E5      !v_mach in km/s
        VSCRATCH(L) = VELO(L) - VMACH(L)        
        IF ((VSCRATCH(L) < 0.) .AND. (Lcand == 0)) THEN
          Lcand = L
        ENDIF
      ENDDO
      IF (Lcand > 1) THEN
        CALL SPLINPOX(Rsonic,0.,RADIUS,VSCRATCH,ND,.FALSE.,Lcand)
        CALL SPLINPOX(Vsonic,Rsonic,VMACH,RADIUS,ND)
      ENDIF
      
      
      IF (XAXISMODE == 'DEPTHINDEX' .OR. XAXISMODE == 'L') THEN 
C***     X-Axis = Depth Index
         XTEXT = CENTER//'DEPTH INDEX L'
         XMIN = 0.
         XMAX = FLOAT(ND)
         XTICK = 5.
         XABST = 10. 
C         DO L=1, ND-2
         DO L=1, ND-1
            X(L) = FLOAT(L) + 0.5
         ENDDO
         IF (Lcand > 1 .AND. Lcand < ND) THEN
           CALL SPLINPOX(Xsonic,Rsonic,X,RI,ND)         
         ENDIF
      ELSEIF (XAXISMODE == 'VELOCITY' .OR. XAXISMODE == 'V') THEN
C***     X-Axis: velocity / v_infty
         XTEXT = CENTER//'v(r) / v' // CHAR(92) // '8'
         XMIN = 0.
         XMAX = 1.
         XTICK = 0.1
         XABST = 0.5 
C         DO L=1, ND-2
         DO L=1, ND-1
           CALL SPLINPOX(VINT, RI(L), VELO, RADIUS, ND)
           X(L) = VINT/VELO(1)
         ENDDO
         IF (Lcand > 1 .AND. Lcand < ND) THEN
           Xsonic = Vsonic/VELO(1)
         ENDIF
      ELSEIF (XAXISMODE == 'RADIUS' .OR. XAXISMODE == 'R') THEN
         XTEXT = CENTER//'LOG (R/R\*-1)'
         XMIN = LOG10 ( 0.9 * (RI(ND-1) - 1.) )
         XMAX = LOG10 ( RADIUS(1) - 1. )
         XTICK = 0.2
         XABST = 1.0 
         DO L=1, ND-1
           X(L) = LOG10( RI(L) - 1. )
         ENDDO
         IF (Lcand > 1 .AND. Lcand < ND) THEN
           Xsonic = LOG10( Rsonic - 1. )
         ENDIF
      ELSE
         WRITE (hCPR,*) 'PLOTACC: Invalid XAXISMODE ************'
         WRITE (hCPR,*) '***** The following plot was aborted:'
         WRITE (hCPR,*) PLOTOPT(:IDX(PLOTOPT))
         RETURN
      ENDIF


      CALL JSYMSET ('G2','TRANSFER')

      HEADLINE = 'ACC:'//MODHEAD(13:)
      WRITE (HEADLINE(90:), '(A8,I7)') ' JOB No.', JOBNUM
C    9 FORMAT (19X,1HJ,I3)

      WRITE (hPLOT, '(A)') '*NEXTPLOT: ' // PLOTOPT
      WRITE (hPLOT, '(A)') 'PLOT: ' // HEADLINE
      WRITE (hPLOT, '(A)') '\FONT=HELVET'      
      WRITE (hPLOT, '(A)') '\DEFINECOLOR 5 1.0 0.58 0.0'
      WRITE (hPLOT, '(A)') '\DEFINECOLOR 7 0.0 0.50 0.0'
      WRITE (hPLOT, '(A)') '\DEFINECOLOR 9 0.6 0.60 0.6'
      WRITE (hPLOT, '(A)') '\PEN=1'     !needed to get rid of definecolor pen bug 

      WRITE (hPLOT, '(A)') '\COLOR=3'
      WRITE (hPLOT, '(A)') '\LINUN XMIN 0 XMAX 0 0 0'

      WRITE (hPLOT, '(A)') '\COLOR=9'
      WRITE (hPLOT, '(A)') '\BGRLUN COLOR=0'
      IF (RCON >= RADIUS(ND) .AND. RCON <= RADIUS(1)) THEN
        CALL SPLINPOX(XICON, RCON, X, RI, ND)
c        IF (XICON >= 0.005) THEN
        IF (XICON >= XMIN + 0.1 * XTICK) THEN
          WRITE (hPLOT, '(A,F6.3,A,F6.3,A)') 
     >      '\LINUN ',XICON,' YMIN ',XICON,' YMAX 0. 0. SYMBOL=10'
          WRITE (hPLOT, '(A,F6.3,A)') 
     >      '\LUNA ',XICON,' YMAX 0. -0.2 0.2 -90  R&Tcon'
        ENDIF
      ENDIF
      IF (Lcand > 1 .AND. Lcand < ND) THEN
        IF (ABS(Xsonic) >= 0.005) THEN
          WRITE (hPLOT, '(A,F12.5,A,F12.5,A)') 
     >      '\LINUN ',Xsonic,' YMIN ',Xsonic,
     >                    ' YMAX 0. 0. SYMBOL=20 SIZE=0.05'
          WRITE (hPLOT, '(A,F12.5,A)') 
     >      '\LUNA ',Xsonic,' YMAX 0. -0.2 0.2 -90  R&Tsonic'
        ENDIF
      ENDIF
      IF (Rcritical >= RI(ND) .AND. Rcritical <= RI(1)) THEN
        CALL SPLINPOX(Xcrit, Rcritical, X, RI, ND)
        IF (Xcrit >= 0.005) THEN
          WRITE (hPLOT, '(A,F6.3,A,F6.3,A)') 
     >      '\LINUN ',Xcrit,' YMIN ',Xcrit,' YMAX 0. 0. SYMBOL=9'
          WRITE (hPLOT, '(A,F6.3,A)') 
     >      '\LUNA ',Xcrit,' YMAX 0. -0.2 0.2 -90  R&Tcrit'
        ENDIF
      ENDIF
      WRITE (hPLOT, '(A)') '\BGRLUN OFF'
      WRITE (hPLOT, '(A)') '\COLOR=1'

      WRITE (hPLOT, '(A,F6.3)') 
     >   'KASDEF LUN XMIN YMIN 1. 1. .5 ' //
     >    '&EWORK Ratio Rad.+Gaspress. / Mech.+Grav. =&N', WORKRATIO  
      WRITE (hPLOT, '(A)') '\INBOX'

      XLENGTH = 2.5
      XOFF = XLENGTH + 14. + 0.5
      WRITE (hPLOT, '(A,F4.1,A)') CHAR(92) // 
     > 'LINREL XMIN YMAX ',XLENGTH,' 0 14. -1. SYMBOL=5'
      WRITE (hPLOT, '(A,F4.1,A)') '\LUN XMIN YMAX ', XOFF, ' M-1.'
     > // ' 0.3 (' // CNORM(:IDX(CNORM)) // '+a&Tmech&M)/' // CNORM

      WRITE (hPLOT, '(A)') CHAR(92) // 'COLOR=2'
      WRITE (hPLOT, '(A,F4.1,A)') CHAR(92) // 
     > 'LINREL XMIN YMAX ',XLENGTH,' 0 14. -1.5 SYMBOL=9'
      WRITE (hPLOT, '(A,F4.1,A)') CHAR(92) // 
     > 'LUN XMIN YMAX ', XOFF, ' M-1.5 0.3 a&Trad&M/' // CNORM

      WRITE (hPLOT, '(A)') CHAR(92) // 'COLOR=4'
      WRITE (hPLOT, '(A,F4.1,A)') CHAR(92) // 
     > 'LINREL XMIN YMAX ',XLENGTH,' 0 14. -2. SYMBOL=5'
      WRITE (hPLOT, '(A,F4.1,A)') CHAR(92) // 
     > 'LUN XMIN YMAX ', XOFF, ' M-2. 0.3 a&Tmech&M/' // CNORM

      WRITE (hPLOT, '(A)') CHAR(92) // 'COLOR=7'
      WRITE (hPLOT, '(A,F4.1,A)') CHAR(92) // 
     > 'LINREL XMIN YMAX ',XLENGTH,' 0 14. -2.5 SYMBOL=9'
      WRITE (hPLOT, '(A,F4.1,A)') CHAR(92) // 
     > 'LUN XMIN YMAX ', XOFF, ' M-2.5 0.3 a&Tcont&M/' // CNORM

      WRITE (hPLOT, '(A)') CHAR(92) // 'COLOR=2'
      WRITE (hPLOT, '(A,F4.1,A)') CHAR(92) // 
     > 'LINREL XMIN YMAX ',XLENGTH,' 0 14. -3 SYMBOL=5'
      WRITE (hPLOT, '(A,F4.1,A)') '\LUN ' //
     > 'XMIN YMAX ', XOFF, ' M-3.0 0.3 (a&Trad&M+a&Tpress&M)/' // CNORM

      WRITE (hPLOT, '(A)') CHAR(92) // 'COLOR=5'
      WRITE (hPLOT, '(A,F4.1,A)') CHAR(92) // 
     > 'LINREL XMIN YMAX ',XLENGTH,' 0 14. -3.5 SYMBOL=9 SIZE=0.05'
      WRITE (hPLOT, '(A,F4.1,A)') CHAR(92) // 
     > 'LUN XMIN YMAX ', XOFF, ' M-3.5 0.3 a&Tpress&M/' // CNORM

      WRITE (hPLOT, '(A)') '\COLOR=6'
      WRITE (hPLOT, '(A,F4.1,A)') 
     > '\LINREL XMIN YMAX ',XLENGTH,' 0 14. -4.0  SYMBOL=5'
      WRITE (hPLOT, '(A,F4.1,A)') '\LUN XMIN YMAX ', XOFF,
     > ' M-4.0 0.3 a&Tthom&M/' // CNORM
      
      WRITE (hPLOT, '(A)') '\COLOR=1'

C***  HTOT: Conversion into radiation temperatures
C***  Note: HTOT may not be calculated if STEAL is used for OUTPUT ONLY
      YMIN = 0.
      YMAX = 0.
cc      YMIN = -1.
cc      YMAX =  1.
      DO L=1, ND-1
         Y1(L) = (ANORM(L) + AMECH(L))/ANORM(L)
         Y2(L) = ARAD(L)/ANORM(L)
         Y3(L) = AMECH(L)/ANORM(L)
         Y4(L) = ACONT(L)/ANORM(L)
         Y5(L) = APRESS(L)/ANORM(L)
cc         IF ((Y1(L) .GT. 0.) .AND. (Y2(L) .GT. 0.)) THEN
cc            YMAX = ALOG10(MAX(10**YMAX, 1.05*Y1(L), 1.05*Y2(L)))
cc            YMIN = ALOG10(MIN(10**YMIN, 0.95*Y1(L), 0.95*Y2(L)))
cc         ENDIF
         IF (Y1(L) > .0) Y1(L) = LOG10(Y1(L))
         IF (Y2(L) > .0) Y2(L) = LOG10(Y2(L))
         IF (Y3(L) > .0) Y3(L) = LOG10(Y3(L))
         IF (Y4(L) > .0) Y4(L) = LOG10(Y4(L))
         IF (Y5(L) > .0) Y5(L) = LOG10(Y5(L))
      ENDDO

C***  Restrict blue curve (AMECH) to inbox area
      YMINVAL = MIN(MINVAL(Y2), MINVAL(Y4))
      YMINVAL = MIN(YMINVAL, MINVAL(Y5))
      DO L = 1, ND-1
        NDIN = L
        IF (Y3(L) < YMINVAL) EXIT
      ENDDO

      CALL PLOTANF (hPLOT,HEADLINE, '&E'//HEADLINE,
     $        XTEXT, CENTER//'log&T10&M(a/g)',
     >        0., XMIN, XMAX, XTICK, XABST, 0.,
     >        0., YMIN, YMAX, .1,  0.2, 0.,
     $        X, Y1, ND-1, 5)

      CALL PLOTCONS (hPLOT,X,Y2,ND-1,'COLOR=2 SYMBOL=9 SIZE=0.1') 
      CALL PLOTCONS (hPLOT,X,Y3,NDIN,'COLOR= 4') 
      CALL PLOTCONS (hPLOT,X,Y4,ND-1,'COLOR=7 SYMBOL=9 SIZE=0.1') 

      DO L=1, ND-1
         Y2(L) = (ARAD(L)+APRESS(L))/ANORM(L)
         IF (Y2(L) > .0) Y2(L) = LOG10(Y2(L))
      ENDDO
      CALL PLOTCONS (hPLOT,X,Y2,ND-1,'COLOR=2') 

      CALL PLOTCONS (hPLOT,X,Y5,ND-1,'COLOR=5 SYMBOL=9 SIZE=0.05') 

      DO L=1, ND-1
         Y2(L) = ATHOM(L)/ANORM(L)
         IF (Y2(L) > .0) Y2(L) = LOG10(Y2(L))
      ENDDO
      CALL PLOTCONS (hPLOT,X,Y2,ND-1,'COLOR=6') 

      RETURN
      END
