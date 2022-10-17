      SUBROUTINE PLOTGAMMA(PLOTOPT, TAUROSS, 
     >                     AGRAV, AMECH, ARAD, APRESS, ACONT, 
     >                     ATHOM, WORKRATIO, VELO, RADIUS, ND, ATMEAN,
     >                     ENTOT, RNE, RCON, T, TEFF, RSTAR, XMU, 
     >                     XMSTAR, Rcritical, bFULLHYDROSTAT, 
     >                     QIONMEAN, GAMMARADMEAN,
     >                     MODHEAD, JOBNUM, hPLOT)
C******************************************************************************
C***  Plot EDDINGTON GAMMA versus DEPTH INDEX, RADIUS, VELOCITY or TAUROSS
C***  for full radiation force, continuum contribution, and Thomson part
C***  (created by ansander in 2014 in analogy to PLOTACC)
C******************************************************************************
 
      IMPLICIT NONE
      INCLUDE 'interfacebib.inc'

      INTEGER, PARAMETER :: NDMAX = 200

      INTEGER, INTENT(IN) :: ND, JOBNUM, hPLOT

      REAL, DIMENSION(ND-1), INTENT(IN) :: AGRAV, AMECH, ARAD, APRESS,
     >                                     ACONT, ATHOM 
      REAL, DIMENSION(ND), INTENT(IN) :: VELO, RADIUS, 
     >                                   ENTOT, RNE, XMU, T, TAUROSS
      REAL, INTENT(IN) :: WORKRATIO, RCON, ATMEAN, XMSTAR, RSTAR, TEFF, 
     >                    Rcritical
     
      LOGICAL, INTENT(IN) :: bFULLHYDROSTAT

      CHARACTER(110) :: MODHEAD, HEADLINE
      CHARACTER(8) :: CENTER, CNORM
      CHARACTER(40) :: XTEXT
      CHARACTER(20) :: CUROPT, XAXISMODE
      CHARACTER PLOTOPT*(*)

      REAL, DIMENSION(NDMAX) :: X, Y1, Y2, Y3, Y4, Y5
      REAL, DIMENSION(ND-1) :: RI, GAMMARAD, GAMMARADCUT, GAMMAEDD,
     >                         GAMMACONT, GAMMATHOM, GAMMAHYDROSTAT
      REAL, DIMENSION(ND) :: VMACH, VSCRATCH, RHO

      INTEGER :: I, L, NPAR, NDIN, Lcand, LCON
      REAL :: XMIN, XMAX, YMIN, YMAX, YMINVAL, GAMMARADMEAN, XOFF,
     >        GEDDL, XLSTAR, XLSTARS, RNEINT, Xcrit, QIONMEAN, XINT,
     >        XABST, XTICK, XLENGTH, XICON, Rsonic, Vsonic, Xsonic
      LOGICAL :: bNormalizeToGEFF

      INTEGER, EXTERNAL :: IDX

      !Physical constants
      REAL, PARAMETER :: PI4 = 12.5663706144    !PI4 = 4*PI
      REAL, PARAMETER :: AMU = 1.66E-24         !Atomic mass unit (gramm)     
      REAL, PARAMETER :: STEBOL = 5.6705E-5     !STEFAN-BOLTZMANN CONSTANT (CGS-UNITS)
      REAL, PARAMETER :: XLSUN = 3.85E33        !Solar Luminosity (CGS-Units)
      REAL, PARAMETER :: XMSUN = 1.989E33       !XMSUN = Solar Mass (g)
      REAL, PARAMETER :: RGAS = 8.3145E7        !Gas Constant (CGS)

      !File and channel handles (=KANAL)
      INTEGER, PARAMETER :: hOUT = 6        !write to wruniqX.out (stdout)
      INTEGER, PARAMETER :: hCPR = 0        !write to wruniqX.cpr (stderr)

      IF (NDMAX < ND) THEN
         WRITE (hCPR,'(A)') 'PLOTACC: NON-FATAL ERROR ******'
         WRITE (hCPR,'(A)') 'PLOTACC: DIMENSION INSUFFICIENT'
         WRITE (hCPR,'(A)') 'PLOTACC: PLOT SKIPPED'
         RETURN
      ENDIF

      CENTER = '\CENTER\'

C***  Decode possible options for x-axis mode
      XAXISMODE = 'DEPTHINDEX'
      bNormalizeToGEFF = .FALSE.
      CALL SARGC (PLOTOPT, NPAR)
      IF (NPAR > 2) THEN
        parloop: DO I=3, NPAR
          CALL SARGV (PLOTOPT, I, CUROPT)
          IF (CUROPT == 'XAXISMODE' .OR. CUROPT == 'X') THEN
            IF (NPAR >= (I+1)) THEN
              CALL SARGV (PLOTOPT, I+1, CUROPT)
            ELSE
              CYCLE parloop
            ENDIF
          ENDIF
          SELECTCASE (CUROPT)
            CASE ('VELOCITY', 'VELO', 'V')
              XAXISMODE = 'VELOCITY'
            CASE ('TAUROSS', 'TAU')
              XAXISMODE = 'TAUROSS'
            CASE ('RADIUS', 'R')
              XAXISMODE = 'RADIUS'
            CASE ('DEPTHINDEX', 'INDEX', 'ND', 'L')
              XAXISMODE = 'DEPTHINDEX'
            CASE ('N', 'NTOT', 'ENTOT', 'PARTDENS')
              XAXISMODE = 'ENTOT'
            CASE ('RHO', 'DENS', 'DENSITY')
              XAXISMODE = 'RHO'
          ENDSELECT
        ENDDO parloop
      ELSE 
C***    Fallback for unique labeling inside steal.plot / wruniq.plot
        PLOTOPT = PLOTOPT(:IDX(PLOTOPT)) // ' L'
      ENDIF
C      IF (NPAR .GT. 2) CALL SARGV (PLOTOPT, 3, XAXISMODE)

      XLSTAR = PI4 * STEBOL * RSTAR*RSTAR * TEFF*TEFF*TEFF*TEFF
      XLSTARS = XLSTAR / XLSUN
      
      DO L=1, ND-1
        RI(L) = 0.5 * ( RADIUS(L) + RADIUS(L+1) )
        GAMMARAD(L) = ARAD(L)/AGRAV(L)
        GAMMARADCUT(L) = MIN(GAMMARAD(L), 0.9)              
        GAMMACONT(L) = ACONT(L)/AGRAV(L)
        GAMMATHOM(L) = ATHOM(L)/AGRAV(L)
        !next should be the same as GAMMATHOM, only calculated to check this
        RNEINT = 0.5 * ( RNE(L) + RNE(L+1) )
        GAMMAEDD(L) = 10.**(-4.51) * (RNEINT/ATMEAN) * XLSTARS / XMSTAR 
        GAMMAHYDROSTAT(L) = (ARAD(L) + APRESS(L)) / AGRAV(L)
      ENDDO
      
      !Find sonic point parameters
      Lcand = 0
      DO L=1, ND
        VMACH(L) = SQRT( RGAS * T(L) / XMU(L) ) / 1.E5      !v_mach in km/s
        VSCRATCH(L) = VELO(L) - VMACH(L)        
        IF ((VSCRATCH(L) < 0) .AND. (Lcand == 0)) THEN
          Lcand = L
        ENDIF
      ENDDO
      IF (Lcand > 1) THEN
        CALL SPLINPOX(Rsonic,0.,RADIUS,VSCRATCH,ND,.FALSE.,Lcand)
        CALL SPLINPOX(Vsonic,Rsonic,VMACH,RADIUS,ND)
      ENDIF
      
      
      IF (XAXISMODE == 'DEPTHINDEX') THEN 
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
           CALL SPLINPOX(Xsonic,Rsonic,X,RI,ND-1)         
         ENDIF
      ELSEIF (XAXISMODE == 'VELOCITY') THEN
C***     X-Axis: velocity / v_infty
         XTEXT = CENTER//'v(r) / v\8'
         XMIN = 0.
         XMAX = 1.
         XTICK = 0.1
         XABST = 0.5 
C         DO L=1, ND-2
         DO L=1, ND-1
            X(L) = 0.5 * (VELO(L) + VELO(L+1)) / VELO(1)
         ENDDO
         IF (Lcand > 1 .AND. Lcand < ND) THEN
           Xsonic = Vsonic/VELO(1)
         ENDIF
      ELSEIF (XAXISMODE == 'TAUROSS') THEN
         XTEXT = CENTER//'LOG Rosseland Optical Depth #t#'
         XMIN = LOG10(0.4 * TAUROSS(2))
         XMAX = LOG10( TAUROSS(ND) )
         XTICK = 0.1
         XABST = 1.0 
         DO L=1, ND-1
           CALL SPLINPOX(XINT, RI(L), TAUROSS, RADIUS, ND)
           X(L) = LOG10(XINT)
         ENDDO
         IF (Lcand > 1 .AND. Lcand < ND) THEN
           CALL SPLINPOX(Xsonic,Rsonic,X,RI,ND-1)         
         ENDIF
      ELSEIF (XAXISMODE == 'RADIUS') THEN
         XTEXT = CENTER//'LOG (R/R\*-1)'
         XMIN = LOG10 ( 0.9 * (RI(ND-1) - 1.) )
         XMAX = LOG10 ( RADIUS(1) - 1. )
         XTICK = 0.1
         XABST = 0.5 
         DO L=1, ND-1
           X(L) = LOG10( RI(L) - 1. )
         ENDDO
         IF (Lcand > 1 .AND. Lcand < ND) THEN
           Xsonic = LOG10( Rsonic - 1. )
         ENDIF
      ELSEIF (XAXISMODE == 'ENTOT') THEN
         XTEXT = CENTER // 'log(&Rn&N&Ttot&M/cm&H-3&M)'
         XMIN = LOG10(ENTOT(1))
         XMAX = LOG10(ENTOT(ND))
         XTICK = 1.
         XABST = 3.
         DO L=1, ND-1
           CALL SPLINPOX(XINT, RI(L), ENTOT, RADIUS, ND)
           X(L) = LOG10(XINT)
         ENDDO      
         IF (Lcand > 1 .AND. Lcand < ND) THEN
           CALL SPLINPOX(Xsonic,Rsonic,X,RI,ND-1)         
         ENDIF      
      ELSEIF (XAXISMODE == 'RHO') THEN
         XTEXT = CENTER // 'log(&R#r#&N/(g cm&H-3&M))'
         DO L=1, ND
           RHO(L) = ENTOT(L) * AMU * ATMEAN
         ENDDO         
         XMIN = LOG10(RHO(1))
         XMAX = LOG10(RHO(ND))
         XTICK = 1.
         XABST = 2.
         DO L=1, ND-1
           CALL SPLINPOX(XINT, RI(L), RHO, RADIUS, ND)
           X(L) = LOG10(XINT)
         ENDDO      
         IF (Lcand > 1 .AND. Lcand < ND) THEN
           CALL SPLINPOX(Xsonic,Rsonic,X,RI,ND-1)         
         ENDIF      
      ELSE
         WRITE (hCPR,*) 'PLOTACC: Invalid XAXISMODE ************'
         WRITE (hCPR,*) '***** The following plot was aborted:'
         WRITE (hCPR,*) PLOTOPT(:IDX(PLOTOPT))
         RETURN
      ENDIF


      CALL JSYMSET ('G2','TRANSFER')

      HEADLINE = 'GEDD:'//MODHEAD(13:)
      WRITE (HEADLINE(90:), '(A8,I7)') ' JOB No.', JOBNUM
C    9 FORMAT (19X,1HJ,I3)

      WRITE (hPLOT, '(A)') '*NEXTPLOT: ' // PLOTOPT
      WRITE (hPLOT, '(A)') 'PLOT: ' // HEADLINE
      WRITE (hPLOT, '(A)') '\FONT=HELVET'      
      WRITE (hPLOT, '(A)') '\DEFINECOLOR 3 0.8 0.90 0.8'
      WRITE (hPLOT, '(A)') '\DEFINECOLOR 5 1.0 0.58 0.0'
      WRITE (hPLOT, '(A)') '\DEFINECOLOR 7 0.0 0.50 0.0'
      WRITE (hPLOT, '(A)') '\DEFINECOLOR 9 0.6 0.60 0.6'
      WRITE (hPLOT, '(A)') '\PEN=1'     !needed to get rid of definecolor pen bug 

      WRITE (hPLOT, '(A)') '\COLOR=1'
      WRITE (hPLOT, '(A)') '\LINUN XMIN 1. XMAX 1. 0 0'

      WRITE (hPLOT, '(A)') '\COLOR=9'
      WRITE (hPLOT, '(A)') '\BGRLUN COLOR=0'
      IF (RCON >= RI(ND-1) .AND. RCON <= RI(1)) THEN
        CALL SPLINPOX(XICON, RCON, X, RI, ND-1)
        IF (XICON >= XMIN + 0.1 * XTICK) THEN
          WRITE (hPLOT, '(A,F8.3,A,F8.3,A)') 
     >      '\LINUN ',XICON,' YMIN ',XICON,' YMAX 0. 0. SYMBOL=10'
          WRITE (hPLOT, '(A,F8.3,A)') 
     >      '\LUNA ',XICON,' YMAX 0. -0.2 0.2 -90  R&Tcon'
        ENDIF
      ENDIF
      IF (Lcand > 0) THEN
        IF (Xsonic >= XMIN + 0.005) THEN
          WRITE (hPLOT, '(A,F8.3,A,F8.3,A)') 
     >      '\LINUN ',Xsonic,' YMIN ',Xsonic,
     >                    ' YMAX 0. 0. SYMBOL=20 SIZE=0.05'
          WRITE (hPLOT, '(A,F8.3,A)') 
     >      '\LUNA ',Xsonic,' YMAX 0. -0.2 0.2 -90  R&Tsonic'
        ENDIF
      ENDIF
      IF (Rcritical > RI(ND-1) .AND. Rcritical <= RI(1)) THEN
        CALL SPLINPOX(Xcrit, Rcritical, X, RI, ND-1)
        IF (Xcrit >= XMIN + 0.005) THEN
          WRITE (hPLOT, '(A,F8.3,A,F8.3,A)') 
     >      '\LINUN ',Xcrit,' YMIN ',Xcrit,' YMAX 0. 0. SYMBOL=9'
          WRITE (hPLOT, '(A,F8.3,A)') 
     >      '\LUNA ',Xcrit,' YMAX 0. -0.2 0.2 -90  R&Tcrit'
        ENDIF
      ENDIF
      WRITE (hPLOT, '(A)') '\COLOR=2'
      WRITE (hPLOT, '(A,F8.3,A,F8.3,A)') 
     >   '\LINUN XMIN ',GAMMARADMEAN,' XMAX ',GAMMARADMEAN,
     >   ' 0. 0. SYMBOL=10 SIZE=0.1'
      WRITE (hPLOT, '(A)') '\BGRLUN OFF'
      WRITE (hPLOT, '(A)') '\COLOR=1'

      
      WRITE (hPLOT, '(A,F6.3)') 
     >   '\LUN XMIN YMAX 1. -1. .3 Q = ', WORKRATIO  
      WRITE (hPLOT, '(A)') '\INBOX'

      XLENGTH = -2.25      

      WRITE (hPLOT, '(A)') '\COLOR=3'
      WRITE (hPLOT, '(A,F4.1,A)') 
     > '\LINREL XMAX YMAX ',XLENGTH,' 0. -4. -1 SYMBOL=9 SIZE=0.1'
      WRITE (hPLOT, '(A)') '\LUN XMAX YMAX ' //
     >     ' L-3.8 M-1.0 0.3 (a&Trad&M+a&Tpress&M)/g'

      WRITE (hPLOT, '(A)') '\COLOR=2'
      XOFF = -4. + 2 * XLENGTH / 3.
      WRITE (hPLOT, '(2(A,F4.1),A)') '\LINREL XMAX YMAX ',
     >   XLENGTH/3.,' 0. ',XOFF,' -1.5 SYMBOL=5'
      XOFF = -4. + 1 * XLENGTH / 3.
      WRITE (hPLOT, '(2(A,F4.1),A)') '\LINREL XMAX YMAX ',
     >   XLENGTH/3.,' 0. ',XOFF,' -1.5 SYMBOL=9 SIZE=0.1'
      XOFF = -4. 
      WRITE (hPLOT, '(2(A,F4.1),A)') '\LINREL XMAX YMAX ',
     >   XLENGTH/3.,' 0. ',XOFF,' -1.5 SYMBOL=10 SIZE=0.1'
      WRITE (hPLOT, '(A)') '\LUN XMAX YMAX ' //
     > ' L-3.8 M-1.5 0.3 #G#&Trad&M' 
      WRITE (hPLOT, '(A)') '\LUN XMAX YMAX ' //
     > ' L-2.8 M-1.55 0.2 &E(used/calc/mean)' 

      WRITE (hPLOT, '(A)') '\COLOR=7'
      WRITE (hPLOT, '(A,F4.1,A)') 
     > '\LINREL XMAX YMAX ',XLENGTH,' 0. -4. -2. SYMBOL=9 SIZE=0.1'
      WRITE (hPLOT, '(A)') '\LUN XMAX YMAX ' //
     > ' L-3.8 M-2. 0.3 #G#&Tcont&M'

      WRITE (hPLOT, '(A)') '\COLOR=5'
      WRITE (hPLOT, '(A,F4.1,A)')
     > '\LINREL XMAX YMAX ',XLENGTH,' 0. -4. -2.5 SYMBOL=9 SIZE=0.05'
      WRITE (hPLOT, '(A)') '\LUN XMAX YMAX ' //
     > ' L-3.8 M-2.5 0.3 #G#&Tthom&M (COLI)'

      WRITE (hPLOT, '(A)') '\COLOR=6'
      WRITE (hPLOT, '(A,F4.1,A)') 
     > '\LINREL XMAX YMAX ',XLENGTH,' 0. -4. -3.  SYMBOL=5'
      WRITE (hPLOT, '(A)') '\LUN XMAX YMAX ' //
     > ' L-3.8 M-3.0 0.3 #G#&Te&M (STEAL)'
      
      WRITE (hPLOT, '(A)') '\COLOR=1'

C***  HTOT: Conversion into radiation temperatures
C***  Note: HTOT may not be calculated if STEAL is used for OUTPUT ONLY
      YMIN = -0.05
      YMAX =  1.5
cc      YMIN = -1.
cc      YMAX =  1.

      CALL PLOTANFS (hPLOT,HEADLINE, '&E'//HEADLINE,
     $        XTEXT, CENTER//'Eddington #G#',
     >        0., XMIN, XMAX, XTICK, XABST, 0.,
     >        0., YMIN, YMAX, .1,  0.2, 0.,
     $        X, GAMMAHYDROSTAT, ND-1, 'COLOR=3 SYMBOL=9 SIZE=0.1')

      CALL PLOTCONS (hPLOT,X,GAMMACONT,ND-1,'COLOR=7 SYMBOL=9 SIZE=0.1') 
      CALL PLOTCONS (hPLOT,X,GAMMAEDD,ND-1,'COLOR=6') 
      CALL PLOTCONS (hPLOT,X,GAMMATHOM,ND-1,'COLOR=5 SYMBOL=9 SIZE=0.05') 
      CALL PLOTCONS (hPLOT,X,GAMMARAD,ND-1,'COLOR=2 SYMBOL=9 SIZE=0.1') 
      CALL PLOTCONS (hPLOT,X,GAMMARADCUT,ND-1,'COLOR=2') 


      RETURN
      END
