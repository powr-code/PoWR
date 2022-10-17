      SUBROUTINE SECONDMODEL_PREP (ZINTER, NPHI, P, NP, NPDIM, NPHIMAX, 
     >     PHIARR, PHIWEIGHT, PHI_VEC, SECONDMODEL_LINE, 
     >     JPFIRST, JPLAST, LPHISTA_ORIG, LPHIEND_ORIG)
C******************************************************************
C***  CALCULATES THE POINTS OF INTERSECTION FOR EITHER:
C***  - THE CONE geometry:
C***     THETA IS THE OPENING ANGLE OF THE CONE
C***     CONEI is its INCLINATION angle
C***  - THE SPHERE geometry:
C***    RSPHERE= radius of the sphere
C***    DSPHERE= radial distance of the sphere's center from the origen
C***    DELTASPHERE= the sphere-center's elevation angle 
C***    ALPHASPHERE= meridian angle to the sphere center (see Manual)
C***    
C***  The pre-existing set of angle points is enlarged to cover the
C***  cone, but not wasting points outside
C******************************************************************

      DIMENSION ZINTER(2, NPDIM, NPHIMAX), P(NP)
      DIMENSION PHIARR(NPHIMAX,NPDIM), PHIWEIGHT(NPHIMAX,NPDIM)
      DIMENSION NPHI(NPDIM)
      DIMENSION PHI_VEC(NPHIMAX)
      CHARACTER SECONDMODEL_LINE*(*), ACTPAR*100, ACTPAR2*100
      CHARACTER SHAPE*10

      PARAMETER (NPHIMAX_TEMP = 5000)
      DIMENSION PHI_VEC_TEMP(NPHIMAX_TEMP)
      LOGICAL BPHI_ORIG(NPHIMAX_TEMP), BPHI_ORIG_LAST

      DATA PI / 3.14159265358979 /

C***  Defaults
      DATA PHI_REFINE / 1. /
      DATA SHAPE / 'UNDEFINED' /
      DATA THETA_DEG / -999. /
      DATA CONEI_DEG / -999. /
      DATA DELTA_DEG / -999. /
      DATA ALPHA_DEG / -999. /
      DATA RSPHERE   / -999. /
      DATA DSPHERE   / -999. /

      RMAX = P(NP)

C***  Decode parameters from input line
      CALL SARGC (SECONDMODEL_LINE, NPAR)
      DO IPAR=1, NPAR
         CALL SARGV (SECONDMODEL_LINE, IPAR, ACTPAR)

         IF (ACTPAR .EQ. 'SHAPE') THEN
C*                        =====
            IF (NPAR .LT. IPAR+1) THEN
               GOTO 90
            ELSE
               CALL SARGV (SECONDMODEL_LINE, IPAR+1,SHAPE)
            ENDIF

         ELSEIF (ACTPAR .EQ. 'THETA') THEN
C*                            =====
            IF (NPAR .LT. IPAR+1) THEN
               GOTO 90
            ELSE
               CALL SARGV (SECONDMODEL_LINE, IPAR+1,ACTPAR2)
               READ (ACTPAR2, '(F10.0)', ERR=98) THETA_DEG
            ENDIF
         
         ELSEIF (ACTPAR .EQ. 'CONEI') THEN
C*                            =====
            IF (NPAR .LT. IPAR+1) THEN
               GOTO 90
            ELSE
               CALL SARGV (SECONDMODEL_LINE, IPAR+1,ACTPAR2)
               READ (ACTPAR2, '(F10.0)', ERR=98) CONEI_DEG
            ENDIF
         
         ELSEIF (ACTPAR .EQ. 'DELTA') THEN
C*                            =====
            IF (NPAR .LT. IPAR+1) THEN
               GOTO 90
            ELSE
               CALL SARGV (SECONDMODEL_LINE, IPAR+1,ACTPAR2)
               READ (ACTPAR2, '(F10.0)', ERR=98) DELTA_DEG
            ENDIF

         ELSEIF (ACTPAR .EQ. 'ALPHA') THEN
C*                            =====
            IF (NPAR .LT. IPAR+1) THEN
               GOTO 90
            ELSE
               CALL SARGV (SECONDMODEL_LINE, IPAR+1,ACTPAR2)
               READ (ACTPAR2, '(F10.0)', ERR=98) ALPHA_DEG
            ENDIF

         ELSEIF (ACTPAR .EQ. 'RSPHERE') THEN
C*                            =====
            IF (NPAR .LT. IPAR+1) THEN
               GOTO 90
            ELSE
               CALL SARGV (SECONDMODEL_LINE, IPAR+1,ACTPAR2)
               READ (ACTPAR2, '(F10.0)', ERR=98) RSPHERE
            ENDIF

         ELSEIF (ACTPAR .EQ. 'DSPHERE') THEN
C*                            =====
            IF (NPAR .LT. IPAR+1) THEN
               GOTO 90
            ELSE
               CALL SARGV (SECONDMODEL_LINE, IPAR+1,ACTPAR2)
               READ (ACTPAR2, '(F10.0)', ERR=98) DSPHERE
            ENDIF

         ELSEIF (ACTPAR .EQ. 'PHI_REFINE') THEN
C*                            =====
            IF (NPAR .LT. IPAR+1) THEN
               GOTO 90
            ELSE
               CALL SARGV (SECONDMODEL_LINE, IPAR+1,ACTPAR2)
               READ (ACTPAR2, '(F10.0)', ERR=98) PHI_REFINE
            ENDIF
         ENDIF
      ENDDO         

C***  Check for mandatory parameters
      IF (SHAPE .EQ. 'UNDEFINED') THEN
         GOTO 92
      ELSEIF  (SHAPE .EQ. 'CONE') THEN
         IF (THETA_DEG .EQ. -999.) GOTO 93
         THETA = THETA_DEG * PI / 180.
         IF (CONEI_DEG .EQ. -999.) GOTO 94
         CONEI = CONEI_DEG * PI / 180.

         IF (THETA_DEG .GE. 89.) THEN
            WRITE (0,'(A,F6.1,A)') '***** ERROR: ' //
     >        'Cone opening angle given: THETA=', THETA_DEG, ' degrees'
            WRITE (0,'(A)') '***** ERROR: ' //
     >        'Cone must have opening angle THETA .le. 89 degrees'
            GOTO 99
         ENDIF

         IF (CONEI_DEG - THETA_DEG .LE. .0) THEN
            WRITE (0,'(A, /, A, F6.1, A, F6.1, A)') '***** ERROR: ' //
     >       ' Invalid choice of geometry:', 'THETA=', THETA_DEG,
     >       'deg,  CONE-Inclination=', CONEI_DEG, 'deg'
            WRITE (0,'(A)') '***** ERROR: ' //
     >       'Cone MUST be seen from the side'
            GOTO 99
         ENDIF

      ELSEIF  (SHAPE .EQ. 'SPHERE') THEN
         IF (DELTA_DEG .EQ. -999.) GOTO 95
         DELTA = DELTA_DEG * PI / 180.
         IF (ALPHA_DEG .EQ. -999.) GOTO 96
         ALPHA = ALPHA_DEG * PI / 180.
         IF (RSPHERE .EQ. -999.) GOTO 97
         IF (DSPHERE .EQ. -999.) GOTO 971
         IF (RSPHERE .GE. DSPHERE+RMAX) GOTO 972

      ELSE
         GOTO 91
      ENDIF

C***  Preparation for the CONE case
      IF (SHAPE .EQ. 'CONE') THEN
         SINI = SIN(CONEI)
         COSI = COS(CONEI)
         COTMINUS = 1. / TAN(CONEI-THETA) 
         COTPLUS  = 1. / TAN(CONEI+THETA) 

C***  Preparation for the SPHERE case
      ELSEIF (SHAPE .EQ. 'SPHERE') THEN
         SIND = SIN (DELTA)
         COSD = COS (DELTA)
         SINA = SIN (ALPHA)
         COSA = COS (ALPHA)
         RSPHERE2 = RSPHERE**2
      ENDIF

      RMAX2 = RMAX * RMAX
      MAXNPHI = 0
      NPHISUM = 0

C***  Loop over impact parameters
      DO JP=1, NP-1

C***  The angle-points will be established  

C***  Only 1 angle point for JP=1 (center)
         IF (JP .EQ. 1) THEN
            PHI_VEC_TEMP(1) = .0
            NPHI_TEMP = 1
            GOTO 10
         ENDIF         

C***  Preparation of angle points
         IF (NPHI(JP) .EQ. 1) THEN
            PHI_VEC(1) = .0
            PHI_VEC(2) = 2.*PI
            NPHI_JP = 2
         ELSE

C**     Copy pre-existing phi points (from wind rotation)
            NPHI_JP = NPHI(JP)
            DO LPHI=1, NPHI_JP
               PHI_VEC(LPHI) = PHIARR (LPHI,JP)
            ENDDO

C**         mirror the list to the southern hemisphere 
            IF ((2 * NPHI_JP - 1) .GT. NPHIMAX) THEN
               WRITE (0,200) NPHIMAX
               WRITE (*,200) NPHIMAX
  200    FORMAT ('*** ERROR: More PHI points needed than dimensioned',
     >          /, 'NPHIMAX= ', I6)
               STOP '*** FATAL ERROR in subr. SECONDMODEL_PREP'
            ENDIF
            DO LPHI=1, NPHI_JP-1         
               PHI_VEC(2*NPHI_JP-LPHI) = 2.*PI -  PHI_VEC(LPHI)
            ENDDO
            NPHI_JP = 2 * NPHI_JP - 1
         ENDIF

C***     Now we add a fine grid of angle points, spaced by the requested
C***     angular resolution DELTA_PHI

C***     Resolution in cone
         IF (SHAPE .EQ. 'CONE') THEN
            NPHI_PER_CONE = NINT(10 * PHI_REFINE)
            DELTA_PHI = 2. * THETA / NPHI_PER_CONE

C***     Resolution in sphere
         ELSEIF (SHAPE .EQ. 'SPHERE') THEN
            NPHI_ACROSS_SPHERE = NINT(10 * PHI_REFINE)
            DELTA_PHI = MIN(2.*PI, RSPHERE / P(JP)) / NPHI_ACROSS_SPHERE

         ELSE
            STOP '*** ERROR: Invalid SHAPE in SECONDMODEL_PREP'
         ENDIF

         PHI_VEC_TEMP(1) = PHI_VEC(1)
         LPHI = 1
         LPHI_TEMP = 1
         DO  
            PHI_NEXT = PHI_VEC_TEMP(LPHI_TEMP) + DELTA_PHI
            LPHI_TEMP = LPHI_TEMP + 1

C***        Error stop if vector PHI_VEC_TEMP too short
            IF (LPHI_TEMP .GT. NPHIMAX_TEMP) THEN
              WRITE (0,201) NPHIMAX_TEMP
              WRITE (*,201) NPHIMAX_TEMP
  201         FORMAT ('*** ERROR: More PHI points needed than' // 
     >                ' dimensioned', /, 'NPHIMAX_TEMP= ', I6)
              STOP '*** FATAL INTERNAL ERROR in subr. SECONDMODEL_PREP'
            ENDIF

            IF (PHI_VEC(LPHI+1) .LT. PHI_NEXT .OR.
     >          PHI_NEXT .GE. 2*PI) THEN
               LPHI = LPHI + 1
               PHI_VEC_TEMP(LPHI_TEMP) = PHI_VEC(LPHI) 
               BPHI_ORIG(LPHI_TEMP) = .TRUE.
            ELSE
               PHI_VEC_TEMP(LPHI_TEMP) = PHI_NEXT
               BPHI_ORIG(LPHI_TEMP) = .FALSE.
            ENDIF
            IF (LPHI .GE. NPHI_JP) EXIT
         ENDDO

         NPHI_TEMP = LPHI_TEMP

   10    CONTINUE
C***     Angle points are now in PHI_VEC_TEMP

C********Loop over all angle points at current JP to find intersections **
         DO LPHI_TEMP=1, NPHI_TEMP
            PHI = PHI_VEC_TEMP (LPHI_TEMP)
            X0  = P(JP) * COS(PHI)
            X02 = X0 * X0
            Y0  = P(JP) * SIN(PHI)
            Y02 = Y0 * Y0 

C***        Intersection points with cone
            IF (SHAPE .EQ. 'CONE') THEN
              ZM = 0.5 * Y0 * (COTPLUS + COTMINUS)
              AAXIS = 0.5 * Y0 * (COTPLUS - COTMINUS)
              AAXIS2 = AAXIS * AAXIS
              BAXIS2 = (Y02 + ZM**2) * (TAN(THETA))**2
              IF (X02 .GE. BAXIS2) THEN
                 Z1 = .0
                 Z2 = .0
              ELSE
                 TERM2 = AAXIS2 * (1- X02/BAXIS2)
                 TERM2 = SQRT(TERM2)
                 Z1 = ZM + TERM2
                 Z2 = ZM - TERM2
              ENDIF

C***        Intersection points with sphere
            ELSEIF (SHAPE .EQ. 'SPHERE') THEN
              XM = DSPHERE * COSD * SINA
              YM = DSPHERE * SIND
              ZM = DSPHERE * COSD * COSA
              DX2 = (X0-XM)**2
              DY2 = (Y0-YM)**2
              TERM2 = RSPHERE2 - DX2 - DY2
              IF (TERM2 .LE. .0) THEN
                 Z1 = .0
                 Z2 = .0
              ELSE
                 TERM2 = SQRT(TERM2)
                 Z1 = ZM + TERM2
                 Z2 = ZM - TERM2
              ENDIF
      
            ELSE
              STOP '*** Internal ERROR in SECONDMODEL_PREP: unknown SHAPE'
            ENDIF

C**         Clipping the intersection line at the RMAX sphere
            ZMAX2 = RMAX2 - X02 - Y02
            ZMAX2 = MAX (.0, ZMAX2)
            ZMAX = SQRT (ZMAX2)
            IF (Z2 .GE. ZMAX .OR. Z1 .LE. -ZMAX) THEN
C**            Intersection interval is entirely outside RMAX sphere
               Z1 = .0
               Z2 = .0
            ELSE
               Z1 = MIN (Z1,  ZMAX)
               Z2 = MAX (Z2, -ZMAX)
            ENDIF

C**         Core rays:
            IF (P(JP) .LT. 1.) THEN
C**            Intersection lines behind the stellar disc are obscured
               ZMIN = SQRT (1. - P(JP)*P(JP))
               IF (Z1 .LE. ZMIN) THEN
                  Z1 = .0
                  Z2 = .0
               ELSE
                  Z2 = ZMIN
               ENDIF
            ENDIF

C***        First angle point is always kept
            IF (LPHI_TEMP .LE. 1) THEN
               ZINTER(1, JP, 1) = Z1
               ZINTER(2, JP, 1) = Z2
               PHIARR(1,JP) = PHI_VEC_TEMP(1)
               LPHI = 1
               BPHI_ORIG_LAST = .TRUE.

            ELSE
C***        Only if last OR current phi intersects -> increase counter; 
C***        else overwrite last point
C***        "original" angle points are also kept and not overwritten
C***        additionally, first non-intersection points are not overwritten
              DZ = Z2 - Z1
              DZLAST = ZINTER(2, JP, LPHI) - ZINTER(1, JP, LPHI)
              IF (DZ .NE. .0 .OR. DZLAST .NE. .0 .OR. BPHI_ORIG_LAST) 
     >           LPHI = LPHI + 1
              ZINTER(1, JP, LPHI) = Z1
              ZINTER(2, JP, LPHI) = Z2
              PHIARR(LPHI,JP) = PHI_VEC_TEMP(LPHI_TEMP)
C**           If last point was inside cone, do not overwrite it
              BPHI_ORIG_LAST= BPHI_ORIG(LPHI_TEMP) .OR. (DZLAST .NE. .0)
            ENDIF

         ENDDO ! phi-loop, index LPHI_TEMP -----------------------------

         NPHI(JP) = LPHI
         IF (NPHI(JP) .GT. MAXNPHI) THEN
           MAXNPHI = NPHI(JP)
           MAXJP   = JP
         ENDIF
         NPHISUM = NPHISUM + NPHI(JP)

      ENDDO ! p-loop --------------------------------------------------

C***  Output of data cube for visualization with gnuplot
      OPEN (20, FILE='secondmodel.dat',
     >      STATUS='UNKNOWN')
      DO J=1, NP
         DO LPHI = 1, NPHI(J)
            X0 = P(J)/P(NP) * COS(PHIARR(LPHI,J))
            Y0 = P(J)/P(NP) * SIN(PHIARR(LPHI,J))
            IF (ZINTER(1,J,LPHI) .NE. .0 .OR.
     >          ZINTER(1,J,LPHI) .NE. .0) THEN
                  WRITE (20, *) X0, Y0, ZINTER(1,J,LPHI)/P(NP)
                  WRITE (20, *) X0, Y0, ZINTER(2,J,LPHI)/P(NP)
                  WRITE (20, *) '   '
                  WRITE (20, *) '   '
            ENDIF
         ENDDO
      ENDDO
      CLOSE (20)

C***  plot of the secondmodel's geometrical (phi,p) grid
      OPEN (66, FILE='secondmodel.plot', STATUS='UNKNOWN')

      CALL PLOT_SECONDMODEL_GRID (P, NP, NPDIM, NPHI, PHIARR,
     >      NPHIMAX, JPFIRST, JPLAST, LPHISTA_ORIG, LPHIEND_ORIG, 
     >      ZINTER)

      CLOSE (66)

C***  Statistical output

      WRITE (*,'(A)') 'SECONDMODEL invoked; specifications:'
      WRITE (*,'(A,$)') 'SHAPE=', SHAPE
      IF (SHAPE .EQ. 'CONE') WRITE (*, '(A,F6.1,A,F6.1)') 
     >         '  THETA=', THETA_DEG, '  CONEI=', CONEI_DEG
      IF (SHAPE .EQ. 'SPHERE') WRITE (*, '(4(A,F6.1))') 
     >         '  RSPHERE=', RSPHERE, '  DSPHERE=', DSPHERE,
     >         '  DELTASPHERE=', DELTASPHERE_DEG, 
     >         '  ALPHASPHERE=', ALPHASPHERE_DEG

      WRITE (0,'(/,A,A)') 
     >       'SECONDMODEL invoked; parameters: SHAPE=', SHAPE
      WRITE (0,'(A,I4,A,I4,A,F10.3)')
     > 'Maximum number of phi angles:', MAXNPHI,
     >     ' at P(', MAXJP, ') =', P(MAXJP)

      NPHI_MEAN = NINT(FLOAT(NPHISUM)/FLOAT(NP))
      WRITE(0,'(A,I4,/)')
     > 'Average number of phi angles per impact parameter:', NPHI_MEAN

C***  Calculate angle integration weights
C***    Note: Normalization of PHIWEIGHT will be done in program FORMAL
cccc    note: this part is identical in ROTATION_PREP 
      DO JP=1, NP-1
        IF (NPHI(JP) .EQ. 1) THEN
            PHIWEIGHT(1,JP) = 1.
        ELSE 
C***       Calculate the angle integration weights by trapezoidal rule
C***       but omitting the factors 0.5 at all weights
           PHIWEIGHT(1,      JP) = PHIARR(2,JP) - PHIARR(1,JP)
           DO LPHI=2, NPHI(JP)-1
             PHIWEIGHT(LPHI,JP) = PHIARR(LPHI+1,JP) - PHIARR(LPHI-1,JP)
           ENDDO
           PHIWEIGHT(NPHI(JP),JP) = PHIARR(NPHI(JP),JP)
     >                           - PHIARR(NPHI(JP)-1,JP)
        ENDIF

      ENDDO


      RETURN

C*********************************************************************
C***  ERROR branches  ************************************************
C*********************************************************************

   90 WRITE (0,*) '*** ERROR: Option ', ACTPAR(:IDX(ACTPAR)), 
     >            ' needs a value (keyword)'
      GOTO 99

   91 WRITE (0,*) '*** ERROR: Invalid keyword after option SHAPE: ', 
     >            SHAPE(:IDX(SHAPE))
      GOTO 99

   92 WRITE (0,*) '*** ERROR: mandatory parameter SHAPE is missing'
      GOTO 99

   93 WRITE (0,*) '*** ERROR: mandatory parameter THETA is missing'
      GOTO 99

   94 WRITE (0,*) '*** ERROR: mandatory parameter CONEI is missing'
      GOTO 99

   95 WRITE (0,*) '*** ERROR: mandatory parameter DELTA is missing'
      GOTO 99

   96 WRITE (0,*) '*** ERROR: mandatory parameter ALPHA is missing'
      GOTO 99

   97 WRITE (0,*) '*** ERROR: mandatory parameter RSPHERE is missing'
      GOTO 99

  971 WRITE (0,*) '*** ERROR: mandatory parameter DSPHERE is missing'
      GOTO 99

  972 WRITE (0,*) '*** ERROR: SPHERE covers the whole atmosphere'
      GOTO 99

   98 WRITE (0,*) '*** ERROR: Parameter cannot be decoded as number:'
      WRITE (0,*) '*** ERROR: ', ACTPAR(:IDX(ACTPAR)), '=', 
     >                           ACTPAR2(:IDX(ACTPAR2))
      GOTO 99


   99 STOP '*** FATAL ERROR in subroutine SECONDMODEL_PREP'

      END
