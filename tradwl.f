      SUBROUTINE TRADWL (KANAL,DELW,ADELW,VELO,ENTOT,ND,LINE,MODHEAD,
     $                   JOBNUM,LEVEL,XLAM, XPLOT, TAUROSS, RADIUS,
     $                   TRANSDWLLINE)
C***********************************************************************
C***  DIRECT TRANSFER OF PLOT DATA: DEPTH OF LINE FORMATION ("DWL-PLOT")
C***  Formal Cards Syntax: TRANS DWL <VELO>,<VELOLOG>,<TAUROSS>,<RSTAR>,
C***  <RSTARLOG>. If empty, default = <DENSITY> 
C***********************************************************************
 
      DIMENSION DELW(ND),ADELW(ND),VELO(ND),ENTOT(ND), XPLOT(ND), TAUROSS(ND), RADIUS(ND)
      CHARACTER*60 HEADER, MODHEAD, XTEXT, YTEXT, XAXIS, TRANSDWLLINE*(*), ACTPAR
      INTEGER NPAR
      CHARACTER*10 LEVEL

C**********************************************************************
C***  DEFINE HEADER LINE                                            ***
C**********************************************************************

      CALL SARGC(TRANSDWLLINE, NPAR)
      IF (NPAR .EQ. 2) THEN
        ACTPAR = 'DENSITY'
      ELSE 
        CALL SARGV (TRANSDWLLINE, 3, ACTPAR)
      ENDIF      
      HEADER         = ' '
      WRITE (HEADER (:8), '(4HLINE,I4)') LINE
      HEADER (11:12) = LEVEL (1:2)
      LAM=NINT(XLAM)
      WRITE (HEADER(13:17),'(I5)') LAM
      HEADER (20:25) = 'MODEL'
C***  DATE AND TIME
      HEADER (27:44) = MODHEAD (15:32)
      HEADER (47:49) = 'J >'
      WRITE (HEADER(51:54),'(I4)') JOBNUM
      
C* Y-AXIS SCALING AND PARAMETERS 
C***  SCALING OF Y-AXIS
      YMIN=.0
      YMAX=.0
      DO 4 L=1,ND
      IF (ADELW(L) .LT. YMIN) YMIN=ADELW(L)
    4 IF (ADELW(L) .GT. YMAX) YMAX=ADELW(L)
      YMIN=AINT(YMIN-0.9999999999)
      YMAX=AINT(YMAX+0.9999999999)
      YSCALE=15./(YMAX-YMIN)
      YTICK=YMAX/4.
      YABST=YTICK*2.0
      YTEXT = CHAR(92) // 'CENTER' // CHAR(92) // 
     >        'NORMALIZED  XI - HILLIER'
C* Plot over radial velocity     
      IF (ACTPAR .EQ. 'VELO') THEN
        XMIN = 0.
        XMAX = VELO(1)
        XSCALE = 0.
        XTICK = 100.
        XABST = 300.
        XTEXT = CHAR(92) // 'CENTER' // CHAR(92) // 
     >        'Radial Velocity (km/s)'        
        DO L=1, ND
          XPLOT(L) = VELO(L)
        ENDDO       
C* Plot over log(velocity)        
      ELSE IF (ACTPAR .EQ. 'VELOLOG') THEN
        XMIN = 0.
        XMAX = ALOG10(VELO(1))
        XSCALE = 0.
        XTICK = 0.25
        XABST = 0.5
        XTEXT = CHAR(92) // 'CENTER' // CHAR(92) // 
     >        'Logarithmic Radial Velocity (Log(v/km/s))'          
        DO L=1, ND
          XPLOT(L) = ALOG10(VELO(L))
        ENDDO
C** OD = plot over optical depth         
      ELSE IF (ACTPAR .EQ. 'TAUROSS') THEN
        XMIN = 0.
        XMAX = TAUROSS(1)
        XSCALE = 0.
        XTICK = 1.
        XABST = 3.
        XTEXT = CHAR(92) // 'CENTER' // CHAR(92) // 
     >        'optical depth'          
        DO L=1, ND
          XPLOT(L) = TAUROSS(L)
        ENDDO
C* Plot over radial distance        
      ELSE IF (ACTPAR .EQ. 'RSTAR') THEN
        XMIN = 0.
        XMAX = RADIUS(1)
        XSCALE = 0.
        XTICK = 100.
        XABST = 300.
        XTEXT = CHAR(92) // 'CENTER' // CHAR(92) // 
     >        'Radial Distance (R*)'        
        DO L=1, ND
          XPLOT(L) = RADIUS(L)
        ENDDO   
C* Plot over log(radius)
      ELSE IF (ACTPAR .EQ. 'RSTARLOG') THEN
        XMIN = 0.
        XMAX = ALOG10(RADIUS(1))
        XSCALE = 0.
        XTICK = 0.25
        XABST = 0.5
        XTEXT = CHAR(92) // 'CENTER' // CHAR(92) // 
     >        'logarithmic Radial Distance (Log(R/R*))'        
        DO L=1, ND
          XPLOT(L) = ALOG10(RADIUS(L))
        ENDDO        
C***  Standard: Plot over log of number density        
      ELSE IF (ACTPAR .EQ. 'DENSITY') THEN
        XMIN = 7.
        XMAX = 16.
        XSCALE = 2.5
        XTICK = 1.0
        XABST = 3.0
        XTEXT = CHAR(92) // 'CENTER' // CHAR(92) // 
     >        'LOG OF NUMBER DENSITY / (CM**-3)'    
        DO L=1, ND
          XPLOT(L) = ALOG10(ENTOT(L))
        ENDDO
      ELSE
        WRITE (0,'(A)') 'DWL x-axis not known'
        WRITE (0,'(A)') 'Syntax: TRANS DWL <VELO>/<VELOLOG>/<TAUROSS>/<RSTAR>/<RSTARLOG>. <> = <DENSITY>' 
        STOP 'ERROR IN SUBROUTINE TRADWL'
      ENDIF
      CALL PLOTANF (KANAL,HEADER,HEADER
     $ ,XTEXT
     $ ,YTEXT
     $ ,XSCALE,XMIN,XMAX,XTICK,XABST,.0
     $ ,YSCALE,YMIN,YMAX,YTICK,YABST,.0
     $ ,XPLOT,ADELW,ND,5)


C*** weiterer Plot
ccc  zur Zeit disabled wegen ungueltiger Ausgabe (suffix U in KASDEF ARC)
      return

      XTEXT = CHAR(92) // 'CENTER' // CHAR(92) // 
     >        'Depth Point Number'
C***  Plot over depth index
      DO L=1, ND
         XPLOT(L) = FLOAT(L)
      ENDDO 
      CALL PLOTANF (KANAL,HEADER,HEADER
     $ ,XTEXT
     $ ,YTEXT
     $ ,0., 0., 80., 5., 10., 0.
     $ ,YSCALE,YMIN,YMAX,YTICK,YABST,.0
     $ ,XPLOT,ADELW,ND,5)
 
      WRITE (KANAL,*) 'PLOT   : DWL-Plot'
      WRITE (KANAL,*) 'KASDEF INBOX'
      WRITE (KANAL,*) 'KASDEF NOBOX'
      WRITE (KANAL,*) 'KASDEF LUN XMIN YMAX 0. D0.1 0.3 ',
     >      'MARKED: 0.1, 3/4*MAX, MAX'
      WRITE (KANAL,*) 'KASDEF DEFINECOLOR 2 0.7 0.7 0.7'
      WRITE (KANAL,*) 'KASDEF COLOR=2'
      WRITE (KANAL,*) 'KASDEF ARC 20. 0.  0. 0.  2.0U 0. 180. FILLED'
      WRITE (KANAL,*) 'KASDEF COLOR=1'
C***  FIND MAXIMUM
      LMAX = 1
      DO L=2, ND
        IF (ADELW(L) .GT. ADELW(LMAX)) LMAX = L
      ENDDO
      R = 2. + FLOAT(LMAX) / 5.
      ADMAX = ADELW(LMAX)
      ADMAX2 = ADMAX * 3. / 4.
      WRITE (KANAL,*) 'KASDEF PEN = 8'
      WRITE (KANAL,'(A,F4.1,A)') 
     >      'KASDEF ARC 20. 0.  0. 0. ',R,'U 0. 180.'
      WRITE (KANAL,*) 'KASDEF PEN = 1'
C***  FIND 0.1
      DO L=2, ND
        IF ((ADELW(L-1)-0.1)*(ADELW(L)-0.1) .LE. 0.) THEN
          R = 2. + FLOAT(L) / 5.
          WRITE (KANAL,*) 'KASDEF PEN = 8'
          WRITE (KANAL,*) 'KASDEF DEFINECOLOR 2 0.7 0.7 0.7'
          WRITE (KANAL,*) 'KASDEF COLOR=2'
          WRITE (KANAL,'(A,F4.1,A)') 
     >          'KASDEF ARC 20. 0.  0. 0. ',R,'U 0. 180.'
          WRITE (KANAL,*) 'KASDEF COLOR=1'
          WRITE (KANAL,*) 'KASDEF PEN = 1'
        ENDIF
        IF ((ADELW(L-1)-ADMAX2)*(ADELW(L)-ADMAX2) .LE. 0.) THEN
          R = 2. + FLOAT(L) / 5.
          WRITE (KANAL,*) 'KASDEF PEN = 8'
          WRITE (KANAL,*) 'KASDEF DEFINECOLOR 2 0.5 0.5 0.5'
          WRITE (KANAL,*) 'KASDEF COLOR=2'
          WRITE (KANAL,'(A,F4.1,A)') 
     >          'KASDEF ARC 20. 0.  0. 0. ',R,'U 0. 180.'
          WRITE (KANAL,*) 'KASDEF PEN = 1'
        ENDIF
      ENDDO      

      WRITE (KANAL,*) 'KASDEF ARC 20. 0.  0. 0. 3.0U 0. 180.'
      WRITE (KANAL,*) 'KASDEF LUN 23. 0. M0. U-0.1 0.2 65'
      WRITE (KANAL,*) 'KASDEF ARC 20. 0.  0. 0. 4.0U 0. 180.'
      WRITE (KANAL,*) 'KASDEF LUN 24. 0. M0. U-0.1 0.2 60'
      WRITE (KANAL,*) 'KASDEF ARC 20. 0.  0. 0. 5.0U 0. 180.'
      WRITE (KANAL,*) 'KASDEF LUN 25. 0. M0. U-0.1 0.2 55'
      WRITE (KANAL,*) 'KASDEF ARC 20. 0.  0. 0. 6.0U 0. 180.'
      WRITE (KANAL,*) 'KASDEF LUN 26. 0. M0. U-0.1 0.2 50'
      WRITE (KANAL,*) 'KASDEF ARC 20. 0.  0. 0. 7.0U 0. 180.'
      WRITE (KANAL,*) 'KASDEF LUN 27. 0. M0. U-0.1 0.2 45'
      WRITE (KANAL,*) 'KASDEF ARC 20. 0.  0. 0. 8.0U 0. 180.'
      WRITE (KANAL,*) 'KASDEF LUN 28. 0. M0. U-0.1 0.2 40'
      WRITE (KANAL,*) 'KASDEF ARC 20. 0.  0. 0. 9.0U 0. 180.'
      WRITE (KANAL,*) 'KASDEF LUN 29. 0. M0. U-0.1 0.2 35'
      WRITE (KANAL,*) 'KASDEF ARC 20. 0.  0. 0. 10.0U 0. 180.'
      WRITE (KANAL,*) 'KASDEF LUN 30. 0. M0. U-0.1 0.2 30'
      WRITE (KANAL,*) 'KASDEF ARC 20. 0.  0. 0. 11.0U 0. 180.'
      WRITE (KANAL,*) 'KASDEF LUN 31. 0. M0. U-0.1 0.2 25'
      WRITE (KANAL,*) 'KASDEF ARC 20. 0.  0. 0. 12.0U 0. 180.'
      WRITE (KANAL,*) 'KASDEF LUN 32. 0. M0. U-0.1 0.2 20'
      WRITE (KANAL,*) 'KASDEF ARC 20. 0.  0. 0. 13.0U 0. 180.'
      WRITE (KANAL,*) 'KASDEF LUN 33. 0. M0. U-0.1 0.2 15'
      WRITE (KANAL,*) 'KASDEF ARC 20. 0.  0. 0. 14.0U 0. 180.'
      WRITE (KANAL,*) 'KASDEF LUN 34. 0. M0. U-0.1 0.2 10'
      WRITE (KANAL,*) 'KASDEF ARC 20. 0.  0. 0. 15.0U 0. 180.'
      WRITE (KANAL,*) 'KASDEF LUN 35. 0. M0. U-0.1 0.2 5'
      WRITE (KANAL,*) 'KASDEF ARC 20. 0.  0. 0. 16.0U 0. 180.'
      WRITE (KANAL,*) 'KASDEF LUN 36. 0. M0. U-0.1 0.2 1'

      XTEXT = CHAR(92) // 'CENTER' // CHAR(92) // 
     >        'Depth Points'
c      YTEXT = CHAR(92) // 'CENTER' // CHAR(92) // 
c     >        'NORMALIZED  XI - HILLIER'
      CALL PLOTANF (KANAL,HEADER,HEADER
     $ ,XTEXT
     $ ,' '
     $ ,0., 0., 40., 25., 25., 0. 
     $ ,0., 0., 30., 20., 20., 0.
     $ ,XDUM, XDUM, 1, 8)

      RETURN
      END
