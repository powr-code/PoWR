      SUBROUTINE PLOTVGRAD (ND ,RADIUS, VELO, MODHEAD, JOBNUM, 
     >                      RMIN, RMAX)
C******************************************************************************
C***  PLOT OF THE VELOCITY GRADIENT STRATIFICATION
C***              DV/DR(R) VERSUS LOG(R/R*-1)
C******************************************************************************

      IMPLICIT NONE
      INCLUDE 'interfacebib.inc'

      INTEGER, PARAMETER :: NDMAX = 100
      INTEGER, INTENT(IN) :: ND, JOBNUM
      CHARACTER(100), INTENT(IN) :: MODHEAD
      REAL, INTENT(IN) :: RMIN, RMAX
      REAL, DIMENSION(ND) :: RADIUS, VELO, VELO2, GRAD1, GRAD2

      !Velocity law common block
      REAL :: VFINAL, VMIN, BETA, VPAR1, VPAR2, RCON, HSCALE,
     >        BETA2, BETA2FRACTION, VPAR1_2, VPAR2_2
      COMMON / VELPAR / VFINAL,VMIN,BETA,VPAR1,VPAR2,RCON,HSCALE,
     >                BETA2, BETA2FRACTION, VPAR1_2, VPAR2_2

      REAL, EXTERNAL :: WRVEL

      REAL, DIMENSION(NDMAX) :: X, Y
      CHARACTER(60) :: HEAD1, HEAD2

      INTEGER :: hPLOT, L, NDBG
      REAL :: XMIN, XMAX, YMIN, YMAX, Vdummy, RP, RM,
     >        XSCALE, XTICK, XABST,
     >        YSCALE, YTICK, YABST
 
      !calculate Gradients
      DO L=1, ND
        IF (VELO(L) > 0.) THEN
          CALL SPLINPOX(Vdummy,RADIUS(L),VELO,RADIUS,ND,DFDX=GRAD1(L))
        ELSE 
          GRAD1(L) = 0.
        ENDIF
      ENDDO

      !GRAD2 = Beta law gradient
      DO L=1, ND
        IF (RADIUS(L) + VPAR2 <= 1.) THEN
          EXIT 
        ENDIF
        RP = VPAR2 + RADIUS(L)
        GRAD2(L) = VPAR1*BETA*(1.-1./RP)**(BETA-1.) /(RP*RP) 
        NDBG = L
      ENDDO

C***  INITIALIZATION
      hPLOT=2
      OPEN (hPLOT, FILE='vgrad.plot', STATUS='UNKNOWN')
C      CALL JSYMSET ('G2','TRANSFER')
C      CALL REMARK ('V-STRATIFICATION TO BE ROUTED')
 
C***  HEADER  ------------------------------------------------------
      HEAD1=' WR VELOCITY GRADIENT STRATIFICATION DV/DR(R)' //
     >          ' VERSUS LOG(R/R*-1)'
      IF (MODHEAD(:5) /= 'DEBUG') THEN
        HEAD2      = MODHEAD(13:32)
        HEAD2(22:) = 'VELOCITY GRADIENT STRATIFICATION dv/dr(r)'
      ELSE
        HEAD2 = 'VELOCITY GRADIENT STRATIFICATION dv/dr(r)'
      ENDIF

C***  X-AXIS: ----------------------------------------------------------
C***  RADIUS RANGE: -3.0 <= LOG(R/R*-1) <= 2.9
      XMAX= LOG10(RADIUS(1)-1.)
      XMIN= LOG10( (RADIUS(ND-1)-1.) * 0.9 ) 
      XSCALE = 22./(XMAX-XMIN)
      XTICK=0.5
      XABST=1.

C***  Y-AXIS:  ---------------------------------------------------------
C***  VELOCITY GRADIENT RANGE
C*** ( 1. + (BETA-1.)/2. - VPAR2 ) * 1.5
      YMAX = 1.5 * MAXVAL(GRAD2(1:NDBG))   !150% of maximum beta law gradient
      YMIN = -9.0
      YSCALE = 15./(YMAX-YMIN)
      YTICK= INT( (YMAX-YMIN) / 20. )
      YABST=5. * YTICK

C***  DATA TABLE ------------------------------------
      DO L=1, ND-1
        X(L) = ALOG10(RADIUS(L)-1.)
C        Y(L) = VELO(L) / 100.
        Y(L) = GRAD1(L)
      ENDDO
 
      WRITE (hPLOT, '(A,A)') 'PLOT: ', HEAD1
      WRITE (hPLOT, '(A)') '\INBOX'

C***  BORDER BETWEEN INNER/OUTER VELOCITY LAW
      WRITE (hPLOT, '(A)') '\BGRLUN COLOR=0'
      IF (RCON > 1.) THEN
        WRITE (hPLOT,'(A)') '\COLOR=9'
        WRITE (hPLOT,'(A,F20.12,A,F20.12,A)') '\LINUN ',LOG10(RCON-1.),
     >    ' YMIN ',LOG10(RCON-1.),' YMAX 0. 0. SYMBOL=10'
        WRITE (hPLOT, '(A,F20.12,A)') 
     >      '\LUNA ',LOG10(RCON-1.),' YMAX 0. -0.2 0.2 -90  R&Tcon&M'
      ENDIF
      IF (RMIN > 1.) THEN
        WRITE (hPLOT,'(A)') '\COLOR=5'
        WRITE (hPLOT,'(A,F20.12,A,F20.12,A)') '\LINUN ',LOG10(RMIN-1.),
     >    ' YMIN ',LOG10(RMIN-1.),' YMAX 0. 0. SYMBOL=9'
        WRITE (hPLOT, '(A,F20.12,A)') 
     >      '\LUNA ',LOG10(RMIN-1.),' YMAX 0. -0.2 0.2 -90  min.'
      ENDIF
      IF (RMAX > 1.) THEN
        WRITE (hPLOT,'(A)') '\COLOR=6'
        WRITE (hPLOT,'(A,F20.12,A,F20.12,A)') '\LINUN ',LOG10(RMAX-1.),
     >    ' YMIN ',LOG10(RMAX-1.),' YMAX 0. 0. SYMBOL=9'
        WRITE (hPLOT, '(A,F20.12,A)') 
     >      '\LUNA ',LOG10(RMAX-1.),' YMAX 0. -0.2 0.2 -90  max.'
      ENDIF
      WRITE (hPLOT, '(A)') '\BGRLUN OFF'
      WRITE (hPLOT, '(A)') '\COLOR=1'


      CALL PLOTANF (hPLOT,HEAD1,HEAD2
     >        ,'\CENTER\log (r/R&T*&M - 1)'
     >        ,'\CENTER\dv/dr(r) / [km/s/Rstar]'
     >        ,XSCALE,XMIN,XMAX,XTICK,XABST,.0
     >        ,YSCALE,YMIN,YMAX,YTICK,YABST,.0
     >        ,X,Y,ND-1, 5)
 
      NDBG = MIN(NDBG, ND-1)
      CALL PLOTCONS (hPLOT,X,GRAD2,NDBG,'COLOR=2 SYMBOL=9 SIZE=0.1') 

  
      CLOSE(hPLOT)

      RETURN
      END
