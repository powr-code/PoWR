      SUBROUTINE PLOTV (ND ,RADIUS, VELO, MODHEAD, JOBNUM)
C******************************************************************************
C***  DIRECT PLOT TRANSFER OF THE VELOCITY STRATIFICATION
C***              V(R) VERSUS LOG(R/R*-1)
C******************************************************************************

      IMPLICIT NONE
      INCLUDE 'interfacebib.inc'

      INTEGER, PARAMETER :: NDMAX = 100
      INTEGER, INTENT(IN) :: ND, JOBNUM
      CHARACTER(100), INTENT(IN) :: MODHEAD
      REAL, DIMENSION(ND) :: RADIUS, VELO

      !Velocity law common block
      REAL :: VFINAL, VMIN, BETA, VPAR1, VPAR2, RCON, HSCALE, 
     >        BETA2, BETA2FRACTION, VPAR1_2, VPAR2_2
      COMMON /VELPAR/ VFINAL,VMIN,BETA,VPAR1,VPAR2,RCON,HSCALE,
     >                BETA2, BETA2FRACTION, VPAR1_2, VPAR2_2

      REAL, DIMENSION(NDMAX) :: X, Y
      CHARACTER(60) :: HEAD1, HEAD2

      INTEGER :: KANAL, L
      REAL :: XMIN, XMAX, YMIN, YMAX, VCON,
     >        XSCALE, XTICK, XABST,
     >        YSCALE, YTICK, YABST
 
C***  INITIALIZATION
      KANAL=2
      OPEN (KANAL, FILE='PLOT', STATUS='UNKNOWN')
      CALL JSYMSET ('G2','TRANSFER')
      CALL REMARK ('V-STRATIFICATION TO BE ROUTED')
 
C***  HEADER  ------------------------------------------------------
      HEAD1=' WR VELOCITY STRATIFICATION V(R) VERSUS LOG(R/R*-1)'
      HEAD2      = MODHEAD(13:32)
      HEAD2(22:) = 'VELOCITY STRATIFICATION v(r)'

C***  X-AXIS: ----------------------------------------------------------
C***  RADIUS RANGE: -3.0 <= LOG(R/R*-1) <= 2.9
      XMAX=2.9
      XMIN=-3.0
      XSCALE = 22./(XMAX-XMIN)
      XTICK=0.5
      XABST=1.

C***  Y-AXIS:  ---------------------------------------------------------
C***  VELOCITY RANGE [IN 100 KM/S]: -1. <= V(R) <= VFINAL+1.
      YMAX = ANINT(VFINAL/100.+1.)
      YMIN = .0
      YSCALE = 15./(YMAX-YMIN)
      YTICK=2.5
      YABST=5.

C***  DATA TABLE ------------------------------------
      DO L=1, ND-1
        X(L) = ALOG10(RADIUS(L)-1.)
        Y(L) = VELO(L) / 100.
      ENDDO
 
      WRITE (KANAL, '(A,A)') 'PLOT: ', HEAD1
      WRITE (KANAL, '(A)') '\INBOX'
C***  BORDER BETWEEN INNER/OUTER VELOCITY LAW
      IF ((RCON > RADIUS(ND)) .AND. (RCON <= RADIUS(1))) THEN
        CALL SPLINPOX(VCON, RCON, VELO, RADIUS, ND)
        WRITE (KANAL, '(A,2F10.3,A)')
     >    '\SYM ', LOG10(RCON-1.), VCON / 100., ' 0. 0. 0.3 4'
      ENDIF
      WRITE (KANAL, '(A,F10.3)') '*RCON=', RCON
      DO L=1, ND-1
        IF (MOD(L, 10) == 0) THEN
          WRITE (KANAL,41) X(L), X(L), X(L), L
   41     FORMAT ('\LINREL ', F7.3, ' YMAX 0. -0.5', /,
     >            '\LINREL ', F7.3, ' YMIN 0.  0.5', /,
     >            '\LUN    ', F7.3, ' YMIN -0.2 0.7 0.3 ', I3)
        ENDIF
      ENDDO

      CALL PLOTANF (KANAL,HEAD1,HEAD2
     >        ,'\CENTER\log (r/R&T*&M - 1)'
     >        ,'\CENTER\v(r) / [100 km/s]'
     >        ,XSCALE,0.,0.,XTICK,XABST,.0
     >        ,YSCALE,YMIN,YMAX,YTICK,YABST,.0
     >        ,X,Y,ND-1, 5)
      CALL PLOTCONS (KANAL,X,Y,ND-1,'SYMBOL=8 SIZE=0.1') 
 
      RETURN
      END
