      SUBROUTINE PLOTGRADI (ND ,RADIUS, GRADI, MODHEAD, JOBNUM)
C******************************************************************************
C***  DIRECT PLOT TRANSFER OF THE VELOCITY GRADIENT
C***              Vector GRADI VERSUS LOG(R/R*-1)
C******************************************************************************

      IMPLICIT NONE
      INCLUDE 'interfacebib.inc'

      INTEGER, PARAMETER :: NDMAX = 100
      INTEGER, INTENT(IN) :: ND, JOBNUM
      CHARACTER(100), INTENT(IN) :: MODHEAD
      REAL, DIMENSION(ND) :: RADIUS, GRADI

      !Velocity law common block
      REAL :: VFINAL, VMIN, BETA, VPAR1, VPAR2, RCON, HSCALE, 
     >        BETA2, BETA2FRACTION, VPAR1_2, VPAR2_2
      COMMON /VELPAR/ VFINAL,VMIN,BETA,VPAR1,VPAR2,RCON,HSCALE,
     >                BETA2, BETA2FRACTION, VPAR1_2, VPAR2_2

      REAL, DIMENSION(NDMAX) :: X, Y
      CHARACTER(60) :: HEAD1, HEAD2

      INTEGER :: KANAL, L
      REAL :: XSCALE, YSCALE, RCONLOG
 
C***  INITIALIZATION
      KANAL=2
      OPEN (KANAL, FILE='PLOT', STATUS='UNKNOWN')
      CALL JSYMSET ('G2','TRANSFER')
      CALL REMARK ('GRADI-Plot TO BE ROUTED')
 
C***  HEADER  ------------------------------------------------------
      HEAD1=' VELOCITY GRADIENT VERSUS LOG(R/R*-1)'
      HEAD2      = MODHEAD(13:32)
      HEAD2(22:) = 'VELOCITY GRADIENT'

C***  AUTO option for both axes
      XSCALE = 0.
      YSCALE = 0.

C***  DATA TABLE ------------------------------------
      DO L=1, ND-1
        X(L) = ALOG10(RADIUS(L)-1.)
        Y(L) = GRADI(L) / 100.
      ENDDO
 
      WRITE (KANAL, '(A,A)') 'PLOT: ', HEAD1
      WRITE (KANAL, '(A)') '\INBOX'
C***  BORDER BETWEEN INNER/OUTER VELOCITY LAW
      IF ((RCON > RADIUS(ND)) .AND. (RCON <= RADIUS(1))) THEN
        RCONLOG = ALOG10(RCON-1.)
        WRITE (KANAL, '(A,F10.6,A,F10.6,A)')
     >    '\LINUNLAB ', RCONLOG, ' YMIN ', RCONLOG,
     >    ' YMAX 0. 0. 0. .4 R&Tcon&M'
      ENDIF

      CALL PLOTANF (KANAL,HEAD1,HEAD2
     >        ,'\CENTER\log (r/R&T*&M - 1)'
     >        ,'\CENTER\dv/dr / [100 km/s /R\*]'
     >        ,XSCALE, 0., 0., 0., 0.,.0
     >        ,YSCALE, 0., 0., 0., 0.,.0
     >        ,X,Y,ND-1, 5)
      CALL PLOTCONS (KANAL,X,Y,ND-1,'SYMBOL=8 SIZE=0.1') 
 
      RETURN
      END
