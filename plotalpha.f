      SUBROUTINE PLOTALPHA (ND ,RADIUS, ALPHAF, MODHEAD, JOBNUM, bOwn)
C******************************************************************************
C***  PLOT FORCE MULTIPLIER PARAMETER ALPHA VS DEPTH POINT OR LOG(R/R*-1)
C******************************************************************************

      IMPLICIT NONE

      INTEGER, PARAMETER :: NDMAX = 100
      INTEGER, INTENT(IN) :: ND, JOBNUM
      CHARACTER(100), INTENT(IN) :: MODHEAD
      REAL, DIMENSION(ND) :: RADIUS, ALPHAF
      LOGICAL, INTENT(IN) :: bOwn           !own plot file or part of wruniq.plot?

      REAL, DIMENSION(NDMAX) :: X, Y
      CHARACTER(60) :: HEAD1, HEAD2
      CHARACTER(10) :: CTIME
      CHARACTER(8) :: CDATE

      INTEGER :: L
      REAL :: XMIN, XMAX, YMIN, YMAX, 
     >        XSCALE, XTICK, XABST,
     >        YSCALE, YTICK, YABST

      !File and channel handles (=KANAL)
      INTEGER, PARAMETER :: hOUT = 6        !write to wruniqX.out (stdout)
      INTEGER, PARAMETER :: hCPR = 0        !write to wruniqX.cpr (stderr)
      INTEGER, PARAMETER :: hPLOT = 2       !write to plot file
 
C***  INITIALIZATION
      IF (bOwn) THEN
        OPEN (hPLOT, FILE='alpha.plot', STATUS='UNKNOWN')
        WRITE (hCPR, '(A)') 'FORCE MULTIPLIERS PLOTTED IN alpha.plot'
      ELSE
        OPEN (hPLOT, FILE='PLOT', STATUS='UNKNOWN')
        CALL JSYMSET ('G2','TRANSFER')
        CALL REMARK ('FORCE MULTIPLIERS TO BE ROUTED')
      ENDIF
        
 
C***  HEADER  ------------------------------------------------------
      HEAD1=' FORCE MULTIPLIER PARAMETER ALPHA VERSUS DEPTH INDEX'
      HEAD2      = MODHEAD(13:32)
      HEAD2(22:) = 'FORCE MULTIPLIER ALPHA(L)'

C***  X-AXIS: ----------------------------------------------------------
C***  DEPTH POINT INDEX
C***  @TODO: IMPLEMENT RADIUS OPTION
      XMAX=FLOAT(ND)
      XMIN=0.
      XSCALE = 22./(XMAX-XMIN)
      XTICK=5.
      XABST=10.

C***  Y-AXIS:  ---------------------------------------------------------
C***  ALPHA RANGE (AUTO)
      YMAX = 0. 
      YMIN = .0
      YSCALE = 0.
      YTICK=2.5 !unwichtig, da AUTO-Option gesetzt
      YABST=5.  !unwichtig, da AUTO-Option gesetzt

C***  DATA TABLE ------------------------------------
      DO L=1, ND-1
        X(L) = FLOAT(L)
        Y(L) = ALPHAF(L)
      ENDDO
 
      WRITE (hPLOT, '(A,A)') 'PLOT: ', HEAD1
      WRITE (hPLOT, '(A)') '\INBOX'
      WRITE (hPLOT, '(A)') '\FONT=HELVET'
      WRITE (hPLOT, '(A)') '\DEFINECOLOR=9 0.5 0.5 0.5'
      WRITE (hPLOT, '(A)') '\COLOR=9'
      WRITE (hPLOT, '(A)') '\LINUN XMIN 0. XMAX 0. 0. 0. SYMBOL=9'
      WRITE (hPLOT, '(A)') '\COLOR=1'

      CALL DATE_AND_TIME(CDATE, CTIME)

      WRITE (hPLOT, '(A,I7)') 
     >   '\LUNA XMAX YMAX 0.7 0. 0.4 -90. JOB No.', JOBNUM

      WRITE (hPLOT, '(11A)') 
     >   '\LUNA XMAX YMIN R0.7 0. 0.3 -90. (calculated at ',
     >   CDATE(1:4), '/', CDATE(5:6), '/', CDATE(7:8), ' ', 
     >   CTIME(1:2), ':', CTIME(3:4), ')'

      CALL PLOTANF (hPLOT,HEAD1,HEAD2
     >        ,'\CENTER\depth index'
     >        ,'\CENTER\force multiplier parameter #a#'
     >        ,XSCALE,XMIN,XMAX,XTICK,XABST,.0
     >        ,YSCALE,0.,0.,YTICK,YABST,.0
     >        ,X,Y,ND-1, 5)
 
C      CLOSE(hPLOT)
      
      RETURN
      END
