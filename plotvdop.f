      SUBROUTINE PLOTVDOP (RADIUS, TAUROSS, VELO, DD_VDOP, ND, NATOM, 
     >                     BPLOTVDOP, SYMBOL, BMICROTURB, DD_VMIC,
     >                     DD_VDOP_LINE, VDOPPLOT_LINE, MODHEAD, 
     >                     BIRONLINES)
C******************************************************************************
C***  Plot of the Dopper-broadening velocity versus v(r) or tauross
C***  Called from: FORMAL
C******************************************************************************
      IMPLICIT NONE

      INTEGER, INTENT(IN) :: ND, NATOM
      REAL, DIMENSION(ND,NATOM) :: DD_VDOP 
      REAL, DIMENSION(ND) :: RADIUS, DD_VMIC, TAUROSS, VELO
      REAL, DIMENSION(ND) :: X, Y
      CHARACTER*60 HEAD1, HEAD2, STYLE
      CHARACTER*(*) DD_VDOP_LINE, VDOPPLOT_LINE, MODHEAD
      CHARACTER*50 ACTPAR1, ACTPAR2, XAXSTR, YAXSTR
      CHARACTER*2 SYMBOL(NATOM)
      INTEGER :: KANAL, L, NA, NPAR1, NPAR2
      REAL :: XMIN, XMAX, YMIN, YMAX, VCON,
     >        XSCALE, XTICK, XABST,
     >        YSCALE, YTICK, YABST
     
      LOGICAL :: BPLOTVDOP, BMICROTURB, BIRONLINES
 
C***  INITIALIZATION
      KANAL=1
      OPEN (KANAL, FILE='PLOT', STATUS='UNKNOWN')
      CALL JSYMSET ('G2','TRANSFER')
 
C***  HEADER  ------------------------------------------------------
      HEAD1='VDOP'
      HEAD2= 'M' //  MODHEAD(12:32) // ' Depth-dependent VDOP'

! C***  X-AXIS: ----------------------------------------------------------
! C***  RADIUS RANGE: -3.0 <= LOG(R/R*-1) <= 2.9
!       XMAX=2.9
!       XMIN=-3.0
!       XSCALE = 22./(XMAX-XMIN)
!       XTICK=0.5
!       XABST=1.
!       
!       
! C***  Y-AXIS:  ---------------------------------------------------------
! C***  VELOCITY RANGE [IN 100 KM/S]: -1. <= V(R) <= VFINAL+1.
!       YMAX = MAXVAL(DD_VDOP)
!       YMAX = YMAX + 0.05*YMAX
!       YMIN = .0
!       YSCALE = 15./(YMAX-YMIN)
!       YTICK=2.5
!       YABST=5.
!         
C*** First data set
      CALL SARGC(VDOPPLOT_LINE, NPAR1) 
      IF (NPAR1 .EQ. 3) THEN
        CALL SARGV (VDOPPLOT_LINE, 3, ACTPAR1)
      ELSE
        ACTPAR1 = ''
      ENDIF
      CALL SARGC(DD_VDOP_LINE, NPAR2)
      IF (NPAR2 .GE. 5) THEN
        CALL SARGV (DD_VDOP_LINE, 5, ACTPAR2)
      ELSE
        ACTPAR2 = ''
      ENDIF
      IF (ACTPAR1 .EQ. '') THEN
        IF ((ACTPAR2 .EQ. '') .OR. (ACTPAR2(:4) .EQ. 'VELO')) THEN
            DO L=1, ND-1
                X(L) = VELO(L)   
            ENDDO
            XAXSTR = '\CENTER\Wind velocity [km/s]'
        ELSE IF (ACTPAR2(:3) .EQ. 'TAU') THEN
            DO L=1, ND-1
                X(L) = TAUROSS(L)
            ENDDO
            XAXSTR = '\CENTER\#t#'
        ELSE IF (ACTPAR2(:1) .EQ. 'R') THEN
            DO L=1, ND-1
            X(L) = ALOG10(RADIUS(L)-1.) 
            ENDDO
            XAXSTR = '\CENTER\log (r/R&T*&M - 1)'
        ELSE 
            WRITE(0,*) "VDOP Plot version not known"
            STOP "Fatal error in subroutine PLOTVDOP"
        ENDIF
      ELSE IF (ACTPAR1 .EQ. 'R') THEN
        DO L=1, ND-1
        X(L) = ALOG10(RADIUS(L)-1.) 
        ENDDO
        XAXSTR = '\CENTER\log (r/R&T*&M - 1)'
      ELSE IF (ACTPAR1 .EQ. 'VELO') THEN
        DO L=1, ND-1
            X(L) = VELO(L)   
        ENDDO
        XAXSTR = '\CENTER\Wind velocity [km/s]'
      ELSE IF (ACTPAR1 .EQ. 'TAU') THEN
        DO L=1, ND-1
            X(L) = TAUROSS(L)
        ENDDO
        XAXSTR = '\CENTER\#t#'
      ELSE
        WRITE(0,*) "VDOP Plot version not known"
        STOP "Fatal error in subroutine PLOTVDOP"
      ENDIF
      DO L=1, ND-1
        Y(L) = DD_VDOP(L,1)
      ENDDO
      YAXSTR = '\CENTER\Doppler broadening velocity VDOP [km/s]'

C***  X-AXIS: ----------------------------------------------------------
C***  RADIUS RANGE: -3.0 <= LOG(R/R*-1) <= 2.9
      XMAX=MAXVAL(X(:ND-1))
      XMIN=MINVAL(X(:ND-1))
      XSCALE = 22./(XMAX-XMIN)
      XTICK=0.5
      XABST=1.

C***  Y-AXIS:  ---------------------------------------------------------
C***  VELOCITY RANGE [IN 100 KM/S]: -1. <= V(R) <= VFINAL+1.
      YMAX = MAXVAL(Y(:ND-1))
      YMIN = 0.
      YSCALE = 15./(YMAX-YMIN)
      YTICK=2.5
      YABST=5.
      
      WRITE (KANAL, '(A,A)') 'PLOT: ', HEAD1
      WRITE (KANAL, '(A)') '\INBOX'
      WRITE (KANAL, '(A)') '\FONT=HELVET'
      WRITE (KANAL, '(A)') '\PENDEF=3'
      WRITE (KANAL, '(A)') '\DEFINECOLOR 3 = .0 .5 .0'

      
C*** Legend
      IF (BMICROTURB) THEN
       WRITE (KANAL, '(A)') '\COLOR=1'
       WRITE (KANAL, '(A)') "\LINREL XMAX YMAX  1 0 0.5 -1"
       WRITE (KANAL, '(A, A)') "\LUN XMAX YMAX L1.7 M-1 0.3 ", SYMBOL(1)
       DO NA=2, NATOM
         WRITE (KANAL, '(A, I2)') '\COLOR=',NA 
         WRITE (KANAL, '(A, F10.5)') 
     >         "\LINREL XMAX YMAX  1 0 0.5", -1.-0.6*(NA-1)
         WRITE (KANAL, '(A, A)') "\NEXTLUN ", SYMBOL(NA)
       ENDDO
      
       WRITE (KANAL, '(A)') '\COLOR=1'
       WRITE (KANAL, '(A, F10.5, A)') "\LINREL XMAX YMAX  1 0 0.5", 
     >                      -1.-0.6*NATOM, " SYMBOL=9 SIZE=0.07"
        WRITE (KANAL, '(A)') "\NEXTLUN v&Tmic&M"
      ELSE
        IF (BIRONLINES) THEN
         WRITE (KANAL, '(A,I3)') '\COLOR=', NATOM
         WRITE (KANAL, '(A)') "\LINREL XMAX YMAX  1 0 0.5 -1"
         WRITE (KANAL, '(A, A)') "\LUN XMAX YMAX L1.7 M-1 0.3 VDOP-FEDAT"
        ENDIF
      ENDIF

      WRITE (KANAL, '(A)') '\COLOR=1'
      CALL PLOTANF (KANAL,HEAD1,HEAD2
     >        ,XAXSTR
     >        ,YAXSTR
     >        ,XSCALE,0.,0.,XTICK,XABST,.0
     >        ,YSCALE,YMIN,YMAX,YTICK,YABST,.0
     >        ,X,Y,ND-1, 5)
      DO NA=2, NATOM
        IF (.NOT. BMICROTURB .AND. .NOT. SYMBOL(NA) .EQ. 'G ') CYCLE
        DO L=1, ND-1
            Y(L) = DD_VDOP(L,NA)
        ENDDO
        STYLE = 'SYMBOL=5 COLOR='
        WRITE(STYLE(15:16), '(I2)') NA
        CALL PLOTCONS (KANAL,X,Y,ND-1, STYLE) 
      ENDDO
      IF (BMICROTURB) THEN
        DO L=1, ND-1
            Y(L) = DD_VMIC(L)
        ENDDO
        STYLE = 'SYMBOL=9 COLOR=1 SIZE=0.07'
        CALL PLOTCONS (KANAL,X,Y,ND-1, STYLE) 
      ENDIF

C***  Set back to initial values (for next blend-range):
      BPLOTVDOP = .FALSE.
      
      RETURN
      END
