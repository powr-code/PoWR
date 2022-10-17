      SUBROUTINE PLOTLIMB (BLIMB, P, SUMINT, NP, JPFIRST, 
     >                    JPLASTLIMB, LIMB_LINE, LIMB_UNIT)
C******************************************************************************
C***  DIRECT PLOT TRANSFER OF THE VELOCITY STRATIFICATION
C***              V(R) VERSUS LOG(R/R*-1)
C******************************************************************************
***      IMPLICIT NONE

      INTEGER, INTENT(IN) :: NP, JPFIRST, JPLASTLIMB, LIMB_UNIT
      REAL, DIMENSION(NP) :: P, SUMINT
      LOGICAL :: BLIMB
      REAL, DIMENSION(NP) :: X, Y
      CHARACTER :: LIMB_LINE*(*)
      CHARACTER(100) :: HEAD1, HEAD2, STYLE, ACTPAR,  
     >                 XAXISSTR, YAXISSTR
      CHARACTER(7) :: PLASTSTR
      INTEGER :: KANAL, JP, NPAR, I
      REAL :: XMIN, XMAX, YMIN, YMAX, VCON,
     >        XSCALE, XTICK, XABST,
     >        YSCALE, YTICK, YABST
     

      CALL SARGV (LIMB_LINE, 1, ACTPAR)
C***  intensities are to be written ==> skip this routine
      IF (ACTPAR .EQ. 'WRITE') THEN
         CLOSE(LIMB_UNIT)
         RETURN
      ENDIF
C***  INITIALIZATION
      KANAL=1
      OPEN (KANAL, FILE='PLOT', STATUS='UNKNOWN')
      CALL JSYMSET ('G2','TRANSFER')
C***  HEADER  ------------------------------------------------------
      HEAD1='LIMB'
      CALL SARGV (LIMB_LINE, 3, ACTPAR)
      WRITE(PLASTSTR, '(F7.2)') P(JPLASTLIMB)
      IF (ACTPAR(:3) .EQ. 'LAM') THEN
         CALL SARGV (LIMB_LINE, 4, ACTPAR)
         HEAD2 = 'Emergent intensities at #l# = ' // 
     >           ACTPAR(:IDX(ACTPAR)) // '\A'
      ELSE IF (ACTPAR .EQ. 'BAND') THEN
         CALL SARGV (LIMB_LINE, 4, ACTPAR)
         HEAD2 = "Emergent intensities in filter " // ACTPAR 
      ENDIF
      HEAD2 = HEAD2(:IDX(HEAD2)) // " till p = " 
      HEAD2 = HEAD2(:IDX(HEAD2)) // PLASTSTR
      DO JP=JPFIRST, JPLASTLIMB
        X(JP) = P(JP)
        Y(JP) = SUMINT(JP)
      ENDDO  
      XAXISSTR = '\CENTER\p [R&T*&M]'
      YAXISSTR = 
     >   '\CENTER\I&T#n#&M [erg cm&H-2&M s&H-1&M #n#&H-1&M sterad&H-1&M]'
      CALL SARGC(LIMB_LINE, NPAR)
      DO I=1, NPAR
            CALL SARGV(LIMB_LINE, I, ACTPAR)
            IF (ACTPAR .EQ. 'NORMALIZED') THEN
                DO JP=JPFIRST, JPLASTLIMB
                    Y(JP) = SUMINT(JP) / SUMINT(JPFIRST)
                ENDDO
                YAXISSTR = '\CENTER\I&T#n#&M / I&T#n#&M(p=0)'
            ELSE IF (ACTPAR .EQ. 'MU') THEN
                DO JP=JPFIRST, JPLASTLIMB               
                   X(JP) = P(JP) / P(JPLASTLIMB)        
                ENDDO
                XAXISSTR = '\CENTER\#m#'
                YAXISSTR = '\CENTER\I&T#n#&M / I&T#n#&M(#m#=0)'
            ENDIF
      ENDDO
      XMAX=MAXVAL(X(JPFIRST:JPLASTLIMB))
      XMIN=MINVAL(X(JPFIRST:JPLASTLIMB))
      XSCALE = 22./(XMAX-XMIN)
      XTICK=0.5
      XABST=1.
            
      YMAX = MAXVAL(Y(JPFIRST:JPLASTLIMB)) 
      YMIN = MINVAL(Y(JPFIRST:JPLASTLIMB))
      YSCALE = 15./(YMAX-YMIN)
      YTICK=2.5
      YABST=5.
      WRITE (KANAL, '(A,A)') 'PLOT: ', HEAD1
      WRITE (KANAL, '(A)') '\INBOX'

      WRITE (KANAL, '(A)') '\COLOR=1'
      CALL PLOTANF (KANAL,HEAD1,HEAD2
     >        ,XAXISSTR
     >        ,YAXISSTR
     >        ,XSCALE,0.,0.,XTICK,XABST,.0
     >        ,YSCALE,YMIN,YMAX,YTICK,YABST,.0
     >        ,X,Y,JPLASTLIMB-JPFIRST, 5)
      

      BLIMB = .FALSE.
      
      RETURN
      
  100 WRITE(0,'(A)') "Wrong syntax"
      WRITE(0,'(A)') "The error occured in the following line"
      WRITE(0,'(A)') LIMB_LINE
      WRITE(0,'(A)') "Possible syntax:"
      WRITE(0,'(A)') "PLOT LIMB BAND UJ|BJ|VJ|US|VS|BS|YS [[MU-]NORMALIZED]"
      WRITE(0,'(A)') "PLOT LIMB LAM x.x [[MU] NORMALIZED]"
      STOP 'Fatal error in subroutine PLOTLIMB'
      
      
      END
