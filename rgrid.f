      SUBROUTINE RGRID (NDDIM,ND,R     ,INCRIT,RMAX, 
     >                  RadiusGridParameters)

C***********************************************************************
C***  SPECIFICATION OF THE GRID OF DEPTH-POINTS   ******************************
C***********************************************************************

      IMPLICIT NONE
      INCLUDE 'interfacebib.inc'

      INTEGER, INTENT(IN) :: NDDIM
      INTEGER, INTENT(INOUT) :: ND 
      REAL, INTENT(IN) :: RMAX
      REAL, DIMENSION(NDDIM), INTENT(INOUT) :: R
      CHARACTER(8), DIMENSION(NDDIM), INTENT(INOUT) :: INCRIT
      CHARACTER(80), DIMENSION(3), INTENT(IN) :: RadiusGridParameters


      INTEGER :: ND2, ND3, NSPEC
      INTEGER :: N, NH, NH2, NDACT
      INTEGER :: N1, N2, ND1, ND1M, NH1
      INTEGER :: I, K, L, M, IERR
      REAL :: QSPEC, QSPEC_OUT
      REAL :: SUM, RI, DELR, Q, RL
      REAL :: RHO, RHOOLD, RHONEW, RHOLAST
      REAL :: V, VL, VLAST, TAUL
      REAL :: XND, XND2, XND3, DLOGTAU, DLOGR, DLOGR2

      !COMMON /COMRPAR/ RPAR
      REAL, DIMENSION(1000) :: RHELP, TAU
      !CHARACTER RPAR*80
      REAL :: SPEC_IN_SPACING, SPEC_OUT_SPACING
      INTEGER :: NPAR, NSPEC_IN, NSPEC_OUT
      CHARACTER(80) :: RPAR, ACTPAR
      CHARACTER(40), DIMENSION(20) :: CURPAR

      LOGICAL :: bOldDecode
      LOGICAL, DIMENSION(4) :: bParamFound

      !File and channel handles (=KANAL)
      INTEGER, PARAMETER :: hOUT = 6        !write to wruniqX.out (stdout)
      INTEGER, PARAMETER :: hCPR = 0        !write to wruniqX.cpr (stderr)

      REAL, EXTERNAL :: WRVEL                 !own function for velocity field

 
C***  DECODE THE INPUT OPTION CARD WHICH SPECIFIES THE RADIUS GRID
      !RadiusGridParameters contains up to 3 lines from CARDS:
      ! RGRID, SPECIAL_OUTER_POINTS and SPECIAL_INNER_POINTS
      !RPAR contains the current CARDS line
      RPAR = RadiusGridParameters(1)
      !DECODE (80,7,RPAR ) XND,XND2,XND3,DLOGTAU
      
      bOldDecode = .FALSE.
      CALL SARGC (RPAR, NPAR)
      IF (NPAR < 5) THEN
        WRITE (hCPR,'(A)') '*** RGRID: NOT ENOUGH PARAMETERS'
        STOP '*** FATAL ERROR WHILE DECODING RGRID CARDS-LINE'
      ENDIF

      !New decoding => allows flexible format and modern syntax
      bParamFound = .FALSE.
      DO i=1, NPAR
        CALL SARGV(RPAR,i,CURPAR(I))
      ENDDO
      IF (NPAR > 2) THEN
        DO i=2, NPAR 
         SELECTCASE (CURPAR(i))
          CASE ('ND','NDTOTAL','ND_TOTAL','NDTOT','ND_TOT')
            IF (NPAR >= (i+1)) THEN
              READ (CURPAR(i+1), '(F5.0)', IOSTAT=IERR) XND
              IF (IERR == 0) THEN
                ND=IFIX(XND)
                bParamFound(1) = .TRUE.
              ENDIF      
            ENDIF
          CASE ('ND2','NDDENS','ND_DENS','DENS','NDRHO','ND_RHO','RHO')
            IF (NPAR >= (i+1)) THEN
              READ (CURPAR(i+1), '(F5.0)', IOSTAT=IERR) XND2
              IF (IERR == 0) THEN
                ND2=IFIX(XND2)
                bParamFound(2) = .TRUE.
              ENDIF                  
            ENDIF
          CASE ('ND3','NDVELO','ND_VELO','NDV','VELO','V')
            IF (NPAR >= (i+1)) THEN
              READ (CURPAR(i+1), '(F5.0)', IOSTAT=IERR) XND3
              IF (IERR == 0) THEN
                ND3=IFIX(XND3)
                bParamFound(3) = .TRUE.
              ENDIF                  
            ENDIF
          CASE ('DLOGTAU')
            IF (NPAR >= (i+1)) THEN
              READ (CURPAR(i+1), '(F10.0)', IOSTAT=IERR) DLOGTAU
              IF (IERR /= 0) THEN
                WRITE(hCPR,'(A)') '*** RGRID: CANNOT READ DLOGTAU'
                STOP 'ERROR'
              ELSE
                bParamFound(4) = .TRUE.
              ENDIF                  
            ENDIF
         ENDSELECT
        ENDDO
      ENDIF

      DO i=1, 4 
        !One or more parameters have not been found switch to old decoding
        IF (.NOT. bParamFound(i)) THEN
          WRITE (hCPR,*) '*** RGRID: Old RGRID decoding used'
          bOldDecode = .TRUE.
        ENDIF
      ENDDO

      IF (bOldDecode) THEN
        !Failsafe: old decoding with fixed format
        READ (UNIT=RPAR, FMT=7, ERR=92) XND,XND2,XND3,DLOGTAU
    7   FORMAT (10X,F5.0,4X,F5.0,4X,F5.0,8X,F10.0 )
        ND=IFIX(XND)
        ND2=IFIX(XND2)
        ND3=IFIX(XND3)
      ENDIF
      
C***  Decode input options for special points
      
      !default values for Special Outer and Inner Points
      NSPEC_IN         = 4
      SPEC_IN_SPACING  = 2.
      NSPEC_OUT        = 0
      SPEC_OUT_SPACING = 2.
      
      RPAR = RadiusGridParameters(2)
      CALL SARGC (RPAR, NPAR)
      IF (NPAR >= 2) THEN         
         CALL SARGV (RPAR, 2, ACTPAR)
         READ (ACTPAR, '(I10)', ERR=92) NSPEC_OUT
         IF (NPAR > 2) THEN
           CALL SARGV (RPAR, 3, ACTPAR)
           IF (ACTPAR == 'SPACING') THEN
             CALL SARGV (RPAR, 4, ACTPAR)
             READ (ACTPAR, '(F10.0)', ERR=92) SPEC_OUT_SPACING
           ENDIF
         ENDIF
      ENDIF
      RPAR = RadiusGridParameters(3)
      CALL SARGC (RPAR, NPAR)
      IF (NPAR >= 2) THEN
         CALL SARGV (RPAR, 2, ACTPAR)
         READ (ACTPAR, '(I10)', ERR=92) NSPEC_IN
         IF (NPAR > 2) THEN
           CALL SARGV (RPAR, 3, ACTPAR)
           IF (ACTPAR == 'SPACING') THEN
             CALL SARGV (RPAR, 4, ACTPAR)
             READ (ACTPAR, '(F10.0)', ERR=92) SPEC_IN_SPACING
           ENDIF
         ENDIF
      ENDIF      
      
      
      NSPEC=NSPEC_IN
      QSPEC=1. / SPEC_IN_SPACING
      QSPEC_OUT=1. / SPEC_OUT_SPACING
      
      IF (ND > NDDIM) THEN
          CALL REMARK ('ND .GT. NDIM')
          STOP 'ERROR'
      ENDIF
      ND1=ND-ND2-ND3-NSPEC - NSPEC_OUT
      IF (ND1 < 2) THEN
          CALL REMARK ('ND1 .LT. 2')
          STOP 'ERROR'
      ENDIF

 
C***  ND1 POINTS EQUALLY SPACED IN LOG TAU
 
C*** FIRST, A TABLE OF 1000 ENTRIES: RHELP(I), TAU(I) , IS ESTABLISHED
      N=1000
      NH1=800
      NH2=200
C***  THE RADIUS POINTS RHELP ARE SPACED LOGARITHMICALLY
      DLOGR=ALOG(RMAX)/FLOAT(NH1-1)
      SUM=.0
      RHOOLD=1./RMAX/RMAX/WRVEL(RMAX)
      TAU(1)=.0
      RHELP(1)=RMAX
      !Hilfsfeingitter aus 1000 (800+200) Punkten, feineres Spacing innen
      DO I=2,NH1-1
        RI=EXP((NH1-I)*DLOGR)
        RHELP(I)=RI
        RHONEW=1./RI/RI/WRVEL(RI)
        DELR=RI-RHELP(I-1)
        SUM=SUM+.5*DELR*(RHOOLD+RHONEW)
        RHOOLD=RHONEW
        TAU(I)=SUM
      ENDDO
      DLOGR2=LOG(RI)/FLOAT(NH2+1)
      DO I=NH1,N
        RI=EXP((N-I)*DLOGR2)
        IF (I == N) RI=1.
        RHELP(I)=RI
        RHONEW=1./RI/RI/WRVEL(RI)
        DELR=RI-RHELP(I-1)
        SUM=SUM+.5*DELR*(RHOOLD+RHONEW)
        RHOOLD=RHONEW
        TAU(I)=SUM
      ENDDO
 
      !Gitterpunkte nach TAU-Kriterium vergeben (gleichm. Spacing in DLOGTAU)
      ND1M=ND1-1
      DO L=2,ND1M
        TAUL=TAU(N)*10**(DLOGTAU*(FLOAT(L-2)/FLOAT(ND1-2)-1.))
        CALL LIPO (R(L),TAUL,RHELP,TAU,N)
        INCRIT(L)='TAU     '
      ENDDO
      R(1)=RMAX
      INCRIT(1)='RMAX    '
      R(ND1)=1.
 
C***  INSERT ND2 POINTS BY DENSITY CRITERION
      N2=ND1+ND2-1
      DO K=ND1,N2
        Q=1.
        RHOLAST=1./R(1)/R(1)/WRVEL(R(1))
        DO L=2,K
          RL=R(L)
          RHO=1./RL/RL/WRVEL(RL)
          IF (RHO/RHOLAST >= Q) THEN 
            Q=RHO/RHOLAST
            M=L
          ENDIF
          RHOLAST=RHO
        ENDDO
        CALL SHIFT (R,M,K)
        R(M)=.5*(R(M-1)+R(M+1))
        CALL SHIFTSTRING (INCRIT,M,K)
        INCRIT(M)='DENSITY '
      ENDDO
 
C***  INSERT ND3 POINTS BY VELOCITY CRITERION
      N1=ND1+ND2
      N2=ND1+ND2+ND3-1
      DO K=N1,N2
        Q=.0
        VLAST=WRVEL(R(1))
        DO L=2,K
          VL=WRVEL(R(L))
          IF (VLAST-VL >= Q) THEN
            Q=VLAST-VL
            M=L
          ENDIF
          VLAST=VL
        ENDDO
        CALL SHIFT (R,M,K)
        R(M)=.5*(R(M-1)+R(M+1))
        CALL SHIFTSTRING (INCRIT,M,K)
        INCRIT(M)='VELOCITY'
      ENDDO
 
C***  Insert NSPEC-OUT Points near the Outer Boundary
      DO I=1, NSPEC_OUT
        NDACT = ND - NSPEC_OUT + I - 2
        CALL SHIFT (R, 2, NDACT)
        CALL SHIFTSTRING (INCRIT, 2, NDACT)
C!!!        R(2) = 0.5 * (R(1) + R(3))
        R(2) = QSPEC_OUT * (R(3) - R(1)) + R(1)
        INCRIT(2)='SPECIAL '
      ENDDO
 
C***  INSERT NSPEC POINTS NEAR THE INNER BOUNDARY
      DO L=ND-NSPEC,ND-1
        R(L)=1.+(R(L-1)-1.)* QSPEC
        INCRIT(L)='SPECIAL '
      ENDDO
 
C***  INNER BOUNDARY
      R(ND)=1.
      INCRIT(ND)='INNER-B.'
 
C***  A special smoothing in radius is applied in order to get
C***    equally spaced depth points. This is done on the logarithmic 
C***    radius grid.         wrh+lars 17-Dec-1997 14:39:34
      DO L=2+NSPEC_OUT, ND-NSPEC_IN-1
        RHELP(L) = 0.25 * 
     >    (ALOG10(R(L-1)) + ALOG10(R(L+1)) + 2.*ALOG10(R(L)))
      ENDDO
      DO L=2+NSPEC_OUT, ND-NSPEC_IN-1
        R(L) = 10.**RHELP(L)
      ENDDO

      RETURN

      !Error message in case of CARDS line decoding fail
   92 WRITE (hCPR,*)
     >   'RGRID: ERROR WHILE DECODING THE FOLLOWING CARDS-LINE:'
      WRITE (hCPR,*) RPAR
      STOP 'ERROR'
      
      END
