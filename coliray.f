      SUBROUTINE COLIRAY (K, U, NDDIM, 
     >        Z, RADIUS, ND, NP, JP, P,
     >        TA, TB, TC, UB, GA, H, QQ, S, V, 
     >        VA, VB, VELO, GRADI, PP,
     >        BCORE, DBDR, XLAMK, ENTOT, DELTAX, W0, 
     >        XJNUED, LJPLOT, OPAK, ETAK, WNUE, 
     >        FWEIGHTL, WP1, WP1LAST, 
     >        XJL, XHL, XKL, XNL, 
     >        CWM0, CWM1, CWM2, CWM3, 
     >        XHO, XHI, XNO, XNI, 
     >        CWM1O, CWM1I, CWM3O, CWM3I,
     >        U_OLD, V_OLD, BSTATIC, IPLOT, BPLOT2, IW_COLIRAY_IPLUS, 
     >        IVERS)

C***********************************************************************
C***  LINE RADIATION TRANSFER IN THE COMOVING FRAME FROM A GIVEN SOURCE FUNCTION
C***  FOR A GIVEN IMPACT-PARAMETER, THE INTEGRATION IS CARRIED OUT IN SPACE
C***  AND FREQUENCY AND RESULTS ARE INTEGRATED UP TO YIELD THE MEAN INTENSITY
C***  OF THE LINE, XJLMEAN = J-NUE-BAR
C***  Note, added 26-Jul-1996  9:57:46:
C***  The differencing scheme with respect to frequencies is 
C***  "fully implicit", i.e. all quantities except the frequency derivative 
C***  are taken at the NEW frequency endpoint of the interval. This 
C***  corresponds to MKH '75, but differs from earlier versions and the 
C***  description in the thesis of W.-R.H.
C***********************************************************************

      DIMENSION U(ND), OPA(ND), RADIUS(ND)
      DIMENSION V(ND), Z(ND,2), WNUE(ND)
      DIMENSION VELO(ND), GRADI(ND), PP(ND)
      DIMENSION W0(ND)
      DIMENSION ETAK(NDDIM), OPAK(NDDIM)
      DIMENSION S(ND), QQ(ND), P(ND)
      DIMENSION XJTOTL(ND), XKTOTL(ND), XNTOTL(ND)
      DIMENSION XJL(ND), XHL(ND), XKL(ND), XNL(ND)
      DIMENSION WP1(NP), WP1LAST(NP)
      DIMENSION CWM0(ND,NP), CWM1(ND,NP), CWM2(ND,NP), CWM3(ND,NP)
      DIMENSION CWM1O(NP), CWM1I(NP), CWM3O(NP), CWM3I(NP)
      DIMENSION U_OLD(ND), V_OLD(ND)

      LOGICAL PLOT, BSTATIC, BPLOT2
C!!!      STATIC  XILASTK, XIPLASTK

      COMMON / COMFUN / DELTAV,XMIN

      LMAX=MIN0(NP+1-JP,ND)
      LZ=LMAX-1

      DO L=1,ND
         U(L) = U_OLD(L)
         V(L) = V_OLD(L)
      ENDDO

      IF (IVERS .EQ. 0) THEN
        XIMINUS=0.
      ELSE IF (IVERS .EQ. 4) THEN
        IF (JP .EQ. 1) THEN
C***        TESTPLOT FACILITY
          PLOT = .FALSE.
C          PLOT = IND .EQ. 1
          CALL SFIT (ND, OPAK, ETAK,
     >               ENTOT, XLAMK, SBOUND, PLOT)
          TAUBOUND = OPAK(1) * RADIUS(1)
          XIMINUS = SBOUND * (1. -EXP(-TAUBOUND))
        ENDIF
      ELSE
        WRITE (0,*) 'Boundary version not valid'
        STOP 'ERROR in Subr. COLIRAY'
      ENDIF      

C***  Special treatment of tangent ray (JP=NP): Only XIMINUS
      IF (JP .EQ. NP) THEN 
         U(1) = XIMINUS
         V(1) = .0
         GOTO 100
      ENDIF

C***  Inner Boundary, NOTEMP-Version
      XIPLUS = BCORE + DBDR*Z(ND,JP)/OPAK(ND)
      IF (XIPLUS .LT. 0.) THEN
         XIPLUS = 0.
         IW_COLIRAY_IPLUS = IW_COLIRAY_IPLUS + 1
      ENDIF

cC***  BLUE-WING BOUNDARY CONDITION if U comes from the BACKJCU file:
cC***  THE FEAUTRIER-FLUX V (AT DEPTH INTERSTICES) IS DERIVED FROM U VIA THE DGL.
cC***   (see below for the alternative version CMFUV=.TRUE.)
c      IF (NEWV .AND. .NOT. CMFUV) THEN
c       DO 5 L=1,LZ
c        X=0.5*(OPA(L)+OPA(L+1))
c        V(L)=(U(L+1)-U(L))/X/(Z(L,JP)-Z(L+1,JP))
c    5  CONTINUE
c      ENDIF

C***  PP(L) = VELOCITY GRADIENT, PROJECTED ON THE PRESENT RAY J, /DELTAX
      IF (.NOT. BSTATIC) THEN
         DO L=1,LMAX
            RL=RADIUS(L)
            Y=Z(L,JP)/RL
            YY=Y*Y
            PP(L)=(YY*GRADI(L)+(1.-YY)*VELO(L)/RL)/DELTAX
         ENDDO
      ELSE
         DO L=1,LMAX
            PP(L)=0.
         ENDDO
      ENDIF

C***  GENERATING INTEGRATION WEIGHTS APPROPRIATE FOR THE CONSIDERED RAY
C***  WHICH WILL BE USED TO INTEGRATE THE 0. TO 3. MOMENTS
      CALL GENW0 (JP,LMAX,ND,NP,Z,RADIUS,P,W0 )
 
C***  At the beginning of the frequency grid, U and V are needed as blue-wing
C***  boundary condition.  
C***  U and V are calculated here by solving the equations 
C***  for velocity gradient = 0 (PP)
      IF (K .EQ. 0) THEN 
        DO L=1, LMAX
          QQ(L) = 0.
        ENDDO
        CALL COLISET (Z(1,JP),
     >         ND,LMAX,TA,TB,TC,UB,VA,VB,GA,H,S,OPAK,ETAK,
     >         QQ,BCORE,DBDR,RADIUS,XIMINUS,XIPLUS,JP)
        DO L=1, LMAX
          U(L) = S(L)
        ENDDO
        CALL INVTRI (TA,TB,TC,U,LMAX)
        CALL GMALU (GA,U,V,LMAX)
      ELSE
        CALL COLISET (Z(1,JP),
     >         ND,LMAX,TA,TB,TC,UB,VA,VB,GA,H,S,OPAK,ETAK,
     >         PP,BCORE,DBDR,RADIUS,XIMINUS,XIPLUS,JP)
        CALL VMALV (VA,VB,V,QQ,LMAX)
        CALL VADD (QQ,S,LMAX)
        CALL MDV (UB,U,LMAX)
        CALL VADD (U,QQ,LMAX)
        CALL INVTRI (TA,TB,TC,U,LMAX)
C***  NOW U IS THE FIELD AT THE NEW INDEX K
        CALL MDV (H,V,LZ)
        CALL GMALU (GA,U,S,LMAX)
        CALL VADD (V,S,LZ)
      ENDIF 

C***  Entry for tangent ray
  100 CONTINUE

C***  ADD THE RADIATION FIELD TO XJL and XKL
      DO L=1, LMAX
        XJL(L) = XJL(L) +   U(L) * CWM0(L,JP)
        XKL(L) = XKL(L) +   U(L) * CWM2(L,JP)
      ENDDO

C***  Add the flux-like intensity to XHL and XNL
      DO L=1, LMAX-1
        XHL(L) = XHL(L) +   V(L) * CWM1(L,JP)
        XNL(L) = XNL(L) +   V(L) * CWM3(L,JP)
      ENDDO

      IF (BPLOT2) THEN
        IF (JP. EQ. 1) THEN
          WRITE (91,'(I8,32(E15.8,1X))') K, XLAMK, 
     >          XJL(IPLOT), XKL(IPLOT), XHL(IPLOT), XNL(IPLOT)
        ENDIF
      ENDIF

C***  Add the radiation field to XHO, XHI, XNO and XNI
      IF (LMAX .EQ. ND) THEN
        XHI = XHI + (XIPLUS - U(ND)) * CWM1I(JP)
C***  Special integration of XHI
C        XHI = XHI - U(ND) * CWM1I(JP)
        XNI = XNI + (XIPLUS - U(ND)) * CWM3I(JP)
      ENDIF
      XHO = XHO + (U(1) - XIMINUS) * CWM1O(JP)
      XNO = XNO + (U(1) - XIMINUS) * CWM3O(JP)

C***  ADD THE RADIATION FIELD TO THE EQUIVALENT WIDTH (BLANKETING FACILITY)
      IF (K .EQ. 1  .OR.  K .EQ. NFL) THEN
         WK = 0.5
         ELSE
         WK = 1.
         ENDIF
      DO 6 L=1, LMAX
      WNUE(L) = WNUE(L) + U(L) * W0(L) * WK
    6 CONTINUE

C***  KEEP THE INTENSITY J-NUE AT ONE DEPTH POINT IF IT IS TO BE PLOTTED
      IF (LJPLOT.NE.0) THEN
         IF (LJPLOT .LE. LMAX) 
     >       XJNUED = XJNUED + U(LJPLOT) * W0(LJPLOT)
      ENDIF

      RETURN
      END
