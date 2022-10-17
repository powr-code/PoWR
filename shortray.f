      SUBROUTINE SHORTRAY(K, XIPLUS, XIPLUS_OLD,  
     >        RADIUS, ND, NP, JP, 
     >        S, S_OLD, XIMINUS, XIMINUS_OLD, 
     >        BCORE, DBDR, XHID, XLAMK, ENTOT, DELTAX, 
     >        OPAK, ETAK, OPAK_OLD, 
     >        XJL, XHL, XKL, XNL, 
     >        CWM0, CWM1, CWM2, CWM3, 
     >        XHO, XHI, XNO, XNI, 
     >        Z, PPP, 
     >        IPLOT, BPLOT2, IW_COLIRAY_IPLUS, IVERS, 
     >        XHOM, XHOP, XNOM, XNOP, OPA)

C***********************************************************************
C***  LINE RADIATION TRANSFER IN THE COMOVING FRAME WITH A GIVEN SOURCE 
C***  FUNCTION (= FORMAL SOLUTION)
C***  SHORT-CHARACTERISTIC METHOD
C***  FOR A GIVEN IMPACT-PARAMETER, THE INTEGRATION IS CARRIED OUT IN SPACE
C***  AND FREQUENCY AND RESULTS ARE INTEGRATED UP TO YIELD THE MEAN INTENSITY
C***  OF THE LINE, XJLMEAN = J-NUE-BAR
C***  Output: The Intensities XIPLUS and XIMINUS 
C***********************************************************************

      INCLUDE 'interfacebib.inc'

      DIMENSION RADIUS(ND), Z(ND), PPP(ND)
      DIMENSION ETAK(ND), OPAK(ND), OPAK_OLD(ND), OPA(ND)
      DIMENSION S(ND), S_OLD(ND)
      DIMENSION XJL(ND), XHL(ND), XKL(ND), XNL(ND)
      DIMENSION CWM0(ND,NP), CWM1(ND,NP), CWM2(ND,NP), CWM3(ND,NP)
      DIMENSION XIPLUS(ND), XIPLUS_OLD(ND), XIMINUS(ND), XIMINUS_OLD(ND)

      LOGICAL PLOT, BPLOT2

      SAVE XIMINUS_OUT

      LMAX=MIN0(NP+1-JP,ND)

C***  Outer Boundary
C***  All versions are independent of angle, 
C***     i.e. I-minus must be set only once per frequency
      IF (JP .EQ. 1) THEN

         IF (IVERS .EQ. 0) THEN
            XIMINUS_OUT=0.

         ELSE IF (IVERS .EQ. 1) THEN
            TAUBOUND = OPAK(1) * RADIUS(1)
            XIMINUS_OUT = S(1) * (1. -EXP(-TAUBOUND))

         ELSE IF (IVERS .EQ. 2) THEN
C***  New Version from Goetz
            PLOT = .FALSE.
            SBOUND = ETAK(1)/OPAK(1)
            TAUBOUND = OPAK(1) * RADIUS(1)
            IF (TAUBOUND .LT. 1.E-3) THEN
               XIMINUS_OUT = SBOUND * (TAUBOUND / 3.)
            ELSE
               EXPTB = EXP(-TAUBOUND)
               FAK1 = 2./TAUBOUND
               FAK2 = FAK1/TAUBOUND
               WPROB = (1.-FAK1+FAK2-FAK2*EXPTB)
               IF (WPROB .GT. 0.9999) THEN
                  WPROB = 0.9999
               ENDIF
               XIMINUS_OUT = SBOUND*WPROB
            ENDIF

         ELSE IF (IVERS .EQ. 4) THEN
C***        TESTPLOT FACILITY -- set valid frequency index for activating!
            PLOT = .FALSE.
c            PLOT = K .EQ. 12482
            CALL SFIT (ND, OPAK, ETAK, ENTOT, XLAMK, SBOUND, PLOT)
            TAUBOUND = OPAK(1) * RADIUS(1)
            IF (TAUBOUND .LT. 1.E-3) THEN
               XIMINUS_OUT = SBOUND * TAUBOUND / 3. 
            ELSE
               EXPTB = EXP(-TAUBOUND)
               FAK1 = 2./TAUBOUND
               FAK2 = FAK1/TAUBOUND
               XIMINUS_OUT = SBOUND * (1.-FAK1+FAK2-FAK2*EXPTB)
            ENDIF

         ELSE IF (IVERS .EQ. 5) THEN
C ***    use only continuum opacity for estimation TBOUND
            PLOT = .FALSE.
c            PLOT = K .EQ. 12482
            CALL SFIT (ND, OPAK, ETAK, ENTOT, XLAMK, SBOUND, PLOT)
            TAUBOUND = OPA(1) * RADIUS(1)
            IF (TAUBOUND .LT. 1.E-3) THEN
               XIMINUS_OUT = SBOUND * TAUBOUND / 3.
            ELSE
               EXPTB = EXP(-TAUBOUND)
               FAK1 = 2./TAUBOUND
               FAK2 = FAK1/TAUBOUND
               XIMINUS_OUT = SBOUND * (1.-FAK1+FAK2-FAK2*EXPTB)
            ENDIF

         ELSE
            WRITE (0,*) 'Boundary version not valid'
            STOP 'ERROR in Subr. SHORT'
         ENDIF      
      ENDIF

      XIMINUS(1) = XIMINUS_OUT

C***  Inward Integration of XIMINUS
      DO L=2, LMAX
        DZ = Z(L-1) - Z(L)
        PPDZ = PPP(L) * DZ / DELTAX
        IF (PPDZ .GT. 1.) THEN
          Q = 1. / PPDZ
          P = 1. - Q
          Z_HAT       = P*Z(L) + Q*Z(L-1)
          CALL SPLINPO_FAST 
     >         (XIMINUS_HAT, Z_HAT, XIMINUS_OLD, Z, LMAX, L, .FALSE.)
          S_HAT      = P*S_OLD(L)      + Q*S_OLD(L-1)
ccc          CALL SPLINPO_FAST_SAME_X 
ccc     >         (S_HAT, Z_HAT, S_OLD, Z, LMAX, L, .FALSE.)

          OPAK_HAT   = P*OPAK_OLD(L)   + Q*OPAK_OLD(L-1)
          TAU1       = 0.5*(OPAK(L) + OPAK_HAT) * Q * DZ
        ELSE          
          Q = PPDZ
          P = 1. - Q
          XIMINUS_HAT = P * XIMINUS(L-1) + Q * XIMINUS_OLD(L-1)
          S_HAT       = P * S      (L-1) + Q * S_OLD      (L-1)
          OPAK_HAT    = P * OPAK   (L-1) + Q * OPAK_OLD   (L-1)
          TAU1        = 0.5 * (OPAK(L) + OPAK_HAT) * DZ
        ENDIF
        EXPTAU1 = EXP(-TAU1)
        IF (ABS(TAU1) .GT. 1.0E-8) THEN
           EDIVTAU = (EXPTAU1-1.) / TAU1
        ELSE 
           EDIVTAU = - 1.0 + TAU1/2.0
        ENDIF

        W0 = 1. + EDIVTAU
        W1 = -EXPTAU1 - EDIVTAU
        XIMINUS(L) = W0*S(L) + W1*S_HAT + EXPTAU1*XIMINUS_HAT
      ENDDO

      IF (LMAX .LT. ND) THEN
C***  Non-Core-Rays
        XIPLUS(LMAX) = XIMINUS(LMAX)
      ELSE
C***  Core-Rays
C        XIPLUS(LMAX) = BCORE + DBDR*Z(ND)/OPAK(ND) 
        XIPLUS(LMAX) = BCORE + 3.*XHID*Z(ND) 

        IF (XIPLUS(LMAX) .LT. 0.) THEN
           XIPLUS_IN = 0.
           IW_COLIRAY_IPLUS = IW_COLIRAY_IPLUS + 1
        ENDIF
      ENDIF

C***  Outward Integration of XIPLUS
      DO L=LMAX-1, 1, -1
        DZ = Z(L) - Z(L+1)
        PPDZ = PPP(L) * DZ / DELTAX
        IF (PPDZ .GT. 1.) THEN
          Q = 1. / PPDZ
          P = 1. - Q
          Z_HAT       = P*Z(L) + Q*Z(L+1)
          CALL SPLINPO_FAST
     >        (XIPLUS_HAT, Z_HAT, XIPLUS_OLD, Z, LMAX, L+1, .FALSE.)
          S_HAT      = P*S_OLD(L)      + Q*S_OLD(L+1)
ccc          CALL SPLINPO_FAST_SAME_X
ccc     >        (S_HAT, Z_HAT, S_OLD, Z, LMAX, L+1, .FALSE.)
          OPAK_HAT   = P*OPAK_OLD(L)   + Q*OPAK_OLD(L+1)
          TAU1       = 0.5*(OPAK(L) + OPAK_HAT) * Q * DZ
        ELSE          
          Q = PPDZ
          P = 1. - Q
          XIPLUS_HAT = P*XIPLUS(L+1) + Q*XIPLUS_OLD(L+1)
          S_HAT      = P*S(L+1)      + Q*S_OLD(L+1)
          OPAK_HAT   = P*OPAK(L+1)   + Q*OPAK_OLD(L+1)
          TAU1       = 0.5*(OPAK(L) + OPAK_HAT) * DZ
        ENDIF
        EXPTAU1 = EXP(-TAU1)
        IF (ABS(TAU1) .GT. 1.0E-8) THEN
           EDIVTAU = (EXPTAU1-1.) / TAU1
        ELSE 
           EDIVTAU = - 1.0 + TAU1/2.0
        ENDIF
        W0 = 1. + EDIVTAU
        W1 = -EXPTAU1 - EDIVTAU
        XIPLUS(L) = W0*S(L) + W1*S_HAT + EXPTAU1*XIPLUS_HAT
      ENDDO



C***  ADD THE RADIATION FIELD TO XJL, XHL, XKL and XNL
      DO L=1, LMAX
        XJL(L) = XJL(L) + 0.5*(XIPLUS(L)+XIMINUS(L)) * CWM0(L,JP)
        XKL(L) = XKL(L) + 0.5*(XIPLUS(L)+XIMINUS(L)) * CWM2(L,JP)
        XHL(L) = XHL(L) + 0.5*(XIPLUS(L)-XIMINUS(L)) * CWM1(L,JP)
        XNL(L) = XNL(L) + 0.5*(XIPLUS(L)-XIMINUS(L)) * CWM3(L,JP)
      ENDDO

      IF (BPLOT2) THEN
        IF (JP. EQ. 1) THEN
          WRITE (91,'(I8,32(E15.8,1X))') K, XLAMK, 
     >          XJL(IPLOT), XKL(IPLOT), XHL(IPLOT), XNL(IPLOT)
        ENDIF
      ENDIF

C***  Add the radiation field to XHO, XHI, XNO and XNI
      IF (LMAX .EQ. ND) THEN
C***    Since 09 Feb 2016: special H (H_spec) for inner boundary
C***    Note the unusual definition of an intensity-like core with
C***    a flux-like weight. This allows to eliminate the backwards
C***    Intensity I- at the inner boundary in the moment equations.
        XHI = XHI + 0.5*(XIPLUS(ND)+XIMINUS(ND)) * CWM1(ND,JP)  
C***    Currently XNI = XNL(ND) is unused
        XNI = XNI + 0.5*(XIPLUS(ND)-XIMINUS(ND)) * CWM3(ND,JP)
      ENDIF
      XHO  = XHO  + 0.5*(XIPLUS(1)-XIMINUS(1)) * CWM1(1,JP)
      XNO  = XNO  + 0.5*(XIPLUS(1)-XIMINUS(1)) * CWM3(1,JP)
      XHOM = XHOM + 0.5*(         -XIMINUS(1)) * CWM1(1,JP)
      XHOP = XHOP + 0.5*(XIPLUS(1)           ) * CWM1(1,JP)
      XNOM = XNOM + 0.5*(         -XIMINUS(1)) * CWM3(1,JP)
      XNOP = XNOP + 0.5*(XIPLUS(1)           ) * CWM3(1,JP)

      RETURN
      END
