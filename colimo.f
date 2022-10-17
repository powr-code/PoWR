      SUBROUTINE COLIMO(K, ND, RADIUS, OPAK, ETAK, ETAKNOTH, 
     >             OPAKNOTH, 
     >             S, XJLMO, XJLMOR2, XHLMO, 
     >             XJLMO_OLD, XHLMO_OLD, 
     >             DLF, DLH, GLF, GLH, VLF, VLH, 
     >             GLF2, GLH2, VLF2, VLH2, 
     >             QLF, QLH, OPAKH, 
     >             EDDIF, EDDIFO, EDDIG, EDDIGO,
     >             EDDIHOUT, EDDIHIN, EDDIHOUTO,
     >             EDDINOUT, EDDININ, EDDINOUTO,
     >             ALH, BLH, CLH, 
     >             A, B, C, W, DX, 
     >             BCORE, DBDR, XIMINUS, BPDONE, XLAMK,
     >             DJDSMO, DJDOMO, 
     >             FULFIL0, FULFIL1, BPLOT, BPLOT2, IPLOT, 
     >             IW_COLIMO_F, IW_COLIMO_G, IW_COLIMO_G2, BSTATIC, 
     >             CLMOETA, CLMOOPA, XHI, 
     >             RADIUS2, EPSG, GEPSB, GEPSBO, 
     >             XHOM, XHOMO, XNOM, XNOMO, 
     >             EDDIHOUTP, EDDINOUTP, EDDIHOUTOP, EDDINOUTOP,
     >             IWARNJMNEG)


C********************************************************
C***  Setup and solution of the tridiagonal equation
C***  -A, B, -C and W are calculated
C********************************************************

      DIMENSION RADIUS(ND), OPAK(ND), ETAK(ND), ETAKNOTH(ND), S(ND)
      DIMENSION OPAKNOTH(ND)
      DIMENSION XJLMO(ND), XJLMOR2(ND), XHLMO(ND)
      DIMENSION XJLMO_OLD(ND), XHLMO_OLD(ND)
      DIMENSION DLF(ND), GLF(ND), VLF(ND)
      DIMENSION DLH(ND-1), GLH(ND-1), VLH(ND-1)
      DIMENSION GLF2(ND), VLF2(ND)
      DIMENSION GLH2(ND-1), VLH2(ND-1)
      DIMENSION QLF(ND), QLH(ND-1), OPAKH(ND)
      DIMENSION EDDIF(ND), EDDIFO(ND)
      DIMENSION EDDIG(ND-1), EDDIGO(ND-1)
      DIMENSION ALH(ND-1), BLH(ND-1), CLH(ND-1)
      DIMENSION A(ND), B(ND), C(ND), W(ND)
      DIMENSION D(ND), E(ND), DINV(ND), DJDSMO(ND), DJDOMO(ND)
      DIMENSION FULFIL0(ND), FULFIL1(ND)
      DIMENSION ETAT(ND), OPAT(ND), CLMOETA(ND), CLMOOPA(ND)
      DIMENSION RADIUS2(ND)
      DIMENSION EPSG(ND), GEPSB(ND), GEPSBO(ND)
      LOGICAL BPDONE, BPLOT, BPLOT2, BSTATIC, BOPAFULL
      CHARACTER*4 REDISMODE

      INTEGER :: IWARNJMNEG
      
C***  DNUEINV is frequency step in Hertz * Doppler velocity in cm/s
      IF (K .EQ. 0 .OR. BSTATIC) THEN
        DNUEINV = 0.
      ELSE
        DNUEINV = 1. / DX
      ENDIF

C***  Calculate CLMOETA nad CLMOOPA
      DO L=1, ND
        CLMOETA(L) = ETAKNOTH(L) 
        CLMOOPA(L) = OPAKNOTH(L) 
      ENDDO

      DO L=1, ND
        GLF(L) = DNUEINV * GLF2(L)
        VLF(L) = DNUEINV * VLF2(L)
        IF (EDDIF(L) .GT. 1.) THEN
          EDDIF(L) = 1.
ccc No reset of EDDIF at last K
ccc          EDDIFO(L) = 1.
          IW_COLIMO_F = IW_COLIMO_F + 1
        ENDIF

        IF (EDDIF(L) .LT. 0.01) THEN
          EDDIF(L) = 0.01
ccc No reset of EDDIF at last K
ccc          EDDIFO(L) = 0.01
          IW_COLIMO_F = IW_COLIMO_F + 1
        ENDIF

        IF (L .EQ. ND) CYCLE
        GLH(L) = DNUEINV * GLH2(L)
        VLH(L) = DNUEINV * VLH2(L)
cc eddig-reset testweise stillgelegt zugunsten von EDDIMIX
        IF (EDDIG(L)*GLH(L)+VLH(L) .LT. -1.E-15) THEN
c          write (0,'(a,i7,1x,f8.3,i3,2(1x,f7.3))') 
c     >      'k, xlamk, l', k, xlamk, l, 
c     >      EDDIG(L), EDDIG(L)*GLH(L)+VLH(L)
          EDDIG(L) = -VLH(L)/GLH(L)
ccc No reset of EDDIG at last K
ccc          EDDIGO(L) = EDDIG(L)
          IW_COLIMO_G = IW_COLIMO_G + 1
          IF (XLAMK .GT. 100. .AND. XLAMK .LT. 2000) THEN
            IW_COLIMO_G2 = IW_COLIMO_G2 + 1
          ENDIF
        ENDIF
      ENDDO

C***  Sphericity-factor (taken from ELIMIN)
      QLF(ND) = 1.
      RRQ = 1.
      RL = RADIUS(ND)
      FL = 3. - 1./EDDIF(ND)
      DO L=ND-1, 1, -1
        RLP = RL
        RL = RADIUS(L)
        FLP = FL

        FL = 3. - 1./EDDIF(L)
ccc   test mit starker Beschraenkung von eddi-f
ccc        if (fl .lt. 0) fl = 0.

        RRQ = RRQ * EXP(FL-FLP) * (RL/RLP)**((FLP*RL-FL*RLP)/(RL-RLP))
        QLF(L) = RRQ / (RL*RL)
      ENDDO

      DO L=1, ND-1
        QLH(L) = 0.5 * (QLF(L)+QLF(L+1))
        OPAKH(L) = 0.5 * (OPAK(L) + OPAK(L+1))
        XNENN = 1. / (OPAKH(L) + GLH(L)*EDDIG(L) + VLH(L))
        ALH(L) = XNENN / DLH(L)
        BLH(L) = GLH(L) * XNENN
        CLH(L) = VLH(L) * XNENN
        GEPSB  (L) = 0.5 * EDDIG (L  ) * EPSG(L  ) * BLH(L  )
        GEPSBO (L) = 0.5 * EDDIGO(L  ) * EPSG(L  ) * BLH(L  )
      ENDDO

      DO L=2, ND-1
C***  Upper Diag. (-C)
        C(L) = (QLF(L+1)*EDDIF(L+1)*ALH(L) / QLH(L) - GEPSB(L))
     >    / DLF(L)
C***  Diag. (B)
        B(L) = (QLF(L)*EDDIF(L)* 
     >         (ALH(L)/QLH(L) + ALH(L-1)/QLH(L-1)) + 
     >         GEPSB(L) - GEPSB(L-1))/DLF(L)  
     >         + GLF(L)*EDDIF(L) + VLF(L)  
     >         + CLMOOPA(L)
C***  Lower Diag. (-A)
        A(L) = (QLF(L-1)*EDDIF(L-1)*ALH(L-1) / QLH(L-1) + GEPSB(L-1)) 
     >    / DLF(L)
C***  Right-Hand-Side
        W(L) = ((BLH(L)*EDDIGO(L)+CLH(L))*XHLMO_OLD(L)-
     >          (BLH(L-1)*EDDIGO(L-1)+CLH(L-1))*XHLMO_OLD(L-1) 
     >          + GEPSBO(L)   * (XJLMO_OLD(L+1) + XJLMO_OLD(L))
     >          - GEPSBO(L-1) * (XJLMO_OLD(L-1) + XJLMO_OLD(L))
     >          ) / DLF(L) 
     >         + (GLF(L)*EDDIFO(L)+VLF(L))*XJLMO_OLD(L) 
     >         + CLMOETA(L)*RADIUS2(L)
      ENDDO


cC***  Outer Boundary (Goetz' second-order version)
c        C(1)  = -2.* QLF(2)*EDDIF(2)*ALH(1) / (QLH(1)*DLF(1))
c        B(1)  = -2.* QLF(1)*EDDIF(1)/DLF(1)*(ALH(1)/QLH(1)) - 
c     >           GLF(1)*EDDIF(1) - VLF(1) -
c     >           2./DLF(1)*EDDIHOUT -
c     >           CLMOOPA(1)
c        W(1)  = -2./DLF(1)*
c     >           (BLH(1)*EDDIGO(1)+CLH(1))*XHLMO_OLD(1) -
c     >           (GLF(1)*EDDIFO(1)+VLF(1))*XJLMO_OLD(1) -
c     >           CLMOETA(1)*RADIUS2(1)

C***  Outer Boundary (Goetz' second-order version), 
C***    new treatment of outer boundary
        C(1)  = 2./DLF(1)* 
     >             (QLF(2)*EDDIF(2)*ALH(1)/QLH(1) - GEPSB(1))
        B(1)  = 2./DLF(1) * 
     >             (QLF(1)*EDDIF(1)*ALH(1)/QLH(1) + EDDIHOUTP
     >              + GEPSB(1) ) 
     >          + GLF(1)*EDDIF(1) + VLF(1) + CLMOOPA(1)
        W(1)  = 2./DLF(1)*
     >              ((BLH(1)*EDDIGO(1)+CLH(1))*XHLMO_OLD(1) 
     >               + GEPSBO(1) * (XJLMO_OLD(1) + XJLMO_OLD(2)) 
     >               - XHOM )
     >           + (GLF(1)*EDDIFO(1)+VLF(1))*XJLMO_OLD(1) 
     >           + CLMOETA(1)*RADIUS2(1) 



C***  Original version (1. Order) 
c      C(1)  = QLF(2)*EDDIF(2) / (QLH(1)*DLF(1))
c      B(1)  = QLF(1)*EDDIF(1) / (QLH(1)*DLF(1)) + 
c     >          GLF(1)*EDDINOUT + (VLF(1)+OPAK(1))*EDDIHOUT
c      W(1)  = (GLF(1)*EDDINOUTO + VLF(1)*EDDIHOUTO) * XJLMO_OLD(1)

C***  Original version (1. Order) 
C***    new treatment of outer boundary
c      C(1)  = QLF(2)*EDDIF(2) / (QLH(1)*DLF(1))
c      B(1)  = QLF(1)*EDDIF(1) / (QLH(1)*DLF(1)) + 
c     >          GLF(1)*EDDINOUTP + (VLF(1)+OPAK(1))*EDDIHOUTP
c       W(1) = -GLF(1)*XNOM - (VLF(1)+OPAK(1))*XHOM
c     >        +(GLF(1)*EDDINOUTOP + VLF(1)*EDDIHOUTOP) * XJLMO_OLD(1)
c     >        -GLF(1)*XNOMO - (VLF(1)+OPAK(1))*XHOMO

C***  Inner Boundary 2. Order
c        B(ND) = -2.* QLF(ND)*EDDIF(ND)/DLF(ND)*(ALH(ND-1)/QLH(ND-1)) -
c     >           GLF(ND)*EDDIF(ND) - VLF(ND) +
c     >           2./DLF(ND)*EDDIHIN -
c     >           CLMOOPA(ND)
c        A(ND) = -2.*QLF(ND-1)*EDDIF(ND-1)*ALH(ND-1)/(QLH(ND-1)*DLF(ND))
c        W(ND) = +2./DLF(ND)*
c     >           (BLH(ND-1)*EDDIGO(ND-1)+CLH(ND-1))*XHLMO_OLD(ND-1) -
c     >           (GLF(ND)*EDDIFO(ND)+VLF(ND))*XJLMO_OLD(ND) -
c     >           CLMOETA(ND)*RADIUS2(ND)

C***  XHI vorgegeben, 2.Ordnung
c        B(ND) = -2.* QLF(ND)*EDDIF(ND)/DLF(ND)*(ALH(ND-1)/QLH(ND-1)) -
c     >           GLF(ND)*EDDIF(ND) - VLF(ND) -
c     >           CLMOOPA(ND)
c        A(ND) = -2.*QLF(ND-1)*EDDIF(ND-1)*ALH(ND-1)/(QLH(ND-1)*DLF(ND))
c        W(ND) = +2./DLF(ND)*
c     >           (BLH(ND-1)*EDDIGO(ND-1)+CLH(ND-1))*XHLMO_OLD(ND-1) -
c     >           (GLF(ND)*EDDIFO(ND)+VLF(ND))*XJLMO_OLD(ND) -
c     >           CLMOETA(ND)*RADIUS2(ND) - 
c     >           2./DLF(ND)*XHI

C***  XHI vorgegeben, 2.Ordnung, neu mit EDDI-Mix
        B(ND) = 2./DLF(ND) * 
     >           (QLF(ND)*EDDIF(ND)*ALH(ND-1)/QLH(ND-1) - GEPSB(ND-1))
     >          + GLF(ND)*EDDIF(ND) + VLF(ND) 
     >          + 2./DLF(ND)*EDDIHIN       !special Eddington factor h_nu from H_spec / J_ray(ND)
     >          + CLMOOPA(ND)              !2nd order additional term (see  Hubeny & Mihalas, P. 393)
        A(ND) = 2./DLF(ND) * 
     >           ( QLF(ND-1)*EDDIF(ND-1)*ALH(ND-1)/QLH(ND-1) 
     >                                                  + GEPSB(ND-1))
        W(ND) = -2./DLF(ND)*
     >           ((BLH(ND-1)*EDDIGO(ND-1)+CLH(ND-1))*XHLMO_OLD(ND-1) 
     >            + GEPSBO(ND-1) * (XJLMO_OLD(ND) + XJLMO_OLD(ND-1)))  
     >          + (GLF(ND)*EDDIFO(ND)+VLF(ND))*XJLMO_OLD(ND) 
     >          + CLMOETA(ND)*RADIUS2(ND)      !2nd order additional term (see  Hubeny & Mihalas, P. 393)
     >          + 2./DLF(ND)*(XHI + BCORE/2.)  !diffusion approximation plus LTE-correction (together w/ EDDIHIN)

C***  XHI vorgegeben, 1. Ordnung
C        B(ND) = QLF(ND)*EDDIF(ND) / (QLH(ND-1)*DLF(ND))
C        A(ND) = QLF(ND-1)*EDDIF(ND-1) / (QLH(ND-1)*DLF(ND))
C        W(ND) = OPAK(ND)*XHI
C***  Alte Version
c        B(ND) = QLF(ND)*EDDIF(ND) / (QLH(ND-1)*DLF(ND)) -
c     >          OPAK(ND)*EDDIHIN
c        A(ND) = QLF(ND-1)*EDDIF(ND-1) / (QLH(ND-1)*DLF(ND))
c        W(ND) = 0.


        IF (BPLOT2) THEN
          WRITE (90,'(I8,32(E15.8,1X))') K, XLAMK, 
     >      DLF(IPLOT), DLH(IPLOT), 
     >      GLF2(IPLOT), VLF2(IPLOT), 
     >      GLH2(IPLOT), VLH2(IPLOT), DNUEINV, QLF(IPLOT), QLH(IPLOT), 
     >      OPAK(IPLOT), OPAKH(IPLOT), XNENN, 
     >      EDDIF(IPLOT), EDDIG(IPLOT), 
     >      ALH(IPLOT), BLH(IPLOT), CLH(IPLOT), 
     >      A(IPLOT), B(IPLOT), C(IPLOT), W(IPLOT), 
     >      B(1), C(1), W(1), 
     >      A(ND), B(ND), W(ND), 
     >      GLF(IPLOT), VLF(IPLOT), GLH(IPLOT), VLH(IPLOT)
        ENDIF

      CALL INVTRI (A, B, C, W, ND)

C**********************************************************************
C***  Calculate Diagonal Elements(^-1) of the Inverse of the T-Matrix
C**********************************************************************
      DO L=1, ND
         A(L) = A(L)/CLMOOPA(L)
         B(L) = B(L)/CLMOOPA(L)
         C(L) = C(L)/CLMOOPA(L)
      ENDDO

C***  B is set negativ due to the form of the matrix -A, B, -C
C***    and Matrix is divided by -OPA
      D(1) = C(1) / (-B(1))
      DO L=2, ND
         D(L) = C(L) / (-B(L) - A(L)*D(L-1))
      ENDDO

      E(ND) = A(ND) / (-B(ND))
      DO L=ND-1, 1, -1
         E(L) = A(L) / (-B(L) - C(L)*E(L+1))
      ENDDO

      DINV(1) = (1. - D(1)*E(2)) * (-B(1))
      DO L=2, ND-1
ccc         DINV(L) = (1. - D(L)*E(L+1)) * (-B(L) - A(L)*D(L-1))
         DINV(L) = (1. - D(L)*E(L+1)) * (B(L) - A(L)*D(L-1))
      ENDDO
ccc      DINV(ND) = -B(ND) - A(ND)*D(ND-1)
      DINV(ND) = B(ND) - A(ND)*D(ND-1)

      DO L=1, ND
         DJDSMO(L) = 1.  / DINV(L)
         DJDOMO(L) = (GLF(L)*EDDIFO(L)+VLF(L)) / (DINV(L)*CLMOOPA(L))
      ENDDO
      DJDSMO(1)    = 0.
      DJDOMO(1)    = 0.
      DJDSMO(ND)   = 0.
      DJDOMO(ND)   = 0.
C**********************************************************************

      DO L=1, ND
C***    WARNING counter increase if negative J encountered    
        IF (W(L) < 0.) IWARNJMNEG = IWARNJMNEG + 1

C***    Set J to zero, if small or neg. values
C***    as recommended by Andreas - wrh  5-Mar-2019
        IF (W(L) < EXP(-499.)) THEN
          XJLMO(L) = 0.
          XJLMOR2(L) = 0.
        ELSE
          XJLMO(L) = W(L)
          XJLMOR2(L) = W(L) / (RADIUS(L)*RADIUS(L))
        ENDIF
      ENDDO


C***  Calculation of the flux from the moment equation (XHLMO)
C***  again suppressing neg. values as recommended by Andreas 
C***  - wrh 5-Mar-2019
      DO L=1, ND-1
        IF (W(L) < EXP(-499.)) THEN
          XHLMO(L) = 0.
        ELSE
          XHLMO(L) = ALH(L)/QLH(L)*
     >                (QLF(L+1)*EDDIF(L+1)*XJLMO(L+1)-
     >                 QLF(L)*EDDIF(L)*XJLMO(L)) + 
     >               (BLH(L)*EDDIGO(L)+CLH(L))*XHLMO_OLD(L)
     >            + GEPSBO(L) * (XJLMO_OLD(L+1) + XJLMO_OLD(L))
     >            - GEPSB (L) * (XJLMO    (L+1) + XJLMO    (L))
        ENDIF
      ENDDO

C***  Check Both Moment equations
      IF (BPLOT) THEN
        DO L=2, ND-1
          FULFIL0(L) = ( 1./DLF(L) * (XHLMO(L)-XHLMO(L-1)) + 
     >      GLF(L) * (-EDDIF(L)*XJLMO(L) + EDDIFO(L)*XJLMO_OLD(L)) + 
     >      VLF(L) * (-XJLMO(L) + XJLMO_OLD(L)) -
     >      OPAKNOTH(L)*XJLMO(L) + ETAKNOTH(L) ) / 
     >      1.
c     >      XJLMO(L) 
          FULFIL1(L) = ( 1./DLH(L) * 
     >      (QLF(L+1)*EDDIF(L+1)*XJLMO(L+1) - QLF(L)*EDDIF(L)*XJLMO(L))
     >      + GLH(L) * (EDDIGO(L)*XHLMO_OLD(L)) 
     >      + VLH(L) * (XHLMO_OLD(L)) 
     >      - XHLMO(L) * (OPAKH(L) + GLH(L)*EDDIG(L) + VLH(L)) ) / 
     >      1.
c     >      XHLMO(L) 
        ENDDO
      ENDIF

      RETURN
      END
