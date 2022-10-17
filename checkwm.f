      SUBROUTINE CHECKWM (CWM0, CWM1, CWM2, CWM3, 
     >              CWM1O, CWM1I, CWM3O, CWM3I,
     >              WP1, WP1LAST, W0, R, ND, NP, Z, P)
 
      DIMENSION CWM0(ND,NP), CWM2(ND,NP)
      DIMENSION CWM1(ND,NP), CWM3(ND,NP)
      DIMENSION CWM1O(NP), CWM1I(NP), CWM3O(NP), CWM3I(NP)
      DIMENSION W0(ND), WP1(NP), WP1LAST(NP)
      DIMENSION R(ND), Z(ND,NP), P(NP)
      DIMENSION HTOTA(ND), HTOTB(ND)
      DIMENSION NTOTA(ND), NTOTB(ND)

      JP = 1
      LMAX=MIN0(NP+1-JP,ND)
      CALL GENW0 (JP,LMAX,ND,NP,Z,R,P,W0 )

      WRITE (0,*)
      WRITE (0,*) 'Check other moments'
      WRITE (0,*) 'JP=', JP
      WRITE (0,'(A2,9A11)')
     >  'L', 'RADIUS', 'W0', 'W0_OLD', 'W1', 'W1_old', 
     >  'W1L_old', 'W2', 'W3'
      DO L=1, LMAX
        RL = R(L)
        R2 = RL*RL
        WRITE (0,'(I2,8(1X,E10.4))')
     >    L, RL, CWM0(L,JP), W0(L)*R2, 
     >           CWM1(L,JP), WP1(JP), WP1LAST(JP), 
     >           CWM2(L,JP), 
     >           CWM3(L,JP)
      ENDDO

      L = 39
      R2 = R(L)*R(L)
      JMAX = NP+1-L
      WRITE (0,*)
      WRITE (0,'(A,2I3,1X,E15.8)') 
     >  'L, JMAX, R(L)^2 =', L, JMAX, R2
      XJA = 0.
      XJB = 0.
      XJC = 0.
      XHA = 0.
      XHB = 0.
      XHC = 0.
      XKA = 0.
      XKB = 0.
      XKC = 0.
      XNA = 0.
      XNB = 0.
      XNC = 0.
C***  O. AND 2. MOMENT
      WRITE (0,*) '0. AND 2. MOMENT'
      DO JP=1, JMAX
c      DO JP=2, JMAX-1
        XIA = 1.
        XIB = Z(L,JP)/R(L)
        XIC = XIB*XIB
        XJA = XJA + XIA*CWM0(L,JP)
        XJB = XJB + XIB*CWM0(L,JP)
        XJC = XJC + XIC*CWM0(L,JP)
        XKA = XKA + XIA*CWM2(L,JP)
        XKB = XKB + XIB*CWM2(L,JP)
        XKC = XKC + XIC*CWM2(L,JP)
        write (0,*) 
     >    jp,cwm0(l,jp),cwm2(l,jp)
      ENDDO
      WRITE (0,'(A,3(1X,E15.8))') 
     >  '0. XJA/R^2, XJB/R^2, XJC/R^2 =', XJA/R2, XJB/R2, XJC/R2
      WRITE (0,'(A,3(1X,E15.8))') 
     >  '2. XKA/R^2, XKB/R^2, XKC/R^2 =', XKA/R2, XKB/R2, XKC/R2
      WRITE (0,*)

C***  1. AND 3. MOMENT at outer and inner boundaries
      WRITE (0,*) '1. AND 3. MOMENT at boundaries'
      R2 = R(1)*R(1)
      XHOA = 0.
      XHOB = 0.
      XHIA = 0.
      XHIB = 0.
      XNOA = 0.
      XNOB = 0.
      XNIA = 0.
      XNIB = 0.
      write (0,*) '####jmax=', jmax, np, nd
      DO JP=1, NP-ND
        write (0,*) 'in  ', cwm1i(jp), cwm3i(jp)
      ENDDO
        write (0,*)
      DO JP=1, NP
        write (0,*) 'out ', cwm1o(jp), cwm3o(jp)
      ENDDO
      DO JP=1, NP
        XIA = 1.
        XIB = Z(1,JP)/R(1)
        XHOA = XHOA + XIA*CWM1O(JP)
        XHOB = XHOB + XIB*CWM1O(JP)
        XNOA = XNOA + XIA*CWM3O(JP)
        XNOB = XNOB + XIB*CWM3O(JP)
        IF (JP .GT. NP-ND) CYCLE
        XIA = 1.
        XIB = Z(ND,JP)/R(ND)
        XHIA = XHIA + XIA*CWM1I(JP)
        XHIB = XHIB + XIB*CWM1I(JP)
        XNIA = XNIA + XIA*CWM3I(JP)
        XNIB = XNIB + XIB*CWM3I(JP)
      ENDDO
      WRITE (0,'(A,3(1X,E15.8))') 
     >  '1. XHOA/R^2, XHOB/R^2 =', XHOA/R2, XHOB/R2
      WRITE (0,'(A,3(1X,E15.8))') 
     >  '1. XHIA/R^2, XHIB/R^2 =', XHIA, XHIB
      WRITE (0,'(A,3(1X,E15.8))') 
     >  '3. XNOA/R^2, XNOB/R^2 =', XNOA/R2, XNOB/R2
      WRITE (0,'(A,3(1X,E15.8))') 
     >  '3. XNIA/R^2, XNIB/R^2 =', XNIA, XNIB

C***  1. UND 3. MOMENT
      R1M = 0.5 * (R(L+1) + R(L))
      R2M = R1M*R1M
      WRITE (0,*) '1. MOMENT L=', L
      WRITE (0,'(A3,4(1X,A15))') 
     >  'JP', 'XIB', 'CWM1', 'XHA/R^2', 'XHB/R^2'
      WRITE (0,'(A3,4(1X,A15))') 
     >  'JP', 'XIB', 'CWM3', 'XNA/R^2', 'XNB/R^2'
      DO JP=1, JMAX-1
        XIA = 1.
C***  Linear in Z
        XIB = SQRT( R2M-P(JP)**2. ) / R1M
C***  Linear in P
C        XIB = (R1M - P(JP)) / R1M
        XHA = XHA + XIA*CWM1(L,JP)
        XHB = XHB + XIB*CWM1(L,JP)
        XNA = XNA + XIA*CWM3(L,JP)
        XNB = XNB + XIB*CWM3(L,JP)
c        write (0,'(A,i3,4(1x,e15.8))') 
c     >    '1. ', jp, xib, CWM1(L,JP)/r2m, xha/r2m ,xhb/r2m
c        write (0,'(A,i3,4(1x,e15.8))') 
c     >    '3. ', jp, xib, CWM3(L,JP), xNa/r2m ,xNb/r2m
      ENDDO
      WRITE (0,*) '1. und 3. MOMENT'
      write (0,'(i3,3(1x,e15.8))')
     >  jp, xib, xha/r2m ,xhb/r2m
      WRITE (0,*) 'XHA, XHB=', XHA/R2M, XHB/R2M
      WRITE (0,*) 'XNA, XNB=', XNA/R2M, XNB/R2M

      STOP 'Test in Subr. CHECKWM'
      RETURN

C***  1. and 3. Moment as in COLIRAY
      WRITE (0,*) '1. MOMENT, second try'
      DO L=1, ND
        HTOTA(L) = 0.
        HTOTB(L) = 0.
        NTOTA(L) = 0.
        NTOTB(L) = 0.
      ENDDO
      DO JP=1, NP-1
        LMAX=MIN0(NP+1-JP,ND)
        DO L=1, LMAX-1
          R1M = 0.5 * (R(L+1) + R(L))
          R2M = R1M*R1M
          XIA = 1.
C***  Linear in Z
          XIB = SQRT( R2M-P(JP)**2. ) / R1M
C***  Linear in P
c          XIB = (R1M - P(JP)) / R1M
          HTOTA(L) = HTOTA(L) + XIA*CWM1(L,JP)
          HTOTB(L) = HTOTB(L) + XIB*CWM1(L,JP)
          NTOTA(L) = NTOTA(L) + XIA*CWM3(L,JP)
          NTOTB(L) = NTOTB(L) + XIB*CWM3(L,JP)
c          HTOTA(L) = HTOTA(L) + XIA*WP1(JP)
c          HTOTB(L) = HTOTB(L) + XIB*WP1(JP)
        ENDDO
      ENDDO

      DO L=1, ND
        R1M = 0.5 * (R(L+1) + R(L))
        R2M = R1M*R1M
        WRITE (0,'(A,I3,4(1X,E15.8))') 
     >    'HTOT, NTOT=', L, HTOTA(L)/R2M, HTOTB(L)/R2M, 
     >                NTOTA(L)/R2M, NTOTB(L)/R2M
      ENDDO

      stop 'test in checkmw'

      RETURN
      END








