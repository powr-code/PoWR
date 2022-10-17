      SUBROUTINE PREPK(
C***  Input
     >    XLAMBDA, NF, XLAMBDA2, NF2, VDOP, CLIGHT, DXMAX, 
     >    VINF, ND, CMFBAND, NLINE, MAXIND, XLAMSOR, 
     >    BPLOT, BXJCE, LINE, 
C***  Output
     >    XLAM0, XLAM0LN, ALN, CMFBANDR, 
     >    XKMID, XKMIN, XKRED, XKMAX, XKC, XKC2, XKCMAX, 
     >    KSPACE, BANDP, BANDM, KOUT, K, KLAST, NK, 
     >    KONCHECK, KONCHECK2, KONTACT, KONTAUP, NLACT, 
     >    LINECHECK, ILINECHECK, 
     >    KCL, KCU, KCDMAX, KCDELTA, KCCHECK, BFF_ACT, IFF_N,
     >    IFF_DK, IFF_MAX)
C****************************************************************
C***  Several quantities concerning frequency grid, linewidth
C***    and Continuum interpolation are preset
C***    Called by COLI
C****************************************************************

      INTEGER, PARAMETER :: TINYINT = SELECTED_INT_KIND(2)

      INTEGER, INTENT(IN) :: IFF_MAX

      DIMENSION XLAMBDA(NF), XLAMBDA2(NF2), LINE(NLINE)
      DIMENSION XKMID(NLINE), XKMIN(NLINE), XLAMSOR(NLINE)
      DIMENSION XKRED(NLINE), XKMAX(NLINE)
      DIMENSION XKC(NF), XKC2(NF2)

      LOGICAL BPLOT, BXJCE, BFF_ACT
      INTEGER (KIND=TINYINT), DIMENSION(IFF_MAX) :: IFF_DK

      REAL, INTENT(IN) :: VINF

C***  Preparing harmonic wavelengths
      XLAM0 = XLAMBDA(1)
      XLAM0LN = ALOG(XLAM0)
      ALN = ALOG(1. + VDOP/CLIGHT*DXMAX)

C***  Preparing other quantities
      CMFBANDR = CMFBAND + 2.*VINF

C***  Preparing line indices (XKMIN, XKMAX)
      DO NL=1, NLINE
        XK = (ALOG(XLAMSOR(NL)) - XLAM0LN) / ALN
        XKMID(NL) = XK
        XKMIN(NL) = XK - CMFBAND/DXMAX - 1.
C***   Opacity
        XKRED(NL) = XK + CMFBAND/DXMAX + 1.
C***   Red Line Wing
        XKMAX(NL) = XK + CMFBANDR/DXMAX + 1.
      ENDDO

C***  Prepare indices of the continua
C***  XKC  contains the Kontinua inserted for a correct flux integration
C***    Only two k-points are necessary
C***  XKC2 contains the true continuum points. Fine spacing is necessary 
C***    a wide k-range
      DO KON=1, NF
        XKC(KON) = (ALOG(XLAMBDA(KON)) - XLAM0LN) / ALN
        IF (BPLOT) THEN
          WRITE (38,'(A,F9.2,A)')
     >      'KASDEF SYM ', XKC(KON), ' -2. 0. 0.2 -0.2 8'
        ENDIF
      ENDDO
      DO KON=1, NF2
        XKC2(KON) = (ALOG(XLAMBDA2(KON)) - XLAM0LN) / ALN
        IF (BPLOT) THEN
          WRITE (38,'(A,F9.2,A)')
     >      'KASDEF SYM ', XKC2(KON), ' -1. 0. 0.2 -0.2 8'
        ENDIF
      ENDDO
C!!!      XKCMAX = 1.  + XKC(NF) + CMFBANDR / DXMAX
C!!!      xkcmax = xkc(nf)
C***  !!!!!!!!!!!!!!!!!!!!!!!!!
      XKCMAX = (ALOG(XLAMBDA(NF)) - XLAM0LN) / ALN
C***  !!!!!!!!!!!!!!!!!!!!!!!!!

C***  Spacing when no lines and no continua are active
C***  At the moment 500 Lambda-points per decade are default
      KSPACE = INT((LOG(10.)/ALN) / 500.)

C***  Do not allow for KSPACE > 200
      KSPACE = MIN(200, KSPACE)

C***  Extended spacing at Continuum points which are not due to edges?
      IF (BXJCE) THEN
        BANDP = CMFBANDR/DXMAX
        BANDM = CMFBAND/DXMAX
      ELSE
        BANDP = 1.
        BANDM = 1.
      ENDIF

      WRITE (0,'(A,I3,A,L1)')
     >  'Frequency grid:  KSPACE=', KSPACE, 
     >  '   Extended non-edge continuum points= ', BXJCE
      WRITE (0,*)
                                 
C***  K is frequency-index
      KOUT = -1
      K = 0
      KLAST = -KSPACE
      NK = 0
      KONCHECK  = 0
      KONCHECK2 = 1
      KONTACT  = 0
      KONTAUP  = 1
      NLACT = 0
      LINECHECK = 1
      ILINECHECK = LINE(LINECHECK)

C***  Start Indices for Continuum Interpolation (KCL = KCU)
      KCL = 0
      KCU = 0
C***  Width of the Interval used for Continuum Interpolation
      KCDMAX = 10
      KCDELTA = KCDMAX
C***  Initialize Counter for Continuum-Edge-Test
      KCCHECK = 1

      BFF_ACT = .FALSE.
      IFF_N = 0
      IFF_DK = 0

      RETURN
      END
