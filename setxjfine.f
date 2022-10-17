      SUBROUTINE SETXJFINE (SFINE_OLD, SFINE_NEW, MAXFINE, ITNEL, L, 
     >                 ND, NDDIM, NDIM, N, VDOP, VELO1, 
     >                 ENTOTL, EN, RSTAR, TL, RNEL, NCHARG, 
     >                 WEIGHT, ELEVEL, EION, EINST, NATOM, 
     >                 KONTHLP, MAXIND, BXJLAPPCORE,
C***    CONTINUA
     >                 XLAMBDA, XLAMBDA2, NF, NF2, 
     >                 XKC, XKC2, ALPHA, 
     >                 SEXPO, ADDCON1, ADDCON2, ADDCON3, IGAUNT, 
     >                 SIGMA1I, KONTLOW, KONTNUP, LASTKON,
     >                 FWTEST, 
C***    LINES
     >                 NLINE, LINE, INDLOW, INDNUP, XLAMSOR, 
     >                 XLAMMIN, XLAMMAX, 
     >                 NUP, LOW, BLASERL, OPAL, ETAL, 
     >                 LIND, LINDS, MAXLIN, 
     >                 XKMIN, XKMAX, XKMID, XKRED, 
     >                 XRED, XBLUE,
     >                 XKRED_CORE, XKBLUE_CORE, 
     >                 
     >                 
C***    IRON
     >                 INDRB, IFRBSTA, IFRBEND, LASTFE, 
     >                 VDOPFE, DXFE, XLAM0FE,
     >                 INDFEACT, MAXFEACT, BFECHECK, BFEWING,
     >                 DFEINDR, SIGMAFE, OPAFE, ETAFE, IFENUP, IFELOW,
     >                 INDEXMAX, BFEMODEL, BNUEFE,
C***    Approximate radiation fields
     >                 LASTIND, XJLAPP, XJL, 
     >                 BXJLAPPNEW, BXJCAPPNEW, BNEWOPER, 
     >                 XJL_PLOTDATA, XJC_PLOTDATA_I, XJC_PLOTDATA_L, 
     >                 IPLOT_XJLAPP, IPLOT_XJCAPP, 
     >                 LPLOT_XJCAPP, NITER_PLOT_JAPP, 
     >                 PWEIGHT, WS, XJC, XJCAPP, XJCAPPNEW, WCHARM, 
     >                 IWARN_NEG_XJCAPP, IWARN_NEG_XJLAPP, 
     >                 XLAMSTART, XLAMEND, 
     >                 OPAC, ETAC, 
C***  New Fine-spaced WCHARM handling
     >                 IFF_MAX, IFF_MAX_MS, FF_INFO, IFF_DK, 
     >                 IFF_WCHARM, WCHARM_FINE, IFF_N_MS, GAMMAC)

C***********************************************************************
C***  Calculates the approximate radiation fields XJCAPP and XJLAPP
C***  accounting for all Opacities and Emissivities
C***  
C***  The same Loop over all frequencies as in COLI is established and
C***  SFINE_OLD and SFINE_NEW, respectively, are calculated 
C***  
C***  The approximate radiation field are then calculated by
C***  Lines :    XJLAPP = XJLfs + INT(CORE) [Snew - Sold]
C***  Continua : XJCAPP = XJCfs + INT(all)  [WCHARM * (Snew - Sold)]
C***  
C***  This Routine replaces SETXJC and SETXJL called by COMA
C***********************************************************************

      INTEGER, PARAMETER :: TINYINT = SELECTED_INT_KIND(2)

      DIMENSION SFINE_OLD(MAXFINE), SFINE_NEW(MAXFINE)
      DIMENSION XLAMBDA(NF), XLAMBDA2(NF2), XKC(NF), XKC2(NF)
      DIMENSION EINST(NDIM,NDIM)
      CHARACTER CMODE*1, NAME*8

C***  Lines: 
      DIMENSION LINE(MAXIND), XLAMSOR(MAXIND)
      DIMENSION XLAMMIN(MAXIND), XLAMMAX(MAXIND)
      DIMENSION XKMIN(MAXIND), XKMAX(MAXIND), XKMID(MAXIND)
      DIMENSION XKRED(MAXIND)
      DIMENSION XRED(MAXIND), XBLUE(MAXIND)
      DIMENSION NUP(MAXIND), LOW(MAXIND), XLAM(MAXIND)
      DIMENSION LIND(MAXLIN), LINDS(MAXLIN)
      DIMENSION WEIGHT(NDIM), OPAL(MAXIND), ETAL(MAXIND)
      DIMENSION XKRED_CORE(MAXLIN), XKBLUE_CORE(MAXLIN)
      DIMENSION XJLAPP(LASTIND), XJL(ND,2)
      DIMENSION PWEIGHT(MAXLIN), WS(MAXLIN)

C***  BNEWOPER switches to new Fine Operator: Default is false
      LOGICAL BXJLAPPNEW, BXJCAPPNEW, BNEWOPER, BXJLAPPCORE
      LOGICAL BKONCHECK, BKONCHECK2
      LOGICAL BLASERL(MAXIND)
      LOGICAL BPLOT, BFIRSTITER

C***  Variables to handle Fine Integration in Steal
C***  FF_INFO : Information on the Frequency grid
C***         (1) : XLAM0
C***         (2) : ALN
C***         (3) : XLAM_FINE_START
C***         (4) : XLAM_FINE_END
C***               (3) + (4) as demanded by the Input Cards
C***         (5) : Start index of Integration
C***         (6) : End index of Integration
C***         (7) : IFF_N Total number of Fine Frequencies
C***         (8) : XLAM_FINE_START
C***         (9) : XLAM_FINE_END
C***               (8) + (9) wavelength of first and last fine point
C***        (10) : KSPACE
C***  BFF_ACT : True, if actual K is in Fine Interval
C***  IFF_DK stores in a compressed way the Delta-K values
C***    IFF_DK can have values between -100 and 100
      REAL, DIMENSION(10) :: FF_INFO
C      DIMENSION IFF_DK(IFF_MAX), 
C      DIMENSION IFF_WCHARM(IFF_MAX)
      INTEGER (KIND=TINYINT), DIMENSION(IFF_MAX) :: IFF_DK
      INTEGER (KIND=TINYINT), DIMENSION(IFF_MAX) :: IFF_WCHARM

      DIMENSION WCHARM_FINE(IFF_MAX)

C***  Iron: Common Block for Iron-Specific Data
      DIMENSION INDRB(MAXIND),
     >              SIGMAFE(INDEXMAX), 
     >              IFRBSTA(MAXIND), IFRBEND(MAXIND),
     >              IFENUP(MAXIND), IFELOW(MAXIND),
     >              OPAFE(NDDIM), ETAFE(NDDIM),
     >              INDFEACT(MAXIND)
      LOGICAL BFECHECK, BFEWING, BFEMODEL, BNUEFE

C***  WPI = SQUARE ROOT OF PI
      DATA WPI / 1.772453851 /

      DATA LASTL / 0 /

      DATA WJMAX / 0.9999999999 /

C***  CLIGHT = VELOCITY OF LIGHT IN KM/S
      DATA CLIGHT /2.99792458E5/


      SAVE XLAM0, XLAM0LN, ALN, CMFBAND, CMFBANDR
      SAVE DXMAX, EXLAM1, EXLAM2, XKCMAX
      SAVE LASTL, BANDM, BANDP, NLINE_NOFE, BSIMPLE
      SAVE XLP1, XLP2, KSPACE, OPA, ETA

C***  Treatment of Continua
      SAVE BNEWINTERVAL, BNEWINTERVALOLD, KCONT, BEMPTY
      SAVE XNUEK, XNUEKOLD, XNUEKOLDOLD, XNUE2, XNUE3
      SAVE FWKOLD, FWKOLDOLD, FWNUEKOLD, FWNUEKOLDOLD
      DIMENSION FWTEST(NF), XJC(ND,NF), XJCAPP(NF), XJCAPPNEW(NF)
      DIMENSION WCHARM(ND,NF)
      LOGICAL BSIMPLE, BNEWINTERVAL, BNEWINTERVALOLD, BEMPTY
      LOGICAL BNEWINT

      LOGICAL BKACTIVE, B_ACTIVE_CORE

C***  BNEWINT handles the fine frequency grid. If true, the same grid as in 
C***    COLI is applied. This is, at least for small Model Atoms, a bit faster
      BNEWINT = .TRUE.

c      if (l .eq. 10) then
c      write (0,'(A,2I2)') 
c     >  'XJAPP Start, L, ITNEL=', L,ITNEL
c      endif

C***  BPLOT switches plot files (same as done in COLI)
      BPLOT = .FALSE.

C***  Check if current Depthpoint is encountered for the first time
      BFIRSTITER = LASTL .NE. L
      LASTL = L

C***  Initialization for the first Depth and first Iteration
      IF (L .EQ. 1 .AND. BFIRSTITER) THEN

        WRITE (0,*) 'XJAPP: NEW=',    BXJLAPPNEW, BXJCAPPNEW
        WRITE (0,*) 'XJAPP: New_Grid, New_Oper=', BNEWINT, BNEWOPER
        WRITE (0,*) 'BNUEFE:    ', BNUEFE
        BSIMPLE = .TRUE.
C***    Plot-Range and Depth
        XLP1 = 0.
        XLP2 = 1.E20
        IPLOT = 60
        IPLOT_ITER = 1
C***    Maximum Spacing of Line-Frequency Points
        DXMAX = 0.3 
C***    Number of Fine frequencies around non-edge Continuum points
        BANDM = 1
        BANDP = 1
C***    Preparing harmonic wavelengths
        IF (BNEWINT) THEN
          XLAM0 = FF_INFO(1)
          ALN   = FF_INFO(2)
        ELSE
          XLAM0 = XLAMBDA(1)
          ALN   = ALOG(1. + VDOP/CLIGHT*DXMAX)
        ENDIF
        XLAM0LN = ALOG(XLAM0)
C***    Spacing when no lines and no continua are active
C***    At the moment 500 Lambda-points per decade are default
        KSPACE = INT((LOG(10.)/ALN) / 500.)
        WRITE (0,'(A,I4)') 'KSPACE=', KSPACE
C***    Maximum Half-Bandwidth for CMF Line Transfer in Doppler Units
        CMFBAND = 4.5
        CMFBANDR = CMFBAND + 2.*VELO1/VDOP
C***    No Opacity and Emissivity of Continuum
        OPA = 0.
        ETA = 0.
C***    Reorder lines in sequence of increasing wavelengths
        NLINE_NOFE = NLINE - LASTFE
        CALL SEQLINECL(NLINE_NOFE, LINE, EINST, INDLOW,
     >         INDNUP, XLAMSOR, ELEVEL,
     $         NDIM, VDOP, CMFBAND, CLIGHT, XLAMMIN, XLAMMAX,
     $         VELO1, EXLAM1, EXLAM2, MAXEXT )
C***    Range of Extended Lines : All
        EXLAM1 = XLAMBDA(1)
        EXLAM2 = XLAMBDA(NF)
C***    Preparing line indices (XKMIN, XKMAX)
        DO NL=1, NLINE_NOFE
          XK = (ALOG(XLAMSOR(NL)) - XLAM0LN) / ALN
          XKMID(NL) = XK
          XKMIN(NL) = XK - CMFBAND/DXMAX - 1.00000001
C***    Opacity
          XKRED(NL) = XK + CMFBAND/DXMAX + 1.00000001
C***    Red Line Wing
C***    NOTE : This could be set to XKRED, because fine frequencies
C***           are only needed for opacities and not for radiation fields
          IF (BXJLAPPCORE) THEN
             XKMAX(NL) = XKRED(NL)
          ELSE
             XKMAX(NL) = XK + CMFBANDR/DXMAX + 1.
          ENDIF
        ENDDO
C***    Prepare indices of the continua
c        write (0,*) 'Prepare indices of the continua'
C***    XKC  contains the Kontinua inserted for a correct flux integration
C***      Only two k-points are necessary
C***    XKC2 contains the true continuum points. Fine spacing is necessary 
C***      a wide k-range
        DO KON=1, NF
          XKC(KON) = (ALOG(XLAMBDA(KON)) - XLAM0LN) / ALN
          IF (BPLOT) THEN
            WRITE (38,'(A,F9.2,A)')
     >        'KASDEF SYM ', XKC(KON), ' -2. 0. 0.2 -0.2 8'
          ENDIF
        ENDDO
        DO KON=1, NF2
          XKC2(KON) = (ALOG(XLAMBDA2(KON)) - XLAM0LN) / ALN
          IF (BPLOT) THEN
            WRITE (38,'(A,F9.2,A)')
     >        'KASDEF SYM ', XKC2(KON), ' -1. 0. 0.2 -0.2 8'
          ENDIF
        ENDDO
        XKCMAX = (ALOG(XLAMBDA(NF)) - XLAM0LN) / ALN
C***    Open Plot-Files
c        write (0,*) 'Open Plot-Files'
        IF (BPLOT) THEN
          WRITE (0,'(A,I3,A)') 'Plot in depth ', IPLOT, ' prepared;  ', 
     >      'Iteration = ', IPLOT_ITER
C          WRITE (0,'(A,I6,A)') 'Plot at frequency ', LPLOT, ' prepared'
          DO NL=1, MAXLIN
            IF (NL .LT. 10) THEN 
              WRITE (NAME,'(A3,I1,A4)') 'cl_', nl, '    '
            ELSE
              WRITE (NAME,'(A3,I2,A3)') 'cl_', nl, '   '
            ENDIF
            OPEN(UNIT=40+NL, FILE=NAME)
            WRITE (40+NL,'(A3)') 'N=?'
          ENDDO
          OPEN (UNIT=40, FILE='coli.dat')
          OPEN (UNIT=39, FILE='coli_cont.dat')
          OPEN (UNIT=38, FILE='coli.kasdef')
          OPEN (UNIT=37, FILE='coli_opa.dat')
          WRITE (37,'(A6,45A16)')
     >    ' * IND', 'LAMBDA', 'OPA-KONT',
     >    'ETA-K', 'ETANOTH-K.', 'S-KONT',
     >    'S-KONT(NOTH)', 'OPA-L', 'ETA-L', 'S-L', 'J-KONT', 'J',
     >    'J-SCHLANGE', 'H-S', 'K-S', 'N-S', 'EDDIF', 'EDDIG',
     >    'XJLMO', 'XJLMOR2', 'XHLMO',
     >    'FULFIL0', 'FULFIL1', 'DJDSMO', 'DJDS', 'SNOTH'
          WRITE (40,'(A3)') 'N=?'
          WRITE (39,'(A3)') 'N=?'
          WRITE (38,'(A3)') 'N=?'
          WRITE (37,'(A3)') 'N=?'
        ENDIF
C***  Set Counters for Warnings
        IWARN_NEG_XJCAPP = 0
        IWARN_NEG_XJLAPP = 0
      ENDIF

C***  Consistency check between Fine grid, prepared in COLI and
C***    fine grid demanded by input cards
      IF (XLAMSTART .NE. FF_INFO(3) .OR. XLAMEND .NE. FF_INFO(4)) THEN
        WRITE (0,*) 'Inconsistency between fine grids of ', 
     >              'COLI and SETXJFINE'
        WRITE (0,'(A,4(1X,F12.5))') 
     >              'XLAMSTART, XLAMEND, FF_INFO(3), FF_INFO(4)=', 
     >               XLAMSTART, XLAMEND, FF_INFO(3), FF_INFO(4)
        STOP 'ERROR in Subr. SETXJFINE'
      ENDIF

C***  Initialize at first Iteration of each Depth Point
C***    Defintion of Line-Cores
c      IF (BFIRSTITER) THEN
c        DO NL=1, NLINE
c          INDEX = LINE(NL)
c          LACT = LIND(NL)
c          LACTS = LINDS(NL)
c          XKBLUE_CORE(NL) = XKMID(NL) - XBLUE(NL) / DXMAX
c          XKRED_CORE(NL)  = XKMID(NL) - XRED(NL) / DXMAX
c          write (0,'(A,4i3,4F20.10)') 'XJAPP: LINECORES: ', 
c     >      NL, index, lacts, lact, XBLUE(NL), XRED(NL), 
c     >      XKBLUE_CORE(NL), XKRED_CORE(NL)
c       ENDDO
cc      test_xjapp = 0.
cc      test = 1. / test_xjapp
cc        stop 'test'
c      ENDIF

      IF (BFIRSTITER) THEN
        CALL LOADWC(IFF_WCHARM, IFF_MAX_MS, L)
        CALL CALCWC(IFF_WCHARM, INT(FF_INFO(7)), L, WCHARM_FINE, 
     >              INT(FF_INFO(5)), GAMMAC)
      ENDIF

C***  Initialize at each Iteration
c      write (0,*) 'Initialize at each Iteration'

      IF (BNEWINT) THEN
        XKCMIN = FF_INFO(5)
        XKCMAX = FF_INFO(6)
      ELSE
        XKCMIN = (ALOG(XLAMSTART) - XLAM0LN) / ALN
        XKCMAX = (ALOG(XLAMEND)   - XLAM0LN) / ALN
      ENDIF
      XK = XKCMIN
      K = INT(XK)
      KSTART = K
      KLAST = K - KSPACE

      LINECHECK = 1
      DO 
        IF (XKMIN(LINECHECK) .LT. XK) THEN
          LINECHECK = LINECHECK + 1
        ELSE
          EXIT
        ENDIF
      ENDDO
      ILINECHECK = LINE(LINECHECK)

C***  Determine Index of first Kontinuum
      KONCHECK = 1
      DO
        IF (XKC(KONCHECK) .LT. XK) THEN
          KONCHECK = KONCHECK + 1
        ELSE
          EXIT
        ENDIF
      ENDDO
      KCFIRST = KONCHECK
      KONCHECK2 = 1

C***  Determine Index last Kontinuum. This is used for the 
C***    approximate radiation fields
      KCFIRST = KONCHECK
      KCLAST = KCFIRST
      DO
        IF (XKC(KCLAST) .LT. XKCMAX) THEN
          KCLAST = KCLAST + 1
        ELSE
          EXIT
        ENDIF
      ENDDO
      KCLAST = KCLAST - 1

C***  Determine Index of first Kontinuum for non-edge points
      DO
        IF (XKC2(KONCHECK2) .LT. XK) THEN
          KONCHECK2 = KONCHECK2 + 1
        ELSE
          EXIT
        ENDIF
      ENDDO

C***  Determine Index of last Line. This is used for the 
C***    approximate radiation fields
      KLLAST = 1
      DO
        IF (XKRED(KLLAST) .LT. XKCMAX) THEN
          KLLAST = KLLAST + 1
          IF (KLLAST .GT. NLINE_NOFE) EXIT
        ELSE
          EXIT
        ENDIF
      ENDDO
      KLLAST = KLLAST - 1
c      write (0,*) 'klllast=', kllast, xkred(kllast-1), xkred(kllast), 
c     >  xkred(kllast+1), xkcmax

c      write (0,*) 'Bereich eingeschraenkt:', 
c     >            xlamstart, xlamend, xkcmin, xkcmax, linecheck, 
c     >            line(linecheck), koncheck, koncheck2

      KOUT = -1
      NK = 0
      NLACT = 0
      DO NL=1, MAXLIN
        LIND(NL) = 0
      ENDDO
      CMODE = 'T'

C***  Start Indices for Continuum Interpolation (KCL = KCU)
      KCL = 0
      KCU = 0
C***  Width of the Interval used for Continuum Interpolation
      KCDMAX = 100
      KCDELTA = KCDMAX
C***  Initialize Counter for Continuum-Edge-Test
      KCCHECK = 1

C***  IRON: SET LOWEST LINE TO '1', NO ACTIVE FE-LINES
      INDFEACT(1) = 1
      MAXFEACT    = 0
      BFEWING = .FALSE.
      BFECHECK = .FALSE.

      KINDEX = 0

C*****************************************************
C***  Main-loop over all frequency-indices
C*****************************************************

c      write (0,*) 'Main-Loop START'
      mainfreqloop: DO !- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 

        KINDEX = KINDEX + 1
        IF (BNEWINT) THEN
          K = K + IFF_DK(KINDEX) + 100
        ENDIF
        XK = FLOAT(K)
        IF (K .EQ. 0) THEN
          XLAMK = XLAM0
        ELSE
          XLAMK = EXP(XLAM0LN + XK*ALN)
        ENDIF
        XNK = FLOAT(NK)

C***  IRON: Test for active FE-Transitions
        IF (BFEMODEL) THEN
           CALL FECHECK (XLAMK, INDRB, IFRBSTA, IFRBEND, LASTFE,
     >                    CLIGHT, VDOPFE, DXFE, XLAM0FE,
     >                    INDFEACT, MAXFEACT, BFECHECK, BFEWING,
     >                    DFEINDR)
        ENDIF

c        if (bfecheck) then
c          write (0,'(a,f10.4,i6,100i4)') 
c     >      'xjapp: k=', xlamk, k, (indfeact(ii), ii=1, maxfeact)
c        endif

C***  Continuum point near actual frequency?
        BKONCHECK = .FALSE.
        IF (KONCHECK .LE. NF) THEN
          DO
            IF (XK .LT. XKC(KONCHECK)-BANDM) THEN
              EXIT
            ENDIF
            IF (XK .GE. XKC(KONCHECK)-BANDM .AND.
     >          XK .LE. XKC(KONCHECK)+BANDP) THEN
              BKONCHECK = .TRUE.
              EXIT
            ENDIF
            IF (XK .GT. XKC(KONCHECK)+BANDP) THEN
              KONCHECK = KONCHECK + 1
              IF (KONCHECK .GT. NF) EXIT
            ENDIF
          ENDDO
        ENDIF
        BKONCHECK2 = .FALSE.
        IF (KONCHECK2 .LE. NF2) THEN
          DO
            IF (XK .LT. XKC2(KONCHECK2)-CMFBAND/DXMAX) THEN
              EXIT
            ENDIF
            IF (XK .GE. XKC2(KONCHECK2)-CMFBAND/DXMAX .AND.
     >          XK .LE. XKC2(KONCHECK2)+CMFBANDR/DXMAX) THEN
              BKONCHECK2 = .TRUE.
              EXIT
            ENDIF
            IF (XK .GT. XKC2(KONCHECK2)+CMFBANDR/DXMAX) THEN
              KONCHECK2 = KONCHECK2 + 1
              IF (KONCHECK2 .GT. NF2) EXIT
            ENDIF
          ENDDO
        ENDIF

C***    Test for Continuum Point in Interpolation Interval
        IF (KCCHECK .LE. NF) THEN
           DO
C***         First check, if Continuum is just passed
              IF (XK .GE. XKC(KCCHECK)) THEN
C***            Cont. Interval Finished
C***            ->  new Continuum, newstart (KCL = KCU)
                 KCCHECK = KCCHECK+1
                 KCL = K
                 KCU = K
                 IF (KCCHECK .GT. NF) EXIT
C***         Then find the Length of the Interpolation Interval
              ELSE
                 KCONTCHECK = INT(XKC(KCCHECK))
                 IF (KCU .LE. (KCONTCHECK-KCDMAX)) THEN
C***            Normal Case -> Length = Maximum Length
                    KCDELTA = KCDMAX
                    EXIT
                 ELSE
C***            Continuum detected -> Shorter Interval is cosen
                    KCDELTA = KCONTCHECK-KCU
                    EXIT
                 ENDIF
              ENDIF
           ENDDO
        ENDIF
 
                   
C***  Explaining Variables :
C***    XLAMSOR      : Array of sorted wavelengths (1..NLINE)
C***    XKMIN, XKMAX : Indices where the starts and ends
C***    LINECHECK    : Actual line which is checked
C***    ILINECHECK   : Index of actual line in original list
C***    LINDS        : Array of linecheck (1..MAXLIN)
C***    LIND         : Array of ILINECHECK (1..MAXLIN)
C***    BLASERL      : Is true for Laser Lines

c      GOTO 999

C***  Which lines are active?
C***    Add new lines (Only until last line KLLAST is reached)
        DO
          IF (XK .LT. XKMIN(LINECHECK) .OR. 
     >        LINECHECK .GT. NLINE_NOFE .OR. 
     >        LINECHECK .GT. KLLAST) EXIT
          IF (XK .GE. XKMIN(LINECHECK) .AND.
     >        XK .LE. XKMAX(LINECHECK)) THEN
C***  Search first free index
            DO NLNEW=1, MAXLIN
              IF (LIND(NLNEW) .EQ. 0) THEN
                LIND(NLNEW) = ILINECHECK
                LINDS(NLNEW) = LINECHECK
                IF (BXJLAPPNEW) THEN
                  XJLAPP(ILINECHECK) = XJL(L,ILINECHECK)
                ENDIF
c                if (bfirstiter) then
c                  prof(ILINECHECK) = 0.
c                endif
                  XKBLUE_CORE(NLNEW) = XKMID(LINECHECK) - 
     >              XBLUE(ILINECHECK) / DXMAX
                  XKRED_CORE(NLNEW)  = XKMID(LINECHECK) - 
     >              XRED(ILINECHECK) / DXMAX
c     >      LINECHECK, index, lacts, lact, 
c     >      XBLUE(ILINECHECK), XRED(ILINECHECK), 
c     >      XKBLUE_CORE(LINECHECK), XKRED_CORE(LINECHECK)

c          write (0,'(a,2i3,2i5,5f12.5)') '!!!!!!!', 
c     >      l, itnel, 
c     >      linecheck, ilinecheck, XKBLUE_CORE(LINECHECK), 
c     >      XKMID(LINECHECK), XKRED_CORE(LINECHECK), 
c     >      XBLUE(ILINECHECK), XRED(ILINECHECK)

C!!!                WS(NLNEW) = 0.
C***  Prepare line quantities
                CALL PRELINE (NUP(NLNEW), LOW(NLNEW), ILINECHECK, N, 
     >                 XLAM(NLNEW), 1, DUMMYA, 
     >                 ELEVEL, NDIM, INDNUP, INDLOW)
                CALL LIOP (EINST(NUP(NLNEW),LOW(NLNEW)), 
     >                 WEIGHT(LOW(NLNEW)), 
     >                 WEIGHT(NUP(NLNEW)), LOW(NLNEW), NUP(NLNEW), 1, 
     >                 XLAMSOR(LINECHECK), ENTOTL, EN, RSTAR, 
     >                 OPAL(NLNEW), ETAL(NLNEW), VDOP)
                IF (BPLOT .AND. L .EQ. IPLOT .AND. 
     >              ITNEL .EQ. IPLOT_ITER) THEN
                  WRITE (38,'(A,F9.2,1X,F3.0,1X,A,I4,1X,I4,F11.3)') 
     >              'KASDEF LUN ', XK, FLOAT(NLNEW), '0. 0.1 0.2 &E', 
     >              LINECHECK, LINE(LINECHECK), 
     >              XLAMSOR(LINECHECK)
                  WRITE (38,'(A,F9.2,A)')
     >              'KASDEF SYM ', XKMID(LINECHECK), 
     >              ' 0.0 0. 0.2 -0.2 8'
                  WRITE (38,'(A,F9.2,1X,F5.1,A)')
     >              'KASDEF SYM ', XKMID(LINECHECK), 
     >              FLOAT(NLNEW), ' 0. -0.2 -0.2 8'
                ENDIF
                NLACT = NLACT + 1
                IF (NLACT .GT. MAXLIN) THEN
                  WRITE (0,*) 'Capacity of MAXLIN exceeded'
                  STOP 'ERROR in Subr. COLI'
                ENDIF
                EXIT
              ENDIF
            ENDDO
          ENDIF
          LINECHECK = LINECHECK + 1
          ILINECHECK = LINE(LINECHECK)
        ENDDO

C***    Check the old lines
        B_ACTIVE_CORE = .FALSE.
        DO NL=1, MAXLIN
          LACT = LIND(NL)
          LACTS = LINDS(NL)
          IF (LACT .EQ. 0) CYCLE
C***     Check, if there are any active Line-Cores
          IF ((FLOAT(K+1) .GE. XKBLUE_CORE(NL)) .AND.
     >        (FLOAT(K-1) .LE.  XKRED_CORE(NL))) THEN
             B_ACTIVE_CORE = .TRUE.
          ENDIF
          IF (XKMAX(LACTS) .LT. XK) THEN
C***  Line has been finished 
            LIND(NL) = 0
            NLACT = NLACT - 1
C***  Set Negative XJLAPP to Zero
            IF (XJLAPP(LACT) .LT. 0.) THEN
              XJLAPP(LACT) = XJL(L,LACT)
              IWARN_NEG_XJLAPP = IWARN_NEG_XJLAPP + 1
c              WRITE (0,'(A,I3,1X,I5)') 
c     >        'Neg XJLAPP set to Zero: L, LINE=', L, LINE(LACT)
            ENDIF
          ELSE
            IF (BPLOT .AND. L .EQ. IPLOT .AND. 
     >          ITNEL .EQ. IPLOT_ITER) THEN
              WRITE (40+NL,'(I8,1X,F3.0)') K, FLOAT(NL)
            ENDIF
          ENDIF
        ENDDO

999     CONTINUE

        IF (BPLOT .AND. L .EQ. IPLOT .AND. 
     >      ITNEL .EQ. IPLOT_ITER .AND. BKONCHECK) THEN
          WRITE (39,'(I6,1X,F3.0)') K, -2.
        ENDIF
        IF (BPLOT .AND. L .EQ. IPLOT .AND. 
     >      ITNEL .EQ. IPLOT_ITER .AND. BKONCHECK2) THEN
          WRITE (39,'(I6,1X,F3.0)') K, -1.
        ENDIF

C***  Output of K
        IF (K .GT. KOUT) THEN
          IF (BFIRSTITER .AND. L .EQ. 1) THEN
            IF (K .EQ. KSTART) THEN
              WRITE (0,'(5X,A)') 
     >          'Fine Grid used in XJAPP'
              WRITE (0,'(5X,A6,1X,A6,1X,A13)') 
     >          'K', 'K-DONE', 'LAMBDA'
            ENDIF
            WRITE (0,'(5X,I6,1X,I6,1X,F13.3)') 
     >        K, NK, XLAMK
          ENDIF
          KOUT = KOUT + 20000
       ENDIF

C***  Calculate XJCAPP from XJC and XJCAPPNEW
      IF ( (FLOAT(K+1) .GT. XKCMAX      .AND. .NOT. BNEWINT) .OR. 
     >     (KINDEX .EQ. INT(FF_INFO(7)) .AND.       BNEWINT) ) THEN
        CMODE = 'F'
c        write (0,*) 'FLOAT(K) .EQ. XKCMAX -> CMODE=', CMODE
        IF (BXJCAPPNEW) THEN
          IF (BFIRSTITER) THEN
            DO KK=KCFIRST + 1, KCLAST - 1
              XJCAPP(KK)    = XJC(L,KK)
            ENDDO
          ELSE
C***  Set Negative XJCAPP to Zero
c          IF (.NOT. BFIRSTITER) THEN
            DO KK=KCFIRST + 1, KCLAST - 1
              IF (KK .EQ. 1) THEN
C***  Do not Scharmer first Continuums frequency
                XJCAPP(1)     = XJC(L,1)
              ELSE
                XJCAPP(KK)    = XJCAPPNEW(KK)/FWTEST(KK) + XJC(L,KK)
              ENDIF
              IF (XJCAPP(KK) .LT. 0.) THEN
                XJCAPP(KK) = XJC(L,KK)
                IWARN_NEG_XJCAPP = IWARN_NEG_XJCAPP + 1
c                WRITE (0,'(A,I3,1X,I5)') 
c     >            'Neg XJCAPP set to XJC: L, K=', L, KK
              ENDIF
            ENDDO
          ENDIF
        ENDIF
      ENDIF

        IF (.NOT. BNEWINT) THEN
C***   IRON: TEST FOR ACTIVE IRON LINE(WING)
          IF (.NOT. BKONCHECK .AND. .NOT. BKONCHECK2 .AND.
     >        .NOT. BFECHECK
     >        .AND. NLACT .EQ. 0 .AND. 
     >        K-KLAST .LT. KSPACE .AND. K .GT. 0 .AND. .NOT. 
     >        CMODE .EQ. 'F') THEN 
            K = K + 1
            IF (K .GT. KCU) THEN
C***          New Start Indices for Continuum Interpolation (KCL = KCU)
               KCL = K
               KCU = K
            ENDIF
            CYCLE
          ENDIF
          IF (BPLOT .AND. L .EQ. IPLOT .AND. 
     >        ITNEL .EQ. IPLOT_ITER) THEN
            WRITE (40,'(I8,1X,F3.0)') K, 0.
          ENDIF
        ENDIF

C***  BKACTIVE = .TRUE. -> Frequency Point is activated
        BKACTIVE = .NOT.BXJLAPPCORE .OR. B_ACTIVE_CORE

C***  Main loop is now restricted to necessary points
C***  Now the radiation transfer is calculated for index K
C***  Prepared quantities :
C***    K, XLAMK                   actual index and its wavelength
C***    BKONCHECK,  KONCHECK       near kontinuum and index of next kontinuum
C***    BKONCHECK2, KONCHECK2      near kontinuum and index of next true kontinuum
C***    NLACT, LIND(NL(1..MAXLIN)) Number of active lines and indices of the lines
C***    

cc        CALL CMFCOOP (XLAMK,1,TL,RNEL,EN,ENTOTL,RSTAR,
cc     $             OPA,ETANOTH,THOMSON,NDIM,N,NCHARG,WEIGHT,ELEVEL,EION,
cc     $             EINST,ALPHA,SEXPO,ADDCON1,ADDCON2,ADDCON3,IGAUNT,
cc     $             SIGMA1I,KONTLOW,KONTNUP,LASTKON,NATOM,KONTHLP,
cc     >             1. ,.FALSE.,.FALSE.,IPLOT,K,KCL,KCU,KCDELTA,
cc     >             OPACL,OPACU,ETACL,ETACU,XLAM0LN,ALN)

        IF (BKACTIVE) THEN
          IF (XLAMK < EXLAM1) THEN
            WRITE (0,*) XLAMK, EXLAM1, EXLAM2
            WRITE (0,*) XLAMK, XLAMBDA(1), XLAMBDA(NF)
            STOP 'ERROR in SETXJFINE: XLAMK too small'
          ELSEIF (XLAMK > EXLAM2) THEN
            WRITE (0,*) XLAMK, EXLAM1, EXLAM2
            WRITE (0,*) XLAMK, XLAMBDA(1), XLAMBDA(NF)
            WRITE (0,*) INT(FF_INFO(7)), NF
            WRITE (0,*) INT(FF_INFO(7)), NF, KINDEX
            STOP 'ERROR in SETXJFINE: XLAMK too large'
          ENDIF
          CALL LIPO (OPA,     XLAMK, OPAC, XLAMBDA, NF)
          CALL LIPO (ETANOTH, XLAMK, ETAC, XLAMBDA, NF)
          THOMSON = .0
c***      damit macht ADDOPA aus OPA -> OPANOTH, was es schon ist
       
        IF (BFECHECK) THEN
           CALL  FEOP_STEAL (XLAMK, N, INDFEACT, MAXFEACT, LASTFE,
     >                    SIGMAFE, OPAFE, ETAFE, INDRB,
     >                    IFRBSTA, IFRBEND, IFENUP, IFELOW,
     >                    VDOPFE, DXFE, XLAM0FE,
     >                    ELEVEL, WEIGHT, RSTAR, EN, ENTOTL, TL,
     >                    BNUEFE)
        ENDIF

C***  Calculate frequency Step in X
C        DELTAX = FLOAT(K-KLAST) * DXMAX
        NDK = K - KLAST
        IF (NDK .NE. 1) THEN
           DELTAX = (EXP(FLOAT(NDK)*ALN) - 1.)*CLIGHT/VDOP
        ELSE
           DELTAX = DXMAX
        ENDIF

        FWEIGHTOLD = FWEIGHTL
        FWEIGHTL = DELTAX * VDOP * 1.E13 / XLAMK

C***  Sum up Opacities
        PARALAS = 0.01
        CALL ADDOPA (1, 1, MAXLIN, MAXIND, LIND, LINDS,
     >         XK, XKMID, XKRED, DELTAX,
     >         PARALAS, LASER,
     >         WS, ETAL, OPAL, ETA, ETANOTH, OPA, ETAK, ETAKNOTH,
     >         OPAK, OPAKNOTH, THOMSON,
     >         PWEIGHT, OPAFE, ETAFE, BFECHECK, BLASERL, NUP, LOW,
     >         LEVEL)

        IF (ETAKNOTH .LT. 0.) THEN
          ETAKNOTH = 0.
        ENDIF

C***  Store SFINE
        IF (BFIRSTITER) THEN
          IF (OPAKNOTH .LE. 0.) THEN
            SFINE_OLD(K) = 0.
          ELSE
            SFINE_OLD(K) = ETAKNOTH / OPAKNOTH
          ENDIF
        ELSE
          IF (OPAKNOTH .LE. 0.) THEN
            SFINE_NEW(K) = 0.
            WCORE = 0.
          ELSE
            SFINE_NEW(K) = ETAKNOTH / OPAKNOTH
            WCORE = WJMAX
          ENDIF
        ENDIF

C-------------------------------------------------------------------------
C---  Now integration of XJCAPPNEW and XJLAPP
C---    Will later be placed in a separate Subroutine!!!
C-------------------------------------------------------------------------

C***  Integration of Lines
        IF (.NOT. BFIRSTITER) THEN
          DSOLD = DS
          DSOLD_WCHARM = DS_WCHARM
          DS = SFINE_NEW(K) - SFINE_OLD(K)

          IF (BNEWINT .AND. BNEWOPER) THEN
            DS_WCHARM = DS * WCHARM_FINE(KINDEX)
          ELSE
C***  Calculate DS multiplied with the Interpolated WCHARM
            DS_WCHARM = DS * 
     >        (WCHARM(L,KCONT)   * DLEFT_WCHARM + 
     >         WCHARM(L,KCONT+1) * DRIGHT_WCHARM)
          ENDIF

C***    Calculate XJLAPP
          IF (BXJLAPPNEW) THEN
            DO NL=1, MAXLIN
              LACT = LIND(NL)
              LACTS = LINDS(NL)
              IF (LACT .EQ. 0) CYCLE
              IF (XKRED_CORE(NL) .LE. XKBLUE_CORE(NL)) THEN
                CYCLE
              ENDIF
              IF (BNEWINT .AND. BNEWOPER) THEN
                XREL = DXMAX * (XKMID(LACTS) - FLOAT(K))
                PHIL = EXP(-XREL * XREL) / WPI
                XJLAPP(LACT) = XJLAPP(LACT) + DS_WCHARM * PHIL * DXMAX
              ELSE
                DKLEFT  = AMAX1(XKBLUE_CORE(NL), FLOAT(K-1))
                DKRIGHT = AMIN1(XKRED_CORE(NL),  FLOAT(K+1))
                DK = 0.5 * AMAX1(0., DKRIGHT-DKLEFT)
                IF (DK .EQ. 0.) CYCLE
                XREL = DXMAX * (XKMID(LACTS) - FLOAT(K))
                PHIL = EXP(-XREL * XREL) / WPI
                XJLAPP(LACT) = XJLAPP(LACT) + DS*PHIL*DK*DXMAX*WCORE
              ENDIF
c           if (lact .eq. 2) then
c             write (0,'(4(F11.5,1x), F8.5, (1x,F11.5), 1x, e20.10, 
c     >1x, F3.1, e20.10)') 
c     >         XKBLUE_CORE(LACTS), XKMID(LACTS), 
c     >         XKRED_CORE(LACTS), XK, dk, xrel, phil, DS, XJLAPP(LACT) 
c           endif
            ENDDO
          ENDIF          


C***    Calculate XJCAPPNEW
          IF (BXJCAPPNEW) THEN
C***  Integration of Continua
C********************************************************************
C***  Averaging Delta S from fine to coarse frequency grid
C********************************************************************

C***  Note: The current (fine) frequency point bears the index K; the 
C***        integration step trails and hits the frequency index K-1 = OLD
C***  KOLD defines the current coarse interval (KCONT, KCONT+1)

C***  Check whether the current frequency point is the first which 
C***  has entered a new coarse interval
            IF (K .EQ. KSTART) THEN 
              BNEWINTERVALOLD = .FALSE.
              KCONT = 1
              BEMPTY = .FALSE.
            ELSE
              BNEWINTERVALOLD =  BNEWINTERVAL
            ENDIF
            BNEWINTERVAL = XLAMK .GE. XLAMBDA(KCONT+1) 

C***  Frequency Integrations Weight for old finegrid point
C***  Integration weight nue*dnue for test function 
C***   Note: First and last interval get additional terms - see below!
            IF (K .EQ. KSTART) THEN
C!!!  ACHTUNG: DER Vektor fwtest ist ueberfluessig; wird in WMODCOLI verwendet. 
              IF (BSIMPLE) THEN
                DO KC = 1, NF
                  FWTEST(KC) = 0.
                ENDDO
              ELSE
                DO KC = 1, NF
                  FWTEST(KC) = 1.
                ENDDO
              ENDIF
              SUMS     = .0
              SUMSW    = .0
              XNUEK = CLIGHT / XLAMK
C***  For K=0 nothing is integrated; Jump to end of Integration
              GOTO 100

            ELSE IF (NK .EQ. 1) THEN
              XNUEKOLD  = XNUEK
              XNUEK     = CLIGHT / XLAMK
              XNUE2 = CLIGHT / XLAMBDA(1)
              XNUE3 = CLIGHT / XLAMBDA(2)
              IF (BSIMPLE) THEN
                DLEFT = (XNUEKOLD - XNUE2) / (XNUE3 - XNUE2)
                DRIGHT = 1. - DLEFT
                DLEFT_WCHARM = DLEFT
                DRIGHT_WCHARM = DRIGHT
                FWKOLD      = 0.5 * (XNUEKOLD - XNUEK)
                FWNUEKOLD   = FWKOLD
              ELSE
                FWKOLD      = 0.5 * (XNUEKOLD - XNUEK)
                FWNUEKOLD   = (2. * XNUEKOLD + XNUEK) *
     >                        (XNUEKOLD-XNUEK) / 6.        
                DLEFT = 1.
                DRIGHT = 1.
                DLEFT_WCHARM = (XNUEKOLD - XNUE2) / (XNUE3 - XNUE2)
                DRIGHT_WCHARM = 1. - DLEFT_WCHARM
              ENDIF
            ELSE
              XNUEKOLDOLD = XNUEKOLD
              XNUEKOLD  = XNUEK
              XNUEK     = CLIGHT / XLAMK
              IF (BNEWINTERVALOLD) THEN
                XNUELEFT = CLIGHT / XLAMBDA(KCONT)
              ELSE
                XNUELEFT = XNUEKOLDOLD
              ENDIF
              IF (BNEWINTERVAL) THEN
                XNUERIGHT = CLIGHT / XLAMBDA(KCONT+1)
              ELSE
                XNUERIGHT = XNUEK
              ENDIF
              IF (BSIMPLE) THEN
                DLEFT = (XNUEKOLD - XNUE2) / (XNUE3 - XNUE2)
                DRIGHT = 1. - DLEFT
                DLEFT_WCHARM = DLEFT
                DRIGHT_WCHARM = DRIGHT
                FWKOLD      = 0.5 * (XNUELEFT-XNUERIGHT)
                FWNUEKOLD   = FWKOLD
              ELSE
                FWKOLD      = 0.5 * (XNUELEFT-XNUERIGHT)
                FWNUEKOLD   = (XNUELEFT + XNUEKOLD + XNUERIGHT) *
     >                        (XNUELEFT-XNUERIGHT) / 6.        
                DLEFT = 1.
                DRIGHT = 1.
                DLEFT_WCHARM = (XNUEKOLD - XNUE2) / (XNUE3 - XNUE2)
                DRIGHT_WCHARM = 1. - DLEFT_WCHARM
              ENDIF
            ENDIF

C***  Here is the entry for an empty interval:
   10       CONTINUE

            IF (KCONT .LT. 1) THEN
              WRITE (0,*) '== WARNING from FREQUINT: KCONT .LT. 1'
              KCONT = 1
            ENDIF
            IF (KCONT .GT. NF-1) THEN
              WRITE (0,*) '== WARNING from FREQUINT: KCONT .GT. NF-1'
              write (0,*) '!!!', kcont, nf-1
              KCONT = NF-1
            ENDIF

            IF (.NOT. BEMPTY) THEN
              IF (BSIMPLE) THEN
                FWTEST(KCONT)  = FWTEST(KCONT)  + DRIGHT * FWKOLD
                FWTEST(KCONT+1)= FWTEST(KCONT+1)+ DLEFT  * FWKOLD
              ENDIF
C***        Average the Jnue's from the fine frequency step for the coarse
C***        continuum frequency mesh used by STEAL
              SUMS   = SUMS  + DSOLD_WCHARM * FWKOLD * DLEFT_WCHARM
C***        The same is done after multiplication with the testfunction NUE
              SUMSW  = SUMSW + DSOLD_WCHARM * FWNUEKOLD * DRIGHT_WCHARM
c marke
c         write (0,'(a,i6,3e20.10)') 'k=', 
c     >     k, djdsc(ltestdjds), sumdjdsc(ltestdjds), 
c     >     sumdjdscw(ltestdjds)
            ENDIF

            IF (BNEWINTERVAL) THEN 
C***  Now the interval (KCONT, KCONT+1) is completed

C***     Complete the integration with the part of the step from 
C***     XLAMKOLD to XLAMBDA(KCONT+1), i.e. XNUEKOLD ... XNUE2
C***     and subtract (!) the excess part of the last integration step
              XNUE1 = CLIGHT / XLAMBDA(KCONT  )
              XNUE2 = CLIGHT / XLAMBDA(KCONT+1)
              XNUE3 = CLIGHT / XLAMBDA(KCONT+2)

              QSPECIAL = (XNUEKOLD - XNUE2) / (XNUEKOLD - XNUEK)
              QSPECIAL2 = 1. - QSPECIAL

              IF (BSIMPLE) THEN
                FWSPECIAL    = 0.5 * (XNUEKOLD - XNUE2)
                FWNUESPECIAL = FWSPECIAL
                DLEFT = 1.
                DRIGHT = 0.
                FWTEST(KCONT)  = FWTEST(KCONT  )
     >                           + DRIGHT * FWSPECIAL
                FWTEST(KCONT+1)= FWTEST(KCONT+1)
     >                           + DLEFT * FWSPECIAL
              ELSE
                FWSPECIAL    = 0.5 * (XNUEKOLD - XNUE2)
                FWNUESPECIAL = (XNUEKOLD - XNUE2) * 
     >                         (2.*XNUE2 + XNUEKOLD) / 6.
                DLEFT = 1.
                DRIGHT = 1.
              ENDIF

              DSSPECIAL_WCHARM = QSPECIAL2 * DSOLD_WCHARM + 
     >                           QSPECIAL * DS_WCHARM
              SUMS   = SUMS  + 
     >                 DSSPECIAL_WCHARM * FWSPECIAL * DLEFT_WCHARM
              SUMSW  = SUMSW + 
     >                 DSSPECIAL_WCHARM * FWNUESPECIAL * DRIGHT_WCHARM

C***     Generate intergration weights at coarse frequency points
C***     which incorporate the test function, i.e. nue*dnue
              FWNUE1 = (2.*XNUE1 + XNUE2) * (XNUE1-XNUE2) / 6.
              FWNUE2 = (2.*XNUE2 + XNUE1) * (XNUE1-XNUE2) / 6.

              FW1 = 0.5 * (XNUE1 - XNUE2)         
              FW2 = FW1
              IF (KCONT .EQ. 1) THEN 
                FP = 1.
                XJCAPPNEW(1) = 0.
              ELSE
                XNUE0 = CLIGHT / XLAMBDA(KCONT-1) 
                FP = (XNUE1 - XNUE2) / (XNUE0 - XNUE2) 
              ENDIF
              FQ = 1. - FP

              if (kcont+1 .gt. nf) write (0,*) '!!!!!!!!1', kcont
              if (kcont .lt. 1) write (0,*) '!!!!!!!!2', kcont
              IF (BSIMPLE) THEN
                XJCAPPNEW(KCONT+1) = SUMS
                XJCAPPNEW(KCONT)   = XJCAPPNEW(KCONT) + SUMSW
              ELSE
C***       Note: At the left point, J has already the contribution 
C***                from the right end of the last interval
C***       Note: The contributions are now weighted 
                XJCAPPNEW(KCONT+1) =
     >               ( FWNUE1 * SUMS - FW1 * SUMSW ) / 
     >               ( FW2 * FWNUE1 - FW1 * FWNUE2)
                XJCAPPNEW(KCONT) = FQ * XJCAPPNEW(KCONT) + FP * 
     >               ( FWNUE2 * SUMS - FW2 * SUMSW)/ 
     >               (FW1 * FWNUE2 - FW2 * FWNUE1)
              ENDIF

C***     Start the next integration with the part of the step from
C***     XLAMBDA(KCONT+1) to XLAMK, i.e. XNUE2 ... XNUEK

              IF (BSIMPLE) THEN
C!!!              FWSPECIAL    = 0.5 * (XNUE2 - XNUEK)
                FWSPECIAL    = 0.
                FWNUESPECIAL = FWSPECIAL
                DLEFT = 0.
                DRIGHT = 0.
              ELSE
                FWSPECIAL    = 0.5 * (XNUE2 - XNUEK)
                FWNUESPECIAL = (XNUE2 - XNUEK) * (2.*XNUE2 + XNUEK) / 6.
                DLEFT = 1.
                DRIGHT = 1.
              ENDIF

              DSSPECIAL_WCHARM = QSPECIAL2 * DSOLD_WCHARM + 
     >                           QSPECIAL * DS_WCHARM
              SUMS   = DSSPECIAL_WCHARM * FWSPECIAL * DLEFT_WCHARM
              SUMSW  = DSSPECIAL_WCHARM * FWNUESPECIAL * DRIGHT_WCHARM

C***     Now increment the KCONT index and check for empty coarse interval 
C***     (i.e. which contains no fine frequency point)
              KCONT = KCONT + 1
              IF (XLAMK .GT. XLAMBDA(KCONT+1)) THEN
                BEMPTY = .TRUE.
                XNUELEFT  = CLIGHT / XLAMBDA(KCONT)
                XNUERIGHT = CLIGHT / XLAMBDA(KCONT+1)
                GOTO 10
              ELSE
                BEMPTY = .FALSE.
              ENDIF

            ENDIF

C***  XJC
c          if (l .eq. 66) then
c          write (0,'(4e20.10,i7,2e20.10)')
c     >      sfine_old(k), sfine_new(k), 
c     >      dsold, dsold_wcharm, kcont, xjc(l,7), xjcappnew(7)
c          endif

          ENDIF          
C***  BFIRSTITER
        ENDIF
C***  BKACTIVE:
        ENDIF
C-------------------------------------------------------------------------
C---  End of Integration
C-------------------------------------------------------------------------

  100   CONTINUE

C***  Test-output at given depth for all wavelengths between XLP1 and XLP2
        IF (BPLOT .AND. L .EQ. IPLOT .AND. ITNEL .EQ. IPLOT_ITER) THEN
          IF (XLAMK .GE. XLP1 .AND. XLAMK .LE. XLP2) THEN
            SK1 = ETA/OPA
            SK2 = ETANOTH/OPA
            SL = ETAK/OPAK
            SLNOTH = ETAKNOTH/OPAKNOTH
            WRITE (37,'(I8,27(1X,E15.7))')
     >        K, XLAMK, OPA, ETA, ETANOTH,
     >        SK1, SK2,
     >        OPAK, ETAK, SL,
     >        -1., -1., 
     >        -1., 
     >        -1., -1., -1., 
     >        -1., -1., 
     >        -1., -1., 
     >        -1., 
     >        -1., -1., -1., 
     >        -1., SLNOTH
C     >        XHI, XHID
          ENDIF
        ENDIF

        KLAST = K
        IF (.NOT. BNEWINT) THEN
          K = K + 1
        ENDIF
        NK = NK + 1

C***  Last frequency has been finished : EXIT Main Loop
        IF (CMODE .EQ. 'F') EXIT mainfreqloop

C***  End of Main-Loop over all frequencies
      ENDDO mainfreqloop !- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

      IF (BFIRSTITER .AND. L .EQ. 1) THEN
        WRITE (0,'(5X,I6,1X,I6,1X,F13.3)') 
     >    K, NK, XLAMK
      ENDIF

c      stop 'Test in XJAPP'

      RETURN

C***  Some other statements
cC***  Renormalization of Line-Profiles (at first ITER)
c      dimension prof(5000)
c          DO NL=1, MAXLIN
c            LACT = LIND(NL)
c            LACTS = LINDS(NL)
c            IF (LACT .EQ. 0) CYCLE
c            XREL = DXMAX * (XKMID(LACTS) - FLOAT(K))
c            PHIL = EXP(-XREL * XREL) / WPI
ccc        write (0,*) '-------------', 
ccc     >    xrel, phil
c            PROF(LACT) = PROF(LACT) + 0.3 * PHIL
c          ENDDO


      END
