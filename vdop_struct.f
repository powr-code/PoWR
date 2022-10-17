      SUBROUTINE VDOP_STRUCT (BDD_VDOP, DD_VDOP_LINE, DD_VDOP, VDOP, 
     >                        VELO, T, 
     >                        ND, NDDIM, NATOM, MAXATOM, DD_VDOPDU, 
     >                        VMICFRAC_DEFAULT, ATMASS, 
     >                        XMAX, XMAXMIN, 
     >                        SYMBOL, VDOPFE, 
     >                        BMICROTURB, BIRONLINES,
     >                        DD_VMIC, TAUROSS, RADIUS, EL_MIN, 
     >                        DD_VMICDU, bDDFECONVOL, IVDOPSTATUS)

C********************************************************************************
C***  Called from formal, 
C***  This routine prepares the depth-dependent DD_VDOP and further,
C***  related, depth-dependent arrays ("DD_..."). 
C***
C***  DD_VDOP(L,NA) returns the Dopplerbroadening -velocity 
C***  at depth point L of atom NA. 
C***  If no depth-dependence is specified, neither for VDOP not VMIC,
C***  then DD_VDOP(L,NA) = VDOP = constant.
C***
C***  Otherwise, 
C***  DD_VDOP(L,NA) = MAX(VDOP, SQRT (Vmic(L)^2 + Vtherm(L,NA)^2))
C***
C***  The squared thermal velocity is given by 2 * k_B * T[L] / m[NA]
C***  
C***  The stratification of Vmic(L) can be specified by the VMIC 
C***  line in FORMAL_CARDS.
C***  Current options: VDOP X.X [VELOFRAC|OUTERVMIC|VMICCONST [X.X]].
C***  In the first two versions, Vmic[L] = f*VELO[L], where the user either
C***  gives f explicity (VELOFRAC) or implicitly 
C***  (OUTERVMIC - f = OUTERVMIC/VELO(1))
C***  In the last version, VMIC is assumed to be constant.
C***  If no value is  stated: 

C***  VDOP is the MINIMUM value of DD_VDOP (relevant for the frequency
C***  spacing to resolve the narrowest lines).  
C***
C***  In case of a SECONDMODEL, this subroutine is caaled a second time, 
C***   with the corresponding parts of the arrays addressed by the CALL.
C***
C********************************************************************************
      IMPLICIT NONE
      
      INTEGER, INTENT(IN) :: ND, NDDIM, NATOM, MAXATOM
      INTEGER NA_VDOPMIN
      INTEGER  :: L, NA, NPAR, I, NA_MIN, ISRCHEQ, IVDOPSTATUS
      REAL, DIMENSION(NATOM), INTENT(IN) :: ATMASS

      REAL, DIMENSION(ND), INTENT(IN) :: VELO, T, TAUROSS, RADIUS
      REAL, DIMENSION(ND) :: DD_VMIC, DD_VMICDU

      REAL, INTENT(IN) :: VMICFRAC_DEFAULT, XMAXMIN 
      REAL :: VDOP, VMIC_FRAC, VDOP_FRAC, DD_VMIC_SQRD, VDOPMIN, MASS 
      REAL :: THERMVELO_SQRD, XMAX, Q, DUMMY
      REAL :: VMIC_MIN, VMIC_MAX, VDOP_MAX, DIST_IN, DIST_OUT, DX, X
      CHARACTER*2 SYMBOL(NATOM)
      REAL, INTENT(IN) :: VDOPFE
      REAL, DIMENSION(NDDIM, MAXATOM):: 
     >       DD_VDOP, DD_VDOPDU
      LOGICAL :: BDD_VDOP, BMICROTURB, BIRONLINES, bDDFECONVOL
      CHARACTER ACTPAR*20, DD_VDOP_LINE*(*), EL_MIN*2
      INTEGER, EXTERNAL :: IDX
C***  WPIINV = 1. / SQRT(PI)
      REAL, PARAMETER :: WPIINV = 0.564189583549 
      REAL, PARAMETER :: PI = 3.141592654 
      REAL VDOPNOIRONMAX
C*** atomic mass unit in g
      REAL, PARAMETER :: AMU = 1.6605387E-24 
C*** IMPORTANT! K_B is given here in units km^2 g s^-2 K^1 to obtain 
C*** a velocity in [km/s]!
      REAL, PARAMETER :: BOLTZK = 1.3807E-26
      

C***  VMIC VERSION => VDOP^2 = VMIC^2 + VTHERM^2 ************
C***  Note: DD_VDOP(L,NA) depends on element (NA)!
      IF (BMICROTURB) THEN
        WRITE(0,'(A)') "VDOP is depth-dependent"
        BDD_VDOP = .TRUE.
        WRITE(0,'(A)') "VDOP^2 = vmic^2 + vtherm^2"
        CALL SARGC (DD_VDOP_LINE, NPAR)
        IF ((NPAR .LT. 2) .OR. (NPAR .GT. 8)) GOTO 102
        CALL SARGV (DD_VDOP_LINE, 2, ACTPAR)
        READ (ACTPAR, '(F20.0)', ERR = 102) VMIC_MIN
C***    VMIC(L) is never smaller than specified by user
        WRITE(0,'(A,F10.5)') "Inner microturbulence:", VMIC_MIN
        IF (NPAR .EQ. 2) THEN
           WRITE(0,'(A,F10.2)') 'VMIC is constant:', VMIC_MIN
           DO L=1, ND
            DD_VMIC(L) = VMIC_MIN
           ENDDO
        ELSE 
            CALL SARGV (DD_VDOP_LINE, 3, ACTPAR)
C***        In this branch: vmic(L) = VMIC_FRAC * VELO(L)
            IF (ACTPAR .EQ. 'VELOFRAC') THEN
                WRITE(0,'(A)') "VMIC is fraction of wind velocity"
                IF (NPAR .EQ. 3) THEN
                    WRITE(0,'(A,F10.5)') 
     >               "fraction set to default value:", VMICFRAC_DEFAULT
                    VMIC_FRAC = VMICFRAC_DEFAULT
                ELSE IF (NPAR .EQ. 4) THEN
                    CALL SARGV (DD_VDOP_LINE, 4, ACTPAR)
                    READ (ACTPAR, '(F20.0)', ERR = 102) VMIC_FRAC
                    WRITE(0,'(A,F10.5)') 
     >               "fraction set by user to", VMIC_FRAC
                ELSE
                    GOTO 102
                ENDIF
                DO L=1, ND
                    DD_VMIC(L) = AMAX1(VMIC_MIN, VMIC_FRAC*VELO(L))
                ENDDO
C***        In this branch, the outer VMIC is specified by user
            ELSE IF (ACTPAR .EQ. 'MAX') THEN
                CALL SARGV (DD_VDOP_LINE, 4, ACTPAR)
                READ (ACTPAR, '(F20.0)', ERR = 102) VMIC_MAX
                WRITE(0,'(A,F6.1)') "Outer microturbulence:", VMIC_MAX
C***            If no further arguments are given, 
C***              VMIC is calculated as VMIC(L) = VMIC_FRAC * VELO(L)
                IF (NPAR .EQ. 4) THEN
                    WRITE(0,'(A)') 
     >              "microturbulence = fraction of wind velocity"
                    VMIC_FRAC = VMIC_MAX/VELO(1)
                    WRITE(0,'(A,F10.5)') 
     >               "fraction implied from maximum microturbulence:", 
     >               VMIC_FRAC
                    DO L=1, ND
                        DD_VMIC(L) = AMAX1(VMIC_MIN, VMIC_FRAC*VELO(L))
                    ENDDO
C***            Otherwise, inner and outer interpolation points 
C***            and method specified by user:
                ELSE IF (NPAR .EQ. 8) THEN
                    CALL SARGV (DD_VDOP_LINE, 5, ACTPAR)
C***                1) Interpolation on velocity
                    IF (ACTPAR .EQ. 'VELO1') THEN
                        CALL SARGV (DD_VDOP_LINE, 7, ACTPAR) 
                        IF (ACTPAR .NE. 'VELO2') GOTO 102
                        CALL SARGV (DD_VDOP_LINE, 6, ACTPAR)
                        READ (ACTPAR, '(F20.0)', ERR = 102) DIST_IN
                        CALL SARGV (DD_VDOP_LINE, 8, ACTPAR)
                        READ (ACTPAR, '(F20.0)', ERR = 102) DIST_OUT
                        IF (DIST_IN > DIST_OUT) THEN
                            DUMMY = DIST_IN
                            DIST_IN = DIST_OUT
                            DIST_OUT = DUMMY
                        ENDIF 
C***                    interpolation boundaries may not exceed model boundaries
                        DIST_IN = AMAX1(DIST_IN, VELO(ND))
                        DIST_OUT = AMIN1(DIST_OUT, VELO(1))
                        WRITE(0,'(A,F10.5, A, F10.5)') 
     >                   "Vmic interpolated between wind velocity", 
     >                   DIST_IN, " and ", DIST_OUT
                        DX = DIST_OUT - DIST_IN
C***                    Interpolation takes place here
                        DO L=1, ND
                            IF (VELO(L) .LE. DIST_IN) THEN
                                DD_VMIC(L) = VMIC_MIN
                            ELSE IF  (VELO(L) .GE. DIST_OUT) THEN
                                DD_VMIC(L) = VMIC_MAX
                            ELSE
                                X = PI * (VELO(L) - DIST_IN ) / DX
                                Q = 0.5 + 0.5 * COS(X)
                                DD_VMIC(L) = Q * VMIC_MIN  + (1.-Q) * VMIC_MAX
                            ENDIF
                        ENDDO
C***                2) Interpolation on Rosseland tau (analog to last block)
                    ELSE IF (ACTPAR .EQ. 'TAU1') THEN
                        CALL SARGV (DD_VDOP_LINE, 7, ACTPAR) 
                        IF (ACTPAR .NE. 'TAU2') GOTO 102
                        CALL SARGV (DD_VDOP_LINE, 6, ACTPAR)
                        READ (ACTPAR, '(F20.0)', ERR = 102) DIST_IN
                        CALL SARGV (DD_VDOP_LINE, 8, ACTPAR)
                        READ (ACTPAR, '(F20.0)', ERR = 102) DIST_OUT
                        IF (DIST_IN < DIST_OUT) THEN
                            DUMMY = DIST_IN
                            DIST_IN = DIST_OUT
                            DIST_OUT = DUMMY
                        ENDIF 
                        DIST_IN = AMIN1(DIST_IN, TAUROSS(ND))
                        DIST_OUT = AMAX1(DIST_OUT, TAUROSS(1))
                        WRITE(0,'(A,F10.5, A, F10.5)') 
     >                   "Vmic interpolated between Rosseland tau", 
     >                   DIST_IN, " and ", DIST_OUT
                        DX = DIST_IN - DIST_OUT
                        DO L=1, ND
                            IF (TAUROSS(L) .GE. DIST_IN) THEN
                                DD_VMIC(L) = VMIC_MIN
                            ELSE IF  (TAUROSS(L) .LE. DIST_OUT) THEN
                                DD_VMIC(L) = VMIC_MAX
                            ELSE
                                X = PI * (DIST_IN - TAUROSS(L)) / DX
                                Q = 0.5 + 0.5 * COS(X)
                                DD_VMIC(L) = 
     >                            Q * VMIC_MIN  + (1.-Q) * VMIC_MAX
                            ENDIF
                        ENDDO
C***                3) Interpolation on Radius (analog to last block)
                    ELSE IF (ACTPAR .EQ. 'R1') THEN
                        CALL SARGV (DD_VDOP_LINE, 7, ACTPAR) 
                        IF (ACTPAR .NE. 'R2') GOTO 102
                        CALL SARGV (DD_VDOP_LINE, 6, ACTPAR)
                        READ (ACTPAR, '(F20.0)', ERR = 102) DIST_IN
                        CALL SARGV (DD_VDOP_LINE, 8, ACTPAR)
                        READ (ACTPAR, '(F20.0)', ERR = 102) DIST_OUT
                        IF (DIST_IN > DIST_OUT) THEN
                            DUMMY = DIST_IN
                            DIST_IN = DIST_OUT
                            DIST_OUT = DUMMY
                        ENDIF                    
                        DIST_IN = AMAX1(DIST_IN, RADIUS(ND))
                        DIST_OUT = AMIN1(DIST_OUT, RADIUS(1))
                        WRITE(0,'(A,F10.5, A, F10.5)') 
     >                   "Vmic interpolated between radii", 
     >                   DIST_IN, " and ", DIST_OUT
                        DX = DIST_OUT - DIST_IN
                        DO L=1, ND
                            IF (RADIUS(L) .LE. DIST_IN) THEN
                                DD_VMIC(L) = VMIC_MIN
                            ELSE IF  (RADIUS(L) .GE. DIST_OUT) THEN
                                 DD_VMIC(L) = VMIC_MAX
                            ELSE
                                X = PI * (RADIUS(L) - DIST_IN ) / DX
                                Q = 0.5 + 0.5 * COS(X)
                                DD_VMIC(L) = 
     >                            Q * VMIC_MIN  + (1.-Q) * VMIC_MAX
                            ENDIF
                        ENDDO
                    ELSE
                        GOTO 102
                    ENDIF
C***                End of interpolation branches 
                ELSE
                    GOTO 102
                ENDIF
C***            End of VMICMAX branch
            ELSE
                GOTO 102
            ENDIF
C***        End of depth-dependent VMIC
        ENDIF
C***    VMIC(L) is now defined for every case
C***    VDOP(L,NA) may now be filled via 
C***      VDOP^2(L,NA) = VMIC^2(L) + VTHERM^2(L,NA)
        DO NA=1, NATOM  
           MASS = ATMASS(NA)* AMU
           DO L=1, ND               
C***        v-thermal squared, in km^2/s^2
            THERMVELO_SQRD = 2 * BOLTZK * T(L) / MASS  
            DD_VMIC_SQRD = DD_VMIC(L) * DD_VMIC(L)
            DD_VDOP(L,NA) = SQRT(THERMVELO_SQRD + DD_VMIC_SQRD)
           ENDDO
        ENDDO
*********** END OF VMIC VERSION ************
      ELSE 
C***    VDOP VERSION => No microturbulence or thermal motion 
C***    Note: DD_VDOP(L,NA) does not depend on the element here 
C***          (except for iron), i.e. VDOP(L,1) = VDOP(L,2) = ... 
C***    Note2: This branch is almost completely analog to the 
C***           VMIC branch - (see detailed comments there)

        CALL SARGC (DD_VDOP_LINE, NPAR)
        IF (NPAR .GT. 8) GOTO 101   
        IF (NPAR .LE. 2) THEN
            WRITE(0,'(A)') "VDOP not depth-dependent"
            DO L = 1,ND
                DO NA=1, NATOM
                    DD_VDOP(L, NA) = VDOP
                ENDDO
            ENDDO
        ELSE 
            BDD_VDOP = .TRUE.
            WRITE(0,'(A)') "VDOP is depth-dependent"
            CALL SARGV (DD_VDOP_LINE, 3, ACTPAR)
            IF (ACTPAR .EQ. 'VELOFRAC') THEN
                WRITE(0,'(A)') "VDOP = fraction of wind velocity"
                IF (NPAR .EQ. 3) THEN
                    WRITE(0,'(A,F10.5)') 
     >                 "fraction set to default value:", VMICFRAC_DEFAULT
                    VDOP_FRAC = VMICFRAC_DEFAULT
                ELSE IF (NPAR .EQ. 4) THEN 
                    CALL SARGV (DD_VDOP_LINE, 4, ACTPAR)
                    READ (ACTPAR, '(F20.0)', ERR = 101) VDOP_FRAC
                    WRITE(0,'(A,F10.5)') 
     >                  "fraction set by user to", VDOP_FRAC
                ELSE
                    GOTO 101
                ENDIF
                DO L=1, ND
                    DO NA = 1, NATOM
                        DD_VDOP(L, NA) = AMAX1(VDOP, VDOP_FRAC*VELO(L))
                    ENDDO
                ENDDO
            ELSE IF (ACTPAR .EQ. 'MAX') THEN
                CALL SARGV (DD_VDOP_LINE, 4, ACTPAR)
                READ (ACTPAR, '(F20.0)', ERR = 101) VDOP_MAX
                WRITE(0,'(A,F6.1)') "Maximum VDOP:", VDOP_MAX
                IF (NPAR .EQ. 4) THEN
                    WRITE(0,'(A)') "VDOP = fraction of wind velocity"
                    VDOP_FRAC = VDOP_MAX/VELO(1)
                    WRITE(0,'(A,F10.5)') 
     >                 "fraction implied from maximum VDOP:", VDOP_FRAC
                    DO L=1, ND
                        DO NA = 1, NATOM
                          DD_VDOP(L,NA) = MAX(VDOP, VDOP_FRAC*VELO(L))
                          IF (SYMBOL(NA) .EQ. 'G ') THEN
                            IF (bDDFECONVOL) THEN
                              DD_VDOP(L,NA) = MAX(VDOPFE, DD_VDOP(L,NA))
                            ELSE 
                              DD_VDOP(L,NA) = VDOPFE
                            ENDIF
                          ENDIF
                        ENDDO
                    ENDDO
                ELSE IF (NPAR .EQ. 8) THEN 
                    CALL SARGV (DD_VDOP_LINE, 5, ACTPAR)
                    IF (ACTPAR .EQ. 'VELO1') THEN
                        CALL SARGV (DD_VDOP_LINE, 7, ACTPAR) 
                        IF (ACTPAR .NE. 'VELO2') GOTO 101
                        CALL SARGV (DD_VDOP_LINE, 6, ACTPAR)
                        READ (ACTPAR, '(F20.0)', ERR = 101) DIST_IN
                        CALL SARGV (DD_VDOP_LINE, 8, ACTPAR)
                        READ (ACTPAR, '(F20.0)', ERR = 101) DIST_OUT
                        IF (DIST_IN > DIST_OUT) THEN
                            DUMMY = DIST_IN
                            DIST_IN = DIST_OUT
                            DIST_OUT = DUMMY
                        ENDIF
                        DIST_IN = AMAX1(DIST_IN, VELO(ND))
                        DIST_OUT = AMIN1(DIST_OUT, VELO(1))
                        WRITE(0,'(A,F10.5, A, F10.5)') 
     >                   "VDOP interpolated between wind velocity", 
     >                   DIST_IN, " and ", DIST_OUT
                        DX = DIST_OUT - DIST_IN
                        DO L=1, ND
                            DO NA = 1, NATOM
                                IF (VELO(L) .LE. DIST_IN) THEN
                                    DD_VDOP(L,NA) = VDOP
                                ELSE IF  (VELO(L) .GE. DIST_OUT) THEN
                                    DD_VDOP(L,NA) = VDOP_MAX
                                ELSE
                                    X = PI * (VELO(L) - DIST_IN ) / DX
                                    Q = 0.5 + 0.5 * COS(X)
                                    DD_VDOP(L,NA) = 
     >                                Q * VDOP  + (1.-Q) * VDOP_MAX
                                ENDIF
                            ENDDO
                       ENDDO
                    ELSE IF (ACTPAR .EQ. 'TAU1') THEN
                        CALL SARGV (DD_VDOP_LINE, 7, ACTPAR) 
                        IF (ACTPAR .NE. 'TAU2') GOTO 101
                        CALL SARGV (DD_VDOP_LINE, 6, ACTPAR)
                        READ (ACTPAR, '(F20.0)', ERR = 101) DIST_IN
                        CALL SARGV (DD_VDOP_LINE, 8, ACTPAR)
                        READ (ACTPAR, '(F20.0)', ERR = 101) DIST_OUT
                        IF (DIST_IN < DIST_OUT) THEN
                            DUMMY = DIST_IN
                            DIST_IN = DIST_OUT
                            DIST_OUT = DUMMY
                        ENDIF
                        DIST_IN = AMIN1(DIST_IN, TAUROSS(ND))
                        DIST_OUT = AMAX1(DIST_OUT, TAUROSS(1))
                        WRITE(0,'(A,F10.5, A, F10.5)')
     >                   "VDOP interpolated between Rosseland tau", 
     >                   DIST_IN, " and ", DIST_OUT
                        DX = DIST_IN - DIST_OUT
                        DO L=1, ND
                            DO NA=1, NATOM
                                IF (TAUROSS(L) .GE. DIST_IN) THEN
                                    DD_VDOP(L,NA) = VDOP
                                ELSE IF  (TAUROSS(L) .LE. DIST_OUT) THEN
                                    DD_VDOP(L,NA) = VDOP_MAX
                                ELSE
                                    X = PI * (DIST_IN - TAUROSS(L)) / DX
                                    Q = 0.5 + 0.5 * COS(X)
                                DD_VDOP(L,NA) = 
     >                              Q * VDOP  + (1.-Q) * VDOP_MAX
                                ENDIF
                            ENDDO
                        ENDDO
                    ELSE IF (ACTPAR .EQ. 'R1') THEN
                        CALL SARGV (DD_VDOP_LINE, 7, ACTPAR) 
                        IF (ACTPAR .NE. 'R2') GOTO 101
                        CALL SARGV (DD_VDOP_LINE, 6, ACTPAR)
                        READ (ACTPAR, '(F20.0)', ERR = 101) DIST_IN
                        CALL SARGV (DD_VDOP_LINE, 8, ACTPAR)
                        READ (ACTPAR, '(F20.0)', ERR = 101) DIST_OUT
                        IF (DIST_IN > DIST_OUT) THEN
                            DUMMY = DIST_IN
                            DIST_IN = DIST_OUT
                            DIST_OUT = DUMMY
                        ENDIF
                        DIST_IN = AMAX1(DIST_IN, RADIUS(ND))
                        DIST_OUT = AMIN1(DIST_OUT, RADIUS(1))
                        WRITE(0,'(A,F10.5, A, F10.5)') 
     >                   "VDOP interpolated between radii", 
     >                   DIST_IN, " and ", DIST_OUT
                        DX = DIST_OUT - DIST_IN
                        DO L=1, ND
                            DO NA=1, NATOM
                                IF (RADIUS(L) .LE. DIST_IN) THEN
                                    DD_VDOP(L,NA) = VDOP
                                ELSE IF  (RADIUS(L) .GE. DIST_OUT) THEN
                                    DD_VDOP(L,NA) = VDOP_MAX
                                ELSE
                                    X = PI * (RADIUS(L) - DIST_IN ) / DX
                                    Q = 0.5 + 0.5 * COS(X)
                                    DD_VDOP(L,NA) = 
     >                                 Q * VDOP  + (1.-Q) * VDOP_MAX
                                ENDIF
                            ENDDO
                        ENDDO
                    ELSE
                        GOTO 101
                    ENDIF
                ELSE
                    GOTO 101
                ENDIF
            ELSE
                GOTO 101
            ENDIF
        ENDIF
*********** END OF VDOP VERSION ************
      ENDIF
C*** Generic (Symbol 'G') has its own VDOP velocity VDOPFE from the iron file
C*** if NO-IRONLINES is given in FORMAL_CARDS, the Doppler velocity of FE is set 
C*** to the largest value in the array to avoid redundent increase of resolution
      VDOPNOIRONMAX = MAXVAL(DD_VDOP)
      DO NA=1, NATOM 
        IF (SYMBOL(NA) .EQ. 'G') THEN
          DO L=1, ND               
            IF (BIRONLINES .AND. bDDFECONVOL) THEN
               DD_VDOP(L,NA) = MAX(VDOPFE, DD_VDOP(L,NA))
            ELSEIF (BIRONLINES) THEN
               DD_VDOP(L,NA) = VDOPFE
            ELSE
               DD_VDOP(L,NA) = VDOPNOIRONMAX
            ENDIF
          ENDDO
        ENDIF
      ENDDO

C*** Final VDOP: minimum of matrix DD_VDOP(L,NA) (**Including iron!!)
C*** NOTE1: It could be that the final VDOP is therefore different 
C***        than stated by user!
C*** NOTE2: It could be that the final VDOP is *SMALLER* than stated 
C***        by user if VDOPFE < VDOP_USER!
      VDOP = MINVAL(DD_VDOP(:ND,:NATOM))
      WRITE(0,'(A,F6.1)') "Minimum VDOP:", VDOP
      IF (VDOP .LE. 0) THEN
        WRITE (0,'(A)') 'Invalid VDOP=', VDOP
        STOP '*** FATAL ERROR in subr. VDOP_STRUCT'
      ENDIF

C***  EL_MIN = element with narrowest lines (smallest DD_VDOP)
C***  -> Handed to PRIPRO for output 
      VDOPMIN = MINVAL(DD_VDOP(:ND,1))  
      NA_VDOPMIN = 1
      DO NA=2, NATOM
        IF (MINVAL(DD_VDOP(:ND,NA)) .LT. VDOPMIN) THEN 
           NA_VDOPMIN = NA  
           VDOPMIN = MINVAL(DD_VDOP(:ND,NA))  
        ENDIF
      ENDDO
      EL_MIN = SYMBOL(NA_VDOPMIN)

      WRITE(0,'(A,F6.1,A)') "Resolution corresponds to VDOP= ", VDOP,
     >                      ' (see output file for reason)'

C***  Give a warning if the minimum VDOP is smaller than VDOPFE      
C***  (In this case the iron lines might be broader than intended)
      IF (VDOP < VDOPFE) THEN
        WRITE (0,1) VDOP, VDOPFE
        WRITE (6,1) VDOP, VDOPFE
    1   FORMAT ('*** WARNING: FEDAT resolution is not fine enough.', /,
     >   '*** WARNING: min(VDOP) = 'F7.2, ', VDOPFE = ', F7.2, /, 
     >   '*** WARNING: Iron lines are therefore more Doppler-broadened',
     >   ' than others' )
      ENDIF

C***  Analog arrays in Doppler units
      DO L=1,ND
          DO NA=1,NATOM
                DD_VDOPDU(L,NA) = DD_VDOP(L,NA) / VDOP
                IF (.NOT. BMICROTURB) THEN
                  DD_VMIC(L) = DD_VDOP(L,1)
                ENDIF
                DD_VMICDU = DD_VMIC(L) / VDOP
          ENDDO
      ENDDO

C***  XMAX = minimum integration interval of lines - determined by maximum Doppler velocity!
      XMAX = AMAX1(XMAX, XMAXMIN*MAXVAL(DD_VDOPDU))   
         
C***  IVDOPSTATUS contains a code number, which criterion has caused 
C***  the VDOP that finally defines the wavelength resolution.
C***  This code number will lead to a corresponding message by PRIPRO 
C***  1: VDOP from MODEL file (default)
C***  2: minimum of all Doppler linewidths (accounting for Vmic) 
C***  3: specified as VDOP in FORMAL_CARDS
C***  4: VDOPFE from FEDAT file, when NO-FECONVOL requested
      IF (BMICROTURB) THEN
         IVDOPSTATUS = 2
      ELSEIF (IDX(DD_VDOP_LINE) .GT. 0) THEN
         IVDOPSTATUS = 3
      ENDIF
C***  minimum VDOP is not from ion, if NOIRINLINES or VDOPFE not minimum
      IF (BIRONLINES .AND. VDOPFE .LE. VDOP 
     >          .AND. .NOT. BDDFECONVOL) THEN      
         IVDOPSTATUS = 4
      ENDIF 

      RETURN
      
C***  ERROR branches *********************************************


 101  WRITE(0,'(A)') '*** ERROR when decoding parameter ' // ACTPAR
      WRITE(0,'(A)') '*** The error occured in the following line:'
      WRITE(0,'(A)') DD_VDOP_LINE
      WRITE(0,'(A)') 'Version not known! possible versions are:' 
      WRITE(0,'(A)') 'VDOP X.X'
      WRITE(0,'(A)') 'VMIC X.X VELOFRAC [X.X]' 
      WRITE(0,'(A)') 'VDOP X.X MAX X.X ' // 
     >          '[TAU1 | R1 | VELO1 X.X TAU2 | R2 | VELO2 X.X]'
      STOP '*** Fatal error in vdop_struct'
      
 102  WRITE (0,'(A)') '*** ERROR when decoding parameter ' // ACTPAR
      WRITE (0,'(A)') '*** The error occured in the following line:'
      WRITE (0,'(A)') DD_VDOP_LINE
      WRITE(0,'(A)') 'Version not known! possible versions are:' 
      WRITE(0,'(A)') 'VMIC X.X'    
      WRITE(0,'(A)') 'VMIC X.X VELOFRAC [X.X]' 
      WRITE(0,'(A)') 'VMIC X.X MAX X.X ' // 
     >       '[TAU1 | R1 |VELO1 X.X TAU2 | R2 | VELO2 X.X]'
      STOP '*** Fatal Error in vdop_struct'

      END


