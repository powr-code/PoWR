      FUNCTION STARKHEII (WAVE, WAVEH, T, XNE)
C**********************************************************************
C***  Called from: SUBROUTINE STARKBROAD
C***  Stark-Broadening for helium (He II) lines
C***  The data table for the current line must have been loaded before 
C***     to the common block VCSSBDAT by Subr. STARKHEIIPREP
C***  Normalization is arbitrary
C**********************************************************************

      PARAMETER (MPDIM     = 50)     ! Max. # OF PROFILE POINTS IN vcssb TABLE
      PARAMETER (MTEMPDIM  =  6)     ! Max. # OF TEMP POINTS IN vcssb TABLE
      PARAMETER (MNEDIM    = 11)     ! Max. # OF DENSITY POINTS IN vcssb TABLE

      COMMON /VCSSBDAT/ MP, MNE, MTEMP, ALOG_NEL(MNEDIM),
     >                   ALOG_T0,     ALOG_T_INC,
     >                   ALOG_ALPHA0, ALOG_ALPHA_INC,
     >                   SVCS(MPDIM,MTEMPDIM,MNEDIM)

C***  Definition of DEX function
      EXP10( X ) = EXP( 2.30258509299405 * X )

      PHI = .0

      ALOG_T = ALOG10( T )
      XNE0    = EXP10( ALOG_NEL(1) )

C***  Interpolation interval, weights in temperature grid
      BTEMP = ( ALOG_T - ALOG_T0 ) / ALOG_T_INC + 1.
      ITEMP = MAX0 ( MIN0(INT(BTEMP), MTEMP-1), 1 )
      WTTEMP = BTEMP - FLOAT(ITEMP)

C***  Interpolation interval, weights in electron-density grid
      ALOG_NE = ALOG10( XNE )
      DO INE = 1, MNE
         IF (ALOG_NEL(INE) .GT. ALOG_NE) EXIT
         INELAST = INE
      END DO

      INE = INELAST - 1
      INE = MAX0 (MIN0(INE, MNE-1), 1)
      WTNE = (ALOG_NE - ALOG_NEL(INE))/(ALOG_NEL(INE+1) - ALOG_NEL(INE))

      IF( WTNE .LT. 0. ) WTNE = 0.    ! NO EXTRAPOLATION
      IF( WTNE .GT. 1. ) WTNE = 1.

C***  Find wavelength interval index: (ialp, ialp+1)
      DL = ABS(WAVE - WAVEH)
C***  Center point
      IF (DL .LT. 1.E-10) THEN
         IALP = 1
         WT   = .0
      ELSE
         F0 = 1.25E-9 * AMAX1(XNE,XNE0)**(2./3.)   ! FIX AT TABLE LIMIT
         ALOG_ALPHA = ALOG10( DL/F0 )
         BALP = (ALOG_ALPHA - ALOG_ALPHA0) / ALOG_ALPHA_INC + 1.
         IALP = INT(BALP)
         IALP = MAX0 (IALP,1)
         WT = BALP - FLOAT(IALP)
      ENDIF

C***  Set profile to zero outside the tabulated wavelength range
      IF (IALP .GE. MP) GOTO 10

      PROFI = (1. - WTNE) * (1. - WTTEMP) * SVCS(IALP,ITEMP  ,INE  ) +
     +        (1. - WTNE) * WTTEMP        * SVCS(IALP,ITEMP+1,INE  ) +
     +           WTNE     * (1. - WTTEMP) * SVCS(IALP,ITEMP  ,INE+1) +
     +           WTNE     * WTTEMP        * SVCS(IALP,ITEMP+1,INE+1)

      PROFIP= (1. - WTNE) * (1. - WTTEMP) * SVCS(IALP+1,ITEMP  ,INE  ) +
     +        (1. - WTNE) * WTTEMP        * SVCS(IALP+1,ITEMP+1,INE  ) +       
     +           WTNE     * (1. - WTTEMP) * SVCS(IALP+1,ITEMP  ,INE+1) +
     +           WTNE     * WTTEMP        * SVCS(IALP+1,ITEMP+1,INE+1)

      PROF = (1. - WT) * PROFI + WT * PROFIP


      PHI = EXP10(PROF)


   10 STARKHEII = PHI
      RETURN
      END

