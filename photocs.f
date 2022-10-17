      SUBROUTINE PHOTOCS (SIGMA,SIGMATH,EDGE,WAVENUM,ALPHA,SEXPO,
     $                    ADDCON1, ADDCON2, ADDCON3, 
     $                    IGAUNT,KON)
C***********************************************************************
C***  CALCULATES SIGMA(NUE), THE FREQUENCY DEPENDENT PHOTO CROSS SECTION
C***  THIS ROUTINE IS ONLY CALLED FROM BFCROSS, COOP AND CMFCOOP
C***********************************************************************

      DIMENSION ALPHA(1),SEXPO(1)
      CHARACTER*8 IGAUNT(1)
      DIMENSION ADDCON1(1), ADDCON2(1), ADDCON3(1)
C***  THE FOLLOWING DATA ARE FOR MIHALAS' GAUNT FACTOR FIT ( HY AND HE II, N=1) 
      DATA A0,A1,A2,A3,AM1,AM2 /
     > 1.2302628,  -3.1927214E-2, 8.9105122E-4, -1.1544111E-5,
     > -0.50812150, 0.10631895 /
      DATA IWARN_PIKB12 / 0 / 
 
      X=EDGE/WAVENUM
      XINV=1./X

C***  The calling programms should not call PHOTOCS beyond the ionisation edge,
C***  althogh this is not a catastrophy. 
C***  The followong warning might be switched off if annoying.
      IF (WAVENUM .LT. EDGE) THEN
         WRITE (0,'(a)') '*** WARNING: PHOTOCS CALLED OUTSIDE EDGE FREQUENCY'
         SIGMA = .0
         RETURN
      ENDIF


C***  VARIABLE IGAUNT IS MISUSED TO CARRY THE KEYWORD FOR NEW
C***  PHOTOIONIZATION CROSS SECTIONS: 
C***  'KOESTER': KOESTER ET AL. 1985, A+A 149, 423
C***   =======
C***             FIT COEFFICIENTS FOR THE MODIFIED (]) FORMULA: 
C***                     SIGMATH = F(A0)
C***                  ALPHA(KON) = A1
C***                  SEXPO(KON) = A2
      IF (IGAUNT(KON) .EQ. 'KOESTER') THEN
          XLN=ALOG(1.E8/WAVENUM)
          XLN2=XLN*XLN
          X0LN=ALOG(1.E8/EDGE)
          X0LN2=X0LN*X0LN
          SIGMA=SIGMATH*X**ALPHA(KON)*EXP(SEXPO(KON)*(XLN2-X0LN2))
C***  'BUTLER..': K. BUTLER (MUNICH), PRIVATE COMMUNICATION
C***   ========
C***  '......12': FORMULA 12 FROM PROGRAM DETAIL
C***              MODIFIED FORM OF THE SEATON FORMULA
      ELSE IF (IGAUNT(KON) .EQ. 'BUTLER12') THEN
          SIGMA=SIGMATH*X**(ALPHA(KON)+SEXPO(KON)*ALOG(X))
C***  'DETAILN3': EXTENDED VERSION (6 COEFFICIENTS!!!) OF FORMULA 12
C***   ========   (SEE "BUTLER12")
      ELSE IF (IGAUNT(KON) .EQ. 'DETAILN3') THEN
          NIII4NO=INT(SEXPO(KON))
C***      PREVENT LARGE PHOTOIONIZATION CROSS SECTIONS FOR SPECIFIED 
C***      N III LEVELS AT X-RAY FREQUENCIES (XLAMBDA <= 30A)
          IF ((NIII4NO .EQ. 4 .OR. NIII4NO .EQ. 5 .OR. NIII4NO .EQ. 6)
     $        .AND. (X .LT. 0.05)) THEN
              SIGMA=0.0
          ELSE
              SUMI36=PHOTON3(X,NIII4NO)
              SIGMA=SIGMATH*X**(ALPHA(KON)+SUMI36)
          ENDIF

      ELSE IF (IGAUNT(KON) .EQ. 'PIKB12  ') THEN
C***          XLN = ALOG(X)
C***          SIGMA = SIGMATH*X**(ALPHA(KON)+XLN*(SEXPO(KON)+XLN*(
C***     >            ADDCON1(KON)+XLN*(ADDCON2(KON)+XLN*ADDCON3(KON)))))
C***      changed by wrh 14-Mar-2005 11:43:35
C***      Corresponding data are in our files for C II and C III
C***      The origin of this fit-formula is unclear. It is rubbish:
C***      The exponent of X:=nue_o/nue is a polynomial of ln(x). 
C***      Taking ln of the whole formula, gives:
C***      ln(SIGMA) = ln (SIGMATH) + ALPHA*ln(X) + SEXPO*(ln(X))**2
C***                               + ADDCON1*(ln(X))**3 ...
C***      As X approaches zero for high frequencies, 
C***      ln(x) becomes -infinity, and the exponent goes towards plus 
C***      infinity when the highest non-zero coefficient (mostly ADDCON1) 
C***      is NEGATIVE - which is the case in all (but one) data!
C***      Therefore, I replace this formula by a hydrogenic slope. 
C***      Corresponding data should be avoided. Therefore, a warning is 
C***      issued when PIKB12 is used for the first time. 
              SIGMA = SIGMATH * X**3
              IF (IWARN_PIKB12 .EQ. 0) THEN
                 IWARN_PIKB12= IWARN_PIKB12 + 1
                 WRITE (0, *) '*** WARNING issued from PHOTOCS:'
                 WRITE (0, *) '*** Obsolete Formula PIKB12 ' 
     >            //'replaced by hydrogenic slope ***'
              ENDIF

C***  'OPAPROIX': ATOMIC DATA FOR OPACITY CALCULATIONS (OPACITY PROJECT):
C***  ==========  IX. The lithium isoelectronic sequence
C***              (Peach, Saraph & Seaton 1988, 
C***               J. Phys. B: At. Mol. Opt. Phys. 21, 3669)
      ELSE IF (IGAUNT(KON) .EQ. 'OPAPROIX') THEN
          XLOG = ALOG10(XINV)
          SIGMA = SIGMATH*XINV**(ALPHA(KON)+XLOG*(SEXPO(KON)+
     +                                      XLOG*ADDCON1(KON)))
      ELSE
C***  DEFAULT - OLD VERSION: SEATON / HYDROGENIC
C***  =======
          IF (ALPHA(KON) .NE. .0) THEN
C***      ALPHA(KON) DEFINED: SEATON FORMULA
            SIGMA=SIGMATH * X**SEXPO(KON)*(ALPHA(KON)+(1.-ALPHA(KON))*X)
          ELSE
C***      ALPHA(KON) NOT DEFINED: HYDROGENIC EXPONENT NUE**(-3)
            SIGMA=SIGMATH*X*X*X
          ENDIF
      ENDIF
 
C***  BOUND-FREE GAUNT FACTORS ARE CALCULATED DEPENDING ON KEYWORD IGAUNT
 
C***  POLYNOMIAL FIT FOR GII(N=1) FROM MIHALAS 1967 APJ 149,P187
C***  THIS MIHALAS GAUNT FACTOR IS ONLY VALID FOR GROUND STATE N=1 ]
      IF (IGAUNT(KON) .EQ. 'MIHALAS' ) THEN
            IF (X .GT. 0.055) THEN
C***           GOOD ACCURACY ONLY FOR WAVELENGTHS (OF HYDROGEN!) FROM
C***           THRESHOLD TO 50 A:
               GAUNT=A0+X*(AM1+X*AM2)+(A1+(A2+A3*XINV)*XINV)*XINV
            ELSE
               GAUNT=1.
            ENDIF
            SIGMA=SIGMA*GAUNT
      ENDIF
 
C***  GII FIT FROM SEATON (1960), REP. PROG. PHYS. 23, P. 313
C***  THIS FORMULA IS VALID FOR ALL HYDROGENIC LEVELS,
C***  I.E.  H I, HE II, ..., C IV, N V, O VI
C***  VARIABLE SEXPO IS MISUSED TO CARRY THE MAIN QUANTUM NUMBER
      IF (IGAUNT(KON) .EQ. 'SEATON' ) THEN
            IF (SEXPO(KON) .LE. .0) THEN
               CALL REMARK ('MAIN QUANTUM NUMBER UNDEFINED')
               STOP 'ERROR IN SUBR. PHOTOCS : SEXPO0'
               ENDIF
            U=XINV - 1.
            DEN=(X/SEXPO(KON))**0.666666666666
            GAUNT=1. + 0.1728 * (U-1.) * DEN -
     -            0.0496 * (U*(U+1.333333333333)+1.) * DEN * DEN
            SIGMA=SIGMA*GAUNT
            ENDIF
 
C***  PREVENT NEGATIVE PHOTOIONIZATION CROSS SECTIONS:
      IF (SIGMA .LT. 0.0)  SIGMA = 0.0

C***  PREVENT LARGE PHOTOIONIZATION CROSS SECTIONS:
      IF (SIGMA .GT. 10.*SIGMATH)  THEN
         WRITE (0,*) 'WARNING : VERY HIGH PHOTOIONISATION CROSS ',
     >               'SECTION DETECTED; SET TO ZERO'

cc        write (0,'(a,2(e15.8,2x))') 'sigma, sigmath=', sigma, sigmath
cc        write (0,'(a,a8)') 'igaunt=',igaunt(kon)
cc        write (0,*) 'kon, wavenum=',kon, wavenum
cc        write (0,*) 'alpha, sexpo...=', alpha(KON), sexpo(KON), 
cc     >              addcon1(KON), addcon2(KON), addcon3(KON)
cc        write (0,*) 'ERROR IN SUBR. PHOTOCS : SIGMA'
cc        STOP 'ERROR IN SUBR. PHOTOCS : SIGMA'

         SIGMA = 0.

      ENDIF

      RETURN
      END
