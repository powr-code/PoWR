      SUBROUTINE PREP_GAMMARAD (bOLDRAD, bFULLHYDROSTAT, GEDD, 
     >      IGEddFix, GAMMARAD, NDold, ARAD, GLOG, RSTAR,
     >      RADIUSold, RSTARold, XMGold, TAURold, RCONold, GEFFLOG,
     >      STAPEL, ATMEAN, XLOGL, XMSTAR, RADGAMMASTART, 
     >      GEDDRAD, bOldStratification, bGAMMARADMEAN, bSaveGEFF)
C***********************************************************************
C***  This subroutine prepares GAMMARAD and related quentities for
C***  the hydrostatic equation
C***  Called from: WRSTART
C***********************************************************************

      IMPLICIT NONE
      
      INTEGER, INTENT(IN) :: NDold, IGEddFix
      REAL, INTENT(IN) :: XMGold, XLOGL, STAPEL, ATMEAN,
     >                    RSTARold, RCONold, RADGAMMASTART
      
      REAL, DIMENSION(NDold), INTENT(IN) :: RADIUSold, TAURold, ARAD
      REAL, DIMENSION(NDold), INTENT(OUT) :: GAMMARAD
      REAL, DIMENSION(NDold) :: ARADL 

      REAL, INTENT(INOUT) :: GEDD, GEFFLOG, GLOG, XMSTAR
      REAL, INTENT(OUT) :: GEDDRAD

      LOGICAL, INTENT(IN) :: bOLDRAD, bFULLHYDROSTAT, bSaveGEFF,
     >                       bOldStratification, bGAMMARADMEAN

      INTEGER :: L
      REAL :: GEDDRADMEAN, TAUNORM, WTAU, RSTAR

      !Physical constants
      REAL, PARAMETER :: GCONST = 6.670E-8  !GRAVITATION CONSTANT (CGS UNITS)
      REAL, PARAMETER :: XMSUN = 1.989E33   !XMSUN = Solar Mass (g)
      
      !File and channel handles (=KANAL)
      INTEGER, PARAMETER :: hOUT = 6        !write to wruniqX.out (stdout)
      INTEGER, PARAMETER :: hCPR = 0        !write to wruniqX.cpr (stderr)

      GEDDRAD = 0.      
      
      IF (bOLDRAD .AND. bFULLHYDROSTAT) THEN
        !calculate effective g via ARAD from OLD MODEL

        IF (GEDD > .0 .AND. IGEddFix == 2) THEN
          !Gamma value has been fixed via CARDS 
          DO L=1, NDold
            GAMMARAD(L) = GEDD
          ENDDO

        ELSE

C***      Interpolate ARAD from interstices to ARADL on full depth points  
          DO L=1, NDold
            IF (L==1) THEN
              ARADL(L) = ARAD(1)
            ELSEIF (L==NDold) THEN
              ARADL(L) = ARAD(NDold-1)
            ELSE
              ARADL(L) = 0.5 * (ARAD(L-1) + ARAD(L))
            ENDIF
            GAMMARAD(L) = ARADL(L)*(RADIUSold(L)*RSTARold)**2/XMGold
C***        restrict GAMMARAD 
            GAMMARAD(L) = MIN(GAMMARAD(L), 0.9)
          ENDDO

          !Calculate mean value, needed for g-g_eff-Relation
          ! if RADGAMMASTART: OLD is used (cannot be done in DECSTAR)
          GEDDRADMEAN = GAMMARAD(NDold) * EXP(-TAURold(NDold))
          TAUNORM = 0.
          DO L=NDold-1, 1, -1
            IF (RADIUSold(L) > RCONold) EXIT  
            WTAU = (TAURold(L+1) - TAURold(L)) * EXP(-TAURold(L))
            GEDDRADMEAN = GEDDRADMEAN + GAMMARAD(L) * WTAU
            TAUNORM = TAUNORM + WTAU
          ENDDO
          IF (TAUNORM > 0.) THEN
            GEDDRADMEAN = GEDDRADMEAN / TAUNORM
          ENDIF
          IF (bSaveGEFF) THEN
            GLOG = LOG10((10.**GEFFLOG)/(1.-GEDDRADMEAN))
            XMSTAR = 10.**GLOG * RSTAR * RSTAR / GCONST / XMSUN
          ELSE
            GEFFLOG = LOG10((10.**GLOG)*(1.-GEDDRADMEAN))
          ENDIF
          
          IF (bGAMMARADMEAN) THEN
            !Use mean value of the inner part instead of individual values (CARDS option)
            DO L=1, NDold
              GAMMARAD(L) = GEDDRADMEAN
            ENDDO
          ENDIF
        ENDIF

      ELSEIF (GEDD > 0. .OR. RADGAMMASTART >= 0.) THEN
C***    for new models only g_THOMSON is available, 
C***    unless RADGAMMA-START is specified
        IF (IGEddFix == 0) THEN        
          IF (RADGAMMASTART > 0) THEN
            GEDD = RADGAMMASTART
          ELSE
            GEDD = 10.**(-4.51) * (STAPEL/ATMEAN) *(10.**XLOGL) /XMSTAR
          ENDIF
          IF (bSaveGEFF) THEN
            GLOG = LOG10((10.**GEFFLOG)/(1.-GEDD))
            XMSTAR = 10.**GLOG * RSTAR * RSTAR / GCONST / XMSUN
          ELSE
            GEFFLOG = ALOG10( (10**GLOG) * (1. - GEDD) ) 
          ENDIF
        ENDIF
        IF (RADGAMMASTART >= 0.) THEN
          DO L=1, NDold
            GAMMARAD(L) = RADGAMMASTART
          ENDDO
        ENDIF
      ENDIF

      IF (bFULLHYDROSTAT) THEN
        IF (bOLDRAD) THEN
          GEDDRAD = GEDDRADMEAN
        ELSEIF (RADGAMMASTART >= 0.) THEN
          GEDDRAD = RADGAMMASTART
        ELSE
          GEDDRAD = GEDD      !poor start with Thomson only
        ENDIF        
      ENDIF
      
C***  Only if not copying the old stratification
      IF (.NOT. bOldStratification) THEN
        IF (bFULLHYDROSTAT) THEN
          WRITE (hCPR,FMT='(A,/,A,F6.3,2(A,F10.5))')
     >      'Values used for hydrostatic domain: ',
     >      '   Full Gamma (ND) = ', GEDDRAD,
     >      '   log g_grav = ', GLOG,
     >      '   log g_eff = ', GEFFLOG

        ELSE
          WRITE (hCPR,FMT='(A,/,A,F6.3,2(A,F10.5))')
     >      'Values used for hydrostatic domain: ',
     >      '   Eddington Gamma = ', GEDD,
     >      '   log g_grav = ', GLOG,
     >      '   log g_eff = ', GEFFLOG
        ENDIF
      ENDIF
      
      RETURN 
      END
