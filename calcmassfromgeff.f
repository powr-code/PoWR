      SUBROUTINE CALCMASSFROMGEFF(bFixGEFF, bFULLHYDROSTAT,
     >                            XMSTAR, GLOG, GEFFLOG, 
     >                            ARAD, APRESS, AGRAV, RADIUS, ND,
     >                            RSTAR, RCON, RI, TAUROSS,
     >                            GAMMARADMEAN, GEDD,
     >                            XLOGL, QIONMEAN, fGAMMACOR)
C***********************************************************************
C***  updates GLOG and XMSTAR if GEFF has been fixed in the CARDS file
C***    and corrects GAMMARADMEAN in the case of full hydrostatic int.
C***  called from ENSURETAUMAX
C***********************************************************************
      IMPLICIT NONE
      
      INTEGER, INTENT(IN) :: ND
      
      REAL, INTENT(IN) :: GEFFLOG, RCON, RSTAR, XLOGL, QIONMEAN
      REAL, INTENT(INOUT) :: XMSTAR, GLOG, GEDD, GAMMARADMEAN
C***  fGAMMACOR is the calculated correction factor for depth-dependend GAMMARAD values
      REAL, INTENT(OUT) :: fGAMMACOR  
      
      LOGICAL, INTENT(IN) :: bFixGEFF, bFULLHYDROSTAT
      
      REAL, DIMENSION(ND) :: RADIUS, TAUROSS
      REAL, DIMENSION(ND-1) :: ARAD, AGRAV, APRESS, RI
            
      
      REAL :: DIFF1, DIFF2, XLSTARS, facMASS, GTRUE
      INTEGER :: L, LCON, Lfm

      !Physical constants:
      REAL, PARAMETER :: RSUN = 6.96E10        !Solar Radius (cm)
      REAL, PARAMETER :: GCONST = 6.670E-8     !GRAVITATION CONSTANT (CGS UNITS)
      REAL, PARAMETER :: XMSUN = 1.989E33      !XMSUN = Solar Mass (g)
      
      !File and channel handles (=KANAL)
      INTEGER, PARAMETER :: hOUT = 6        !write to wruniqX.out (stdout)
      INTEGER, PARAMETER :: hCPR = 0        !write to wruniqX.cpr (stderr)
      
      IF (.NOT. bFixGEFF) THEN
        WRITE (hCPR,'(A)')
     >    'WARNING: MASS UPDATE CALLED EVEN THOUGH MASS IS FIXED '
     >    // ' --  UPDATE HAS BEEN IGNORED!'
        fGAMMACOR = 1.
        RETURN
      ENDIF
                  
      XLSTARS = 10**(XLOGL)

      Lfm = ND
      IF (bFixGEFF .AND. bFULLHYDROSTAT) THEN
        LCON = ND
        DO WHILE (RADIUS(LCON) < RCON) 
          LCON = LCON - 1
        ENDDO
        
        !Calculate GAMMARADMEAN without limitation here:
        CALL CALCGAMMARADMEAN(ARAD, AGRAV, RADIUS, TAUROSS, 
     >                        ND, RCON, RI, GAMMARADMEAN)
        
        !We need to update the mass, so we need to find a point where
        ! GAMMARADMEAN is realized
        DIFF1 = ARAD(LCON)/AGRAV(LCON) - GAMMARADMEAN
        DO L=LCON+1, ND-1
          DIFF2 = ARAD(L)/AGRAV(L) - GAMMARADMEAN
          IF (DIFF1 * DIFF2 <= 0.) THEN
            !point found => calculate the factor that must be
            ! multiplied to AGRAV in order to fullfill the hydro equation
            ! (As Rstar is fixed, this is a correction factor to the mass)
            facMass = (ARAD(L)+APRESS(L))/AGRAV(L)
            LfM = L
            EXIT
          ENDIF
          DIFF1 = DIFF2
        ENDDO
      ENDIF
      GAMMARADMEAN = MIN(GAMMARADMEAN, 0.9) !now limit the actual value

      IF (bFULLHYDROSTAT) THEN
          !For the calculation of the mass we need to select one GAMMA => use mean value         
          !The Mass factor is needed to have a proper value for XMSTAR and GLOG
          ! as for those cases where GEFF is fixed the main iteration does
          ! not really react on changed GAMMA values and therefore GAMMARAD in
          ! the next iteration would be either way too large or too small.
C          GEDD = 1. - 1. / facMass + GAMMARADMEAN / facMass
      
          GTRUE = 10**GEFFLOG + ARAD(LfM)
          facMass = GTRUE / 10**(GLOG)
          GEDD = ARAD(LfM) / GTRUE
          GLOG = LOG10(GTRUE)
          GEDD = MIN(GEDD, 0.9)     !now limit the actual value
          fGAMMACOR = GEDD / GAMMARADMEAN
C          WRITE (hCPR,*) 'LfM, fMass, fGAMMACOR ', LfM, facMass, fGAMMACOR
          GAMMARADMEAN = GEDD
          XMSTAR = 10**GLOG * RSTAR**2 / XMSUN / GCONST          
      ELSE
          XMSTAR = 10.**( GEFFLOG - 4.4371 ) * (RSTAR/RSUN)**(2.)
     >                 + 10.**(-4.51) * QIONMEAN * XLSTARS
          GLOG = LOG10( GCONST * XMSTAR * XMSUN / (RSTAR * RSTAR) )
          GEDD = 1. - 10.**( GEFFLOG - GLOG )
      ENDIF        
      
      
      END
