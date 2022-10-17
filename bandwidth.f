      SUBROUTINE BANDWIDTH (ND, RADIUS, OPAC, OPAL, LINPRO, AVOIGT, 
     >                      NDDIM,NBL, MAXLAP, XMAX, XMAXBROAD, XMAXLIN,
     >                      PHITAB, NFDIMPHITAB, 
     >                      DXMAX, LEVELNUP, LEVELLOW,
     >                      BDD_VDOP, GRIEMPAR, VDOP, ALN, bHYDROGEN,
     >                      TAUMAX, TAUMINBROAD)
C***********************************************************************
C***  Called from: STARKBROAD
C***  This subroutine estimates the necessary line bandwidth
C***  such that the *static* optical depth drops under a specified limit 
C***********************************************************************
      INTEGER, INTENT(IN) :: NFDIMPHITAB, ND, NBL, MAXLAP
      REAL, INTENT(IN) :: TAUMINBROAD

      REAL, DIMENSION(-NFDIMPHITAB:NFDIMPHITAB, ND) :: PHITAB
      REAL, DIMENSION(MAXLAP,NDDIM) :: AVOIGT, GRIEMPAR
      REAL, DIMENSION(ND) :: RADIUS, OPAC, OPAL
      CHARACTER(8) :: LINPRO
      CHARACTER(10) :: LEVELNUP, LEVELLOW
      LOGICAL :: BDD_VDOP, bHYDROGEN

      REAL, EXTERNAL :: PHIHOLTSMARK
      
      REAL, PARAMETER :: CLIGHT  = 2.99792458E5    ! c in km / s
      
C***  Optical Depth limit: the smaller this value, the larger the bandwidth 
C***  (Jan 2017: Default now in formal.f, can be adjusted via FORMAL_CARDS)
      TAUCRIT = TAUMINBROAD

C***  Branch for pressure-broadening  (depth-dependent VOIGT profiles) 
      IF (     LINPRO .EQ. 'BRD-HeI ' .OR.
     >         LINPRO .EQ. 'VOIGT   ' .OR.
     >         LINPRO .EQ. 'Q-STARK ') THEN     

C***     TAU = integral OPAL * PHI(X) * dr
C***     the Voigt profile is radius-dependent (i.e. inside the sum)

         X = XMAX

    3    CONTINUE
         PHI = VOIGTH(AVOIGT(NBL,1), X)
         TAU  = 0.5 * (RADIUS(1) - RADIUS(2)) * OPAL(1) * PHI
         TAUC = 0.5 * (RADIUS(1) - RADIUS(2)) * OPAC(1)

         DO L=2, ND-1
            PHI = VOIGTH(AVOIGT(NBL,L), X)
            TAU = TAU + 0.5 * 
     >            (RADIUS(L-1) - RADIUS(L+1)) * OPAL(L) * PHI
            TAUC = TAUC + 0.5 * 
     >            (RADIUS(L-1) - RADIUS(L+1)) * OPAC(L)
C***        Integrate only up to a point 
C***        where the continuum opacity exceeds the specified TAUMAX     
            IF (TAUC > TAUMAX) EXIT
         ENDDO

         IF (TAUC < TAUMAX) THEN
           PHI = VOIGTH(AVOIGT(NBL,ND), X) 
           TND = 0.5 * (RADIUS(ND-1) - RADIUS(ND)) * OPAL(ND) * PHI
           TAU = TAU + TND
         ENDIF
         
         IF (TAU .GT. TAUCRIT) THEN
            X = X + DXMAX
            GOTO 3
         ENDIF

         XMAXLIN = MAX(XMAX,X)
         XMAXBROAD = MAX(XMAXBROAD,X)

C***  Branch for pressure-broadening  (depth-dependent Holtsmark profiles) 
      ELSE IF (LINPRO == 'L-STARK ') THEN     
         
         X = XMAX
         DLAMDOPREL = (VDOP / CLIGHT)
         
    4    CONTINUE
         
C***     Approximative treatment for linear stark broadening of hydrogenic ions
C***     using a joined profile function of a Doppler and a Holtsmark profile
C***     (Note: This is inferior to precise tables!)
C***     Details:
C***     - The parameter GRIEMPAR has been precalculated in LINSTARK
C***     - Normalized profile function is calculated in PHIHOLTSMARK
C***     - BETA and BETADOP need to be provided in Angstroem/cm, which is
C***         guaranteed by the definition of GRIEMPAR in LINSTARK
         BETA = GRIEMPAR(NBL,1) * (EXP(ALN*X) - 1.)
         BETADOP = GRIEMPAR(NBL,1) * DLAMDOPREL
         PHI = PHIHOLTSMARK(BETA, BETADOP, bHYDROGEN)
                  
         TAU  = 0.5 * (RADIUS(1)   - RADIUS( 2)) * OPAL( 1) * PHI
         TAUC = 0.5 * (RADIUS(1)   - RADIUS( 2)) * OPAC( 1)

         DO L=2, ND-1
            BETA = GRIEMPAR(NBL,L) * (EXP(ALN*X) - 1.)
            BETADOP = GRIEMPAR(NBL,L) * DLAMDOPREL
            PHI = PHIHOLTSMARK(BETA, BETADOP, bHYDROGEN)
            TAU = TAU + 0.5 * 
     >            (RADIUS(L-1) - RADIUS(L+1)) * OPAL(L) * PHI
C***        Integrate only up to a point 
C***        where the continuum opacity exceeds the specified TAUMAX     
            TAUC = TAUC + 0.5 * 
     >            (RADIUS(L-1) - RADIUS(L+1)) * OPAL(L) * PHI
            IF (TAUC > TAUMAX) EXIT
         ENDDO

         IF (TAUC < TAUMAX) THEN
            BETA = GRIEMPAR(NBL,ND) * (EXP(ALN*X) - 1.)
            BETADOP = GRIEMPAR(NBL,ND) * DLAMDOPREL
            PHI = PHIHOLTSMARK(BETA, BETADOP, bHYDROGEN)
            TND = 0.5 * (RADIUS(ND-1) - RADIUS(ND)) * OPAL(ND) * PHI
            TAU = TAU + TND
         ENDIF
         
         IF (TAU > TAUCRIT) THEN
            X = X + DXMAX
            IF (X < 1./DLAMDOPREL) THEN
C***          Increase BANDWIDTH if X < C/VDOP
              GOTO 4
            ELSE 
C***          Cut profile, but write warning into CPR and OUT            
              WRITE (0,'(4A)') 'WARNING: BANDWIDTH CANNOT'
     >          //  ' FULLY COVER L-STARK PROFILE FOR ',              
     >                    LEVELNUP, ' - ', LEVELLOW
              WRITE (*,'(4A)') 'WARNING: BANDWIDTH CANNOT'
     >          //  ' FULLY COVER L-STARK PROFILE FOR ',              
     >                    LEVELNUP, ' - ', LEVELLOW
            ENDIF
         ENDIF
         
         XMAXLIN = MAX(XMAX,X)
         XMAXBROAD = MAX(XMAXBROAD,X)
      
C***  Branch for tabulated STARK profiles (H I and He II)
      ELSE IF (LINPRO(:3) .EQ. 'BRD') THEN

         KMIN = INT(XMAX/DXMAX)

C***     Note: only the bue wing is tested
         DO K=KMIN, NFDIMPHITAB
C***        TAU = integral OPAL * PHI(X,r) * dr
            X = K * DXMAX

            TAU = .5 * (RADIUS(1) - RADIUS(2)) * OPAL(1) * PHITAB(K,1)
            TAUC= .5 * (RADIUS(1) - RADIUS(2)) * OPAC(1)
            DO L=2, ND-1
               TAU = TAU + 0.5 * (RADIUS(L-1) - RADIUS(L+1)) * OPAL(L)
     >             * PHITAB(K,L)
               TAUC= TAUC + 0.5 * (RADIUS(L-1) - RADIUS(L+1)) * OPAC(L)
               IF (TAUC > TAUMAX) EXIT
            ENDDO
            IF (TAUC < TAUMAX)
     >        TAU = TAU +
     >        .5 * (RADIUS(ND-1) - RADIUS(ND)) * OPAL(ND) * PHITAB(K,ND)
            IF (TAU .LT. TAUCRIT) GOTO 2

         ENDDO

C***     Print a warning if TAU is still larger than TAUCRIT at the end 
C***       of the tabulated profile
         WRITE (0,'(5A)') 'WARNING: PROFILE NOT FULLY COVERED BY ',
     >                    'TABLE DIM. NFDIMPHITAB: ', 
     >                    LEVELNUP, ' - ', LEVELLOW

C***     XMAX is raised to the frequecy which was required 
C***       to reach TAU < TAUCRIT (at max NFDIMPHITAB*DXMAX) 
    2    XMAXLIN = MAX(XMAX,X)
         XMAXBROAD = MAX(XMAXBROAD,X)

C***  Error branch
      ELSE
         WRITE (0,*) '*** UNKNOWN PROFILE TYPE: ', LINPRO
         STOP ' *** FATAL ERROR IN SUBROUTINE BANDWIDTH'
      ENDIF


      RETURN
      END

