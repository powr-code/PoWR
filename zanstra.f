      SUBROUTINE ZANSTRA (NF, EMFLUX, XLAMBDA, FWEIGHT, TEFF, RSTAR )
C***********************************************************************
C***   Printout of UV Continuum Photons and Zanstra Temperatures
C***********************************************************************
 
      DIMENSION XLAMBDA(NF),EMFLUX(NF),FWEIGHT(NF)

      PARAMETER (NPH = 8)
      DIMENSION PHLAM(NPH)
      CHARACTER*30 CPHLAM(NPH)
 
      DATA  (PHLAM(I),I=1,NPH) / 911.55, 911.33, 504.30, 227.85, 353.02,
     >                           225.69, 302.67, 195.49 /
      DATA (CPHLAM(I),I=1,NPH) /
     >                            'Ly edge for Hydrogen models  ',
     >                            'Ly edge for Helium models    ',
     >                            'He I edge for Helium models  ',
     >                            'He II edge for Helium models ',
     >                            'O II edge                    ',
     >                            'O III edge                   ',
     >                            'Ne II edge                   ',
     >                            'Ne III edge                  '/
      
C***  HC = H * C  ( ERG * ANGSTROEM )
      DATA HC / 1.98648E-8/
      DATA PI / 3.141592654 /
 
 
C***  CALCULATION OF H-LYMAN PHOTONS
      PRINT *, 'NUMBER OF PHOTONS BLUEWARDS OF A SPECIFIED EDGE'
      PRINT *, '-----------------------------------------------'
      PRINT 15, 'EDGE', 'INTEGRATED UP TO', 
     >         'LOG OF NUMBER OF PHOTONS PER SECOND'
   15 FORMAT (/, 1X, A4, 40X, A16, 3X, A35)
      DO 17 I=1, NPH
      SNNUE=.0
        DO 6 K=1,NF
          IF (XLAMBDA(K) .GT. PHLAM(I) ) GOTO 7
          KMAX=K
          ANNUE=EMFLUX(K)*XLAMBDA(K)/HC
          SNNUE=SNNUE+ANNUE*FWEIGHT(K)
    6   CONTINUE
    7   CONTINUE
        IF (SNNUE.LE.0.) THEN 
           SNNUE=UNDEF
           ELSE
           SNNUE=ALOG10(4.*PI*PI*RSTAR*RSTAR*SNNUE)
           ENDIF
        XLAMMAX=(XLAMBDA(KMAX)+XLAMBDA(KMAX+1))*.5
        PRINT 18, CPHLAM(I), PHLAM(I), XLAMMAX, SNNUE
   18   FORMAT (1X, A30, 3X, F7.2, 3X, F7.2, 12X, F6.2)
   17 CONTINUE
 

C***  Zanstra Temperatures using visual monochromatic flux at V for reference

C**   Loop: H I and He II Zanstra temperatures, respectively
C***   Note: The Zanstra temperature for He I is not calculated. 
C***         In order to include this, one must only change the next 
C***         statement to:
C***      DO 70 I=2, 4
C***      Moreover, the corresponding WRITE statement further below must be 
C***      activated. This is not done by standard, because it need changes 
C***      in the program colofilter.exe that is used on the PoWR homepage
C***      for extracting colors etc. from wruniq.out 
      DO 70 I=2, 4, 2
      XLAMV = 5450.       

C***  Model
      CALL LIPO (VMODEL, XLAMV, EMFLUX, XLAMBDA, NF)
      UVMODEL  = .0 
      DO  K=1,NF
       IF (XLAMBDA(K) .GT. PHLAM(I) ) GOTO 20
      UVMODEL  = UVMODEL  + 
     >           FWEIGHT(K) *   EMFLUX(K)         * XLAMBDA(K) 
C     (Faktor hc kuerzt sich auf beiden Termen raus)
      ENDDO
   20 CONTINUE
      QMODEL = VMODEL / UVMODEL

C***  Blackbody: Temperature loop
      TZAFIRST = 5000.
      TZBFIRST = 500000.
      TZA = TZAFIRST
      TZB = TZBFIRST
      INIT = 0
   10 CONTINUE
      IF (INIT .EQ. 0) TZ = TZA
      IF (INIT .EQ. 1) TZ = TZB
      VPLANCK = BNUE(XLAMV, TZ)

      
      UVPLANCK = .0 
      DO  K=1,NF
       IF (XLAMBDA(K) .GT. PHLAM(I) ) GOTO 21
      UVPLANCK = UVPLANCK + 
     >           FWEIGHT(K) * BNUE(XLAMBDA(K),TZ) * XLAMBDA(K) 
C     (Faktor hc kuerzt sich raus)
      ENDDO
   21 CONTINUE
      QPLANCK = VPLANCK / UVPLANCK
      
      F = QMODEL - QPLANCK 

      IF (INIT .EQ. 0) THEN
         IF (F .GT. .0) GOTO 99
         INIT = 1
      ELSE  IF (INIT .EQ. 1) THEN
         IF (F .LT. .0) GOTO 99
         INIT = 2
      ELSE
       IF (F .GT. .0) THEN
         TZB = TZ 
         ELSE
         TZA = TZ 
       ENDIF
      TZ = 0.5 * (TZA + TZB)
      ENDIF 

      IF (TZB-TZA .GT. 10.) GOTO 10

      ITZ = IFIX(TZ)
      IF (I .EQ. 2) WRITE 
     >  (6, '(/,A,I6,A)') 'Zanstra Temperature   H I : ', ITZ, ' K' 
C*** HeI Zanstra not activated, see comment above!!
ccc      IF (I .EQ. 3) WRITE 
ccc     >  (6, '(A,I6,A)') 'Zanstra Temperature  He I : ', ITZ, ' K' 
      IF (I .EQ. 4) WRITE 
     >  (6, '(A,I6,A)')   'Zanstra Temperature  He II: ', ITZ, ' K' 

   70 Continue 
      WRITE (6, '(A, I5, A)') 
     >       'Note: Reference Flux is the continuum at', 
     >       IFIX(XLAMV), ' Ang' 

      RETURN

C***  Error branch
   99 CONTINUE
      PRINT *, '>>>>>>>>> Troubles in SUBROUTINE ZANSTRA:   <<<<<<<<' 
      PRINT *, '>>>>>>>>> No Zanstra temperature determined <<<<<<<<' 
      RETURN 

      END
