      SUBROUTINE COLLI(NDIM,N,ENLTE,TL,ENE,NCHARG,ELEVEL,EINST,CRATE,
     $                 EION,COCO,KEYCBB,WEIGHT,ALTESUM,NATOM,NOM,KODAT,
     $                 INDNUP,INDLOW,LASTIND,
     $                 KONTNUP,KONTLOW,LASTKON,KEYCBF,IONGRND)
C*******************************************************************************
C***  COLLISIONAL TRANSITION RATES STORED IN MATRIX CRATE **********************
C***  BOUND-BOUND: DEPENDING ON THE ELEMENT (HE, H, N, C, O)
C***  OMEGA(UP-LOW) IS CALCULATED DEPENDING ON KEYWORD KEYCBB  *********
C***  BOUND-FREE: NOT DEPENDING ON THE ELEMENT
C*******************************************************************************
 
      DIMENSION EINST(NDIM,NDIM),CRATE(NDIM,NDIM)
      DIMENSION ENLTE(NDIM),NCHARG(NDIM),ELEVEL(NDIM),WEIGHT(NDIM)
      DIMENSION EION(NDIM),ALTESUM(4,NDIM)
      DIMENSION IONGRND(NDIM)
      DIMENSION INDNUP(LASTIND),INDLOW(LASTIND)
      DIMENSION KONTNUP(LASTKON),KONTLOW(LASTKON)
C***  ARRAY "KEYCBF" IS PROVIDED FOR FUTURE TESTS OF DIFFERENT COLLISIONAL
C***  IONIZATION FORMULAES:
      DIMENSION KEYCBF(LASTKON)
      DIMENSION NOM(N)
      DIMENSION KODAT(NATOM)
      CHARACTER*4 KEYCBB(LASTIND)

C***  SUPPRESS REPETITIVE WARNING OF NEGATIVE CROSS SECTIONS
      DATA CBBWARN / 0. /
 
C***  C1 = H * C / K    ( CM * ANGSTROEM )
      DATA C1 / 1.4388 /

      TROOT=SQRT(TL)
      T32=TL*TROOT
 
C***  INITIALIZE ALL ELEMENTS OF MATRIX "CRATE":
      DO 1 J=1,N
      DO 1 I=1,N
      CRATE(I,J)=.0
    1 CONTINUE
 
C***  LINE TRANSITIONS  *******************************************************
      DO 11 IND=1,LASTIND
      NUP=INDNUP(IND)
      LOW=INDLOW(IND)
      WAVENUM=ELEVEL(NUP)-ELEVEL(LOW)
      WN2=WAVENUM*WAVENUM
      WN3=WN2*WAVENUM
 
C***  'NONE': TRANSITIONS WITH UNKNOWN COLLISIONAL COEFFICIENTS
C***          (COLLISIONAL CROSS SECTION SIGMA(LOW,UP) IS SET TO  PI*A0**2)
      IF (KEYCBB(IND) .EQ. 'NONE') THEN
C                           ====
              OMEGA=5.465E-11*TROOT*(1.+C1*WAVENUM/TL)*WEIGHT(LOW)/
     /               WEIGHT(NUP)
              GOTO 99
 
C***  'ZERO': NO COLLISIONAL TRANSITION
      ELSE IF (KEYCBB(IND) .EQ. 'ZERO') THEN
C                                ====
              OMEGA=0.
              GOTO 99
      ENDIF
 
C***  HELIUM  ==========================================================
      IF (NOM(LOW) .EQ. KODAT(2)) THEN
            CALL CBBHE (OMEGA,IND,NUP,LOW,TL,TROOT,T32,NDIM,N,NCHARG,
     $                  EINST,COCO,KEYCBB,WEIGHT,WAVENUM,WN2,WN3)
 
C***  HYDROGEN  ========================================================
      ELSE IF (NOM(LOW) .EQ. KODAT(1)) THEN
            CALL CBBH  (OMEGA,IND,NUP,LOW,TL,TROOT,T32,NDIM,N,NCHARG,
     $                  EINST,COCO,KEYCBB,WEIGHT,WAVENUM,WN2,WN3)
 
C***  NITROGEN  ========================================================
      ELSE IF (NOM(LOW) .EQ. KODAT(7)) THEN
            CALL CBBN  (OMEGA,IND,NUP,LOW,TL,TROOT,T32,NDIM,N,NCHARG,
     $                  EINST,COCO,KEYCBB,WEIGHT,WAVENUM,WN2,WN3)
 
C***  GENERIC ION ======================================================
      ELSE IF (NOM(LOW) .EQ. KODAT(26)) THEN
            CALL CBBFE (OMEGA, NUP, LOW, TL, TROOT, WAVENUM, WN3, 
     >                  EINST, NDIM)
 
      ELSE

C*** ALL OTHER ELEMENTS ================================================
         CALL CBBMORE  (OMEGA,IND,NUP,LOW,TL,TROOT,T32,NDIM,N,NCHARG,
     $                  EINST,COCO,KEYCBB,WEIGHT,WAVENUM,WN2,WN3)
      ENDIF
 
ccc      WRITE (0,*) 'OMEGA', LOW, NUP, OMEGA, OMEGA*ENLTE(NUP)/ENLTE(LOW)

C***  COLLISION RATE COEFFICIENTS CRATE
   99 CONTINUE
      IF (OMEGA .LT. 0.0 .AND. CBBWARN .EQ. 0.) THEN
         CBBWARN = 1.
         PRINT 98, NUP, LOW
   98    FORMAT (' STEAL> WARNING: NEGATIVE BOUND-BOUND COLLISIONAL ',
     $           'CROSS SECTION DETECTED (LEVELS: UP= ',I3,', LOW= ',
     $            I3,')')
      ENDIF
      CRATE(NUP,LOW)=ENE*OMEGA
c        if (enlte(nup) .eq. 0. .or. enlte(low) .eq. 0.) then
c          write (0,'(a,2i4,2g25.15)')
c     >      'COLLI2: warning lu, enlte(low) enlte(nup)=',
c     >      low,nup,enlte(low),enlte(nup)
c        endif

      CRATE(LOW,NUP)=ENE*OMEGA*ENLTE(NUP)/ENLTE(LOW)
   11 CONTINUE

C***  This is a modification added 24-Aug-2010 by helge + wrh
C***  - see WR-Memo 100825.txt
C***  Here we assign a collisonal rate coefficient for those superlines 
C***  which have ZERO radiative cross-section and were therefore
C***  skipped in the index numbering (i.e. no call of CBBFE occured)
      DO LOW=1, N  
C***    Only Iron:
        IF (NOM(LOW) .EQ. KODAT(26)) THEN 
           DO NUP=LOW+1, N
C***          Within the same ion?
              IF (NCHARG(LOW) .EQ. NCHARG(NUP)) THEN
C***             No CRATE given by CBBFE?
                 IF (CRATE(NUP,LOW) .EQ. .0) THEN
                     WAVENUM=ELEVEL(NUP)-ELEVEL(LOW)
C***                 Bohr's cross section 
                     OMEGA = 5.465E-11*TROOT * (1.+C1*WAVENUM/TL)
     >                       * WEIGHT(LOW)/ WEIGHT(NUP)
                     CRATE(NUP,LOW)=ENE*OMEGA
                     CRATE(LOW,NUP)=ENE*OMEGA*ENLTE(NUP)/ENLTE(LOW)
                 ENDIF
              ENDIF              
           ENDDO
        ENDIF
      ENDDO

C***  COLLISIONAL IONIZATION   ************************************************
C***  OMEGA(LOW-UP) ACCORDING TO JEFFERIES P. 121 (EQ. 6.39)
C***  G = FACTOR DEPENDING ON THE CHARGE OF THE ION (=UPPER LEVEL)
      DO 21 KON=1,LASTKON
      NUP=KONTNUP(KON)
      LOW=KONTLOW(KON)
      G=.3
      IF (NCHARG(NUP) .EQ. 2) G=.2
      IF (NCHARG(NUP) .EQ. 1) G=.1
      EDGE = ELEVEL(NUP)+EION(LOW)-ELEVEL(LOW)
      EXPFAC=EXP(-C1*EDGE/TL)
      OMEGA=G*1.08E-5*TROOT*EINST(LOW,NUP)*EXPFAC/EDGE

C***  BOUND-FREE COLLISIONAL RATE COEFFICIENTS: CRATE
      CRATE(LOW,NUP)=ENE*OMEGA
      CRATE(NUP,LOW)=ENE*OMEGA*ENLTE(LOW)/ENLTE(NUP)
   21 CONTINUE 


C***  ADD THE COLLISIONS TO ADDITIONAL UPPER LEVELS WHICH ARE ASSUMED IN LTE.
C***  NOTE: THE SUMMATION HAS BEEN PERFORMED IN ADVANCE, AND IS IMPLICITELY
C***  CONTAINED IN THE PARAMETERS 'ALTESUM' ENTERING NOW FROM THE ATOMIC DATA
C***  FILE.
C***  NOTE:   IN THIS VERSION ALL DERIVATIVES OF THE "ALTESUM" TERMS
C    =======     --> SUBR. DERIV: WITH RESPECT TO ELECTRON DENSITY
C***             --> SUBR. DERIVT:     - || -     TEMPERATURE
C***          ARE NEGLECTED
      DO 31 LOW=1, N
      IF (ALTESUM(1,LOW) .GT. .0) THEN
C***     UPPER LEVEL IS THE GROUND LEVEL OF THE PARENT ION
         NUP=IONGRND(LOW)
         EDGE=ELEVEL(NUP)+EION(LOW)-ELEVEL(LOW)
         X=1000./TL
         FOFT=(ALTESUM(3,LOW)*X+ALTESUM(2,LOW))*X
         FOFT=10.**FOFT
         AOFT=ALTESUM(1,LOW)*FOFT
         OMSUM=4.06/TROOT*EXPFAC/EDGE/EDGE/EDGE * AOFT
C***     BOUND-FREE COLLISIONAL RATE COEFFICIENTS: CRATE
         CRATE(LOW,NUP)=CRATE(LOW,NUP)+ENE*OMSUM
         CRATE(NUP,LOW)=CRATE(NUP,LOW)+ENE*OMSUM*ENLTE(LOW)/ENLTE(NUP)
      ENDIF
   31 CONTINUE
 
      RETURN
      END
