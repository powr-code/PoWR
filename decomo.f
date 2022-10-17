      SUBROUTINE DECOMO (LSOPA, LSINT, NOTEMP, LSHTOT, LPLOHTOT, 
     >            MODHIST, BUNLU, BPLOTRTAU1, 
     >            NPLOTOPA, OPTIONPLOTOPA, NPLOTOPADIM)
C*******************************************************************************
C***  DECODING INPUT OPTIONS FOR PROGRAM "COMO"
C*******************************************************************************
 
      CHARACTER(80) :: KARTE
      CHARACTER(80), DIMENSION(NPLOTOPADIM) :: OPTIONPLOTOPA
      LOGICAL NOTEMP, BUNLU, BPLOTRTAU1
      INTEGER :: NPLOTOPA, NPLOTOPADIM

C***  DEFAULT VALUES
      LSOPA=-1
      LSINT=-1
      LSHTOT=-1
      LPLOHTOT=-1
      NCON=0
      NOTEMP=.FALSE.
      BUNLU = .FALSE.
      BPLOTRTAU1 = .FALSE.
      NPLOTOPA = 0

      OPEN (1, FILE='CARDS', STATUS='UNKNOWN')
      REWIND 1
 
    8 READ (1,4, END=99) KARTE
    4 FORMAT (A)
 
      IF (KARTE(:10) .EQ. 'PRINT INT ') THEN
C                          ==========
            DECODE (80,7,KARTE) XL
    7       FORMAT (10X,F10.0)
   77       FORMAT (11X,F10.0)
            LSINT=IFIX(XL)
            IF (LSINT.EQ.0) LSINT=1
      ELSE IF (KARTE(:10) .EQ. 'PRINT OPA ') THEN
C                               ==========
            DECODE (80,7,KARTE) XL
            LSOPA=IFIX(XL)
            IF (LSOPA.EQ.0) LSOPA=1
      ELSE IF (KARTE(:11) .EQ. 'PRINT HTOTC') THEN
C                               ===========
            DECODE (80,77,KARTE) XL
            LSHTOT=IFIX(XL)
            IF (LSHTOT.EQ.0) LSHTOT=1
      ELSE IF (KARTE(:10) .EQ. 'PLOT HTOTC') THEN
C                               ==========
            LPLOHTOT=1
      ELSE IF (KARTE(:10) .EQ. 'PLOT RTAU1') THEN
C                               ==========
            BPLOTRTAU1 = .TRUE.
      ELSE IF (KARTE(:8) .EQ. 'PLOT OPA') THEN
C                              ==========
            NPLOTOPA = NPLOTOPA + 1
            IF (NPLOTOPA .LE. NPLOTOPADIM) THEN 
               OPTIONPLOTOPA(NPLOTOPA) = KARTE
            ELSE
               WRITE (*,*) 
     >            'WARNING: more PLOT OPA options than dimensioned' 
               WRITE (0,*) 
     >            'WARNING: more PLOT OPA options than dimensioned' 
            ENDIF
      ELSE IF (KARTE(:8) .EQ. 'NO TEMPE') THEN
C                              ========
            NOTEMP=.TRUE.
            IF (KARTE(30:40) .NE. ' ') 
     $         CALL DECNOT (MODHIST,MODHIST,KARTE,NOTEMP,'COMO')
      ELSE IF (KARTE(:7) .EQ. 'UNLUTEC') THEN
C                              =======
            BUNLU=.TRUE.
      ENDIF

      GOTO 8

C***  END-OF-FILE REACHED:
   99 CONTINUE
      CLOSE (1)
      RETURN

      END
