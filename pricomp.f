      SUBROUTINE PRICOMP (NDIM,EINST,N,NCHARG,NOM,NATOM,ABXYZ,ATMASS,
     $                   STAGE,NFIRST,NLAST,ELEMENT,SYMBOL,LASTIND,
     $                   INDLOW,INDNUP,NAUTO,LOWAUTO,
     $                   EAUTO,KONTNUP,KONTLOW,LASTKON,XMASS,KRUDAUT)
C***********************************************************************
C***  PRINTOUT OF THE CHEMICAL COMPOSITION OF THE WR MODEL ATMOSPHERE
C***********************************************************************
 
      DIMENSION EINST(NDIM,NDIM)
      DIMENSION NCHARG(N),NOM(N)
      DIMENSION ABXYZ(NATOM)
      DIMENSION ATMASS(NATOM),STAGE(NATOM),NFIRST(NATOM),NLAST(NATOM)
      DIMENSION INDLOW(LASTIND),INDNUP(LASTIND)
      DIMENSION KONTNUP(LASTKON),KONTLOW(LASTKON)
      DIMENSION LOWAUTO(NAUTO),EAUTO(NAUTO),KRUDAUT(NAUTO)
      CHARACTER*80 KARTE
      CHARACTER*10 ELEMENT(NATOM)

      PARAMETER ( MAXION = 24 )
      CHARACTER*4 ROMION(0:MAXION)
      CHARACTER*2 SYMBOL(NATOM)


      DATA (ROMION(I),I=0,MAXION) /'   I','  II',' III','  IV','   V',
     >                             '  VI',' VII','VIII','  IX','   X',
     >                             '  XI',' XII','XIII',' XIV','  XV',
     >                             ' XVI','XVII',' 18 ',' XIX','  XX',
     >                             ' XXI','XXII',' 23 ','XXIV',' XXV'/

 
      PRINT 1
    1 FORMAT (1H1,//,
     $        10X,'C H E M I C A L  C O M P O S I T I O N',
     $         4X,'(FILE: "DATOM")',/,
     $        10X,38('='),/)

C***  DECODING FIRST INPUT CARDS (COMMENTS) FROM TAPE 4 = DATOM
      OPEN (4, FILE='DATOM', STATUS='UNKNOWN')
    9 READ (4,6, END=99) KARTE
    6 FORMAT(A)
      IF (KARTE(:1) .EQ. '*') THEN
C***     EXPLICIT END ("*END") FOR INPUT CARDS IS GIVEN:
         IF (KARTE(2:4) .EQ. 'END') GOTO 99
         PRINT 2, KARTE
    2    FORMAT (1X,A)
         GOTO 9
      ENDIF
   99 CLOSE (4)
      PRINT 22
   22 FORMAT (//)

C***  CALCULATION OF THE MEAN ATOMIC WEIGHT "ATMEAN";
C***  SEARCH FOR THE INDEX OF MAIN ELEMENT "HELIUM"
      ATMEAN=0.0
      NAHE=0
      DO 10 NA=1, NATOM
      ATMEAN=ATMEAN+ABXYZ(NA)*ATMASS(NA)
      IF (SYMBOL(NA) .EQ. 'HE') NAHE=NA
   10 CONTINUE
      IF (NAHE .GT. 0) THEN
         PRINT 11
   11    FORMAT
     >   (1X,' REL. ABUNDANCES:',12X,' BY NUMBER ',5X,'  BY MASS  ',
     >    11X,' REL. TO HELIUM:',12X,' BY NUMBER ',5X,'  BY MASS  ')
      ELSE
         PRINT 12
   12    FORMAT
     >   (1X,' REL. ABUNDANCES:',12X,' BY NUMBER ',5X,'  BY MASS  ')
      ENDIF

      DO 19 NA=1,NATOM
      ABNA=ABXYZ(NA)
      ATMASNA=ATMASS(NA)
      FRACM=ABNA*ATMASNA/ATMEAN
      IF (NAHE .GT. 0) THEN
         ABHE=ABXYZ(NAHE)
         PRINT 13, SYMBOL(NA),ABNA,FRACM,
     >             SYMBOL(NA),ABNA/ABHE,ABNA*ATMASNA/(ABHE*ATMASS(NAHE))
      ELSE
         PRINT 13, SYMBOL(NA),ABNA,FRACM
      ENDIF
   13 FORMAT (23X,A2,5X,1PE11.4,5X,E11.4,32X,A2,5X,E11.4,5X,E11.4)
   19 CONTINUE
 
      NION=1
      DO 20 NLEV=2,N
      IF ((NOM(NLEV) .NE. NOM(NLEV-1)) .OR.
     $    (NCHARG(NLEV) .NE. NCHARG(NLEV-1))) NION=NION+1
   20 CONTINUE
      NDR=0
      NLTE=0
      DO 26 IND=1,NAUTO
      IF (EAUTO(IND) .GT. 0.) NDR=NDR+1
      IF (EAUTO(IND) .LT. 0.) NLTE=NLTE+1
   26 CONTINUE
      IF ((NDR+NLTE) .NE. NAUTO) STOP 'NAUTO'
      PRINT 21, NATOM,NION,N,LASTKON,LASTIND,NDR,NLTE
   21 FORMAT (///,1X,' STATISTICS:',5X,I3,'  ELEMENTS',/,18X,I3,'  IONS'
     $        ,/,17X,I4,'  LEVELS'
     $        ,/,17X,I4,'  CONTINUUM TRANSITIONS'
     $        ,/,16X,I5,'  LINE TRANSITIONS'
     $        ,/,17X,I4,'  STABIL. TRANSITIONS FROM AUTOIONIZING '
     $                 ,'STATES'
     $        ,/,17X,I4,'  TRANSITIONS FROM LTE LEVELS',///)
 
      PRINT 31
   31 FORMAT (1X,' INDEX',3X,'ELEMENT',6X,'ATOMIC MASS',3X,'IONS',3X,
     $        'MAIN ION',/,1X,52('-'))
      DO 39 NA=1,NATOM
      NAION=1
      DO 38 NLEV=NFIRST(NA)+1,NLAST(NA)
      IF (NCHARG(NLEV) .NE. NCHARG(NLEV-1)) NAION=NAION+1
   38 CONTINUE
      ISTAGE = INT(STAGE(NA))-1
      IF (ISTAGE .GT. MAXION) THEN
         WRITE (0,*) '*** MAXION too small ***************'
         STOP        '*** FATAL ERROR IN SUBR. PRICOMP ***'
      ENDIF
      PRINT 32,NA,ELEMENT(NA),ATMASS(NA),NAION,SYMBOL(NA),
     $         ROMION(ISTAGE)
   32 FORMAT (3X,I2,5X,A10,5X,F6.2,7X,I2,5X,A2,1X,A4)
   39 CONTINUE

      PRINT 50, ATMEAN, XMASS
   50 FORMAT (///, '  MEAN ATOMIC MASS:', F6.2, 5X, 
     >        'MEAN PARTICLE MASS:', F6.2)

 
      PRINT 41
   41 FORMAT (///,1X,' ELEMENT',3X,'ION',3X,'CHARGE',3X,'LEVELS',3X,
     $        'CONTINUA',3X,
     $        'LINES',1X,'(RUD.)',3X,
     $        'LINES(LTE)',1X,'(RUD.)',3X,
     $        'STAB.LINES',1X,'(RUD.)',
     $        /,1X,99('-'))
      DO 49 NA=1,NATOM
      N1=NFIRST(NA)
   44 NCH1=NCHARG(N1)
      IF (NCH1 .EQ. NCHARG(NLAST(NA))) THEN
          NION=NLAST(NA)-N1+1
      ELSE
          NION= ISRCHNE(NLAST(NA)-N1+1,NCHARG(N1+1),1,NCH1)
      ENDIF
      NKONT=0
      NLINE=0
      NRUD1=0
      NLTE=0
      NRUD2=0
      NDR=0
      NRUD3=0
      DO 43 KON=1,LASTKON
      LOW=KONTLOW(KON)
      IF ((NOM(LOW) .EQ. NA) .AND. (NCHARG(LOW) .EQ. NCH1)) 
     $    NKONT=NKONT+1
   43 CONTINUE
      DO 45 IND=1,LASTIND
      LOW=INDLOW(IND)
      IF ((NOM(LOW) .EQ. NA) .AND. (NCHARG(LOW) .EQ. NCH1)) THEN
          NLINE=NLINE+1
          IF (EINST(LOW,INDNUP(IND)) .EQ. -2.) NRUD1=NRUD1+1
      ENDIF
   45 CONTINUE
      DO 46 IND=1,NAUTO
      LOW=LOWAUTO(IND)
      IF ((NOM(LOW) .EQ. NA) .AND. (NCHARG(LOW) .EQ. NCH1)) THEN
         IF (EAUTO(IND) .LT. 0.) THEN
            NLTE=NLTE+1
            IF (KRUDAUT(IND) .EQ. 1) NRUD2=NRUD2+1
         ENDIF
         IF (EAUTO(IND) .GT. 0.) THEN
            NDR=NDR+1
            IF (KRUDAUT(IND) .EQ. 1) NRUD3=NRUD3+1
         ENDIF
      ENDIF
   46 CONTINUE
      PRINT 42, SYMBOL(NA),ROMION(NCH1),NCH1,NION,NKONT,NLINE,NRUD1,
     >          NLTE,NRUD2,NDR,NRUD3
   42 FORMAT (4X,A2,6X,A4,4X,I2,7X,I2,8X,I2,6X,I4,3X,I4,8X,I3,6X,I3,
     >        7X,I4,6X,I3)
      N1=N1+NION
      IF (N1 .LE. NLAST(NA)) GOTO 44
   49 CONTINUE

      RETURN
      END
