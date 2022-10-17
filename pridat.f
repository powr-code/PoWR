      SUBROUTINE PRIDAT (NDIM,N,LEVEL,NCHARG,WEIGHT,ELEVEL,EION,EINST,
     $                   KEY,NF,ALPHA,SEXPO,
     $                   ADDCON1, ADDCON2, ADDCON3, 
     $                   IGAUNT, SIGMATHK, SEXPOK, EDGEK, MAXATOM,
     >                   COCO,KEYCBB, ALTESUM,
     $                   NATOUT,NATOM,ELEMENT,NOM,ABXYZ,ATMASS,
     $                   NAUTO,MAXAUTO,LOWAUTO,WAUTO,EAUTO,AAUTO,
     $                   IONAUTO,
     $                   INDNUP,INDLOW,LASTIND,IONGRND,
     $                   KEYCBF,KONTNUP,KONTLOW,LASTKON,
     $                   MAXNDR,ELEVDR,NCHARGDR,INDNUPDR, 
     $                   LEVELDR,NOMDR,IONDR,WEIGHTDR,KRUDAUT) 
C*******************************************************************************
C***  PRINTOUT OF THE ATOMIC DATA DECODED *************************************
C***  ADDITIONAL PRINTOUT OF RELATIVE ABUNDANCES (READ FROM INPUT CARDS)
C***  AND MASS FRACTIONS
C*******************************************************************************

      DIMENSION NCHARG(NDIM),WEIGHT(NDIM),ELEVEL(NDIM)
      DIMENSION EION(NDIM),EINST(NDIM,NDIM),KEY(NF),IONGRND(NDIM)
      DIMENSION ALPHA(LASTKON),SEXPO(LASTKON),IGAUNT(LASTKON)
      DIMENSION ADDCON1(LASTKON), ADDCON2(LASTKON), ADDCON3(LASTKON)
      DIMENSION KONTNUP(LASTKON),KONTLOW(LASTKON),KEYCBF(LASTKON)
      DIMENSION INDNUP(LASTIND),INDLOW(LASTIND)
      DIMENSION SIGMATHK(MAXATOM,MAXATOM),SEXPOK(MAXATOM,MAXATOM)
      DIMENSION EDGEK(MAXATOM,MAXATOM)
      DIMENSION NOM(N)
      DIMENSION ABXYZ(NATOM),ATMASS(NATOM)
      DIMENSION COCO(4,LASTIND)
      DIMENSION ALTESUM(4,NDIM)
      DIMENSION LOWAUTO(NAUTO),WAUTO(NAUTO),EAUTO(NAUTO)
     $         ,AAUTO(NAUTO),IONAUTO(NAUTO),KRUDAUT(NAUTO)
      DIMENSION NCHARGDR(MAXNDR),ELEVDR(MAXNDR),NOMDR(MAXNDR)
      DIMENSION WEIGHTDR(MAXNDR),INDNUPDR(MAXAUTO)
      CHARACTER*3 KRUDI,DRRUDI
      CHARACTER*4 KEYCBB(LASTIND)
      CHARACTER*10 LEVEL(NDIM), LEVELI, PARION
      CHARACTER*10 ELEMENT(NATOM),NAME
      CHARACTER*10 IONDR(MAXNDR),LEVELDR(MAXNDR)
      CHARACTER CHF*9, CHCOCO(4)*8
 
C***  C4 = 1 / (SIGMA-CLASSIC * 8 PI)     (IN ANGSTROEM**2)
      DATA C4 /1.499E-16/

C***  CALCULATION OF THE MEAN ATOMIC WEIGHT "ATMEAN"
      ATMEAN=0.
      DO 40 NA=1,NATOM
      ATMEAN=ATMEAN+ABXYZ(NA)*ATMASS(NA)
   40 CONTINUE

C***  COLLECTING INFORMATIONS ABOUT UPPER DR-LEVELS
      IF (NAUTO .GT. 0) CALL DRDAT( NDIM, MAXNDR, MAXAUTO, 
     $            LOWAUTO, EION, LEVEL, NCHARG, NOM, IONGRND,
     $            NDR, ELEVDR, NCHARGDR, INDNUPDR, LEVELDR,
     $            NOMDR, IONDR, WEIGHTDR) 

      IND=0
      PRINT 9
    9 FORMAT(1H1,//,20X,'A T O M I C   D A T A   U S E D :',
     $  /,20X,33('='))
 
C***  ELEMENTS ---------------------------------------------------------
C***  OUTPUT ONLY FOR ONE OR ALL ELEMENTS
      IF (NATOUT .NE. 0) THEN
         NASTART=NATOUT
         NASTOP=NATOUT
      ELSE
         NASTART=1
         NASTOP=NATOM
      ENDIF
      DO 29 NA=NASTART,NASTOP
      ABNA=ABXYZ(NA)
      ABLOG=ALOG10(ABNA)
      FRACM=ABNA*ATMASS(NA)/ATMEAN
      PRINT 10, ABNA,ELEMENT(NA),ABLOG,FRACM
   10 FORMAT (//,20X,33('/'),10X,'RELATIVE ABUNDANCE (BY NUMBER):',2X,
     $           G10.3,/,
     $       20X,33('/'),41X,12('='),/,
     $       20X,'/',31X,'/',/,
     $       20X,'/',11X,A10,10X,'/',33X,'LOG(AB)=  ',F6.2,/,
     $       20X,'/',31X,'/',/,
     $       20X,33('/'),27X,'MASS FRACTION:',2X,G10.3,/,
     $       20X,33('/'),41X,12('='),/)
 
C***  LEVELS -----------------------------------------------------------
      PRINT 11
   11 FORMAT (//,10X,'1. LEVELS :',/,10X,11('-'),/,
     $ '  NR      NAME       WEIGHT     ENERGY    CHARGE   IONIZATION'
     $ ,' POT.',/,
     $ 30X,'(KAYSER)',15X,'(KAYSER)',/)
      DO 2 J=1,N
      IF (NOM(J) .NE. NA) GOTO 2
      NW=IFIX(WEIGHT(J))
      PRINT 3, J,LEVEL(J),NW,ELEVEL(J),NCHARG(J),EION(J)
    3 FORMAT (1X,I3,3X,A10,I10,F11.2,I10,5X,F12.2)
    2 CONTINUE
 
C***  LINE TRANSITIONS  ------------------------------------------------
      PRINT 4
    4 FORMAT (//,10X,'2. LINE TRANSITIONS :',/,10X,20('-'),//,
     $ '  IND  UP LOW   UP           LOW       WAVELENGTH   A(UP-LOW)',
     $ 2X,'F(LOW-UP)',12X,'COLLISIONAL  TRANSITIONS',/,
     $ '   (INDICES)    (CONFIGURATION)       (ANGSTROEM) (PER SECOND)'
     $ ,18X,'(KEYWORD)',10X,'(COEFFICIENTS)',/)
      DO 1 IND=1,LASTIND
      I=INDNUP(IND)
      IF (NOM(I) .NE. NA) GOTO 1
      J=INDLOW(IND)
      KRUDI='   '
      IF (EINST(J,I) .EQ. -2.) KRUDI='RUD'
      WLENG=1.E8/(ELEVEL(I)-ELEVEL(J))
      FVALUE = C4 * WLENG * WLENG * EINST(I,J) * WEIGHT(I) / WEIGHT(J)
C***  SHORTENING OF OUTPUT:  VALUES==0.0 ARE OMITTED
      IF (FVALUE .NE. 0.0) THEN
         WRITE (CHF,'(1PG9.3)') FVALUE
      ELSE
         CHF='       '
      ENDIF
      DO 17 ICOCO=1,4
      IF (COCO(ICOCO,IND) .NE. 0.0) THEN
         WRITE (CHCOCO(ICOCO),'(1PE8.1)') COCO(ICOCO,IND)
      ELSE
         CHCOCO(ICOCO)='        '
      ENDIF
   17 CONTINUE
      PRINT 5, IND,I,J,LEVEL(I),LEVEL(J),WLENG,EINST(I,J),CHF,KRUDI,
     $         KEYCBB(IND),(CHCOCO(M),M=1,4)
    5 FORMAT(1X,3I4,3X,A10,3X,A10,F10.2,1PE12.2,2X,A9,3X,A3,5X,A4,
     $       4(2X,A8))
    1 CONTINUE
 
C***  CONTINUA  --------------------------------------------------------
      PRINT 7
    7 FORMAT(//,10X,'3. CONTINUA :',/,10X,13('-'),//,
     $ '  KON LOW  UP   LOW          UP         THRESHOLD     ',
     $ 'PHOTO CROSS-SECTION    SEATON COEFF.         GAUNT FACTOR'/
     $ '  (INDICES)       (CONFIGURATIONS)    (ANGSTROEM)     ',
     $ 17H  (10**-18 CM**2) ,'      ALPHA      S           FORMULA'/)
      LASTI=0
      DO 16 KON=1,LASTKON
      I=KONTLOW(KON)
      IF (NOM(I) .NE. NA) GOTO 16
      J=KONTNUP(KON)
      WEDGE=1.E8/(ELEVEL(J)+EION(I)-ELEVEL(I))
      IF (I .EQ. LASTI) THEN
         LEVELI='    ||    '
      ELSE
         LEVELI=LEVEL(I)
      ENDIF
*       WRITE (*,*) KON,I,J,LEVELI,LEVEL(J),WEDGE,EINST(I,J),
*     >             SEXPO(KON)
*       WRITE (*,'(A8)') IGAUNT(KON)
      IF (ALPHA(KON) .EQ. .0) THEN
         PRINT 8,
     $     KON,I,J,LEVELI,LEVEL(J),WEDGE,EINST(I,J),SEXPO(KON),
     $     IGAUNT(KON),
     $     ADDCON1(KON), ADDCON2(KON), ADDCON3(KON)
    8    FORMAT
     $   (1X,3I4,3X,A10,3X,A10,F10.2,F15.3,11X,'HYDROGENIC',F7.3,8X,A8,
     $    F7.3,1X,F7.3,1X,F7.3)
      ELSE
         PRINT 13,
     $     KON,I,J,LEVELI,LEVEL(J),WEDGE,EINST(I,J),ALPHA(KON),
     $     SEXPO(KON),IGAUNT(KON),
     $     ADDCON1(KON), ADDCON2(KON), ADDCON3(KON)
   13    FORMAT
     $   (1X,3I4,3X,A10,3X,A10,F10.2,F15.3,F18.3,F10.3,8X,A8,
     $    F7.3,1X,F7.3,1X,F7.3)
      ENDIF
      LASTI=I
   16 CONTINUE

C***  K-SHELL IONISTION
        PRINT 50
   50   FORMAT(//,10X,'4. K-SHELL-IONISATION :',/,
     $         10X,20('-'),//,
     $'   ELEMENT       SIGMA         EXPONENT    IONISATION-ENERGY ',/,
     $'            (10**-18 CM**2)                    (KAYSER) ',/)

      DO ISTATE=1, MAXATOM
         IF (EDGEK(NA,ISTATE) .NE. .0) 
     >     PRINT 51, ELEMENT(NA), ISTATE, SIGMATHK(NA,ISTATE), 
     >                                    SEXPOK  (NA,ISTATE),
     >                                    EDGEK   (NA,ISTATE)
   51    FORMAT (3X,A10,1X,I2,3X,F10.3,3X,F10.3,6X,F10.1) 
      ENDDO
 
C***  DIELECTRONIC RECOMBINATIONS  -------------------------------------
      DO 31 I=1,NAUTO
      IF (NOM(LOWAUTO(I)) .EQ. NA) GOTO 32
   31 CONTINUE
      GOTO 20
   32 CONTINUE

      PRINT 90
   90 FORMAT (//,10X,'4. DR-LEVELS :',/,10X,14('-'),/,
     $ '  NR      NAME       WEIGHT     ENERGY    CHARGE   PARENT ION'
     $ ,/,30X,'(KAYSER)',/)

      DO 91 J=1,NDR
      IF (NOMDR(J) .NE. NA) GOTO 91
      NW=IFIX(WEIGHTDR(J))
      PRINT 92, J+N,LEVELDR(J),NW,ELEVDR(J),NCHARGDR(J),IONDR(J)
   92 FORMAT (1X,I3,3X,A10,I10,F11.2,I10,3X,A10)
   91 CONTINUE

      PRINT 33
   33 FORMAT (//,10X,'5. STABILIZING LINE TRANSITIONS',/,
     $ 10X,31('-'),//,
     $ '  IND  UP LOW INDDR  DR-LEVEL     LOW       WAVELENGTH',
     $ '   A(UP-LOW)  F(LOW-UP)  RUD',/,
     $ '   (INDICES)         (CONFIGURATION)       (ANGSTROEM)',
     $ ' (PER SECOND)',/)

      DO 35 INDDR=1,NAUTO
      I = INDNUPDR(INDDR)
      J = LOWAUTO(INDDR)
      WLENG=1.E8/(ELEVDR(I)-ELEVEL(J))
      FVALUE = C4 * WLENG * WLENG * AAUTO(INDDR) * 
     *         WEIGHTDR(I) / WEIGHT(J)
      WRITE (CHF,'(1PG9.3)') FVALUE
      DRRUDI='   '
      IF (KRUDAUT(INDDR) .EQ. 1) DRRUDI=' X '
      IF (NOM(J) .EQ. NA) PRINT 34, 
     >        INDDR+LASTIND,I+N,J,INDDR,LEVELDR(I),LEVEL(J),
     >        WLENG,AAUTO(INDDR),CHF,DRRUDI
   34 FORMAT(1X,3I4,1X,I4,3X,A10,3X,A10,F10.2,1PE12.2,2X,A9,2X,A3)
   35 CONTINUE

C***  SUM OF TRANSITIONS TO LTE LEVELS  --------------------------------
   20 DO 21 I=1,N
      IF ((NOM(I) .EQ. NA) .AND. (ALTESUM(1,I) .GT. .0)) GOTO 22
   21 CONTINUE
      GOTO 29
   22 PRINT 23
   23 FORMAT (//,10X,'5. SUM OF TRANSITIONS TO LTE LEVELS',/,
     $ 10X,35('-'),//,
     $ '     LOW         LOW            UP          A-SUM',
     $      '   TEMPERATURE FUNCTION'/
     $ ' (INDEX)   (CONFIG.)   (SUM RANGE)   (PER SECOND)',
     $      '     COEFF.1   COEFF.2'/)
      DO 25 LOW=1,N
      IF ((NOM(LOW) .EQ. NA) .AND. (ALTESUM(1,LOW) .GT. .0)) PRINT 24,
     $  LOW,LEVEL(LOW),ALTESUM(4,LOW),ALTESUM(1,LOW),ALTESUM(2,LOW),
     $   ALTESUM(3,LOW)
   24 FORMAT (I8,2X,A10,6X,A8,1P,E15.3,0P,F12.4,F10.4)
   25 CONTINUE

   29 CONTINUE
C***  ENDLOOP ----------------------------------------------------------
 
      RETURN
      END
