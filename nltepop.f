      SUBROUTINE NLTEPOP (NDIM,N,ENLTE,TL,ENE,WEIGHT,NCHARG,EION,ELEVEL,
     $                   EN,EINST,XLAMBDA,FWEIGHT,XJC,NF,L,LEVEL,
     $                   XJL,ND,CRATE,RRATE,RATCO,SIGMAKI,ALTESUM,COCO,
     $                   KEYCBB,NOM,NATOM,ABXYZ,KODAT,NFIRST,NLAST,
     $                   NAUTO,MAXAUTO,LOWAUTO,WAUTO,EAUTO,AAUTO,
     $                   RDIEL,RAUTO,IONAUTO,IONGRND,
     $                   INDNUP,INDLOW,LASTIND,
     $                   KONTNUP,KONTLOW,LASTKON,KODRNUP,KODRLOW,
     $                   LASTKDR,KEYCBF,MAXATOM,SIGMATHK,SEXPOK,EDGEK,
     $                   LAINDHI,KRUDAUT, ZERO_RATES, POPMIN)
 
C******************************************************************************
C***  CALCULATION OF NEW NLTE POPULATION NUMBERS
C***  SOLUTION OF THE LINEAR RATE EQUATION SYSTEM
C***  Calling tree: WRSTART - POPZERO - NLTEPOP
C******************************************************************************

      DIMENSION EINST(NDIM,NDIM),CRATE(NDIM,NDIM),RRATE(NDIM,NDIM)
      DIMENSION RATCO(NDIM,NDIM)
      DIMENSION ENLTE(NDIM),EN(NDIM),NCHARG(NDIM),WEIGHT(NDIM)
      DIMENSION EION(NDIM),ELEVEL(NDIM)
      DIMENSION IONGRND(NDIM),RDIEL(NDIM),RAUTO(NDIM)
      DIMENSION KODRLOW(LASTKDR)
      DIMENSION ABXYZ(NATOM),KODAT(NATOM),NFIRST(NATOM),NLAST(NATOM)
      DIMENSION XLAMBDA(NF),FWEIGHT(NF)
      CHARACTER*4 CKEY
      CHARACTER LEVEL(NDIM)*10
      LOGICAL ZERO_RATES(N,ND)
 
C***  CHOOSE ADAEQUATE VERSION FOR MATRIX INVERSION
      CKEY = 'OWN '

C***  SET UP THE COEFFICIENT MATRICES CRATE AND RRATE FOR ALL ELEMENTS
      CALL COLLI (NDIM,N,ENLTE,TL,ENE,NCHARG,ELEVEL,EINST,CRATE,
     $            EION,COCO,KEYCBB,WEIGHT,ALTESUM,NATOM,NOM,KODAT,
     $            INDNUP,INDLOW,LASTIND,KONTNUP,KONTLOW,LASTKON,KEYCBF,
     $            IONGRND)
      CALL RADIO (NDIM,N,ENLTE,TL,WEIGHT,EION,ELEVEL,EINST,
     $            RRATE,XLAMBDA,FWEIGHT,XJC,NF,L,XJL,ND,SIGMAKI,
     $            ENE,NAUTO,MAXAUTO,LOWAUTO,WAUTO,EAUTO,AAUTO,
     $            RDIEL,RAUTO,IONAUTO,IONGRND,
     $            INDNUP,INDLOW,LASTIND,KONTNUP,KONTLOW,LASTKON,
     $            NATOM,MAXATOM,SIGMATHK,SEXPOK,EDGEK,KODAT,NFIRST,
     $            NLAST,NCHARG,LAINDHI,KRUDAUT)

C***  LOOP FOR EACH ELEMENT  -------------------------------------------
      DO 11 NA=1,NATOM
      NFIRNA=NFIRST(NA)
      NLANA=NLAST(NA)
      NDELP1 = NLANA - NFIRNA + 1

C***  Find Index of highest-populated level within current element
C***  from LTE popnumbers for later replacing this column by number 
C***  conservation 

      IMAXPOP = ISMAX(NDELP1, ENLTE(NFIRNA), 1)

C***  ADDING RADIATIVE AND COLLISIONAL TERMS INTO RATCO (=RATE COEFFICIENTS)
      DO 4 I=NFIRNA,NLANA
      ISHIFT=I-NFIRNA+1
      DO 4 J=NFIRNA,NLANA
      JSHIFT=J-NFIRNA+1
      RATCO(ISHIFT,JSHIFT)=CRATE(I,J)+RRATE(I,J)
    4 CONTINUE
 
C***  ADDING ADDITIONAL D-R TERMS INTO RATCO (=RATE COEFFICIENTS)
      DO 14 KDR=1,LASTKDR
      LOW=KODRLOW(KDR)
      IF ((LOW .LT. NFIRNA) .OR. (LOW .GT. NLANA)) GOTO 14
      LOSHIFT=LOW-NFIRNA+1
      NUP=IONGRND(LOW)
      NUSHIFT=NUP-NFIRNA+1
      RATCO(LOSHIFT,NUSHIFT)=RATCO(LOSHIFT,NUSHIFT)+RAUTO(LOW)
      RATCO(NUSHIFT,LOSHIFT)=RATCO(NUSHIFT,LOSHIFT)+RDIEL(LOW)
   14 CONTINUE

C***  DIAGONAL ELEMENTS : - SUM OF THE ROW (I.E. OVER COLUMN INDEX)
      DO 1 I=1,NDELP1
      SUM=.0
      DO 2 J=1,NDELP1
    2 SUM=SUM+RATCO(I,J)
    1 RATCO(I,I)=-SUM
 
C***  Check if all Rate Coefficients in one column are non-zero 
C***  (otherwise: the matrix is singular!)
C***  NOTE! Special, very small POPMIN must be adopted here! 
C***  It is intended that, for the WRSTART job, this POPMIN also 
C***  applies for the subsequent COMA branch

      POPMIN = 1.E-99
      CALL FLAG_ZERORATES (1, NDELP1, RATCO, NDIM, IMAXPOP, 
     >                     EN(NFIRNA), POPMIN, ZERO_RATES(NFIRNA,L))

C***  COLUMN IMAXPOP: NUMBER CONSERVATION
      DO 5 I=1,NDELP1
      RATCO(I,IMAXPOP)=1.
    5 CONTINUE

C***  If ZERO_RATES: Replace diagonal element by 1.0
C***  and the rest of this column by 0.0
C***  Note: BUG in the index of ZERO_RATES removed 14-Apr-2009 by wrh 
      DO J = 1, NDELP1
         IF (.NOT. ZERO_RATES(J+NFIRNA-1,L)) CYCLE
         DO I = 1, NDELP1
            RATCO(I,J) = .0
            RATCO(J,I) = .0
         ENDDO
         RATCO(J,J) = 1.
      ENDDO

ccc      if (l .eq. 20) call primat (ratco, ndelp1, ndim, level(nfirna))

C***  INVERSION OF RATE COEFFICIENT MATRIX RATCO
      CALL INV(NDELP1,NDIM,RATCO,CKEY)
 
C***  POPULATION NUMBERS EN(J) = LAST ROW OF INVERSE MATRIX
      AB=ABXYZ(NA)
      DO J=NFIRNA,NLANA
         EN(J)=AB*RATCO(IMAXPOP,J-NFIRNA+1)
         IF (ZERO_RATES(J-NFIRNA+1,L)) EN(J)=POPMIN
      ENDDO

   11 CONTINUE
C***  ENDLOOP  ---------------------------------------------------------
 
      RETURN
      END
