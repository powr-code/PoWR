      SUBROUTINE RADIO (NDIM,N,ENLTE,TL,WEIGHT,EION,ELEVEL,EINST,
     >                  RRATE,XLAMBDA,FWEIGHT,XJC,NF,L,XJL,ND,SIGMAKI,
     >                  ENE,NAUTO,MAXAUTO,LOWAUTO,WAUTO,EAUTO,AAUTO,
     >                  RDIEL,RAUTO,IONAUTO,IONGRND,
     >                  INDNUP,INDLOW,LASTIND,KONTNUP,KONTLOW,LASTKON,
     >                  NATOM,MAXATOM,SIGMATHK,SEXPOK,EDGEK,KODAT,
     >                  NFIRST,NLAST,NCHARG,LAINDHI,KRUDAUT)
 
C*******************************************************************************
C***  RADIATIVE RATES ARE CALCULATED AT ONE DEPTH POINT L FROM THE GIVEN
C***  RADIATION FIELD AND ATOMIC CROSS SECTIONS
C***  RADIATIVE TRANSITION RATES STORED IN MATRIX RRATE  ***********
C***  This is the LINEAR version of the rate equations (NO radiative 
C***  brackets, no ALO).
C***  Is is used if all GAMMA's are ZERO
C***  Compare the non-linear version in RADNET
C***  Calling tree: STEAL - POPZERO - NLTEPOP - RADIO 
C*******************************************************************************
 
C***  Dimension of the core-charge data locally provided here
      PARAMETER ( MAXATOMDIM = 26)

      DIMENSION KODAT(MAXATOM)
      DIMENSION SIGMATHK(MAXATOM,MAXATOM),SEXPOK(MAXATOM,MAXATOM)
      DIMENSION EDGEK(MAXATOM,MAXATOM)
      DIMENSION NFIRST(NATOM),NLAST(NATOM)
      DIMENSION NCHARG(NDIM)
      DIMENSION RRATE(NDIM,NDIM),EINST(NDIM,NDIM)
      DIMENSION ELEVEL(NDIM),IONGRND(NDIM)
      DIMENSION EION(NDIM),ENLTE(NDIM),WEIGHT(NDIM)
      DIMENSION RDIEL(NDIM),RAUTO(NDIM)
      DIMENSION XJC(ND,2),XJL(ND,2)
      DIMENSION INDNUP(LASTIND),INDLOW(LASTIND)
      DIMENSION KONTNUP(LASTKON),KONTLOW(LASTKON)
      DIMENSION SIGMAKI (NF,LASTKON)
      DIMENSION XLAMBDA(NF),FWEIGHT(NF)
      DIMENSION LOWAUTO(MAXAUTO),WAUTO(MAXAUTO),EAUTO(MAXAUTO)
     >         ,AAUTO(MAXAUTO),IONAUTO(MAXAUTO),KRUDAUT(MAXAUTO)
      DIMENSION KODATIND(MAXATOMDIM)

C***  CI : FACTOR IN SAHA EQUATION (MIHALAS, P. 113)
      DATA CI / 2.07E-16 /
C***  C1 = H * C / K    ( CM * KELVIN )
      DATA C1 / 1.4388 /
C***  C2 = 2 * H * C    ( H AND C IN CGS UNITS )
      DATA C2 / 3.9724E-16 /
C***  C3 = 4 * PI / H / C ( CGS UNITS )
      DATA C3 / 0.06327E18 /
 
C***  INITIALIZE ALL ELEMENTS OF MATRIX "RRATE"
      DO 1 J=1,N
      DO 1 I=1,N
      RRATE(I,J)=0.0
    1 CONTINUE

      IF (NATOM .GT. MAXATOMDIM) THEN
         WRITE (0,*) '*** ERROR: MAXATOMDIM TOO SMALL'
         STOP 'ERROR IN RADIO'
      ENDIF

C***  Now find for each NATOM the corresponding KODAT index
      DO NA=1, NATOM
         KODATIND(NA) = 0
         DO J = 1, MAXATOM
            IF (NA .EQ. KODAT(J)) KODATIND(NA) = J
         ENDDO
         IF (KODATIND(NA) .EQ. 0) THEN
            WRITE (0,*) '*** ERROR: ELEMENT NOT FOUND'
            STOP 'ERROR IN RADIO'
         ENDIF
         IF (KODATIND(NA) .GT. MAXATOMDIM) THEN
            WRITE (0,*) '*** ERROR: NCORECHARGE NOT FOUND'
            STOP 'ERROR IN RADIO'
         ENDIF
      ENDDO


C***  LINE TRANSITIONS  ********************************************************
      DO 11 IND=1,LASTIND
      NUP=INDNUP(IND)
      LOW=INDLOW(IND)
      WAVENUM=ELEVEL(NUP)-ELEVEL(LOW)
      W3=WAVENUM*WAVENUM*WAVENUM
 
C***  CHECK WHETHER THIS TRANSITION IS ONLY RUDIMENTAL
      IF (EINST(LOW,NUP) .EQ. -2.) THEN
C***        LINE IS RUDIMENTAL
            CALL XRUDI (XJ,WAVENUM,XJC,XLAMBDA,ND,NF,L)
            ELSE
C***        LINE IS EXPLICIT
            XJ=XJL(L,IND)
            ENDIF
 
C***  CALCULATION OF LINE RATES
      EMINDU=EINST(NUP,LOW)*XJ/C2/W3
      RRATE(LOW,NUP)=EMINDU*WEIGHT(NUP)/WEIGHT(LOW)
      RRATE(NUP,LOW)=EINST(NUP,LOW)+EMINDU
   11 CONTINUE
 
C***  CONTINUUM TRANSITIONS   **************************************************
C***  SIGMAKI = PRECALCULATED CROSS SECTION IN CM**2
      DO 21 KON=1,LASTKON
      NUP=KONTNUP(KON)
      LOW=KONTLOW(KON)
C***  EDGE = THRESHOLD ENERGY IN KAYSER *******
      EDGE = ELEVEL(NUP)+EION(LOW)-ELEVEL(LOW)
      EDGELAM=1.E8/EDGE
 
C***  RATE INTEGRAL
      SUM=.0
      REC=.0
C***  FIND EDGE FREQUENCY INDEX
      NFEDGE=ISRCHFGT (NF,XLAMBDA,1,EDGELAM) - 1
      DO 20 K=1,NFEDGE
      WAVENUM=1.E8/XLAMBDA(K)
      W3=WAVENUM*WAVENUM*WAVENUM
      XJCLK=XJC(L,K)
      SIGMA=SIGMAKI(K,KON)

C***  IONIZATION ***********************
      SUM=SUM+SIGMA*XJCLK*FWEIGHT(K)/WAVENUM

C***  RECOMBINATION ********************
      REC=REC +
     >    SIGMA*(XJCLK+C2*W3)*FWEIGHT(K)*EXP(-C1*WAVENUM/TL)/WAVENUM

   20 CONTINUE

C***  CALCULATION OF BOUND-FREE RATES
      RRATE(LOW,NUP)=C3*SUM
      RRATE(NUP,LOW)=C3*REC*ENLTE(LOW)/ENLTE(NUP)
   21 CONTINUE


C***  K-SHELL AND AUGER-IONISATION **********************************
      DO 40 NA=1,NATOM

         NLANA  = NLAST(NA)
         NFIRNA = NFIRST(NA)
         LASTISTATE = -1

         DO 44 LOW=NFIRNA,NLANA
C***        Auger-ionisation needs 4 electrons as minimum.
C***        We assume that k-shell absortption is always followed
C***        by auto-ionization (Auger). In reality, there is radiative
C***        decay (emission of an X-ray photon) as competing process.
C***        For light elements (core-charge Z less than 33), Auger ionization
C***        is dominating (cf. Baum, Diplomarbeit).

C***        Is there a higher ionization stage available for this element?
            IF (NCHARG(LOW) .EQ. NCHARG(NLANA)) EXIT

            ISTATE = NCHARG(LOW) + 1

C***        ARE THERE K_SHELL DATA FOR CURRENT ELEMENT/ION?
            IF (SIGMATHK(NA,ISTATE) .LE. 0.) CYCLE

C***        Has the Rate Integral not been calculated for this ion?
            IF (ISTATE .NE. LASTISTATE) THEN
               LASTISTATE = ISTATE
               EDGELAM = 1.E8 / EDGEK(NA,ISTATE)
               NFEDGEK = ISRCHFGT(NF,XLAMBDA,1,EDGELAM)-1
C***           RATE INTEGRAL
               SUM=0.
               DO K=1, NFEDGEK
                  WAVENUM = 1.E8 / XLAMBDA(K)
                  XJCLK = XJC(L,K)
                  CALL KSIGMA(SIGMAK, SIGMATHK(NA,ISTATE), 
     >                     EDGEK(NA,ISTATE), WAVENUM, SEXPOK(NA,ISTATE))
                  SUM = SUM + SIGMAK * XJCLK * FWEIGHT(K) / WAVENUM
               ENDDO
               RATE = C3 * SUM
            ENDIF

            IF (KODATIND(NA)-NCHARG(LOW) .GE. 4) THEN
               IF (NCHARG(LOW) .GT. NCHARG(NLANA)-2) GOTO 47
               NLOW = LOW
***            FIND UPPER LEVELINDEX I: GROUND_STATE WITH NCHARGE(NLOW)+2
               DO I = NLOW+1, NLANA
                 IF (NCHARG(LOW) .EQ. NCHARG(I)-2)  THEN
                    RRATE(LOW,I) = RATE + RRATE(LOW,I)
                    GOTO 44
                 ENDIF
               ENDDO
            ENDIF
   47       CONTINUE

C***  K-shell ionisation without Auger-process happens when the ion
C***  has exactly three electrons. (Two or one electrons are already
C***  considered as normal ionisation.) Note that we go through this branch
C***  also in case of more than 3 electrons (by means of the ".GE."), if
C***  in the previous Auger-Block an upper limit with ncharg+2 was not found
C**   (i.e. does not exist in the model atom).
C***  The upper level of the Auger process should be the first *excited*
C***  level of the helium-like ion. However, as that ion (i.e. C V, N VI,
C***  O VII) is represented by only one level in our standard model atom,
C***  the ground level is taken instead.

            IF (KODATIND(NA)-NCHARG(LOW) .GE. 3) THEN
C***           FIND UPPER LEVELINDEX I: GROUND_STATE WITH NCHARGE(NLOW)+1
               DO I = NLOW+1, NLANA
                  IF (NCHARG(LOW) .EQ. NCHARG(I)-1)  THEN
                     RRATE(LOW,I) = RATE + RRATE(LOW,I)
                    GOTO 44
                  ENDIF
               ENDDO
            ENDIF
   44       CONTINUE

   40 CONTINUE               


C***  DR: DIELECTRONIC RECOMBINATION / AUTOIONIZATION  *************************
      DO 30 LOW=1,N
      RDIEL(LOW) = .0
      RAUTO(LOW) = .0
   30 CONTINUE

C***  SUM OVER ALL AUTO-IONIZATION LEVELS (INDEX: J)
      DO 31 J=1,NAUTO
      LOW=LOWAUTO(J)
      NUP=IONAUTO(J)
C***  WAVENUMBER OF STABILIZING TRANSITION
      WSTABIL=EION(LOW)-ELEVEL(LOW)+EAUTO(J)
      WSTAB3=WSTABIL*WSTABIL*WSTABIL

C***  GETTING LINE RADIATION FIELD AT WAVELENGTH OF STABILIZING TRANSITION
C***  IF NO "DRLINE" CARD IS ACTIVE, LAINDHI IS EQUAL TO LASTIND
      IF ( (LAINDHI.EQ.LASTIND) .OR. (KRUDAUT(J).EQ.1) ) THEN
C***     ALL DR-LINES ARE SET RUDIMENTAL (NO "DRLINE" CARD) OR THE 
C***     SPECIFIED DR-TRANSITION IS A RUDIMENTAL IN THE DATOM-FILE
         CALL XRUDI (XJLSTAB,WSTABIL,XJC,XLAMBDA,ND,NF,L)
      ELSE
C***     "DRLINE" CARD IS ACTIVE AND SPECIFIED DR-LINE IS A NON-RUDIMENTAL
         XJLSTAB = XJL(L,LASTIND+J)
      ENDIF

      DRINDU=AAUTO(J)*XJLSTAB/C2/WSTAB3
C***  AUTOIONIZATION RATE
      DRLOWUP=DRINDU*WAUTO(J)/WEIGHT(LOW)
C***  FUNCTION PHI(T) FROM SAHA'S EQUATION
      SAHAPHI=WAUTO(J)/WEIGHT(NUP)*CI/TL/SQRT(TL)
     *                *EXP(-C1*(EAUTO(J)-ELEVEL(NUP))/TL)
C***  DIELECTRONIC RECOMBINATION RATE
      DRUPLOW=ENE*SAHAPHI*(AAUTO(J)+DRINDU)

C***  ADD D-R TERMS INTO THE ONE-DIMENSIONAL ARRAYS
C***  (NOTE: IN THIS VERSION TH EUPPER LEVLELS OF THE DR-TRANSITIONS
C***         ARE ASSUMED TO BE THE GROUND STATES OF THE PARENT IONS)
      RDIEL(LOW) = RDIEL(LOW) + DRUPLOW
      RAUTO(LOW) = RAUTO(LOW) + DRLOWUP
   31 CONTINUE

      RETURN
      END
