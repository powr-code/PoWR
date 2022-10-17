      SUBROUTINE RADNET(NDIM,N,ENLTE,TL,WEIGHT,NCHARG,EION,ELEVEL,EINST,
     $                  SLNEW,XRED,XBLUE,EN,RRATE,XLAMBDA,FWEIGHT,
     $                  XJC,NF,XJL,SIGMAKI,ENTOTL,NOTEMP,IONGRND,
     $                  NAUTO,MAXAUTO,LOWAUTO,WAUTO,EAUTO,AAUTO,IONAUTO,
     $                  DRRATEN,RDIEL,RAUTO,DRJLW,DRJLWE,DRLJW,
     $                  INDNUP,INDLOW,LASTIND,KONTNUP,KONTLOW,LASTKON,
     $                  NFEDGE,NATOM,MAXATOM,SIGMATHK,SEXPOK,EDGEK,
     $                  KODAT,NFIRST,NLAST,LAINDHI,KRUDAUT,L,ND, 
     >                  NRB_CONT, EXPFAC, WCHARM)
C*******************************************************************************
C***  RADIATIVE RATE COEFFICIENT MATRIX RRATE IS CALCULATED AT DEPTH POINT L
C***  FROM THE GIVEN RADIATION FIELD
C***  NOTE THAT THIS SUBROUTINE CALCULATES NETTO RATE COEFFICIENTS
C***  FOR NON-RUDIMENTAL LINE TRANSITIONS (Net Radiative Brackets)
C***  Calling tree: STEAL - LINPOP - COMA - RADNET 
C*******************************************************************************

C***  Dimension of the core-charge data locally provided here
      PARAMETER ( MAXATOMDIM = 26)
      
c      character*10 level(n) 
      DIMENSION KODAT(MAXATOM)
      DIMENSION SIGMATHK(MAXATOM,MAXATOM),SEXPOK(MAXATOM,MAXATOM)
      DIMENSION EDGEK(MAXATOM,MAXATOM)
      DIMENSION NFIRST(NATOM),NLAST(NATOM)
      DIMENSION RRATE(NDIM,NDIM),EINST(NDIM,NDIM)
      DIMENSION XJC(NF), XJL(LASTIND)
      DIMENSION SLNEW(LASTIND), XRED(LASTIND), XBLUE(LASTIND)
      DIMENSION NCHARG(NDIM),ELEVEL(NDIM)
      DIMENSION EION(NDIM),ENLTE(NDIM),WEIGHT(NDIM),EN(NDIM)
      DIMENSION DRRATEN(NDIM),IONGRND(NDIM)
      DIMENSION SIGMAKI (NF,LASTKON)
      DIMENSION XLAMBDA(NF),FWEIGHT(NF), EXPFAC(NF)
      DIMENSION LOWAUTO(MAXAUTO),WAUTO(MAXAUTO),EAUTO(MAXAUTO)
      DIMENSION AAUTO(MAXAUTO),IONAUTO(MAXAUTO),KRUDAUT(MAXAUTO)
      DIMENSION INDNUP(LASTIND),INDLOW(LASTIND)
      DIMENSION KONTNUP(LASTKON),KONTLOW(LASTKON),NFEDGE(LASTKON)
      DIMENSION RDIEL(N),RAUTO(N),DRJLW(N),DRJLWE(N),DRLJW(N)
      DIMENSION KODATIND(MAXATOMDIM)
      DIMENSION WCHARM(ND,NF)

      LOGICAL NRB_CONT(LASTKON)
      LOGICAL NOTEMP

C***  CI : FACTOR IN SAHA EQUATION (MIHALAS, P. 113)
      DATA CI / 2.07E-16 /
C***  C1 = H * C / K    ( CM * KELVIN )
      DATA C1 / 1.4388 /
C***  C2 = 2 * H * C    ( H AND C IN CGS UNITS )
      DATA C2 / 3.9724E-16 /
C***  C3 = 4 * PI / H / C  (CGS UNITS)
      DATA C3 / 6.32684E16 /
C***  PI8 = 8 * PI
      DATA PI8 / 25.1327412288 /

      SQRTL=SQRT(TL)
      T32 = TL * SQRTL
      ENE = EN(N+1) * ENTOTL

C***  INITIALIZE ALL ELEMENTS OF MATRIX "RRATE"
      DO 1 J=1, N
         DO 1 I=1, N
            RRATE(I,J) = .0
    1 CONTINUE

      IF (MAXATOM .GT. MAXATOMDIM) THEN
         WRITE (0,*) '*** ERROR: MAXATOMDIM TOO SMALL'
         STOP 'ERROR IN RADNET'
      ENDIF

C***  Now find for each NATOM the corresponding KODAT index
      DO NA=1, NATOM
         KODATIND(NA) = 0
         DO J = 1, MAXATOM
            IF (NA .EQ. KODAT(J)) KODATIND(NA) = J
         ENDDO
         IF (KODATIND(NA) .EQ. 0) THEN
            WRITE (0,*) '*** ERROR: ELEMENT NOT FOUND'
            STOP 'ERROR IN RADNET'
         ENDIF
         IF (KODATIND(NA) .GT. MAXATOMDIM) THEN
            WRITE (0,*) '*** ERROR: NCORECHARGE NOT FOUND'
            STOP 'ERROR IN RADNET'
         ENDIF
      ENDDO


C***  LINE TRANSITIONS  ************************************************
      DO 11 IND=1,LASTIND
      NUP=INDNUP(IND)
      LOW=INDLOW(IND)
      WAVENUM=ELEVEL(NUP)-ELEVEL(LOW)
      W3=WAVENUM*WAVENUM*WAVENUM
 
      IF (EINST(LOW,NUP) .EQ. -2.) THEN
C***     TRANSITION IS ONLY RUDIMENTAL
         IF (EINST(NUP,LOW) .EQ. .0) THEN
C***        Zero f-value --> zero rates
            RRATE(LOW,NUP) = .0
            RRATE(NUP,LOW) = .0
         ELSE  
C***        RADIATION FIELD FROM INTERPOLATION OF CONT.
C***        Note: ND and L = 1 because XJC is one-dimensional here
            CALL XRUDI (XJ,WAVENUM,XJC,XLAMBDA,1,NF,1)
            EMINDU=EINST(NUP,LOW)*XJ/C2/W3
            RRATE(LOW,NUP)=EMINDU*WEIGHT(NUP)/WEIGHT(LOW)
            RRATE(NUP,LOW)=EINST(NUP,LOW)+EMINDU
         ENDIF
      ELSE 
C***     NON-RUDIMENTAL TRANSITIONS: 
         IF (XRED(IND) .GE. XBLUE(IND) .OR. SLNEW(IND) .LT. 1.E-100) THEN
C***        NO NETTO BRACKET IN CASE OF NON-EXISTING SCHARMER CORE
C***        Note: Iron lines always have NRB by setting XBLUE=1 in SUBR. LCORE
            XJ=XJL(IND)
            EMINDU=EINST(NUP,LOW)*XJ/C2/W3
            RRATE(LOW,NUP)=EMINDU*WEIGHT(NUP)/WEIGHT(LOW)
            RRATE(NUP,LOW)=EINST(NUP,LOW)+EMINDU
         ELSE
C***        NET RADIATIVE BRACKETS
            RRATE(LOW,NUP)=.0
            RRATE(NUP,LOW)=EINST(NUP,LOW)*(1.-XJL(IND)/SLNEW(IND))
         ENDIF
      ENDIF
   11 CONTINUE
 
C******************************************************************
C***  CALCULATION OF BOUND-FREE RATES (2 Branches: NRBs / no NRBs)
C******************************************************************

C***  SIGMAKI = PRECALCULATED CROSS SECTION IN CM**2
      DO 21 KON=1,LASTKON
      NUP=KONTNUP(KON)
      LOW=KONTLOW(KON)
C***  EDGE = THRESHOLD ENERGY IN KAYSER *******
      EDGE = ELEVEL(NUP)+EION(LOW)-ELEVEL(LOW)
      EDGELAM=1.E8/EDGE
 
      KONEDGE=NFEDGE(KON)
      REC=.0

C***************************************************************
C***  Criterion to decide if this continuum is treated with NRB
cc      NRB_CONT(KON) = EN(LOW) .GT. 1.E-5 
cc     >                .AND. WCHARM(L,KONEDGE) .GT. .0
        NRB_CONT(KON) = .TRUE.

C**************************************************************

C***  Branch for NET RADIATIVE BRACKET
      IF (NRB_CONT(KON)) THEN

C***  RATE INTEGRAL
      DO K = 1, KONEDGE
        WAVENUM=1.E8/XLAMBDA(K)
        W2=WAVENUM*WAVENUM
        W3=W2*WAVENUM
        XJCLK=XJC(K) 
ccc        XJCLK=AMAX1(XJC(K), .0)
        SIGMA=SIGMAKI(K,KON)
        HNUEKT = C1 * (WAVENUM-EDGE) / TL
C***    HNUEKT .GT. 600.  =>   EXP(-HNUEKT)->0.
        IF (HNUEKT .GT. 500.) CYCLE
C***    CALCULATE BOUND-FREE SOURCE FUNCTION FOR TRANSITION LOW-UP ONLY
        EXFAC = EXP(-HNUEKT)
        G = WEIGHT(LOW)/WEIGHT(NUP) * EXFAC * ENE * CI / T32 
        DIVISOR = C2 * W3 * EN(NUP) * G
        IF (DIVISOR .LT. 1.E-100) CYCLE
        SBFINV = (EN(LOW) - EN(NUP) * G) / DIVISOR 
        REC = REC + 
     >    SIGMA * EXFAC * W2 * (1. - XJCLK * SBFINV) * FWEIGHT(K)  
      ENDDO
C***  Net rates:
      RRATE(LOW,NUP) = .0
      RRATE(NUP,LOW) = 
     >   PI8 * ENE / T32 * CI * WEIGHT(LOW) / WEIGHT(NUP) * REC 

      ELSE

C***  Branch without NET RADIATIVE BRACKET
      SUM=.0
C***  RATE INTEGRAL
      DO K = 1, KONEDGE
        WAVENUM=1.E8/XLAMBDA(K)
        W3=WAVENUM*WAVENUM*WAVENUM
        XJCLK=XJC(K)
ccc        XJCLK=AMAX1(XJC(K), .0)
        SIGMA=SIGMAKI(K,KON)
C***    IONIZATION ***********************
           
        FAC = SIGMA * FWEIGHT(K) / WAVENUM
        SUM = SUM + XJCLK * FAC

C***    RECOMBINATION ********************
        REC=REC +
     $    (XJCLK+C2*W3) * EXPFAC(K) * FAC
      ENDDO

C***  Bound-free rates (no NRBs)
      RRATE(LOW,NUP)=C3*SUM
      RRATE(NUP,LOW)=C3*REC*ENLTE(LOW)/ENLTE(NUP)

      ENDIF


   21 CONTINUE
C***  End-of-loop over the bound-free continua

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

C***        Is there a higher ionization stake available for this element?
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
                  XJCLK = XJC(K)
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
                    RRATE(LOW,I) = RATE
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

 
C***  DR: DIELECTRONIC RECOMBINATION / AUTOIONIZATION  *****************
C***  INITIALIZE SOME VECTORS WHICH WILL BE USED LATER
      DO 30 LOW=1, N
      RDIEL(LOW) = .0
      RAUTO(LOW) = .0
      DRRATEN(LOW) = .0
      DRJLW(LOW) = .0
      DRJLWE(LOW) = .0
      DRLJW(LOW) = .0
   30 CONTINUE

C***  SUM OVER ALL AUTO-IONIZATION TRANSITIONS (INDEX: J)
      DO 31 J=1, NAUTO
      LOW=LOWAUTO(J)
      NUP=IONAUTO(J)
C***  WAVENUMBER OF STABILIZING TRANSITION
      WSTABIL = EION(LOW) - ELEVEL(LOW) + EAUTO(J)
      WSTAB3 = WSTABIL * WSTABIL * WSTABIL

C***  GETTING LINE RADIATION FIELD AT WAVELENGTH OF STABILIZING TRANSITION
C***  IF NO "DRLINE" CARD IS ACTIVE, LAINDHI IS EQUAL TO LASTIND
      IF ( (LAINDHI.EQ.LASTIND) .OR. (KRUDAUT(J).EQ.1) ) THEN
C***     ALL DR-LINES ARE SET RUDIMENTAL (NO "DRLINE" CARD) OR THE 
C***     SPECIFIED DR-TRANSITION IS A RUDIMENTAL IN THE DATOM-FILE
C***     Note: ND and L = 1 because XJC is one-dimensional here
         CALL XRUDI (XJLSTAB,WSTABIL,XJC,XLAMBDA,1,NF,1)
      ELSE
C***     "DRLINE" CARD IS ACTIVE AND SPECIFIED DR-LINE IS A NON-RUDIMENTAL
         XJLSTAB = XJL(LASTIND+J)
      ENDIF

      DRINDU = AAUTO(J) *XJLSTAB / C2 / WSTAB3
C***  AUTOIONIZATION RATE
      DRLOWUP = DRINDU * WAUTO(J) / WEIGHT(LOW)
C***  FUNCTION PHI(T) FROM SAHA'S EQUATION
      SAHAPHI=WAUTO(J)/WEIGHT(NUP)*CI/TL/SQRTL
     *                *EXP(-C1*(EAUTO(J)-ELEVEL(NUP))/TL)
C***  DIELECTRONIC RECOMBINATION RATE
      DRUPLOW = ENTOTL * EN(N+1) * SAHAPHI * (AAUTO(J)+DRINDU)

      IF (.NOT. NOTEMP) THEN
C***      STORE "DRRATEN" FOR LATER USE IN SUBR. DERIVT:
          DRRATEN(LOW) = DRRATEN(LOW) + DRUPLOW * EAUTO(J)
C***      STORE SOME DR-RATES FOR LATER USE IN SUBR. TEMPEQ:
          DRJLW(LOW) = DRJLW(LOW) + DRUPLOW * WSTABIL
          DRJLWE(LOW) = DRJLWE(LOW) + DRUPLOW * WSTABIL * EAUTO(J)
          DRLJW(LOW) = DRLJW(LOW) + DRLOWUP * WSTABIL
      ENDIF

C***  ADD D-R TERMS INTO THE ONE-DIMENSIONAL DR-ARRAYS
C***  (NOTE: IN THIS VERSION THE UPPER LEVELS OF THE DR-TRANSITIONS
C***         ARE ASSUMED TO BE THE GROUND STATES OF THE PARENT IONS)
      RDIEL(LOW) = RDIEL(LOW) + DRUPLOW
      RAUTO(LOW) = RAUTO(LOW) + DRLOWUP
   31 CONTINUE

      RETURN
      END
