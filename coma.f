      SUBROUTINE COMA (CRATE,RRATE,RATCO,DM,N,NRANK ,NDIM,V1,ABXYZ,
     $    ENLTE,TL,ENE,NCHARG,ELEVEL,EINST,EION,WEIGHT,ALTESUM,XLAMBDA,
     $    FWEIGHT,XJC,NF,L,XJL,ND,XJLAPP,SLOLD,LASTIND,INDLOW,
     $      INDNUP,NOM,NATOM,KODAT,NFIRST,NLAST,PHI,PWEIGHT,DELTAX,XMAX,
     $      NFL,OPAC,SCNEW,DOPA,DETA,OPAL,SLNEW,DOPAL,DETAL,SIGMAKI,
     $      ETAC,NFEDGE,EXPFAC,SCOLIND,SCNEIND,OPACIND,SIGMAFF,MAXION,
     $      NOTEMP,TLOLD,KONTLOW,KONTNUP,LASTKON,RUDLINE,IONGRND,
     $      XRED,XBLUE,WCHARM,EN,RSTAR,SCOLD,XJCAPP,VDOP,COCO,
     $      KEYCBB, NRB_CONT, ZERO_RATES, POPMIN, 
     $      IONAUTO,NAUTO,MAXAUTO,LOWAUTO,WAUTO,EAUTO,AAUTO,DRRATEN,
     $      RDIEL,RAUTO,DRJLW,DRJLWE,DRLJW,IBLENDS,MAXLAP,XLAMZERO,
     $      KODRNUP,KODRLOW,LASTKDR,KEYCBF, OPALOLD,
     $      BETA,PHIL,NBLENDS,BROYDEN,ATEST,BTEST,MAXATOM,
     $      SIGMATHK,SEXPOK,EDGEK,XDATA,RL,
     $      XJCLP1, OPAC1, RADIUS, ITNEL, TEFF, OPATHOM, 
     $      LAINDHI, KRUDAUT, LEVEL, 
     >   WFELOW, WFENUP, EN1, BDIAG,  
     >   FERATLU, FERATUL, LASTFE, FERATLU0, FERATUL0, 
C*** Quantities for fine-frequency grid 
     >   SFINE_OLD, SFINE_NEW, MAXFINE, KONTHLP, MAXIND, XKC, XKC2,
     >   SIGMA1I, NLINE, LINEINDEX, XLAMSOR, XLAMMIN, XLAMMAX,
     >   NUPACT, LOWACT, BLASERL, ETAL, LIND, LINDS,
C***  IRON
     >   INDRB, IFRBSTA, IFRBEND, 
     >   VDOPFE, DXFE, XLAM0FE,
     >   INDFEACT, MAXFEACT, BFECHECK, BFEWING,
     >   DFEINDR, SIGMAFE, OPAFE, ETAFE, IFENUP, IFELOW,
     >   INDEXMAX, BFEMODEL, BNUEFE, 
C***
     >   XLAMBDA2, NF2, NDDIM, VELO1, XKMIN, XKMAX, XKMID, XKRED, 
     >   ALPHA, SEXPO, ADDCON1, ADDCON2, ADDCON3, IGAUNT, 
     >   XKRED_CORE, XKBLUE_CORE, 
     >   BXJLAPPNEW, BXJCAPPNEW, BNEWOPER, BXJLAPPCORE,
     >   XJL_PLOTDATA, XJC_PLOTDATA_I, XJC_PLOTDATA_L, 
     >   IPLOT_XJLAPP, IPLOT_XJCAPP, LPLOT_XJCAPP, NITER_PLOT_JAPP, 
     >   BPLOTAPP, PWEIGHTCL, WS, FWTEST, 
     >   IWARN_NEG_XJCAPP, IWARN_NEG_XJLAPP, 
     >   XJCAPPNEW, XJLAPPNEW, 
     >   GAMMAC, GAMMAL, GAMMAR, 
     >   XLAM_FINE_START, XLAM_FINE_END, IMAXPOP, bBLOCKINVERSION,
C***  New Fine-spaced WCHARM handling
     >   IFF_MAX, IFF_MAX_MS, FF_INFO, IFF_DK, IFF_WCHARM, WCHARM_FINE, 
     >   IFF_N_MS, bFFASSET,
C***  FERAT correction     
     >   CORRS, DEXFAC, bFeTCORR)
C*******************************************************************************
C***  THIS ROUTINE SETS UP THE RATE COEFFICIENT MATRIX RATCO
C***  AND ITS VECTOR DERIVATIVE DM (LINEARIZED MATRIX)
C***  AND THE RIGHT-HAND SIDE VECTOR V1
C***  CALLED FROM: SUBROUTINE LINPOP
C*******************************************************************************

      INTEGER, INTENT(IN) :: NDIM, NRANK, LASTKDR, LASTFE, LASTIND, 
     >                       ND, NF, N, NATOM, MAXATOM

      REAL, DIMENSION(NDIM,NDIM) :: CRATE, RRATE
      REAL, DIMENSION(NDIM) :: RDIEL, RAUTO, ELEVEL
      INTEGER, DIMENSION(NDIM) :: NCHARG, IONGRND
      REAL, DIMENSION(NRANK,NRANK) :: RATCO, DM
      REAL, DIMENSION(NRANK) :: EN, V1
      INTEGER, DIMENSION(LASTKDR) :: KODRLOW
      REAL, DIMENSION(NATOM) :: ABXYZ
      INTEGER, DIMENSION(NATOM) :: NFIRST, NLAST, IMAXPOP
      INTEGER, DIMENSION(MAXATOM) :: KODAT
      REAL, DIMENSION(ND,NF) :: XJC 
      REAL, DIMENSION(NF) :: XJCAPP
      INTEGER, DIMENSION(LASTIND) :: INDLOW, INDNUP
      REAL, DIMENSION(LASTIND) :: DEXFAC
      REAL, DIMENSION(LASTFE) :: FERATLU, FERATUL, FERATLU0, FERATUL0
      LOGICAL :: NOTEMP, BROYDEN, BFIRSTITER, BFERATE, bBLOCKINVERSION, 
     >           bFFASSET, BXJLAPPNEW, BXJCAPPNEW, BNEWOPER, BPLOTAPP, 
     >           bFeTCORR
      CHARACTER(8) :: NAME
      LOGICAL, DIMENSION(LASTKON) :: NRB_CONT 
      LOGICAL, DIMENSION(N, ND) :: ZERO_RATES
      LOGICAL, DIMENSION(LASTIND) :: BDIAG
      REAL, DIMENSION(ND) :: CORRS

C*** for test output
      CHARACTER(10), DIMENSION(NDIM) :: LEVEL
      
      REAL :: CORRT, WAV0
      
      REAL, PARAMETER :: C1 = 1.4388  !C1 = H * C / K    (CM * KELVIN)

      INTEGER, SAVE :: LASTL
      DATA LASTL / 0 /

C***  CHECK IF CURRENT DEPTHPOINT IS ENCOUNTERED FOR THE FIRST TIME
      BFIRSTITER = LASTL .NE. L
      LASTL = L

      NPLUS1=N+1
C***  CALCULATE CONTINUUM OPACITIES AND EMISSIVITIES FROM CURRENT POPNUMBERS.
C***  ONLY TRUE OPACITIES ARE ACCOUNTED FOR.
      RNEL=EN(NPLUS1)
      ENTOTL=ENE/RNEL

      CALL COOPFRQ (NF,OPAC,ETAC,XLAMBDA,EXPFAC,SIGMAKI,N,NCHARG,
     $              WEIGHT,ELEVEL,EION,NFEDGE,EN,NOM,RSTAR,ENTOTL,
     $              RNEL,TL,SIGMAFF,MAXION,RL,XDATA,
     $              SIGMATHK,SEXPOK,EDGEK,KODAT,MAXATOM,
     $              KONTNUP,KONTLOW,LASTKON,OPATHOM)

C***  CALCULATE NEW CONT. SOURCE FUNCTION AND SCHARMER'S RADIATION FIELD
      CALL SETXJC (XJCAPP,XJC,OPAC,ETAC,SCNEW,SCOLD,WCHARM,
     $             XLAMBDA,L,NF,ND,TL,TLOLD,NOTEMP)

C***  CALCULATE LINE RADIATION FIELD WITH APPROXIMATE LAMBDA OPERATOR TERMS
      CALL SETXJL (LASTIND,INDLOW,INDNUP,XRED,XBLUE,OPACIND,
     $             SCNEIND,SCOLIND,SLNEW,SLOLD,OPAL,XJLAPP,
     $             NF,XLAMBDA,SCNEW,OPAC,OPALOLD,ITNEL,LAINDHI,
     $             NFL,PHI,PWEIGHT,NDIM,EINST,ELEVEL,EN,WEIGHT,ND,XJL,
     $             ENTOTL,RSTAR,VDOP,DELTAX,XMAX,L,TL,TLOLD,NOTEMP,
     $             IBLENDS,MAXLAP,XLAMZERO,BETA,PHIL,ATEST,NBLENDS,
     $             KRUDAUT,MAXAUTO,WFELOW,WFENUP,EN1,BDIAG)

      IF (BPLOTAPP .AND. ITNEL .EQ. NITER_PLOT_JAPP) THEN
        CALL PREPLOTAPP(ND, NF, L,
     >    XJL, XJLAPP, XJL_PLOTDATA, XRED, XBLUE,
     >    XJC, XJCAPP, WCHARM, 
     >    XJC_PLOTDATA_I, XJC_PLOTDATA_L,
     >    IPLOT_XJLAPP, IPLOT_XJCAPP, LPLOT_XJCAPP, 
     >    0)
      ENDIF
      
C***  New Subr. SETXJFINE to calculate accurate Approximate Radiation fields
C***    XJCAPP and XJLAPP
      IF (bFFASSET .AND. 
     >    ((BXJLAPPNEW .AND. (GAMMAL .GT. 0. .OR. GAMMAR .GT. 0.)) .OR.
     >     (BXJCAPPNEW .AND. GAMMAC .GT. 0.)) ) THEN
        CALL SETXJFINE (SFINE_OLD, SFINE_NEW, MAXFINE, ITNEL, L, 
     >                 ND, NDDIM, NDIM, N, VDOP, VELO1, 
     >                 ENTOTL, EN, RSTAR, TL, RNEL, NCHARG, 
     >                 WEIGHT, ELEVEL, EION, EINST, NATOM, 
     >                 KONTHLP, MAXIND, BXJLAPPCORE,
C***    CONTINUA
     >                 XLAMBDA, XLAMBDA2, NF, NF2, 
     >                 XKC, XKC2, ALPHA, 
     >                 SEXPO, ADDCON1, ADDCON2, ADDCON3, IGAUNT, 
     >                 SIGMA1I, KONTLOW, KONTNUP, LASTKON,
     >                 FWTEST, 
C***    LINES
     >                 NLINE, LINEINDEX, INDLOW, INDNUP, XLAMSOR, 
     >                 XLAMMIN, XLAMMAX, 
     >                 NUPACT, LOWACT, BLASERL, OPAL, ETAL, 
     >                 LIND, LINDS, MAXLAP, 
     >                 XKMIN, XKMAX, XKMID, XKRED, 
     >                 XRED, XBLUE, 
     >                 XKRED_CORE, XKBLUE_CORE,
     >                 
     >                 
C***    IRON
     >                 INDRB, IFRBSTA, IFRBEND, LASTFE, 
     >                 VDOPFE, DXFE, XLAM0FE,
     >                 INDFEACT, MAXFEACT, BFECHECK, BFEWING,
     >                 DFEINDR, SIGMAFE, OPAFE, ETAFE, IFENUP, IFELOW,
     >                 INDEXMAX, BFEMODEL, BNUEFE, 
C***    Approximate radiation fields
     >                 LASTIND, XJLAPP, XJL,
     >                 BXJLAPPNEW, BXJCAPPNEW, BNEWOPER, 
     >                 XJL_PLOTDATA, XJC_PLOTDATA_I, XJC_PLOTDATA_L, 
     >                 IPLOT_XJLAPP, IPLOT_XJCAPP, 
     >                 LPLOT_XJCAPP, NITER_PLOT_JAPP, 
     >                 PWEIGHTCL, WS, XJC, XJCAPP, XJCAPPNEW, WCHARM, 
     >                 IWARN_NEG_XJCAPP, IWARN_NEG_XJLAPP, 
     >                 XLAM_FINE_START, XLAM_FINE_END, 
     >                 OPAC, ETAC, 
C***  New Fine-spaced WCHARM handling
     >                 IFF_MAX, IFF_MAX_MS, FF_INFO, IFF_DK, 
     >                 IFF_WCHARM, WCHARM_FINE, IFF_N_MS, GAMMAC)
      ENDIF

C***  Store Old XJCAPP and XJLAPP : MODE = 1
      IF (BPLOTAPP .AND. ITNEL .EQ. NITER_PLOT_JAPP) THEN
        CALL PREPLOTAPP(ND, NF, L,
     >    XJL, XJLAPP, XJL_PLOTDATA, XRED, XBLUE,
     >    XJC, XJCAPP, WCHARM, 
     >    XJC_PLOTDATA_I, XJC_PLOTDATA_L,
     >    IPLOT_XJLAPP, IPLOT_XJCAPP, LPLOT_XJCAPP, 
     >    1)
      ENDIF

C***  SETUP THE COLLISIONAL AND RADIATIVE RATE COEFFICIENTS
      CALL COLLI (NDIM,N,ENLTE,TL,ENE,NCHARG,ELEVEL,EINST,CRATE,
     $            EION,COCO,KEYCBB,WEIGHT,ALTESUM,NATOM,NOM,KODAT,
     $            INDNUP,INDLOW,LASTIND,KONTNUP,KONTLOW,LASTKON,KEYCBF,
     $            IONGRND)
C***  RADIATIVE RATES ARE CALCULATED WITH THE MODIFIED RADIATION FIELD
C***  NOTE THAT THE ARRAYS XJCAPP AND XJLAPP ARE ONE-DIMENSIONAL :
      CALL RADNET(NDIM,N,ENLTE,TL,WEIGHT,NCHARG,EION,ELEVEL,EINST,
     $            SLNEW,XRED,XBLUE,EN,RRATE,XLAMBDA,FWEIGHT,
     $            XJCAPP,NF,XJLAPP,SIGMAKI,ENTOTL,NOTEMP,IONGRND,
     $            NAUTO,MAXAUTO,LOWAUTO,WAUTO,EAUTO,AAUTO,IONAUTO,
     $            DRRATEN,RDIEL,RAUTO,DRJLW,DRJLWE,DRLJW,
     $            INDNUP,INDLOW,LASTIND,KONTNUP,KONTLOW,LASTKON,
     $            NFEDGE,NATOM,MAXATOM,SIGMATHK,SEXPOK,EDGEK,KODAT,
     $            NFIRST,NLAST,LAINDHI,KRUDAUT, L, ND, NRB_CONT, 
     >            EXPFAC, WCHARM)
 
C***  THE RADIATIVE RATES FOR THE IRON SUPERLINES ARE PRE-CALCULATED
C***      IN PROGRAM COLI AND NOW READ FROM THE MODEL FILE.
C***      The rates calculated by RADNET are only approximate, but are 
C***      differentially applied to the pre-calculated rates 
C***      in order to account for the Scharmer pre-estimate. 

C***  FIRST CHECK FOR PRESENCE OF IRON RATES 
C***       (IMPORTANT FOR FIRST STEAL AFTER WRSTART)
      IF (BFIRSTITER .AND. L .EQ. 1) THEN 
          CALL READMS (3, FERATLU, LASTFE, 'FERLU  1', IERR)
          BFERATE = IERR .GE. 0                 
          IF (.NOT.BFERATE .AND. LASTFE .GT. 0) WRITE (0,*) 
     >     'STEAL: *** WARNING: FERATES NOT YET ON MODEL FILE  ***'
      ENDIF

      IF (.NOT. BFERATE) GOTO 20

      IF (BFIRSTITER) THEN
         WRITE (NAME, '(A5, I3)') 'FERLU', L
         CALL READMS (3, FERATLU, LASTFE, NAME, IERR)
         WRITE (NAME, '(A5, I3)') 'FERUL', L
         CALL READMS (3, FERATUL, LASTFE, NAME, IERR)
      ENDIF

      DO INDFE = 1, LASTFE
         IND = LASTIND - LASTFE + INDFE
         NUP = INDNUP(IND)
         LOW = INDLOW(IND)

      if (.false.) then      
C***  Note: RRATE contains Net Rates for iron lines. 
C***        FERAT (from COLI) contains "normal" (not NET) rates.  
         IF (BFIRSTITER) THEN
ccc kann weg wegen Netto-Rate!            FERATLU0(INDFE) = RRATE(LOW,NUP)
            FERATUL0(INDFE) = RRATE(NUP,LOW)
         ENDIF         
ccc kann weg wegen Netto-rate        RRATE(LOW,NUP) = 0.
         RRATE(NUP,LOW) = RRATE(NUP,LOW) +
     +          (FERATUL(INDFE) - EN(LOW)/EN(NUP)*FERATLU(INDFE)) -
     -          FERATUL0(INDFE)
      endif

        IF (BFIRSTITER) THEN
            FERATLU0(INDFE) = RRATE(LOW,NUP)
            FERATUL0(INDFE) = RRATE(NUP,LOW)
            RRATE(LOW,NUP) = FERATLU(INDFE)
            RRATE(NUP,LOW) = FERATUL(INDFE)

C***        Correction factor accounting for temperature changes:
C***        (motivated by Goetz branch)
            IF (bFeTCORR) THEN
               WAV0    = (ELEVEL(NUP) - ELEVEL(LOW))
C***           CORRS reflects (inverse) departure of S from B
               CORRT = MIN(1.,CORRS(L))
C***           CORRT <= 1: damp only, do not amplify for S < B
               DEXFAC(IND) = EXP(CORRT*C1*WAV0*(1./TL - 1./TLOLD))
            ELSE
               DEXFAC(IND) = 1.
            ENDIF
            
         ELSE
            RRATE(LOW,NUP) = FERATLU(INDFE) +
     >                       RRATE(LOW,NUP) - FERATLU0(INDFE)
            RRATE(NUP,LOW) = FERATUL(INDFE) +
     >                       RRATE(NUP,LOW) - FERATUL0(INDFE)
         ENDIF


C***     Account for temperature changes between the
C***     COLI calculation of the rates and now
C***     (This affects only the UL rates and only in the 
C***      case that only one set of cross section is stored in FEDAT)
         RRATE(NUP,LOW) = DEXFAC(IND) * RRATE(NUP,LOW)
C***     Note: Low-up rate does not contain the Boltzmann factor
         
      ENDDO

   20 CONTINUE
C***  END OF BRANCH FOR FE-RATES FROM COLI *****************

C***  ADD RADIATIVE AND COLLISIONAL TERMS INTO RATE COEFFICIENT MATRIX RATCO
      DO 1 J=1,N
        DO 1 I=1,N
           RATCO(I,J)=-RRATE(I,J)-CRATE(I,J)
    1 CONTINUE
 
C***  ADD ADDITIONAL D-R TERMS INTO RATE COEFFICIENT MATRIX RATCO
      DO 13 KDR=1,LASTKDR
      LOW=KODRLOW(KDR)
      NUP=IONGRND(LOW)
      RATCO(LOW,NUP)=RATCO(LOW,NUP)-RAUTO(LOW)
      RATCO(NUP,LOW)=RATCO(NUP,LOW)-RDIEL(LOW)
   13 CONTINUE

C***  DIAGONAL ELEMENTS: -SUM OF THE ROW (I.E. OVER COLUMN INDEX)
      DO 2 I=1,N
      SUM=.0
      DO 3 J=1,N
    3 SUM=SUM+RATCO(I,J)
    2 RATCO(I,I)=-SUM

C***  V1 = RIGHT-HAND SIDE VECTOR (INHOMOGENEITY)
      DO 9 J=1,NRANK
    9 V1(J)=.0
 
C***  COLUMN IMAXPOP (I.E. INDEX OF MAX. POPNUMBER PER ELEMENT):
C***         NUMBER CONSERVATION FOR EACH ELEMENT (NA)
C***  REMARK: TOTAL NUMBER CONSERVATION IS IMPLICITLY ENSURED

      DO 23 NA=1,NATOM
         NFIRNA=NFIRST(NA)
         NLANA=NLAST(NA)

C***     Check if all Rate Coefficients in one column are non-zero
C***     (otherwise: the matrix is singular!)
C***     and store logical flag ZERO_RATES for later use 
         IF (ITNEL .EQ. 1) THEN  
            CALL FLAG_ZERORATES(NFIRNA, NLANA, RATCO, NRANK,
     >                       IMAXPOP(NA), EN, POPMIN, ZERO_RATES(1,L))
C***        new, wrh 25-Feb-2015:
C***        in case of iron, levels MUST be popmin if they are also flagged
C***        at the next-inner depth point 
cc            IF (NA .EQ. KODAT(26) .AND. L .LT. ND) THEN
            IF (L .LT. ND) THEN
               DO J=NFIRNA, NLANA
                  IF (ZERO_RATES(J,L+1)) THEN
                     DO LL=L, 1, -1
                        ZERO_RATES(J,LL)=.TRUE.
                     ENDDO
                  ENDIF
               ENDDO
            ENDIF
         ENDIF

         DO 22 I=NFIRNA,NLANA
   22    RATCO(I,IMAXPOP(NA))=1.

         V1(IMAXPOP(NA))=ABXYZ(NA)


C***     If ZERO_RATES: Replace diagonal element by 1.0
C***     and the rest of this column by 0.0
         DO J = NFIRNA, NLANA
            IF (.NOT. ZERO_RATES(J,L)) CYCLE
            DO I = NFIRNA, NLANA
               RATCO(I,J) = .0
               RATCO(J,I) = .0
            ENDDO
            RATCO(J,J) = 1.
            V1(J) = POPMIN
            EN(J) = POPMIN
         ENDDO
   23 CONTINUE
C***  End-of-loop over all elements

C***  COLUMN N+1 : CHARGE CONSERVATION
      DO 4 I=1,N
    4 RATCO(I,NPLUS1)=NCHARG(I)
      RATCO(NPLUS1,NPLUS1)=-1.
C***  ROW N+1 : ZERO
      DO 5 J=1,N
    5 RATCO(NPLUS1,J)=.0
 
ccc   test output
cc      if (itnel .eq. 1 .and. l .eq. 10)
cc     >     call primat (RATCO, NPLUS1, NRANK, 'RATCO')
cc      if (itnel .eq. 1 .and. l .eq. 10)
cc     >     call primat (RRATE, NPLUS1, NDIM, 'RRATE')



C***  CALCULATION OF DERIVATIVES UNNECCESSARY IF BROYDEN=.TRUE.
      IF (BROYDEN) RETURN

C***  DERIVATIVE MATRIX DM
C***  FIRST TERMS : THE ORIGINAL MATRIX RATCO
      DO 6 I=1,NRANK
      DO 6 J=1,NRANK
    6 DM(I,J)=RATCO(I,J)
 
      DO 10 I=1,NRANK
C***  CONSTRUCT DERIVATIVE VECTORS DOPA, DETA WITH RESPECT TO EN(I)
      CALL       DCOOP (I,DOPA,DETA,XLAMBDA,NF,TL,RNEL,ENTOTL,EN,RSTAR,
     $                  WCHARM,ND,L,NFEDGE,EXPFAC,NDIM,N,NCHARG,WEIGHT,
     $                  ELEVEL,EION,EINST,SIGMAKI, KONTLOW,
     $                  KONTNUP,LASTKON,SIGMAFF,MAXION,NOM,KODAT,
     $                  SIGMATHK,SEXPOK,EDGEK,MAXATOM,bBLOCKINVERSION)
C***  CONSTRUCT DERIVATIVE VECTORS DOPAL, DETAL (LINES) WITH RESPECT TO EN(I)
      CALL       DLIOP (I,ENTOTL,DOPAL,DETAL,XRED,XBLUE,VDOP,RSTAR,N,
     $            NDIM,EINST,WEIGHT,ELEVEL,LASTIND,INDLOW,INDNUP,EN)
 
C***  COMPUTE DERIVATIVE MATRIX DM 
      CALL DERIV (DM,NRANK,I,NPLUS1,EN,CRATE,RRATE,EXPFAC,NFEDGE,
     $      WCHARM,ND,L,ENLTE,PHI,PWEIGHT,NFL,DELTAX,XMAX,
     $      XRED,XBLUE,DETAL,DOPAL,SLNEW,OPAL,XJLAPP,XJCAPP,
     $      FWEIGHT,DOPA,DETA,OPAC,SCNEW,XLAMBDA,NF,SCNEIND,OPACIND,
     $      NDIM,N,EINST,SIGMAKI,RUDLINE,
     $      LASTIND,INDLOW,INDNUP,KONTLOW,KONTNUP,LASTKON,
     $      IBLENDS,MAXLAP,XLAMZERO,BETA,PHIL,NBLENDS,VDOP,ATEST,BTEST,
     $      RDIEL,IONGRND,KODRLOW,LASTKDR,WFELOW,WFENUP,BDIAG,WEIGHT, 
     >      ELEVEL, NRB_CONT, TL, ENE, NOM, bBLOCKINVERSION)

   10 CONTINUE
 
C***  COLUMNS IMAXPOP (I.E. COLUMNS CONTAINING THE EQUATIONS OF NUMBER
C***  CONSERVATION FOR ELEMENT NA)  ARE NOT CHANGED
      DO 99 NA=1,NATOM
      DO 99 I=1,NRANK
      DM(I,IMAXPOP(NA))=RATCO(I,IMAXPOP(NA))
   99 CONTINUE

C***  COLUMNS WITH (almost) ZERO RATES ARE REPLACED BY THE EQUATION 
C***      n(j) * 1 = .0
 
      DO J=1, N
         IF (.NOT. ZERO_RATES(J,L)) CYCLE
         DO I=1,NRANK
            DM(I,J)=RATCO(I,J)
         ENDDO
cc         if (l .eq. 34)
cc     >   write (0,'(A,I3,A,A,A,G12.3)') 
cc     >    'ZERO_RATES at IT=', itnel, '  LEVEL=', LEVEL(J), 
cc     >     ' POPNUM=', EN(J)
      ENDDO

ccc   test output
cc      if (itnel .eq. 1 .and. l .eq. 17) then
cc          call primat (RATCO, NPLUS1, NRANK, 'RATCO')
cc          call primat (DM, NPLUS1, NRANK, 'DM')
cc      endif

      RETURN
      END
