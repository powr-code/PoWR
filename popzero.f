      SUBROUTINE POPZERO (T,RNE,POPNUM,DEPART,ENTOT,ITNE,NDIM,N,ENLTE,
     $                   WEIGHT,NCHARG,EION,ELEVEL,EN,EINST,LEVEL,
     $                   XLAMBDA,FWEIGHT,XJC,NF,XJL,IFRRA,ITORA,ALPHA,
     $                   SEXPO,
     $                   ADDCON1, ADDCON2, ADDCON3, 
     $                   IGAUNT,MODHEAD,JOBNUM,
     $                   ND,LSRAT,CRATE,RRATE,RATCO,
     $                   SIGMAKI,ALTESUM,COCO,KEYCBB,NOM,NATOM,ABXYZ,
     $                   KODAT,NFIRST,NLAST,NATOUT,
     $                   NAUTO,MAXAUTO,LOWAUTO,WAUTO,EAUTO,AAUTO,
     $                   RDIEL,RAUTO,IONAUTO,IONGRND,
     $                   INDNUP,INDLOW,LASTIND,
     $                   KONTNUP,KONTLOW,LASTKON,KODRNUP,KODRLOW,
     $                   LASTKDR,KEYCBF,MAXATOM,SIGMATHK,SEXPOK,EDGEK,
     $                   LAINDHI,KRUDAUT, ZERO_RATES, POPMIN,
     >                   PRILEVRA, SGAIN, SLOSS)
 
C*******************************************************************************
C***  CALCULATION OF THE NLTE POP.NUMBERS BY SOLVING THE RATE EQUATIONS.
C***   --- NO APPROXIMATE LAMBDA OPERATORS ARE INCLUDED HERE -----
C***  THIS ROUTINE IS SIMILAR TO SUBR. NEWPOP
C***  IT IS USED BY MAIN PROGRAM STEAL IF ALL GAMMA'S ARE ZERO
C***  IT HAS THE ADVANTAGE THAT THE RATE COEFFICIENTS ARE PROVIDED
C***  WHICH CAN BE PRINTED WITH SUBR. PRIRAT
C***  CALLED FROM: STEAL
C*******************************************************************************
 
      DIMENSION T(ND),ENTOT(ND),RNE(ND),POPNUM(ND,N),ITNE(ND)
      DIMENSION DEPART(ND,N)
      DIMENSION NCHARG(NDIM),EN(NDIM),ENLTE(NDIM)
      DIMENSION ABXYZ(NATOM),KODAT(NATOM),NFIRST(NATOM),NLAST(NATOM)
      LOGICAL KONVER
      CHARACTER JOB*7, PRILEVRA*10
      LOGICAL ZERO_RATES(N,ND)
      CHARACTER STRING3*3, STRING15*15, LEVEL(N)*(*)
      DIMENSION NZERORATES(NDIM), LMINZERORATES(NDIM), LMAXZERORATES(NDIM)

      COMMON /GIIIERR/  NTUP,NTLOW,NFUP,NFLOW,NHELP
      COMMON / COMNEGI / NEGINTL,INEGMIN,INEGMAX,LNEGMIN,LNEGMAX
      COMMON / COMITWA / ITWARN, ITMAX

      NTUP=0
      NTLOW=0
      NFUP=0
      NFLOW=0
      NEGINTL=0
      ITWARN = 0
 
C***  EPSNE = ACCURACY LIMIT FOR THE ELECTRON DENSITY ITERATION
      EPSNE=0.001
C***  MAX. NUMBER OF ITERATIONS
      ITMAX=10

C***  Initialize counters for ZERO_RATES
      DO J=1, N
        NZERORATES(J) = 0
        LMINZERORATES(J) = ND+1
        LMAXZERORATES(J) = 0
      ENDDO


C***  GENERATE ONCE FOR ALL PHOTOCROSSSECTIONS AT ALL FREQUENCIES
C***  SIGMAKI(K,LOW) IN CM**2
      CALL       BFCROSS (SIGMAKI,NF,N,ELEVEL,EION,EINST,NDIM,
     $                    XLAMBDA,ALPHA,SEXPO,
     $                    ADDCON1, ADDCON2, ADDCON3, 
     $                    IGAUNT,
     $                    KONTNUP,KONTLOW,LASTKON)
 
C***  LOOP OVER ALL DEPTH POINTS  --------------------------------------
      DO 1 L=1,ND
      TL=T(L)
      ENTOTL=ENTOT(L)
      RNEL=RNE(L)
 
C***  ITERATION FOR THE ELECTRON DENSITY
      ITNEL=0
   13 ITNEL=ITNEL+1
      ENE=RNEL*ENTOTL
      CALL LTEPOP (N,ENLTE,TL,ENE,WEIGHT,NCHARG,EION,ELEVEL,NOM,ABXYZ,
     $             NFIRST,NLAST,NATOM)

C***  The ZERORATE flagging needs already some estimate of the LTE 
C***  popnumbers EN. Therefore, we initialize with LTE:
      IF (ITNEL .EQ. 1) THEN
         DO I=1, N
            EN(I) = ENLTE(I)
         ENDDO
      ENDIF

      CALL NLTEPOP (NDIM,N,ENLTE,TL,ENE,WEIGHT,NCHARG,EION,ELEVEL,EN,
     $             EINST,XLAMBDA,FWEIGHT,XJC,NF,L,LEVEL,XJL,ND,
     $             CRATE,RRATE,RATCO,SIGMAKI,ALTESUM,COCO,KEYCBB,NOM,
     $             NATOM,ABXYZ,KODAT,NFIRST,NLAST,
     $             NAUTO,MAXAUTO,LOWAUTO,WAUTO,EAUTO,AAUTO,
     $             RDIEL,RAUTO,IONAUTO,IONGRND,
     $             INDNUP,INDLOW,LASTIND,
     $             KONTNUP,KONTLOW,LASTKON,KODRNUP,KODRLOW,LASTKDR,
     $             KEYCBF,MAXATOM,SIGMATHK,SEXPOK,EDGEK,
     >             LAINDHI,KRUDAUT, ZERO_RATES, POPMIN)
      RNEOLD=RNEL
      RNEL=0.0
      DO 12 J=1,N
   12 RNEL=RNEL+NCHARG(J)*EN(J)
 
      KONVER= ABS(RNEOLD-RNEL) .LT. EPSNE .AND. ITNEL .GT. 1
      IF (.NOT. KONVER .AND. ITNEL .LT. ITMAX) GOTO 13
      IF (.NOT. KONVER) THEN
         ITNEL = -ITMAX
         ITWARN = ITWARN + 1
         ENDIF
 
      ITNE(L)=ITNEL
      RNE(L)=RNEL
      DO 14 J=1,N
      DEPART(L,J)=EN(J)/ENLTE(J)
C***      POP1(L,J)=EN(J)
   14 POPNUM(L,J)=EN(J)

C***  Update counters for ZERO_RATES (based on the final iteration)
      DO J=1, N
         IF (ZERO_RATES(J,L)) THEN
            NZERORATES(J) = NZERORATES(J) + 1
            LMINZERORATES(J) = MIN0 (LMINZERORATES(J),L)
            LMAXZERORATES(J) = MAX0 (LMAXZERORATES(J),L)
         ENDIF
      ENDDO
 
C***  ===  PRINTOUT OF RATE COEFFICIENTS ETC.  ===
      IF (LSRAT.NE.-1) THEN
         IF ((L.GE.IFRRA.AND.L.LE.ITORA).OR.ITORA.EQ.0) THEN
            LM1=L-1
            IF (IFRRA.GT.0) LM1=L-IFRRA
            IF (PRILEVRA .NE. ' ') THEN
              IF  (((LM1)/LSRAT)*LSRAT.EQ.(LM1).OR.L.EQ.ND)
     >        CALL PRI1RAT (N, LEVEL, NDIM, L, CRATE, RRATE, RATCO, EN,
     >        MODHEAD, JOBNUM, NFIRST, NLAST, NATOM,
     >        NAUTO, RDIEL, RAUTO, IONGRND, KODRLOW, LASTKDR,
     >        PRILEVRA, SGAIN, SLOSS)
            ELSE
              NETTO=0
              IF  (((LM1)/LSRAT)*LSRAT.EQ.(LM1).OR.L.EQ.ND)
     $        CALL PRIRAT (DUMMY  ,N,LEVEL,NDIM,L,CRATE,RRATE,RATCO,EN,
     $           IFRRA,MODHEAD,JOBNUM,NETTO,NFIRST,NLAST,NATOM,NATOUT,
     $           NAUTO,RDIEL,RAUTO,IONGRND,KODRLOW,LASTKDR)
            ENDIF
         ENDIF
      ENDIF
 
C***********************************************************************
C***  test output: "image" of specified matrix (NOTE:  ndim <= 125 )
      ldepth = 0
      if (l .eq. ldepth) then
      call primat (rrate,n,ndim,'RRATE')
      call primat (crate,n,ndim,'CRATE')
      endif
C***********************************************************************

    1 CONTINUE
C***  ENDLOOP  ---------------------------------------------------------

C***  Printout of ZERO_RATES 
       WRITE (0,'(/,A)') 'Levels for which the equations have been removed'
      WRITE (0,'(  A)') '------------------------------------------------'
      DO J=1,N
         IF (NZERORATES(J) .EQ. 0) CYCLE
         IF (NZERORATES(J) .EQ. ND) THEN
            STRING15 = ' = all depths!'
         ELSE
            STRING15 = ''
         ENDIF
         IF (NZERORATES(J) .EQ.
     >      LMAXZERORATES(J) + 1 - LMINZERORATES(J)) THEN
            STRING3 = 'all'
         ELSE
            WRITE (STRING3, '(I3)') NZERORATES(J)
         ENDIF

         WRITE (0, '(A,I4,A,I3,A,I3,A )')
     >      'Level', J,': ' // LEVEL(J) // ' at ' // STRING3 //
     >      ' depth points between L=', LMINZERORATES(J),
     >      ' and', LMAXZERORATES(J), STRING15
      ENDDO
 


 
      RETURN
      END
