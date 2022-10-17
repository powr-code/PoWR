      SUBROUTINE DRTRANS (XLAM,LINE,LOW,NUP,INDLAP,XLAMLAP,
     >                    DELXLAP,NBLINE,MAXLAP,INDLOW,INDNUP,LASTIND,
     >                    MAXIND,LEVEL,WEIGHT,ELEVEL,N,EINST,NDIM,
     >                    POPNUM,T,ND,ALN,VDOP,EION,ENTOT,RNE,
     >                    MAXSUBL,NSUBLOW,BROAD,LINPRO,AVOIGT, 
     >                    DENSCON,NMOD, MAXMOD, NDDIM, 
     >                    MAINQN, NCHARG, NOM, IND_ORIGLEV)
C***********************************************************************
C***  SUBROUTINE FOR DRTRANSIT HANDLING OF MAIN PROGRAM "FORMAL"
C***  CALLED FROM SUBR. PREFORM IN CASE OF DECODED OPTION "DRTRANSIT"
C***  INPUT OPTIONS: "/AUTONIVEAU" - DATA FOR UPPER (AUTOIONIZATION) LEVEL
C***                 "/LOWERLEVEL" - DATA FOR LOWER SUBLEVEL
C***                 "/ADDLINE"    - DATA FOR HANDLING OF ADDITIONAL LINE
C***                 "-DRTRANS"    - END OF DRTRANSIT INPUT 
C***  ACTION: 1. READ INPUT DATA (LEVELS, LINES) FOR DRTRANSIT HANDLING
C***          2. CALCULATE POPNUMBER FROM SAHA-BOLTZMANN FORMULA
C***             (I.E. LTE POPNUMBER)
C***          3. SPLIT LOWER LEVEL INTO LTE-LEVELS (OPTIONAL)
C***********************************************************************

      DIMENSION ND(MAXMOD)
      DIMENSION INDLAP(MAXLAP),XLAMLAP(MAXLAP),DELXLAP(MAXLAP)
      DIMENSION AVOIGT(MAXLAP,NDDIM,MAXMOD)
      DIMENSION WEIGHT(N),ELEVEL(N),EION(N), MAINQN(N), NOM(N)
      DIMENSION EINST(NDIM,NDIM), IND_ORIGLEV(NDIM), NCHARG(NDIM)
      DIMENSION INDLOW(LASTIND),INDNUP(LASTIND)
      DIMENSION POPNUM(NDDIM,NDIM,NMOD)
      REAL, DIMENSION (NDDIM,NMOD) :: T, ENTOT, RNE, DENSCON
      DIMENSION NSUBLOW(MAXSUBL)
      CHARACTER KARTE*80
      CHARACTER*10 LEVEL(N), LEV, LEVUP,LEVLOW
      CHARACTER*8 LINPRO(MAXLAP), LINPROBL
      LOGICAL BROAD

C***  CI : FACTOR IN SAHA EQUATION (MIHALAS, P. 113)
      DATA CI / 2.07E-16 /
C***  C1 = H * C / K    ( CM*KELVIN )
      DATA C1 / 1.4388 /

C***  DEFAULTS FOR DRTRANSIT HANDLING:
      NBLSAVE=NBLINE
C***  NLOW: NUMBER OF /LOWERLEVELS from splitting LEVELS
      NLOW=0
C***  NSUBLOW: POINTER TO ADDITIONAL LEVELS IN THE ORIGINAL ARRAYS
      DO 10 I=1,MAXSUBL
   10 NSUBLOW(I)=0

C***  1. READ INPUT FOR DRTRANSIT HANDLING OF SPECIFIED LINE  ----------
    1 READ (2,'(A)') KARTE
      IF (KARTE(1:1) .EQ. '*') GOTO 1
      IF (KARTE      .EQ. '' ) GOTO 1

      IF (KARTE(:8) .EQ. '-DRTRANS') GOTO 20
C                         ========

C************************************************************************
      IF (KARTE(:11) .EQ. '/AUTONIVEAU') THEN
C                          ===========
         READ (KARTE,3) LEV, NW, ELEV
    3    FORMAT (12X,A10,1X,I4,1X,F20.0)

C***     FIND INDEX of AUTOLEVEL (IF ALREADY EXISTING):
         JNUP=0
C***     a) ALREADY EXISTING ?
         DO J=1,N
            IF (LEVEL(J) .EQ. LEV) THEN
               IF ((ELEVEL(J) .NE. ELEV) .OR.
     >             (WEIGHT(J) .NE. FLOAT(NW))) THEN
  202             FORMAT (' >>>>> SUBR. DRTRANS: WARNING: ', //, 
     >                    'INCONSISTENT ASSIGNMENT OF LEVEL ', A)
  203             FORMAT ('ELEVEL=',2(F9.1,1X),F7.4)
                  WRITE (0,202) LEV
                  WRITE (0,203) ELEVEL(J), ELEV, 1.-ELEVEL(J)/ELEV
                  WRITE (*,202) LEV
                  WRITE (*,203) ELEVEL(J), ELEV, 1.-ELEVEL(J)/ELEV
               ENDIF
               JNUP=J
               EXIT
            ENDIF
         ENDDO

C***     b) LEVEL NOT YET DEFINED -> append to the list 
         IF (JNUP .EQ. 0) THEN
            N=N+1
            IF (N .GT. NDIM) THEN
               PRINT *, ' >>>>> SUBR. DRTRANS: ERROR STOP (N .GT. NDIM)'
               CALL REMARK ('DRTRANS: N GREATER THAN NDIM')
               STOP 'NDIM exceeded in DRTRANS'
            ENDIF

            LEVEL(N)=LEV
            WEIGHT(N)=FLOAT(NW)
            ELEVEL(N)=ELEV
C*          The AUTO level belongs to the lower ion - copy from there:
            MAINQN(N)= MAINQN(LOW)
            NCHARG(N)= NCHARG(LOW)
            NOM(N)   = NOM(LOW)
            EION(N)  = EION(LOW)
            IND_ORIGLEV(N) = LOW

C***        POPNUMBER OF THE NEW AUTOIONIZATION LEVELS: RELATIVE TO THE 
C***        GIVEN UPPER STATE by SAHA-BOLTZMANN (LTE ratio)
C***        Note: AUTONIVEAU population is *not* subtracted from 
C***              population of the parent level! --> inconsistent! ллллллллллллллл
            WAUTO= WEIGHT(N)
            WNUP = WEIGHT(NUP)
            EAUTO= ELEVEL(N)
            EILOW= EION(LOW)
            ENUP = ELEVEL(NUP)

            DO IMOD=1, NMOD
               DO L=1, ND(IMOD)
                  TL = T(L,IMOD)
                  SQRTL = SQRT(TL)
                  ENTOTL = ENTOT(L,IMOD) * DENSCON(L,IMOD)
                  RNEL = RNE(L,IMOD)
                  IPOP = L+ (NUP-1)*ND(IMOD)
                  POPNUP = POPNUM(IPOP,1,IMOD)
                  BOLTZ = EXP(-C1*(EAUTO-EILOW-ENUP)/TL)
                  SAHAPHI = WAUTO/WNUP*CI/TL/SQRTL * BOLTZ

ccc wrong versions of BOLTZ:  wrh, 26-Oct-2021
cc???     >           BOLTZ = EXP(-C1*(EAUTO-EILOW)/TL)
cc???     >           BOLTZ = EXP(-C1*(ELOW-EILOW-ENUP)/TL)

                  POPAUTO = POPNUP * ENTOTL * RNEL * SAHAPHI
                  IPOP = L+ (N-1)*ND(IMOD)
                  POPNUM(IPOP,1,IMOD) = POPAUTO

               ENDDO
            ENDDO 
         ENDIF

C************************************************************************
      ELSE IF (KARTE(:11) .EQ. '/LOWERLEVEL') THEN
C                               ===========
         N=N+1
         IF (N .GT. NDIM) THEN
            PRINT *, ' >>>>> SUBR. DRTRANS: ERROR STOP (N .GT. NDIM)'
            CALL REMARK ('DRTRANS: N GREATER THAN NDIM')
            STOP 'NDIM1'
         ENDIF
         NLOW=NLOW+1
         IF (NLOW .GT. MAXSUBL) THEN
            PRINT *,
     >           ' >>>>> SUBR. DRTRANS: ERROR STOP (NLOW .GT. MAXSUBL)'
            CALL REMARK ('DRTRANS: NLOW GREATER THAN MAXSUBL')
            STOP 'NLOW'
         ENDIF
         NSUBLOW(NLOW)=N
         READ (KARTE,3) LEVEL(N),NW,ELEVEL(N)
         WEIGHT(N)=FLOAT(NW)

C*       The sublevel stems from the lower level - copy from there:
         MAINQN(N)= MAINQN(LOW)
         NCHARG(N)= NCHARG(LOW)
         NOM(N)   = NOM(LOW)
         EION(N)  = EION(LOW)
         IND_ORIGLEV(N) = LOW

C************************************************************************
      ELSE IF (KARTE(:8) .EQ. '/ADDLINE') THEN
C                              ========

C***     First: if the LOWERLEVEL is to be split:
C***     BOLTZMANN (LTE)
         IF (NLOW .GT. 0) THEN
            DO IMOD=1, NMOD
               DO 30 L=1, ND(IMOD)
                  IPOP = L+ (NSUBLOW(1)-1)*ND(IMOD)
                  POPNUM(IPOP,1,IMOD) = 1.
                  DO J=2,NLOW
                    IPOP   = L+ (NSUBLOW(J  )-1)*ND(IMOD)
                    IPOPM1 = L+ (NSUBLOW(J-1)-1)*ND(IMOD)
                    POPNUM(IPOP,1,IMOD) = 
     >               EXP(C1*(ELEVEL(NSUBLOW(J-1))
     -               -ELEVEL(NSUBLOW(J)))/T(L,IMOD))*WEIGHT(NSUBLOW(J))
     /               /WEIGHT(NSUBLOW(J-1))*POPNUM(IPOPM1,1,IMOD)
                  ENDDO

C***              NORMALIZATION
                  SUM=0.
                  DO J=1,NLOW     
                     IPOP = L+ (NSUBLOW(J)-1)*ND(IMOD)
                     SUM = SUM + POPNUM(IPOP,1,IMOD)
                  ENDDO

                  IPOP   = L+ (LOW-1)*ND(IMOD)
                  SUM = SUM / POPNUM(IPOP,1,IMOD)

                  DO J=1,NLOW
                     IPOP = L+ (NSUBLOW(J)-1)*ND(IMOD)
                     POPNUM(IPOP,1,IMOD) 
     >                 = POPNUM(IPOP,1,IMOD) / SUM
                  ENDDO

   30          CONTINUE
            ENDDO
         ENDIF

C***     Now executing the ADDLINE actions
         LASTIND=LASTIND+1

         IF (LASTIND .GT. MAXIND) THEN
            PRINT *,
     >        ' >>>>> SUBR. DRTRANS: ERROR STOP (LASTIND .GT. MAXIND)'
            CALL REMARK ('DRTRANS: LASTIND GREATER THAN MAXIND')
            STOP 'MAXIND'
         ENDIF
         IF (NBLINE+1 .GT. MAXLAP) THEN
            PRINT *,
     >          ' >>>>> SUBR. DRTRANS: ERROR STOP (NBLINE .GT. MAXLAP)'
            CALL REMARK ('DRTRANS: NBLINE GREATER THAN MAXLAP')
            STOP 'MAXLAP'
         ENDIF

         READ (KARTE,23) LEVUP,LEVLOW,AUPLOW
   23    FORMAT (9X,A10,2X,A10,G10.0)

C***     FIND UPPER INDEX:
         DO 24 J=1,N
         JNUP=J
         IF (LEVEL(J) .EQ. LEVUP) GOTO 25
   24    CONTINUE

   90    FORMAT ('DRTRANS: UPPER LINE LEVEL NOT FOUND: ', A10)
         WRITE (0, 90)  LEVUP
         STOP 'UPPER'

C***     FIND LOWER INDEX:
   25    DO 26 J=1,N
         JLOW=J
         IF (LEVEL(J) .EQ. LEVLOW) GOTO 27
   26    CONTINUE

   91    FORMAT ('DRTRANS: LOWER LINE LEVEL NOT FOUND: ', A10)
         WRITE (0, 91) LEVLOW
         STOP 'LOWER'

   27    INDNUP(LASTIND)=JNUP
         INDLOW(LASTIND)=JLOW

C***     Check for correct sequence UP - LOW
         IF (ELEVEL(JNUP) .LE. ELEVEL(JLOW)) THEN
            WRITE (0,'(A)') '*** FATAL ERROR in FORMAL_CARDS:'
            WRITE (0,'(A)') '*** E(UP) must be higher than E(LOW) !'
            WRITE (0,'(A)') '*** error occured in the following line:'
            WRITE (0,'(A)') KARTE(:IDX(KARTE))
            STOP '*** ERROR in subroutine DRTRANS'
         ENDIF

C***     default wavelength 
         XLAMSUB=1.E8/(ELEVEL(JNUP)-ELEVEL(JLOW))

C***     read wavelength if given, covert from AIR to VAC
C***     read optional VOIGT parameter only if BROAD=.TRUE.
          CALL READ_LINECARD_PARAMETERS (KARTE(43:), XLAMSUB,
     >            BROAD, LINPROBL, AVOIGTBL)
          IF (BROAD .AND. (LINPROBL .EQ. '')) LINPROBL = 'DRTRANS '

C***     If first line in blend block: it defines reference lambda
         IF (NBLINE .EQ. 0) XLAM = XLAMSUB

C***     In case AUPLOW < 0 it means f-Value --> convert
         IF (AUPLOW .LT. 0.) THEN
            EINST(JNUP,JLOW)=-6.669E15/XLAMSUB/XLAMSUB*WEIGHT(JLOW)
     /                                      /WEIGHT(JNUP)*AUPLOW
         ELSE
            EINST(JNUP,JLOW)=AUPLOW
         ENDIF

C***     Insert current ADDLINE in the list sorted by increasing DELTAX
         CALL INSERT_LINE (LINE, LASTIND, NBLINE, INDLAP,
     >                   XLAMLAP,DELXLAP,
     $                   XLAMSUB, LINPROBL, AVOIGTBL, XLAM, MAXLAP, ALN,
     >                   ND, LINPRO, AVOIGT, NMOD, NDDIM, MAXMOD )

C************************************************************************
      ELSE
         WRITE (0,*) 'UNRECOGNIZED INPUT CARD IN SUBR. DRTRANS!'
         WRITE (0,*) 'KARTE=', KARTE( :IDX(KARTE))

      ENDIF
C************************************************************************

      GOTO 1


   20 CONTINUE
C***  END OF INPUT FOR DRTRANSIT HANDLING  -----------------------------


C***  NO SUBLINES DECODED: RETURN WITHOUT ANY DR-TRANSITION
      IF (NBLINE .EQ. NBLSAVE) THEN
        WRITE (0,*) '*** WARNING: DRTRANSIT BLOCK WITHOUT LINES!'
      ENDIF

   99 CONTINUE


      RETURN
      END
