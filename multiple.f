      SUBROUTINE MULTIPLE (XLAM, LINE, LOW, NUP, INDLAP, XLAMLAP,
     >                    DELXLAP,NBLINE,MAXLAP,INDLOW,INDNUP,LASTIND,
     >                    MAXIND,LEVEL,WEIGHT,ELEVEL,N,EINST,NDIM,
     >                    POPNUM,T,ND,ALN,VDOP, 
     >                    MAXSUBL,NSUBLOW,NSUBNUP,BROAD,LINPRO,AVOIGT,
     >                    NMOD, MAXMOD, NDDIM, 
     >                    MAINQN, NCHARG, EION, NOM, IND_ORIGLEV)
C***********************************************************************
C***  SUBROUTINE FOR MULTIPLET HANDLING OF MAIN PROGRAM "FORMAL"
C***  CALLED FROM SUBR. PREFORM IN CASE OF DECODED OPTION "MULTIPLET"
C***  INPUT OPTIONS: "/LOWER"   - DATA FOR LOWER SUBLEVEL
C***                 "/UPPER"   - DATA FOR UPPER SUBLEVEL
C***                 "/SUBLINE" - DATA FOR SUBLINE HANDLING
C***                 "-MULTI"   - END OF MULTIPLET INPUT 
C***  ACTION: 1. READ INPUT DATA (LEVELS, LINES) FOR MULTIPLET HANDLING
C***          2. ARRANGE SUBLINES IN A SEQUENCE OF INCREASING FREQUENCIES
C***          3. SPLIT ORIGINAL ENERGY LEVELS INTO SPECIFIED SUBLEVELS
C***             WITH POPNUMBERS CALCULATED FROM THE BOLTZMANN FORMULA
C***             (I.E. LTE POPNUMBERS)
C***********************************************************************

C***  DEFINE FORTRAN CHANNEL FOR LOGFILE MESSAGES
      PARAMETER ( LOGCHAN = 0 )
    
      DIMENSION ND(NMOD)
      DIMENSION INDLAP(MAXLAP),XLAMLAP(MAXLAP),DELXLAP(MAXLAP)
      DIMENSION AVOIGT(MAXLAP,NDDIM,NMOD)
      DIMENSION WEIGHT(N), ELEVEL(N), MAINQN(N), NOM(N)
      DIMENSION NCHARG(N), EION(N)
      DIMENSION EINST(NDIM,NDIM), IND_ORIGLEV(NDIM)
      DIMENSION INDLOW(LASTIND),INDNUP(LASTIND)
      DIMENSION POPNUM(NDDIM,NDIM,NMOD),T(NDDIM,NMOD)
      DIMENSION NSUBLOW(MAXSUBL),NSUBNUP(MAXSUBL)
      CHARACTER KARTE*80
      CHARACTER*10 LEVEL(N), LEV, LEVUP, LEVLOW
      CHARACTER*8 LINPRO(MAXLAP)
      LOGICAL BROAD, BAIR, BVAC

C***  C1 = H * C / K    ( CM*KELVIN )
      DATA C1 / 1.4388 /
      
C***  DEFAULTS FOR MULTIPLET HANDLING:
      NBLSAVE=NBLINE
C***  NLOW, NNUP: NUMBER OF LOWER, UPPER SUBLEVELS
      NLOW=0
      NNUP=0
C***  NSUBLOW, NSUBNUP: POINTER TO SUBLEVELS IN THE ORIGINAL ARRAYS
      DO 10 I=1,MAXSUBL
      NSUBLOW(I)=0
   10 NSUBNUP(I)=0

C***  1. READ INPUT FOR MULTIPLET HANDLING OF SPECIFIED LINE  ----------
    1 READ (2,2) KARTE
      IF (KARTE(1:1) .EQ. '*' .OR. KARTE .EQ. ' ') GOTO 1
    2 FORMAT (A)
      IF (KARTE(:6) .EQ. '-MULTI') GOTO 20
C                         ======
      IF (KARTE(:6) .EQ. '/LOWER') THEN
C                         ======
         READ (KARTE,3) LEV,NW,ELEV
    3    FORMAT (12X,A10,1X,I4,1X,F20.0)

C***     FIND LOWER INDEX:
         JLOW=0
         DO 200 J=1,N
            IF (LEVEL(J) .EQ. LEV) THEN
               IF ((ELEVEL(J) .NE. ELEV).OR.
     >             (WEIGHT(J) .NE. FLOAT(NW))) THEN
                  WRITE (0,*) ' >>>>> SUBR. MULTIPLE: WARNING ',
     >                    '(MULTIPLE ASSIGNMENT OF LEVEL "',LEV,'")'
                  WRITE (0,'(A,2(F9.1,1X),F7.4)')
     >                    'ELEVEL=', ELEVEL(J), ELEV, 1.-ELEVEL(J)/ELEV
                  WRITE (*,*) ' >>>>> SUBR. MULTIPLE: WARNING ',
     >                    '(MULTIPLE ASSIGNMENT OF LEVEL "',LEV,'")'
                  WRITE (*,'(A,2(F9.1,1X),F7.4)')
     >                    'ELEVEL=', ELEVEL(J), ELEV, 1.-ELEVEL(J)/ELEV
               ELSE
                 JLOW=J
               ENDIF
            ENDIF
 200     CONTINUE

C***     LEVEL NOT YET DEFINED
         IF (JLOW .EQ. 0) THEN
            N=N+1
            IF (N .GT. NDIM) THEN
               PRINT *,
     >            ' >>>>> SUBR. MULTIPLE: ERROR STOP (N .GT. NDIM)'
               CALL REMARK ('MULTIPLE: N GREATER THAN NDIM')
               STOP 'NDIM1'
            ENDIF
            LEVEL(N)=LEV
            WEIGHT(N)=FLOAT(NW)
            ELEVEL(N)=ELEV
            MAINQN(N)=MAINQN(LOW)
            NCHARG(N)=NCHARG(LOW)
            NOM(N) = NOM(LOW)
            EION(N)  =EION(LOW)
            IND_ORIGLEV(N) = LOW
            JLOW=N
         ENDIF

         NLOW=NLOW+1
         IF (NLOW .GT. MAXSUBL) THEN
            PRINT *,
     >           ' >>>>> SUBR. MULTIPLE: ERROR STOP (NLOW .GT. MAXSUBL)'
            CALL REMARK ('MULTIPLE: NLOW GREATER THAN MAXSUBL')
            STOP 'NLOW'
         ENDIF
         NSUBLOW(NLOW)=JLOW

      ELSE IF (KARTE(:6) .EQ. '/UPPER') THEN
C                              ======
         READ (KARTE,3) LEV, NW, ELEV 

C***     FIND UPPER INDEX:
         JNUP=0
         DO 201 J=1,N
            IF (LEVEL(J) .EQ. LEV) THEN
               IF ((ELEVEL(J) .NE. ELEV).OR.
     >             (WEIGHT(J) .NE. FLOAT(NW))) THEN
                  WRITE (0,*) ' >>>>> SUBR. MULTIPLE: WARNING ',
     >                    '(MULTIPLE ASSIGNMENT OF LEVEL "',LEV,'")'
                  WRITE (0,'(A,2(F9.1,1X),F7.4)')
     >                    'ELEVEL=', ELEVEL(J), ELEV, 1.-ELEVEL(J)/ELEV
                  WRITE (*,*) ' >>>>> SUBR. MULTIPLE: WARNING ',
     >                    '(MULTIPLE ASSIGNMENT OF LEVEL "',LEV,'")'
                  WRITE (*,'(A,2(F9.1,1X),F7.4)')
     >                    'ELEVEL=', ELEVEL(J), ELEV, 1.-ELEVEL(J)/ELEV
               ELSE
                 JNUP=J
               ENDIF
            ENDIF
 201     CONTINUE

C***     LEVEL NOT YET DEFINED
         IF (JNUP .EQ. 0) THEN
            N=N+1
            IF (N .GT. NDIM) THEN
               PRINT *,
     >            ' >>>>> SUBR. MULTIPLE: ERROR STOP (N .GT. NDIM)'
               CALL REMARK ('MULTIPLE: N GREATER THAN NDIM')
               STOP 'NDIM2'
            ENDIF
            LEVEL(N)=LEV
            WEIGHT(N)=FLOAT(NW)
            ELEVEL(N)=ELEV
            MAINQN(N)=MAINQN(NUP)
            NCHARG(N)=NCHARG(NUP)
            NOM(N) = NOM(NUP)
            EION(N)  =EION(NUP)
            IND_ORIGLEV(N) = NUP
            JNUP=N
         ENDIF

         NNUP=NNUP+1
         IF (NNUP .GT. MAXSUBL) THEN
            PRINT *,
     >           ' >>>>> SUBR. MULTIPLE: ERROR STOP (NNUP .GT. MAXSUBL)'
            CALL REMARK ('MULTIPLE: NNUP GREATER THAN MAXSUBL')
            STOP 'NNUP'
         ENDIF
         NSUBNUP(NNUP)=JNUP

      ELSE IF (KARTE(:8) .EQ. '/SUBLINE') THEN
C                              ========
         LASTIND=LASTIND+1
         IF (LASTIND .GT. MAXIND) THEN
            PRINT *,
     >        ' >>>>> SUBR. MULTIPLE: ERROR STOP (LASTIND .GT. MAXIND)'
            CALL REMARK ('MULTIPLE: LASTIND GREATER THAN MAXIND')
            STOP 'MAXIND'
         ENDIF
         IF (NBLINE+1 .GT. MAXLAP) THEN
            PRINT *,
     >          ' >>>>> SUBR. MULTIPLE: ERROR STOP (NBLINE .GT. MAXLAP)'
            CALL REMARK ('MULTIPLE: NBLINE GREATER THAN MAXLAP')
            STOP 'MAXLAP'
         ENDIF

         READ (KARTE,23) LEVUP,LEVLOW,AUPLOW
   23    FORMAT (9X,A10,2X,A10,G10.0)

C***     FIND UPPER INDEX:
         DO 24 J=1,N
         JNUP=J
         IF (LEVEL(J) .EQ. LEVUP) GOTO 25
   24    CONTINUE

   90    FORMAT ('*** ERROR: UPPER LINE LEVEL NOT FOUND: ', A10)
         WRITE (LOGCHAN, 90)  LEVUP
         STOP 'ERROR STOP IN SUBR. MULTIPLE'

C***     FIND LOWER INDEX:
   25    DO 26 J=1,N
         JLOW=J
         IF (LEVEL(J) .EQ. LEVLOW) GOTO 27
   26    CONTINUE

   91    FORMAT ('*** ERROR: LOWER LINE LEVEL NOT FOUND: ', A10)
         WRITE (LOGCHAN, 91) LEVLOW
         STOP 'ERROR STOP IN SUBR. MULTIPLE'

   27    INDNUP(LASTIND)=JNUP
         INDLOW(LASTIND)=JLOW
C***     default wavelength 
         XLAMSUB=1.E8/(ELEVEL(JNUP)-ELEVEL(JLOW))

C***     read wavelength if given, covert from AIR to VAC
C***     VOIGT parameter?
         CALL SARGC (KARTE(43:),NPAR)
         IF (NPAR .GT. 0) 
     >     CALL READ_LINECARD_PARAMETERS (KARTE(43:), XLAMSUB,
     >            BROAD, LINPROBL, AVOIGTBL)

         IF (NBLINE .EQ. 0) XLAM = XLAMSUB 

         IF (AUPLOW .LT. 0.) THEN
            EINST(JNUP,JLOW)=-6.669E15/XLAMSUB/XLAMSUB*WEIGHT(JLOW)
     /                                      /WEIGHT(JNUP)*AUPLOW
         ELSE
            EINST(JNUP,JLOW)=AUPLOW
         ENDIF


C***     Insert current SUBLINE in the list sorted by increasing DELTAX
         CALL INSERT_LINE (LINE, LASTIND, NBLINE, INDLAP, XLAMLAP,DELXLAP,
     $                   XLAMSUB, LINPROBL, AVOIGTBL, XLAM, MAXLAP, ALN,
     >                   ND, LINPRO, AVOIGT, NMOD, NDDIM, MAXMOD )

      ELSE
         WRITE (0,*) 'UNRECOGNIZED INPUT CARD IN SUBR. MULTIPLE!'
         WRITE (0,*) 'KARTE=', KARTE( :IDX(KARTE))
      ENDIF

      GOTO 1

   20 CONTINUE
C***  END OF INPUT FOR MULTIPLET HANDLING  -----------------------------

C***  NO SUBLINES DECODED
      IF (NBLINE .EQ. NBLSAVE) THEN
         WRITE (0,'(A)') 
     >    ' >>>>> SUBR. MULTIPLE: WARNING - NO SUBLINE IN MULTIPLET:'
         WRITE (0,'(1X,A,1X,A)') LEVEL(NUP) , LEVEL(LOW)
         WRITE (*,'(A)') 
     >    ' >>>>> SUBR. MULTIPLE: WARNING - NO SUBLINE IN MULTIPLET:'
         WRITE (*,'(1X,A,1X,A)') LEVEL(NUP) , LEVEL(LOW)
         GOTO 99
      ENDIF

C***  3. SPLITTING OF THE ENERGY LEVELS: BOLTZMANN (LTE POPNUMBERS)
      DO IMOD=1, NMOD
        CALL MULTSPLI(ND(IMOD), N, NSUBLOW, MAXSUBL, NLOW, POPNUM(1,1,IMOD), 
     >                ELEVEL, 
     >                T(1,IMOD), WEIGHT, NNUP, NSUBNUP, LOW, NUP, 
     >                LEVEL)
      ENDDO


C***  Copy the constant AVOIGT values to second_model 
C***  (all lines, entry at L=1)
      DO IMOD=2, NMOD
         DO NBL = 1, NBLINE
          AVOIGT(NBL,1,IMOD) = AVOIGT(NBL,1,1)
         ENDDO
      ENDDO

   99 CONTINUE
      RETURN
      END
