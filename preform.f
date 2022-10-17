      SUBROUTINE PREFORM (KARTE,N,ELEVEL,LINE,INDLOW,INDNUP,LASTIND,
     $                   CLIGHT, VDOP, INDLAP, XLAMLAP, DELXLAP, ALN,
     $                   XLAM, NBLINE, MAXLAP, MAXIND, MAXATOM, 
     $                   LEVEL, WEIGHT, EINST, NDIM, POPNUM,
     >                   T, ND, NOM, NCHARG, EION, ENTOT, RNE, MAXSUBL,
     $                   NSUBLOW, NSUBNUP, BROAD, LINPRO,
     >                   AVOIGT, NMOD, NDDIM, MAXMOD, DENSCON, MAINQN,
     >                   MULTIIND, NMULTI, DD_VDOP, NATOM, IND_ORIGLEV)
C*******************************************************************************
C***  This subroutine reads atomic data from FORMAL_CARDS   
C***  CALLED FROM MAIN PROGRAM "FORMAL" 
C***  after DECFORM has encountered the input-line 
C***      BLEND (i.e. beginning of a BLEND-Block)
C***   or LINE ... (calculation of one isolated line)
C***   This subroutine is left when encountering the input-line
C***      -BLEND (i.e. end of a BLEND-Block)
C*******************************************************************************
 
C***  POPNUM enters, because MULTIPLET is inserting additional (LTE) levels
C***  The same holds for the the DRTRANSIT option
C***  Hence, all entering parmeters must be indexed with IMOD
      DIMENSION POPNUM(NDDIM, NDIM, MAXMOD)
      DIMENSION ND(MAXMOD)
      REAL, DIMENSION (NDDIM,MAXMOD) :: DENSCON, ENTOT, RNE, T
      DIMENSION DD_VDOP(NDDIM, MAXATOM, MAXMOD) 
        
      INTEGER, DIMENSION(NDIM) :: NOM
      DIMENSION ELEVEL(N)
      DIMENSION INDNUP(LASTIND),INDLOW(LASTIND), MULTIIND(LASTIND)
      INTEGER, DIMENSION(MAXSUBL) :: NSUBLOW, NSUBNUP
      DIMENSION INDLAP(MAXLAP),XLAMLAP(MAXLAP),DELXLAP(MAXLAP)
      DIMENSION AVOIGT(MAXLAP,NDDIM,MAXMOD)
      DIMENSION EINST(NDIM,NDIM), IND_ORIGLEV(NDIM)
      CHARACTER KARTE*80, INDSTR*4, KARTE2*80
      CHARACTER*8 LINPRO(MAXLAP), LINPROBL
      LOGICAL BLEND, BROAD, BFINDERR, BAIR, BVAC

C***  4 * PI
      DATA PI4 / 12.56637062 /

      ALN = ALOG (1. - VDOP / CLIGHT)

      BFINDERR = .FALSE.
      NMULTI = 0
      XLAM = -1.
      NBLINE = 0

C***  SAVE NUMBER OF ORIGINAL LEVELS AND LINES:
C***  MULTIPLET HANDLING (IN SUBR. MULTIPL) MAY INCREASE THESE NUMBERS 
      NNEW=N
      INDNEW=LASTIND

C***  levels which aren't the product 
C***   of "splitting" are their own "original levels" 
      DO NLEV=1, N
        IND_ORIGLEV(NLEV) = NLEV
      ENDDO

C***  DEFAULTS:
      BLEND=.FALSE.

      DO NBL = 1, MAXLAP
         LINPRO(NBL) = ''
         AVOIGT(NBL,1,1) = -1.
      ENDDO

C***  Check if the current CARD-line was "BLEND"  -------------------
C***    otherwise, it must have been an isolated LINE card -> re-read 
      IF (KARTE(:5) .EQ. 'BLEND') THEN 
C                         =====
         BLEND = .TRUE.
      ELSE
         BACKSPACE (2)
      ENDIF

C******************************************************************
C***  Begin of loop over all LINE, MULTIPLET or DRTRANSIT blocks
C******************************************************************

C***  Read next line, omit comment lines
    4 READ (2,'(A)',END=100) KARTE
      IF (KARTE(:1) .EQ. '*'  .OR.  KARTE(:1) .EQ. ' '
     >   .OR. KARTE(:1) .EQ. '.' ) GOTO 4
       
C***  Remove leading "+" signs (historical relic)
      IF (KARTE(:5) .EQ. '+LINE') KARTE=KARTE(2:)
      IF (KARTE(:10) .EQ. '+MULTIPLET') KARTE=KARTE(2:)
      IF (KARTE(:10) .EQ. '+DRTRANSIT') KARTE=KARTE(2:)

      IF (KARTE(:6) .EQ. '-BLEND') THEN
         IF (.NOT. BLEND) 
     >       WRITE (0,*) '*** WARNING: "-BLEND" encountered ',
     >                     'albeit no BLEND block was opened before'
         IF (NBLINE .EQ. 0)
     >       WRITE (0,*) '*** WARNING: BLEND-Block contained no lines'
         GOTO 20
       ENDIF

C******************************************************************
      IF (KARTE(:4) .EQ. 'LINE') THEN
C                         ====

         CALL SARGV (KARTE, 2, INDSTR)
         CALL FINDIND (INDBL,INDSTR,LEVEL,N,INDNUP,INDLOW,LASTIND
     >                ,BFINDERR, LOWBL, NUPBL)
         IF (BFINDERR) GOTO 4

C***     default wavelength 
         XLAMBL=1.E8/(ELEVEL(NUPBL)-ELEVEL(LOWBL))

C***     read wavelength if given, covert from AIR to VAC
         CALL SARGREST (KARTE, NPAR, 3, IRESTSTART, IRESTEND)
         IF (NPAR .GT. 0)
     >     CALL READ_LINECARD_PARAMETERS (KARTE(IRESTSTART:), XLAMBL,
     >            BROAD, LINPROBL, AVOIGTBL)

C***     In case that this line is the first one encounetred, 
C***        its wavelength serves as reference value:
         IF (NBLINE .EQ. 0) XLAM = XLAMBL

         CALL INSERT_LINE (LINE, INDBL, NBLINE, INDLAP, XLAMLAP, DELXLAP,
     $                   XLAMBL, LINPROBL, AVOIGTBL, XLAM, MAXLAP, ALN, 
     >                   ND, LINPRO, AVOIGT, NMOD, NDDIM, MAXMOD )

C***     for isolated line: finish
         IF (.NOT. BLEND) GOTO 20

C*****************************************************************
      ELSEIF (KARTE(:9) .EQ. 'MULTIPLET') THEN
C                             ---------

         CALL SARGV (KARTE, 2, INDSTR)
         CALL FINDIND (INDBL,INDSTR,LEVEL,N,INDNUP,INDLOW,LASTIND
     >                ,BFINDERR, LOW, NUP)

C***     If LEVELS were not found: Skip the whole multiplet
         IF (BFINDERR) THEN
   44       READ (2,'(A)') KARTE
            IF (KARTE(:10) .EQ. '-MULTIPLET') GOTO 4
            GOTO 44
         ENDIF

C***    Store LINE index for checking if this transition was split
            NMULTI = NMULTI + 1
            MULTIIND(NMULTI) = INDBL

            CALL MULTIPLE (XLAM, LINE, LOW, NUP, INDLAP, XLAMLAP,
     >                 DELXLAP, NBLINE,MAXLAP,INDLOW,INDNUP,INDNEW,MAXIND,
     >                 LEVEL,WEIGHT,ELEVEL,NNEW,EINST,NDIM,
     >                 POPNUM, T, ND, ALN, VDOP,
     >                 MAXSUBL, NSUBLOW, NSUBNUP, BROAD, LINPRO, 
     >                 AVOIGT, NMOD, MAXMOD, NDDIM, 
     >                 MAINQN, NCHARG, EION, NOM, IND_ORIGLEV)



      ELSE IF (KARTE(:9) .EQ. 'DRTRANSIT') THEN
C                              ---------

         CALL FINDLDR (LOW, NUP, LEVEL, NOM, NCHARG, N, BFINDERR)

C***     If LEVELS were not found: Skip the whole DRTRANSIT block
         IF (BFINDERR) THEN
  55       READ (2,'(A)') KARTE
            IF (KARTE(:10) .EQ. '-DRTRANSIT') GOTO 4
            GOTO 55
         ENDIF

         CALL DRTRANS (XLAM, LINE, LOW, NUP, INDLAP, XLAMLAP, DELXLAP,
     >                 NBLINE, MAXLAP, INDLOW, INDNUP, INDNEW, MAXIND,
     >                 LEVEL, WEIGHT, ELEVEL, NNEW, EINST, NDIM,
     >                 POPNUM, T, ND,
     >                 ALN, VDOP, EION, ENTOT, RNE,
     >                 MAXSUBL, NSUBLOW, BROAD, LINPRO, 
     >                 AVOIGT, DENSCON, NMOD, MAXMOD, NDDIM, 
     >                 MAINQN, NCHARG, NOM, IND_ORIGLEV) 

C***  ERROR branch *******************************
      ELSE
         WRITE (0,*) 'ERROR: INVALID OPTION FOUND'
         WRITE (0,*) KARTE(:IDX(KARTE))
      ENDIF

      GOTO 4
C***  END OF BLEND BLOCK  ****************************************


   20 CONTINUE

      RETURN

  100 WRITE (0,*) '*** UNEXPECTED E-O-F in FORMAL_CARDS ***'
      WRITE (0,*) '-BLEND missing?'
      RETURN

      END
