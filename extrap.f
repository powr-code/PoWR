C***  MAIN PROGRAM EXTRAP  *********************************************
      SUBROUTINE EXTRAP
C***********************************************************************
C***  THIS PROGRAM MAKES A (LOGARITHMIC) EXTRAPOLATION OF THE POPULATION
C***  NUMBERS FROM THE LAST THREE ITERATIONS
C***  AITKEN: EXTRAPOLATION OF POP.NUMBER FROM LAST 3 ITERATIONS AT
C***          SINGLE DEPTH POINT L
C***  NG: EXTRAPOLATION OF POP.NUMBER FROM LAST 3 ITERATIONS AT ALL
C***      DEPTH POINTS L (L=1...ND)
C***      IF (.NOT.NOTEMP): ALSO EXTRAPOLATION OF TEMPERATURE T(R)
C***  DEFAULT EXTRAPOLATION: ---  NG  ---
C***  THE NG-EXTRAPOLATION IS CONTROLLED FOR THE INCREASE OF POPNUMBERS
C***********************************************************************
 
C***  DEFINE ARRAY DIMENSIONS ******************************************
      PARAMETER ( MAXATOM =          26 )
      PARAMETER ( NDIM    =        1560 )
      PARAMETER ( NFDIM   = 2*NDIM + 400 )
      PARAMETER ( MAXKONT =     NFDIM/2 )
      PARAMETER ( MAXKODR =        NDIM )
      PARAMETER ( MAXIND  =       20000 )
      PARAMETER ( MAXFEIND  =       1500 )
      PARAMETER ( NDDIM   =          89 )
      PARAMETER ( MAXHIST =        4000 )
C***  NUMBER OF ENTRYS STORED IN THE GAMMA HISTORY
      PARAMETER (MAXGAHIST = 100)
 
C***  HANDLING OF DIELECTRONIC RECOMBINATION / AUTOIONIZATION (SUBR. DATOM)
      PARAMETER ( MAXAUTO = 2850 )
      COMMON / COMAUTO / LOWAUTO(MAXAUTO),WAUTO(MAXAUTO)
     $                  ,EAUTO(MAXAUTO),AAUTO(MAXAUTO),IONAUTO(MAXAUTO)
     $                  ,KRUDAUT(MAXAUTO)
      COMMON // NCHARG(NDIM),WEIGHT(NDIM),ELEVEL(NDIM),EION(NDIM)
     $ ,IONGRND(NDIM),KODRNUP(MAXKODR),KODRLOW(MAXKODR)
     $ ,MAINQN(NDIM),EINST(NDIM,NDIM),NOM(NDIM)
     $ ,ALTESUM(4,NDIM),COCO(4,MAXIND)
     $ ,ABXYZ(MAXATOM),KODAT(MAXATOM),ATMASS(MAXATOM),STAGE(MAXATOM)
     $ ,NFIRST(MAXATOM),NLAST(MAXATOM)
     $ ,RNE(NDDIM),RADIUS(NDDIM),W(NDDIM)
     $ ,TNEW(NDDIM),T(NDDIM),TOLD1(NDDIM),TOLD2(NDDIM), TOLD3(NDDIM)
     $ ,POPNEW(NDDIM,NDIM),POPNUM(NDDIM,NDIM),POP1(NDDIM,NDIM)
     $ ,POP2(NDDIM,NDIM), POP3(NDDIM,NDIM)
     $ ,ALPHA(MAXKONT),SEXPO(MAXKONT)
     $ ,ADDCON1(MAXKONT), ADDCON2(MAXKONT), ADDCON3(MAXKONT) 
     $ ,KONTNUP(MAXKONT),KONTLOW(MAXKONT)
     $ ,INDLOW(MAXIND),INDNUP(MAXIND)
     $ ,MDUMMY(MAXHIST)
      CHARACTER*8 IGAUNT(MAXKONT), KEYCBF(MAXKONT)
      CHARACTER(MAXHIST*8) :: MODHIST
      DIMENSION GAHIST(26,MAXGAHIST)
C***  PREVENT UNNECESSARY OUTPUT FROM SUBR. PRICORR:
      COMMON / COMNEGI / NEGINTL,INEGMIN,INEGMAX,LNEGMIN,LNEGMAX
      COMMON / COMNEGT / NEGT
      COMMON / COMITWA / ITWARN, ITMAX
      DATA NEGINTL, NEGT, ITWARN / 0, 0, 0 /
      DIMENSION SIGMATHK(MAXATOM,MAXATOM),SEXPOK(MAXATOM,MAXATOM)
      DIMENSION EDGEK(MAXATOM,MAXATOM)
      CHARACTER KARTE*80,MODHEAD*100
      CHARACTER LEVEL(NDIM)*10
      CHARACTER*10 ELEMENT(MAXATOM), ARG(5)
      CHARACTER*4 KEYCBB(MAXIND)
      CHARACTER*2 SYMBOL(MAXATOM)
      LOGICAL NGACC, NOTEMP, NOEXTEMP, NGWEIGHT, BEX4
 
C***  IRON: COMMON BLOCK FOR IRON-SPECIFIC DATA
C***  include "dimblock"
      INTEGER, PARAMETER :: INDEXMAX = 1E7, NFEREADMAX = 3E5    !std
C      INTEGER, PARAMETER :: INDEXMAX = 4E7, NFEREADMAX = 5E5     !vd20
C      INTEGER, PARAMETER :: INDEXMAX = 1E8, NFEREADMAX = 6E5     !xxl

      COMMON /IRON/ FEDUMMY(NFEREADMAX),
     >              INDRB(MAXFEIND),INDRF(MAXFEIND), SIGMAFE(INDEXMAX),
     >              IFRBSTA(MAXFEIND), IFRBEND(MAXFEIND),
     >              IFENUP(MAXFEIND), IFELOW(MAXFEIND),
     >              SIGMAINT(MAXFEIND)
      LOGICAL :: BFEMODEL

C***  Operating system:
      COMMON / COMOS / OPSYS
      CHARACTER*8 OPSYS

      CHARACTER TIM1*10, TIM2*10

      IF (OPSYS .EQ. 'CRAY' .OR. OPSYS .EQ. 'SGI') THEN
        CALL CLOCK(TIM1)
      ELSE
        CALL TIME(TIM1)
      ENDIF

      CALL       DATOM (NDIM,N,LEVEL,NCHARG , WEIGHT,ELEVEL,EION,MAINQN,
     $                  EINST,ALPHA,SEXPO,
     $                  ADDCON1, ADDCON2, ADDCON3, 
     $                  IGAUNT,COCO,KEYCBB,ALTESUM,
     $                  INDNUP,INDLOW,LASTIND,MAXIND,MAXATOM,NATOM,
     $                  ELEMENT,SYMBOL,NOM,KODAT,ATMASS,STAGE,
     $                  SIGMATHK,SEXPOK,EDGEK,NFIRST,
     $                  NLAST,NAUTO,MAXAUTO,LOWAUTO,WAUTO,EAUTO,AAUTO,
     $                  IONAUTO,KRUDAUT,KONTNUP,KONTLOW,LASTKON,MAXKONT,
     $                  IONGRND,KODRNUP,KODRLOW,LASTKDR,MAXKODR,KEYCBF,
C***  IRON: ADDITIONAL PARAMETERS FOR IRON-GROUP LINE BLANKETING
     >            'EXTRAP', INDEXMAX, NFEREADMAX, MAXFEIND,
     >             LASTFE, SIGMAFE, INDRB, INDRF,
     >             IFENUP, IFELOW, IFRBSTA, IFRBEND, FEDUMMY,
     >             VDOPFE, DXFE, XLAM0FE, SIGMAINT, BFEMODEL)
 
C***  READING OF THE MODEL FILE  ---------------------------------------
      CALL OPENMS(3, IDUMMY, IDUMMY, 1, IERR)
      CALL READMS (3,JOBNUM,1,'JOBNUM  ', IERR)
      CALL READMS (3,MODHEAD,13,'MODHEAD ', IERR)
      CALL READMS(3,ND,1,'ND      ', IERR)
      IF (ND.GT.NDDIM) THEN
            CALL REMARK ('TOO MANY DEPTH POINTS')
            STOP 'ERROR'
            ENDIF
      CALL READMS (3,RADIUS,ND,  'R       ', IERR)
      CALL READMS (3,POPNUM,ND*N,'POPNUM  ', IERR)
      CALL READMS (3,POP1  ,ND*N,'POP1    ', IERR)
      CALL READMS (3,POP2  ,ND*N,'POP2    ', IERR)
      CALL READMS (3,POP3  ,ND*N,'POP3    ', IERR)
C***  READING TEMPERATURE ARRAYS FOR EXTRAPOLATION OF T(R)
      CALL READMS (3,T    ,ND,   'T       ', IERR)
      IERR1=1
      CALL READMS (3,TOLD1,ND,   'TOLD1   ',IERR1)
      IERR2=1
      CALL READMS (3,TOLD2,ND,   'TOLD2   ',IERR2)
      IERR3=1
      CALL READMS (3,TOLD3,ND,   'TOLD3   ',IERR2)
      IF (((IERR1 .LT. 0) .AND. (IERR1 .NE. -10)) .OR. 
     $    ((IERR2 .LT. 0) .AND. (IERR2 .NE. -10)) .OR. 
     $    ((IERR3 .LT. 0) .AND. (IERR3 .NE. -10))) THEN
         CALL REMARK ('ERROR WHEN READING TOLD FROM MODEL FILE')
         STOP 'TOLD'
      ENDIF

C***  READ 'ABXYZ' (REL. ABUNDANCES OF ALL ELEMENTS) AND CHECK WHETHER
C***  THE RECORD EXISTS
      IERR=1
      CALL READMS (3,ABXYZ,NATOM,'ABXYZ   ',IERR)
      IF ((IERR .LT. 0) .AND. (IERR .NE. -10)) THEN
         CALL REMARK ('ERROR WHEN READING ABXYZ FROM MODEL FILE')
         STOP 'ERROR'
      ENDIF
C***  NOT EXISTING RECORD 'ABXYZ': DEFAULT IS AN ATOMIC DATA FILE "DATOM"
C***  CONTAINING "HELIUM" AS THE ONLY ELEMENT
      IF (IERR .EQ. -10) THEN
         IF (NATOM .EQ. 1) THEN
            ABXYZ(1)=1.
         ELSE
            CALL REMARK ('NOT EXISTING RECORD ABXYZ')
            STOP 'ERROR'
         ENDIF
      ENDIF
      CALL READMS(3,LAST,1,'MODHIST ', IERR)
      IF (LAST >= MAXHIST) THEN
            CALL REMARK ('MODHIST DIMENSION INSUFFICIENT')
            STOP 'ERROR'
            ENDIF
      CALL READMS (3,MODHIST,LAST,'MODHIST ', IERR)

C***  Read GAMMA History
      CALL READMS (3, GAHIST, 26*MAXGAHIST, 'GAHIST  ',IERR)
C***  Move Entries in GAMMA HISTORY
      DO I=MAXGAHIST-1, 1, -1
        DO J=1, 26
          GAHIST(J,I+1) = GAHIST(J,I)
        ENDDO
      ENDDO
C***  This Job is EXTRAP
      GAHIST(2,1) = 1.
C***  No Broyden Statistics and GAMMAS
        GAHIST( 3,1) = 0.
        GAHIST( 4,1) = 0.
        GAHIST( 5,1) = 0.
        GAHIST( 6,1) = 0.
        GAHIST(22,1) = 0.
        GAHIST(23,1) = 0.
        GAHIST(24,1) = 0.
        GAHIST(25,1) = 0.

C***  DEFAULTS (INPUT OPTIONS):  ***************************************
      LSPOP    = -1
      NGACC    = .TRUE.
      NOTEMP   = .FALSE.
      NOEXTEMP = .FALSE.
      BEX4     = .TRUE.
      NEWWRC   = 6

C***  DECODING INPUT CARDS  --------------------------------------------
      OPEN (1, FILE='CARDS', STATUS='UNKNOWN')
      REWIND 1
    7 READ(1,2, END=6) KARTE
    2 FORMAT (A)
      CALL SARGC (KARTE, NPAR)
      IF (NPAR .LT. 1) GOTO 7
      IF (NPAR .GT. 5) NPAR = 5
      DO 3 I=1, NPAR
    3 CALL SARGV (KARTE, I, ARG(I))

      IF ( KARTE(:8) .EQ. 'NO TEMPE' ) THEN
C                          ========
            NOTEMP=.TRUE.
            IF (KARTE(30:40) .NE. ' ')
     $          CALL DECNOT (LAST,MODHIST,KARTE,NOTEMP,'EXTRAP')
            GOTO 7
            ENDIF
      IF ( KARTE(:9) .EQ. 'PRINT POP' ) THEN
C                          =========
            DECODE (80,4,KARTE) XL
    4       FORMAT (9X,F10.0)
            LSPOP=IFIX(XL)
            IF (LSPOP.EQ.0) LSPOP=1
            GOTO 7
            ENDIF
      IF (ARG(1) .EQ. 'EXTRAP' .AND. ARG(2) .EQ. 'NOTEMP') THEN
C                      ======                     ======
            NOEXTEMP = .TRUE.
            GOTO 7
            ENDIF
      IF (ARG(1) .EQ. 'EXTRAP' .AND. ARG(2) .EQ. 'NGWEIGHT') THEN
C                      ======                     ======
            NGWEIGHT = .TRUE.
            GOTO 7
            ENDIF

      IF ( KARTE(:6) .EQ. 'AITKEN' ) THEN
C                          ======
         NGACC = .FALSE.
         BEX4 = .FALSE.
         GOTO 7
         ENDIF
      IF (ARG(1) .EQ. 'NG3') THEN
C                      ===
            BEX4 = .FALSE.
            GOTO 7
            ENDIF

      IF (ARG(1) .EQ. 'NEWWRC') THEN
C                      ======
            READ (ARG(2),'(F10.0)') XL
            NEWWRC=IFIX(XL)
            GOTO 7
            ENDIF

      GOTO 7
    6 CLOSE (1)

      IF (NOEXTEMP) NOTEMP = .TRUE.

C***  NO TEMPERATURE EXTRAPOLATION IN CASE OF NON-EXISTING RECORDS
C***                               "TOLD."  (OLD MODEL FILE)
      IF ((IERR1 .EQ. -10) .OR. (IERR3 .EQ. -10)) NOTEMP=.TRUE.

      IF (BEX4) THEN
         CALL NG4 (ND, N, RNE, NCHARG, 
     >             POPNEW, POPNUM, POP1, POP2, POP3, 
     >             NATOM, ABXYZ, NFIRST, NLAST, 
     >             TNEW, T, TOLD1, TOLD2, TOLD3, NOTEMP, TRESH, NOUT)
      ELSE
C***  LOGARITHMIC EXTRAPOLATION  ***************************************
        IF (.NOT. NGACC) THEN
C***          THIS VERSION: NO EXTRAPOLATION OF T(R)
          DO 10 L=1,ND
   10     TNEW(L)=T(L)
          CALL AITKEN (ND,N,RNE,NCHARG,POPNEW,POPNUM,POP1,POP2,NATOM,
     $                   ABXYZ,NFIRST,NLAST, TNEW, T, TOLD1, TOLD2, 
     $                   NOTEMP)
        ELSE
          CALL NG3 (ND,N,RNE,NCHARG,POPNEW,POPNUM,POP1,POP2,NATOM,
     $      ABXYZ,NFIRST,NLAST,RADIUS,W,TNEW,T,TOLD1,TOLD2,
     $      NOTEMP, NGWEIGHT, NDONE, NTDONE)
        ENDIF
      ENDIF
 
C***  CYCLIC CHANGING OF THE ARRAY NAMES FOR THE LAST THREE ITERATIONS
      CALL CHANGE ('POP3    ','HELP    ', 3)
      CALL CHANGE ('POP2    ','POP3    ', 3)
      CALL CHANGE ('POP1    ','POP2    ', 3)
      CALL CHANGE ('POPNUM  ','POP1    ', 3)
      CALL CHANGE ('HELP    ','POPNUM  ', 3)
      CALL CHANGE ('TOLD3   ','THELP   ', 3)
      CALL CHANGE ('TOLD2   ','TOLD3   ', 3)
      CALL CHANGE ('TOLD1   ','TOLD2   ', 3)
      CALL CHANGE ('T       ','TOLD1   ', 3)
      CALL CHANGE ('THELP   ','T       ', 3)
 
C***  UPDATING THE MODEL FILE  -----------------------------------------
      JOBNUM=JOBNUM+1
C      IF (JOBNUM .GE. 1000) JOBNUM=JOBNUM-1000
      CALL WRITMS (3,JOBNUM,1,'JOBNUM  ',-1, IDUMMY, IERR)
      CALL WRITMS (3,POPNEW,ND*N,'POPNUM  ',-1, IDUMMY, IERR)
      CALL WRITMS (3,RNE,ND,'RNE     ',-1, IDUMMY, IERR)
      CALL WRITMS (3,TNEW,ND,'T       ',-1, IDUMMY, IERR)
 
C***  PRINTOUT  -------------------------------------------------------
      IF (LSPOP.GT.0)
     $  CALL PRIEX (ND,N,RNE,LEVEL,POPNUM,JOBNUM,MODHEAD,LSPOP,
     $              TNEW,NOTEMP)
      CALL INHIBIT (POPNEW, N, ND, NCHARG, RNE, 
     $              NATOM, ABXYZ, NFIRST, NLAST, 1.E-99)


C***  Commented Lines; full call as in STEAL
      CALL PRICORR (POPNEW, POPNUM, LEVEL, N, ND, MODHEAD, LSPOP,
C*   $              CORMAX, RTCMAX, JOBNUM, REDUCE,
     $              CORMAX, RTCMAX, JOBNUM, 1.,      
C*   >              GAMMAC, GAMMAL, GAMMAR, GAMMAD,
     >               0.,      0.,      0.,     0., 
C*   $              T, TNEW, EPSILON, DELTAC, SMPOP,
     $              T, TNEW, 0.,       1.,    0.,  
C*   $              BUNLU,  DUNLU, DUNLUR, DUNLU2, BTND, TNDCORR,
     >              .NOT.NOTEMP, 0.,    0.,     0.,   .FALSE.,  0., 
C*   >              GAHIST, MAXGAHIST, STHLP,
     >              GAHIST, MAXGAHIST, .FALSE., 
C*   >              IWARN_NEG_XJCAPP, IWARN_NEG_XJLAPP)
     >              0,                 0              )

      IF (NOTEMP) THEN
        GAHIST(26,1) = 0.
      ELSE
        GAHIST(26,1) = 1.
      ENDIF

      WRITE (*,'(1X,3A,1E8.2,A,I2)')
     >  'EXTRAP: POP. Numbers below treshhold and at outer points ', 
     >  'are not taken into ',
     >  'account : TRESH=',TRESH, '   NOUT=', NOUT

      IF (BEX4) THEN
        PRINT *,'This was a new NG4-extrapolation -------------------'
      ELSE
        IF (NGACC) THEN
          PRINT *,
     >      'This was an old NG3-extrapolation -------------------'
          PRINT *,'THE EXTRAPOLATION WAS DONE FOR ',NDONE,' OF ',ND,
     >            ' * ', N ,' POPNUMBERS AND ',NTDONE,' TEMPRERATURES'
        ELSE
          PRINT *,'This was an AITKEN-extrapolation  -------------'
        ENDIF
      ENDIF

C***  UPDATING THE MODEL HISTORY  --------------------------------------
      IF (CORMAX > 1.E-100) CORMAX=ALOG10(CORMAX)
      IF (BEX4) THEN
        WRITE(UNIT=BUFFER40, FMT=22) JOBNUM, CORMAX
C        ENCODE (32,22,MODHIST(LAST+1)) JOBNUM,CORMAX
C***  NOTE : The word Cor. must not be written in capitals in order to 
C***         distinguish from COR. of the Main Program STEAL
   22   FORMAT ('/',I7,'. EXTRAP4 Cor.=',F8.4)
        CALL ADDHISTENTRY(MODHIST,-1,MAXHIST,40,BUFFER40)
        IF (NOTEMP) THEN
          CALL ADDHISTENTRY(MODHIST,-1,MAXHIST,8,' NOTEMP ')
        ENDIF
      ELSE
        IF (.NOT. NGACC) THEN
          WRITE(UNIT=BUFFER40, FMT=20) JOBNUM, CORMAX
   20     FORMAT ('/',I7,'. EXTRAP  COR.=',F8.4)
          CALL ADDHISTENTRY(MODHIST,-1,MAXHIST,40,BUFFER40)
        ELSE
          WRITE(UNIT=BUFFER40, FMT=30) JOBNUM, CORMAX
   30     FORMAT ('/',I7,'. EXTRAP(NG)  COR.=',F8.4)
          CALL ADDHISTENTRY(MODHIST,-1,MAXHIST,40,BUFFER40)
        ENDIF
        IF (NGACC) THEN
          WRITE(UNIT=BUFFER32, FMT=31) NDONE, ND, N
   31     FORMAT (I5,' DPS OF ',I2,'*',I3,' DONE')
          CALL ADDHISTENTRY(MODHIST,-1,MAXHIST,32,BUFFER32)
          IF (.NOT. NOTEMP) THEN
            WRITE(UNIT=BUFFER32, FMT=32) RTCMAX, NTDONE
   32       FORMAT ('   TC=',F8.4, ' DONE FOR ',I2,' DPS')
            CALL ADDHISTENTRY(MODHIST,-1,MAXHIST,32,BUFFER32)
          ENDIF
        ELSE
          IF (NOTEMP) THEN
            CALL ADDHISTENTRY(MODHIST,-1,MAXHIST,8,' NOTEMP ')
          ENDIF
        ENDIF
      ENDIF
      CALL WRITMS (3,MODHIST,MAXHIST,'MODHIST ',-1, IDUMMY, IERR)
      CALL WRITMS (3, GAHIST, 26*MAXGAHIST, 'GAHIST  ',-1, IDUMMY, IERR)
 
   99 CALL CLOSMS (3, IERR)

C***  NEXT JOB: REPEAT-CYCLE if NEWWRC NE 1
      IF (NEWWRC .NE. 1) THEN
        CALL REMARK ('EXTRAP: NEXTJOB=REPEAT')
        PRINT *, 'EXTRAP: NEXTJOB=REPEAT'
        CALL JSYMSET ('G1','REPEAT')
      ELSE
C***  NEXT JOB: WRCONT-CYCLE if NEWWRC EQ 1
        CALL REMARK ('EXTRAP: NEXTJOB=WRCONT')
        PRINT *, 'EXTRAP: NEXTJOB=WRCONT'
        CALL JSYMSET ('G1','WRCONT')
      ENDIF
      CALL JSYMSET ('G3','MOREJOBS')

      CALL JSYMSET ('G0','0')

      CALL STAMP (OPSYS, 'EXTRAP', TIM1)

      STOP 'O.K.'
      END
