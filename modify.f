C***  MAIN PROGRAM MODIFY  *********************************************
      SUBROUTINE MODIFY
C***********************************************************************
C***  THIS PROGRAM INTERPOLATES THE LEVEL POPULATIONS BETWEEN TWO
C***  GIVEM POINTS
C***********************************************************************
 
      IMPLICIT NONE

C***  DEFINE ARRAY DIMENSIONS ******************************************
C***  IRON: ADD GENERIC ION TO MAXATOM
      INTEGER, PARAMETER :: MAXATOM  =          26 
      INTEGER, PARAMETER :: NDIM     =        1560 
      INTEGER, PARAMETER :: NFDIM    = 2*NDIM + 400 
      INTEGER, PARAMETER :: MAXKONT  =     NFDIM/2 
      INTEGER, PARAMETER :: MAXKODR  =        NDIM 
      INTEGER, PARAMETER :: MAXIND   =       20000 
      INTEGER, PARAMETER :: MAXFEIND =        1500 
      INTEGER, PARAMETER :: NDDIM    =          89 
      INTEGER, PARAMETER :: MAXHIST  =        4000 
 
C***  HANDLING OF DIELECTRONIC RECOMBINATION / AUTOIONIZATION (SUBR. DATOM)
      INTEGER, PARAMETER :: MAXAUTO = 2850 
      INTEGER, DIMENSION(MAXAUTO) :: LOWAUTO, IONAUTO, KRUDAUT
      REAL, DIMENSION(MAXAUTO) :: WAUTO, EAUTO, AAUTO
      COMMON / COMAUTO / LOWAUTO, WAUTO, EAUTO, AAUTO, IONAUTO, KRUDAUT

      INTEGER, DIMENSION(MAXKODR) :: KODRNUP, KODRLOW
      INTEGER, DIMENSION(NDIM) :: MAINQN, NOM, IONGRND, NCHARG
      REAL, DIMENSION(NDIM) :: WEIGHT, ELEVEL, EION
      REAL, DIMENSION(NDIM, NDIM) :: EINST
      REAL, DIMENSION(4, NDIM) :: ALTESUM
      REAL, DIMENSION(4, MAXIND) :: COCO
      REAL, DIMENSION(MAXATOM) :: ABXYZ, ATMASS, STAGE
      INTEGER, DIMENSION(MAXATOM) :: NFIRST, NLAST, KODAT
      REAL, DIMENSION(NDDIM) :: RNE, ENTOT, T
      REAL, DIMENSION(NDDIM, NDIM) :: POPNUM
      REAL, DIMENSION(MAXKONT) :: ADDCON1, ADDCON2, ADDCON3,
     >                            ALPHA, SEXPO
      INTEGER, DIMENSION(MAXKONT) :: KONTNUP, KONTLOW
      CHARACTER*8 IGAUNT(MAXKONT), KEYCBF(MAXKONT)
      INTEGER, DIMENSION(MAXIND) :: INDLOW, INDNUP
      REAL, DIMENSION(MAXATOM,MAXATOM) :: SIGMATHK, SEXPOK, EDGEK
      CHARACTER(8*MAXHIST) :: MODHIST

      INTEGER :: I, N, L, IFRO, ITO, JOBNUM, LASTIND, NATOM, NAUTO, ND,
     >           LASTKON, LASTKDR, LASTFE, JOBNUM_SAVE, IDUMMY, IERR,
     >           LAST, LSPOP, IARG, NARG
      REAL :: A, B, X, Q, FEDUMM,
     >        VDOPFE, DXFE, XLAM0FE, XL, FROM, FTO, ASECOND

      CHARACTER(255) :: HISTENTRY
      CHARACTER(100) :: MODHEAD
      CHARACTER(80) :: KARTE
      CHARACTER(10), DIMENSION(NDIM) :: LEVEL
      CHARACTER(10), DIMENSION(MAXATOM) :: ELEMENT
      CHARACTER(7) :: MODE
      CHARACTER(4), DIMENSION(MAXIND) :: KEYCBB
      CHARACTER(2), DIMENSION(MAXATOM) :: SYMBOL
      CHARACTER(10), DIMENSION(10) :: ARGUMENT
      CHARACTER(64) :: BUFFER64
      CHARACTER(72) :: BUFFER72
      LOGICAL NOTEMP, BAUTO, NOPOP
 
C***  IRON: COMMON BLOCK FOR IRON-SPECIFIC DATA
C***  include "dimblock"
      INTEGER, PARAMETER :: INDEXMAX = 1E7, NFEREADMAX = 3E5    !std
C      INTEGER, PARAMETER :: INDEXMAX = 4E7, NFEREADMAX = 5E5     !vd20
C      INTEGER, PARAMETER :: INDEXMAX = 1E8, NFEREADMAX = 6E5     !xxl

      REAL, DIMENSION(NFEREADMAX) :: FEDUMMY
      INTEGER, DIMENSION(MAXFEIND) :: INDRB, INDRF, IFRBSTA, IFRBEND,
     >                                IFENUP, IFELOW
      REAL, DIMENSION(MAXFEIND) :: SIGMAINT
      REAL, DIMENSION(INDEXMAX) :: SIGMAFE
      COMMON /IRON/ FEDUMM, INDRB, INDRF, SIGMAFE, IFRBSTA, IFRBEND,
     >              IFENUP, IFELOW, SIGMAINT
      LOGICAL :: BFEMODEL

C***  Operating system:
      CHARACTER(8) :: OPSYS
      COMMON / COMOS / OPSYS

      CHARACTER(10) :: TIM1, TIM2

      !File and channel handles (=KANAL)
      INTEGER, PARAMETER :: hCPR = 0        !write to wruniqX.cpr (stderr)
      INTEGER, PARAMETER :: hMODINPUT = 1   !MODIFY_INPUT file
      INTEGER, PARAMETER :: hMODEL = 3      !MODEL file (fort.3)
      INTEGER, PARAMETER :: hHIST = 21      !write to MODHIST file
      INTEGER, PARAMETER :: hOUT = 6        !write to wruniqX.out (stdout)

C***  Link data to identify program version
      CHARACTER(30) :: LINK_DATE
      CHARACTER(10) :: LINK_USER
      CHARACTER(60) :: LINK_HOST
      COMMON / COM_LINKINFO / LINK_DATE, LINK_USER, LINK_HOST

      INTEGER, EXTERNAL :: IDX

C***  Write Link Data (Program Version) tp CPR file
      WRITE(hCPR,'(2A)') '>>> MODIFY started: Program Version from '
     >      , LINK_DATE
      WRITE(hCPR,'(4A)') '>>> created by ', LINK_USER(:IDX(LINK_USER)),
     >      ' at host ', LINK_HOST(:IDX(LINK_HOST))

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
     >            'MODIFY', INDEXMAX, NFEREADMAX, MAXFEIND,
     >             LASTFE, SIGMAFE, INDRB, INDRF,
     >             IFENUP, IFELOW, IFRBSTA, IFRBEND, FEDUMMY,
     >             VDOPFE, DXFE, XLAM0FE, SIGMAINT, BFEMODEL)
 
C***  READING OF THE MODEL FILE ----------------------------------------
      CALL OPENMS (hMODEL, IDUMMY, IDUMMY, 1, IERR)
      CALL READMS (hMODEL, ND, 1, 'ND      ', IERR)
      IF(ND > NDDIM) THEN
            CALL REMARK ('TOO MANY DEPTH POINTS')
            STOP 'ERROR'
            ENDIF
      CALL READMS (hMODEL, MODHEAD, 13, 'MODHEAD ' , IERR)
      CALL READMS (hMODEL, JOBNUM ,  1, 'JOBNUM  ' , IERR)
C***  Save JOBNUM
      JOBNUM_SAVE = JOBNUM
      CALL READMS (hMODEL, POPNUM ,ND*N,'POPNUM  ' , IERR)
      CALL READMS (hMODEL, RNE    , ND, 'RNE     ' , IERR)
      CALL READMS (hMODEL, ENTOT  , ND, 'ENTOT   ' , IERR)
      CALL READMS (hMODEL, T      , ND, 'T       ' , IERR)
C***  READ 'ABXYZ' (REL. ABUNDANCES OF ALL ELEMENTS) AND CHECK WHETHER
C***  THE RECORD EXISTS
      IERR=1
      CALL READMS (hMODEL,ABXYZ,NATOM,'ABXYZ   ',IERR)
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
      CALL READMS(hMODEL,LAST,1,'MODHIST ' , IERR)
      IF (LAST+9 >= MAXHIST) THEN
          !not really sufficient check because MODIFY can write more than one line
          ! so the limit might still be too small even if this check does not fail
          CALL REMARK ('MODHIST DIMENSION INSUFFICIENT')
          STOP 'ERROR'
      ENDIF
      CALL READMS (hMODEL,MODHIST,LAST,'MODHIST ', IERR)
C*********************************************************************** 

C***  DECODING INPUT CARDS  ********************************************
      OPEN (hMODINPUT, FILE='MODIFY_INPUT', STATUS='UNKNOWN')
      REWIND(hMODINPUT)
      LSPOP=-1
      MODE='NOTHING'
      NOTEMP = .FALSE.
      BAUTO = .FALSE.
      NOPOP = .FALSE.

C***  IF INPUT-FILE IS FINISHED (END=66) THEN UPDATE THE MODEL
    7 READ(hMODINPUT,'(A)',END=66) KARTE

      DO IARG=1, 10
       ARGUMENT(IARG) = ""
      ENDDO

      CALL SARGC (KARTE,NARG)
      DO IARG=1, NARG
       CALL SARGV (KARTE,IARG,ARGUMENT(IARG))
      ENDDO

      IF ( KARTE(:11) .EQ. 'AUTO MODIFY' ) THEN
C                           ===========
            BAUTO = .TRUE.

      ELSEIF ( KARTE(:9) .EQ. 'PRINT POP' ) THEN
C                              =========
            DECODE (80,4,KARTE) XL
    4       FORMAT (9X,F10.0)
            LSPOP=IFIX(XL)
            IF (LSPOP.EQ.0) LSPOP=1
C***        READ NEXT INPUT CARD
            GOTO 7

      ELSEIF ( ARGUMENT(1)(:9) .EQ. 'INTERPOLA' ) THEN
C                                    =========
            READ (ARGUMENT(4),'(F10.0)') FROM
            READ (ARGUMENT(7),'(F10.0)') FTO
            IFRO=IFIX(FROM)
            ITO=IFIX(FTO)
            MODE='INTERPO'
            JOBNUM = JOBNUM + 1
C***        DO INTERPOLATION
            GOTO 80

      ELSEIF ( ARGUMENT(1)(:5) .EQ. 'INNER' 
C                                    =====
     $    .OR. ARGUMENT(1)(:5) .EQ. 'OUTER') THEN
C                                    =====
            IF ( ARGUMENT(2)(:5) .EQ. 'EXTRA' ) THEN
C                                      =====
              READ (ARGUMENT(5),'(F10.0)') FROM
              READ (ARGUMENT(8),'(F10.0)') ASECOND
              IFRO=IFIX(FROM)
              ITO=IFIX(ASECOND)
              JOBNUM = JOBNUM + 1 
              IF ( ARGUMENT(1)(:5) .EQ. 'INNER' ) THEN
                IF (IFRO .LT. ITO) THEN
                  CALL REMARK('EXTRA IN: SECOND POINT WRONG')
                  STOP 'ERROR'
                ENDIF
                MODE='EXTR IN'
              ELSEIF (ARGUMENT(1)(:5) .EQ. 'OUTER') THEN  
                IF (IFRO .GT. ITO) THEN
                  CALL REMARK('EXTRA OUT: SECOND POINT WRONG')
                  STOP 'ERROR'
                ENDIF
                MODE='EXTROUT'
              ENDIF
C***          DO EXTRAPOLATION
              GOTO 80
         
            ELSE

              CALL REMARK('INNER / OUTER: SECOND ARGUMENT WRONG')
              STOP 'ERROR'              

            ENDIF

      ELSEIF (KARTE(:7) .EQ. 'NO TEMP') THEN
C                             =======
            NOTEMP = .TRUE.
C***        READ NEXT INPUT CARD
            GOTO 7

      ELSEIF (KARTE(:4) .EQ. 'TEMP') THEN
C                             ====
            NOTEMP = .FALSE.
C***        READ NEXT INPUT CARD
            GOTO 7

      ELSEIF (KARTE(:6) .EQ. 'NO POP') THEN
C                             ======
            NOPOP = .TRUE.
C***        READ NEXT INPUT CARD
            GOTO 7

      ELSEIF (KARTE(:3) .EQ. 'POP') THEN 
C                             ===
            NOPOP = .FALSE.
C***        READ NEXT INPUT CARD
            GOTO 7
            
      ENDIF

      GOTO 7
C***  LOOP OF DECODING-OF-INPUT-CARDS  **********************************

C***  INTER- OR EXTRAPOLATION OF THE POPULATIONS  **********************
C***  ENTRY-POINT FOR INTERPOLATION AND EXTRAPOLATION
   80 CONTINUE
          CALL INTEPO (ND,N,RNE,NCHARG,POPNUM,ENTOT,NATOM,
     $                           ABXYZ,NFIRST,NLAST,IFRO,ITO,MODE,NOPOP)

C***  INTERPOLATION OF THE TEMPERATURE STRATIFICATION
      IF (.NOT. NOTEMP) THEN
          A = ALOG10(ENTOT(IFRO))
          B = ALOG10(ENTOT(ITO ))

          IF (MODE .EQ. 'INTERPO') THEN
              DO L=IFRO+1, ITO-1
               X = ALOG10(ENTOT(L))
               Q = (X-A) / (B-A)
               T(L) = (1.-Q) * T(IFRO) + Q * T(ITO)
              ENDDO
          ENDIF

          IF (MODE .EQ. 'EXTR IN') THEN
              IF (IFRO .EQ. ITO) THEN
                DO L=IFRO+1,ND
                 T(L) = T(IFRO)
                ENDDO
              ELSE
                DO L=IFRO+1,ND
                 X = ALOG10(ENTOT(L))
                 Q = (X-B) / (A-B)
                 T(L) = (1.-Q) * T(ITO) + Q * T(IFRO)
                ENDDO
              ENDIF
          ENDIF

          IF (MODE .EQ. 'EXTROUT') THEN
              IF ((IFRO .EQ. ITO)) THEN
                DO L=IFRO-1,1,-1
                 T(L) = T(IFRO)
                ENDDO
              ELSE
                DO L=IFRO-1,1,-1
                 X = ALOG10(ENTOT(L))
                 Q = (X-A) / (B-A)
                 T(L) = (1.-Q) * T(IFRO) + Q * T(ITO)
                ENDDO
              ENDIF
          ENDIF
      ENDIF

C***  TO AVOID NEGATIVE OR SMALL POPNUMBERS: CALL INHIBIT
      CALL INHIBIT(POPNUM, N, ND, NCHARG, RNE, 
     >             NATOM, ABXYZ, NFIRST, NLAST, 1.E-99)

C***  UPDATING THE MODEL HISTORY
C***  JOBNUM IS INCREASED FOR EVERY INTERPOLATION OR EXTRAPOLATION CARD
C***  JOBNUM is now set to JOBNUM_SAVE + 1 (if JOBNUM was increased)
      JOBNUM = MIN(JOBNUM, JOBNUM_SAVE+1)
C***  Now no Jobnumbers greater than 999 are written (now obsolete in i7 versions)
C      IF (JOBNUM <= 999) THEN
C        JOBNUM2 = JOBNUM
C      ELSE
C        JOBNUM2 = 999
C      ENDIF

      IF (MODE == 'INTERPO') THEN
C        ENCODE (55,50,MODHIST(LAST+1)) JOBNUM2,IFRO,ITO
        WRITE(UNIT=BUFFER64, FMT=50) JOBNUM, IFRO, ITO
   50   FORMAT ('/',I7,'. MODIFY  INTERPOLATION  FROM POINT ', 
     >          I3,' TO POINT ',I3,'    ')
        CALL ADDHISTENTRY(MODHIST,-1,MAXHIST,64,BUFFER64)
C        LAST=LAST+7 !56
      ENDIF
      IF (MODE == 'EXTR IN') THEN
C        ENCODE (63,51,MODHIST(LAST+1)) JOBNUM2,IFRO,ITO
        WRITE(UNIT=BUFFER72, FMT=51) JOBNUM, IFRO, ITO
   51   FORMAT ('/',I7,'. MODIFY  INNER EXTRAPOLATION FROM POINT ',I3,
     >  ' SECOND POINT ',I3,'  ') 
        CALL ADDHISTENTRY(MODHIST,-1,MAXHIST,72,BUFFER72)
C        LAST = LAST + 8 !64
      ENDIF
      IF (MODE == 'EXTROUT') THEN
C        ENCODE (63,52,MODHIST(LAST+1)) JOBNUM2,IFRO,ITO
        WRITE(UNIT=BUFFER72, FMT=52) JOBNUM, IFRO, ITO
   52   FORMAT ('/',I7,'. MODIFY  OUTER EXTRAPOLATION FROM POINT ',I3,
     >  ' SECOND POINT ',I3,'  ')
C        LAST = LAST + 8 !64
        CALL ADDHISTENTRY(MODHIST,-1,MAXHIST,72,BUFFER72)
      ENDIF
      IF (NOTEMP) THEN
C        LAST = LAST + 1 !8
C        MODHIST(LAST) = ' NOTEMP'
        CALL ADDHISTENTRY(MODHIST,-1,MAXHIST,8,' NOTEMP ')
      ENDIF
      IF (BAUTO) THEN
C        LAST = LAST + 1 !8
C        MODHIST(LAST) = ' AUTO'
        CALL ADDHISTENTRY(MODHIST,-1,MAXHIST,8,' AUTO   ')
      ENDIF
C      MODHIST(1)=LAST
      CALL WRITMS (hMODEL,MODHIST,MAXHIST,'MODHIST ',-1, IDUMMY, IERR)
C***  LOOK FOR NEXT INPUT-CARD
      GOTO 7
 
C***  ENTRY, IF INPUT FILE IS FINISHED
C***  UPDATING THE MODEL, IF MODE <> 'NOTHING'
   66 CONTINUE
      IF (MODE .EQ. 'NOTHING') THEN
         CALL REMARK ('NO INTERPOLATION OPTION DECODED')
         PRINT *, 'MODIFY: NO INTERPOLATION OPTION DECODED'
         GOTO 99
      ENDIF

C***  CYCLIC CHANGING OF THE ARRAY NAMES FOR THE LAST THREE ITERATIONS
      CALL CHANGE ('POP3    '  ,'HELP    '  , hMODEL)
      CALL CHANGE ('POP2    '  ,'POP3    '  , hMODEL)
      CALL CHANGE ('POP1    '  ,'POP2    '  , hMODEL)
      CALL CHANGE ('POPNUM  '  ,'POP1    '  , hMODEL)
      CALL CHANGE ('HELP    '  ,'POPNUM  '  , hMODEL)
      CALL CHANGE ('TOLD3   '  ,'THELP   '  , hMODEL)
      CALL CHANGE ('TOLD2   '  ,'TOLD3   '  , hMODEL)
      CALL CHANGE ('TOLD1   '  ,'TOLD2   '  , hMODEL)
      CALL CHANGE ('T       '  ,'TOLD1   '  , hMODEL)
      CALL CHANGE ('THELP   '  ,'T       '  , hMODEL)

      CALL WRITMS (hMODEL, POPNUM, ND*N, 'POPNUM   ', -1, IDUMMY, IERR) 
      CALL WRITMS (hMODEL, RNE   ,   ND, 'RNE      ', -1, IDUMMY, IERR)
      CALL WRITMS (hMODEL, T     ,   ND, 'T        ', -1, IDUMMY, IERR)
      CALL WRITMS (hMODEL, JOBNUM ,  1,  'JOBNUM   ', -1, IDUMMY, IERR)

      IF (LSPOP.GT.0)
     $   CALL    PRIMO (ND,N,RNE,LEVEL,POPNUM,JOBNUM,MODHEAD,LSPOP,
     $                  IFRO,ITO)
 
C***  CLOSE FILES
   99 CLOSE (hMODINPUT)
      CALL CLOSMS (hMODEL, IERR)
      CALL JSYMSET ('G0','0')

      !write model history entry into explicit history file
      CALL GETHISTENTRY(HISTENTRY,JOBNUM,MODHIST,MAXHIST)
      OPEN (hHIST, FILE='MODHIST', STATUS='UNKNOWN',
     >             ACTION='READWRITE', POSITION='APPEND')
      WRITE (hHIST,FMT='(A)') TRIM(ADJUSTL(HISTENTRY))
      CLOSE(hHIST)

      CALL STAMP (OPSYS, 'MODIFY', TIM1)

      STOP 'O.K.'
      END
