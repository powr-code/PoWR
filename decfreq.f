      SUBROUTINE DECFREQ (XLAMBDA, NF, NFDIM, TREF, OLDFGRID, KEY, 
     >                    MODOLD, XLAMBLUE)
C*******************************************************************************
C***  PURPOSE: READS THE FREQUENCY GRID XLAMBDA(K) (WAVELENGTHS IN A) 
C***       EITHER: INPUT FROM TAPE6 = FGRID
C***       ALTERNATIVELY: IF (OLDFGRID) : TAKEN FROM OLD MODEL FILE
C***       - NOTE: ONLY THE HAND-MADE POINTS ARE DEFINED HERE, 
C***               FURTHER POINTS (EDGES!) ARE INSERTED BY SUBR. FGRID
C***  CALLING TREE: WRSTART - FGRID - DECFREQ  
C*******************************************************************************
 
      COMMON /COMTEFF/ TEFF,TMIN,TMODIFY,SPHERIC
      CHARACTER KARTE*80, ACTPAR*20
      DIMENSION XLAMBDA(NFDIM), KEY(NFDIM)
      LOGICAL OLDFGRID
 
C***  IF NO REFERENCE TEMPERATURE IS SPECIFIED, TEFF IS DEFAULT
      TREF=TEFF

C***  BRANCH FOR OLDFGRID OPTION:READ FROM CHANNEL 9 = OLD MODEL FILE 
      IF (OLDFGRID) THEN

      CALL OPENMS (9, IDUMMY, IDUMMY, 1, IERR)
      IERR=1
      CALL READMS (9, MODOLD, 13, 'MODHEAD ', IERR)
      IF (IERR .LT. 0) THEN
         CALL REMARK (' OLD MODEL FILE NOT AVAILABLE FOR READING FGRID')
         PRINT *,     ' OLD MODEL FILE NOT AVAILABLE FOR READING FGRID'
         STOP 'ERROR'
         ENDIF
      CALL READMS (9, NF     , 1 , 'NF      ', IERR)
C***  ARRAY BOUND CHECK
      IF (NF .GT. NFDIM) THEN
         CALL REMARK (' OLD MODEL HAS TOO MANY FREQUENCY POINTS')
         PRINT *,      'OLD MODEL HAS TOO MANY FREQUENCY POINTS'
         STOP 'ERROR'
         ENDIF
      CALL READMS (9, XLAMBDA, NF, 'XLAMBDA ', IERR)
      CALL READMS (9, KEY    , NF, 'KEY     ', IERR)
      CALL CLOSMS (9, IERR)

C***  REMOVE ALL ENTRIES WHICH ARE NOT HAND-MADE FREQUENCY POINTS 
C***    (RECOGNIZED BY BLANK ENTRY IN ARRAY "KEY")
      KK = 0
      DO 10 K=1, NF
        IF (KEY(K) .NE. 8H         ) GOTO 10
        KK = KK + 1
        XLAMBDA(KK) = XLAMBDA(K)
   10 CONTINUE
      NF = KK 

      GOTO 7
      ENDIF
C***  END OF THE BRANCH with OLD FGRID ****************************

C*** If BLUEMOST is given in CARDS, use this number instead of reading FGRID 
      IF (XLAMBLUE .GT. .0) THEN
         NF = 1
         XLAMBDA(1) = XLAMBLUE
         RETURN
      ENDIF


C***  ELSE: DECODING INPUT CARDS FROM TAPE 7 = FGRID
      NF=0
      OPEN (7, FILE='FGRID', STATUS='UNKNOWN')
    1 READ (7,6, END=7) KARTE
    6 FORMAT(A)

      IF (KARTE(:1) .EQ.'*' ) THEN
            PRINT 2,KARTE
    2       FORMAT (1X,A)
            GOTO 1
            ENDIF
      IF (KARTE(:4) .EQ. 'TREF' ) THEN
            CALL SARGV (KARTE, 2, ACTPAR)
            READ (ACTPAR,'(F20.0)', ERR=99) TREF
            GOTO 1
            ENDIF
      NF=NF+1
      IF(NF.GT.NFDIM) THEN
            CALL REMARK ('FREQUENCY DIMENSION NFDIM INSUFFICIENT')
            STOP 'ERROR in SUBR. DECFREQ'
            ENDIF
            READ (KARTE,'(F20.0)', ERR=99) XLAMBDA(NF)
    5 FORMAT(F10.0)
      GOTO 1
 
C*********** ERROR STOPS ***********************************************

    7 IF(NF .LE. 0) THEN
            CALL REMARK ('*** NO FREQUENCY SCALE ENCOUNTERED')
            CALL REMARK ('You must either use OLD FGRID')
            CALL REMARK ('or provide a file FGRID')
            CALL REMARK ('or specify BLUEMOST=x.x in the CARDS file')
            STOP '*** FATAL ERROR in SUBR. DECFREQ'
            ENDIF
 
      DO 21 K=2,NF
      IF((XLAMBDA(K-1)-XLAMBDA(K))*(XLAMBDA(1)-XLAMBDA(NF)).LE..0) THEN
            CALL REMARK ('WAVELENGTH SCALE OUT OF SEQUENCE:')
            WRITE (0,'(F20.2)') XLAMBDA(K)
            STOP 'ERROR in SUBR. DECFREQ'
            ENDIF
   21 CONTINUE

      RETURN

C*************
   99 WRITE (0,*) '*** File FGRID: CANNOT DECODE FLOATING POINT NUMBER!'
      WRITE (0,*) 'The error occured in the following line:'
      WRITE (0,*) KARTE
      STOP 'ERROR in SUBR. DECFREQ'

      END
