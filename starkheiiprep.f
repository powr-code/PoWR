      SUBROUTINE STARKHEIIPREP 
     >   (MAINQN, NUP, LOW, LINPRO, LEVEL, PATH_VCSSB)
C*******************************************************************
C***  Called from: SUBROUTINE STARKBROAD 
C***  This Subr. reads the line-broadening tables for He II
C***  The tables are searched for the transition low-nup
C*******************************************************************
      CHARACTER LINPRO*8
      CHARACTER*10 LEVEL(2)
      CHARACTER*256 FILENAME
      CHARACTER*256 PATH_VCSSB
      DIMENSION MAINQN(2)

      PARAMETER (MPDIM     = 50)     ! Max. # OF PROFILE POINTS IN vcssb TABLE
      PARAMETER (MTEMPDIM  =  6)     ! Max. # OF TEMP POINTS IN vcssb TABLE
      PARAMETER (MNEDIM    = 11)     ! Max. # OF DENSITY POINTS IN vcssb TABLE

      COMMON /VCSSBDAT/ MP, MNE, MTEMP, ALOG_NEL(MNEDIM), 
     >                   ALOG_T0,     ALOG_T_INC,
     >                   ALOG_ALPHA0, ALOG_ALPHA_INC, 
     >                   SVCS(MPDIM,MTEMPDIM,MNEDIM)

      IF (PATH_VCSSB .EQ. 'default') THEN  
         FILENAME = '/home/corona/wrh/work/wrdata/VCSSB.DAT'
      ELSE
         FILENAME = PATH_VCSSB(:IDX(PATH_VCSSB)) // '/' // 'VCSSB.DAT'
      ENDIF
      KANAL=13

      MAINQNLOW = MAINQN(LOW)
      MAINQNNUP = MAINQN(NUP)
C***  Check if Principle Quantum numbers are known
      IF (MAINQNLOW .LE. 0) GOTO 93
      IF (MAINQNNUP .LE. 0) GOTO 94

      OPEN(KANAL, FILE=FILENAME(:IDX(FILENAME)), 
     >  STATUS = 'OLD', ERR=99 )

C***  Loop over the lines ***************************************
      DO 
        READ(KANAL, *, ERR=95, END=98) NL, NU
        READ(KANAL, *, ERR=95, END=96) MNE, MTEMP, MP
C***    CHECK DIMENSIONS
	IF( MP     .GT. MPDIM    .OR.
     >      MTEMP  .GT. MTEMPDIM .OR.
     >      MNE    .GT. MNEDIM   ) GOTO 97

        READ (KANAL,*,ERR=95) (ALOG_NEL(J), J = 1, MNE),
     >                   ALOG_T0,     ALOG_T_INC,
     >                   ALOG_ALPHA0, ALOG_ALPHA_INC

        READ (KANAL,*,ERR=95) (( (SVCS(K,M,J), K = 1, MP),
     >                           M = 1, MTEMP), J = 1, MNE)
C***    Leave loop if the line was found, otherwise read next line

cc        write (0,'(A,4I3)') 'NL, LOW, NU, NUP =', NL, LOW, NU, NUP

        IF (NL .EQ. MAINQNLOW .AND. NU .EQ. MAINQNNUP) GOTO 1
      END DO

   1  CONTINUE

C***  Regular exit: line data found
      GOTO 101

C***  Error branches *******************************************

   93 WRITE (0,*) '*** WARNING: Stark data table search needs MAINQN ',
     >   'for ', LEVEL(LOW)
      WRITE (*,*) '*** WARNING: Stark data table search needs MAINQN ',
     >   'for ', LEVEL(LOW)
      LINPRO = 'Q-STARK '
      GOTO 100

   94 WRITE (0,*) '*** WARNING: Stark data table search needs MAINQN ',
     >   'for ', LEVEL(NUP)
      WRITE (*,*) '*** WARNING: Stark data table search needs MAINQN ',
     >   'for ', LEVEL(NUP)
      LINPRO = 'Q-STARK '
      GOTO 100

   95 WRITE (0,*) '*** READ ERROR on file VCSSB.DAT'
      STOP       ' *** FATAL ERROR IN Subr. STARKHEIIPREP'

   96 WRITE (0,*) '*** ERROR: E-O-F reached prematurely on VCSSB table'
      STOP       ' *** FATAL ERROR IN Subr. STARKHEIIPREP'


   97 WRITE (0,*) '*** ERROR: INSUFFICIENT DIMENSION FOR VCSSB TABLE'
      STOP       ' *** FATAL ERROR IN Subr. STARKHEIIPREP'

   98 WRITE (0,'(6A)') '*** WARNING: no Stark broadening table found ',
     >   'for ', LEVEL(LOW), ' - ', LEVEL(NUP),
     >   '  -> using L-STARK approximation instead'
      WRITE (*,*) '*** WARNING: no Stark broadening table found ',
     >   'for ', LEVEL(LOW), ' - ', LEVEL(NUP), 
     >   '  -> using L-STARK approximation instead'
      LINPRO = 'L-STARK '
      GOTO 100

   99 WRITE (0,*) '*** ERROR: No Stark data found for He II '
      WRITE (0,*) '*** Cannot open file: ', FILENAME(:IDX(FILENAME))
      WRITE (0,*) '*** You may want to define a different path by'
      WRITE (0,*) '*** using the FORMAL_CARDS option PATH_VCSSB = ...'
      WRITE (*,*) '*** ERROR: No Stark data found for He II '
      WRITE (*,*) '*** Cannot open file: ', FILENAME(:IDX(FILENAME))
      WRITE (*,*) '*** You may want to define a different path by'
      WRITE (*,*) '*** using the FORMAL_CARDS option PATH_VCSSB = ...'
      STOP       ' *** FATAL ERROR IN Subr. STARKHEIIPREP'


  100 CONTINUE
  101 CLOSE (KANAL)

      RETURN
      END
  
