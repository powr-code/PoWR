      SUBROUTINE MSINFO
C*********************************************************************
C***  Tool for the inspection od mass storage files (e.g. MODEL file)
C***
C***     MAXADR increased from 25600 to 256000 (2000 Records)
C***     This is compatible with Version 3.1 of SUBR. CMSSTORE
C***                wrh 17-Apr-2008 11:40:11
C***     
C*******************************************************************

      IMPLICIT NONE

      INTEGER, PARAMETER :: MAXADR = 2000 * 128
      INTEGER, DIMENSION(MAXADR) :: IADR

      CHARACTER(8) :: LINE, ACTION, CNAME, CNAME2, MODE, CSTAT
      CHARACTER(12) :: FNINPUT
      CHARACTER(40) :: CFORMAT      !NOTE : CFORMAT is longer than 8 characters. 

      INTEGER :: ICHANNEL, IERR, NPAR, IPAR
      REAL :: DUMMY

      LOGICAL :: bFIRST, bINPUTEXIST, bFILEEXIST, bACTION, bINFO

C***  File and channel handles (=KANAL)
      INTEGER, PARAMETER :: hINPUT = 10        !input parameter file msinfo.input

      INTEGER, EXTERNAL :: IDX

C***  Link data to identify program version
      CHARACTER(30) :: LINK_DATE
      CHARACTER(10) :: LINK_USER
      CHARACTER(60) :: LINK_HOST
      COMMON / COM_LINKINFO / LINK_DATE, LINK_USER, LINK_HOST

      ICHANNEL = 33
      MODE = 'CRAY'
      CSTAT = 'AUTO'
      FNINPUT = 'msinfo.input'

      INQUIRE (FILE=FNINPUT, EXIST=bINPUTEXIST)

      IF (.NOT. bINPUTEXIST) THEN

        NPAR = COMMAND_ARGUMENT_COUNT()
        IF (NPAR < 1) THEN
          bINFO = .TRUE.
        ELSE
          CALL GET_COMMAND_ARGUMENT(1, CNAME2)
          !check and copy file
          INQUIRE (FILE=CNAME2, EXIST=bFILEEXIST)
          IF (.NOT. bFILEEXIST) THEN
            WRITE (0,*) 'File ', CNAME2, ' could not be found'
            STOP 'Error in MSINFO'
          ELSE
            CALL STORAGE (ICHANNEL, IADR, MAXADR, CSTAT, CNAME2, 
     >                    DUMMY, 1, 'OPEN', MODE, IERR)
          ENDIF
          WRITE (*,*) '*'
          WRITE (*,*) '* MSINFO: Information ueber MS-File ', CNAME2
          WRITE (*,*) '*'
          IF (NPAR > 1) THEN
            DO IPAR=2, NPAR              
              CALL GET_COMMAND_ARGUMENT(IPAR, ACTION)
              bACTION = .FALSE.
              SELECTCASE(ACTION)
                CASE ('INFO', 'INFO-L')
                  CNAME = ' '
                  CFORMAT = ' '
                  bACTION = .TRUE.
                CASE ('INFO-D')
                  IF (NPAR > IPAR) THEN
                    CALL GET_COMMAND_ARGUMENT(IPAR+1, CNAME)
                    IF (CNAME(1:4) == 'INFO') THEN
                      GOTO 120
                    ELSE
                      bACTION = .TRUE.
                    ENDIF
                  ELSE
                    GOTO 120
                  ENDIF
                  IF (NPAR > IPAR+1) THEN      
                    CALL GET_COMMAND_ARGUMENT(IPAR+2, CFORMAT)
                    IF (CFORMAT(1:4) == 'INFO') THEN 
                      !This is no format string but already the next action
                      CFORMAT = 'AUTO'
                    ENDIF
                  ELSE
                    CFORMAT = 'AUTO'
                  ENDIF
              ENDSELECT
              IF (bACTION) THEN
                CALL STORAGE (ICHANNEL, IADR, MAXADR, CNAME, CFORMAT, 
     >                        DUMMY, 1, ACTION, MODE, IERR)
                WRITE (*,*) '*'
                WRITE (*,*) '*'
              ENDIF
            ENDDO
          ENDIF
        ENDIF
      ENDIF

      IF (bINFO) THEN
        WRITE (*,*) 'Bedienungsanleitung f??r msinfo:'
        WRITE (*,*) '================================='
        WRITE (*,*) ''
C        WRITE (*,*) 'Das Kommand benoetigt den File-Namen und das Kommando.'
C        WRITE (*,*) 'An beliebiger Stelle hinter dem File-Namen kann das keyword red stehen.'
C        WRITE (*,*) 'Der Output wird dann in ein gesondertes File geschrieben.'
        WRITE (*,*) 'Syntax:'
        WRITE (*,*) 'msinfo <file> <Befehl1> <Befehl2> <...>'
        WRITE (*,*) ''
        WRITE (*,*) 'Folgende Kommandos sind m??glich:'
        WRITE (*,*) 'INFO      :  Kurzinfo (default)'
        WRITE (*,*) 'INFO-L    :  Langinfo'
        WRITE (*,*) 'INFO-D    :  Ausgabe des Inhalts einer Variablen'
        WRITE (*,*) '               (Zusaetzliche Parameter werden benoetigt)'
        WRITE (*,*) 'Beispiel  :  msinfo <file> INFO-D <Name der Variablen> <Ausgabeformat>'
        WRITE (*,*) ''
        WRITE (*,*) 'msinfo MODEL INFO INFO-D ND (I6)'
        GOTO 29
      ENDIF

      IF (bINPUTEXIST) THEN
        OPEN (UNIT=hINPUT, FILE=FNINPUT, ERR=100)
        READ (hINPUT, '(A)', END=110) LINE

        WRITE (*,*) '*'
        WRITE (*,*) '* MSINFO: Information ueber MS-File ', LINE
        WRITE (*,*) '*'

        CNAME2 = '        '

        CALL STORAGE (ICHANNEL, IADR, MAXADR, CSTAT, CNAME2, DUMMY, 1, 
     >              'OPEN', MODE, IERR)

        bFIRST = .TRUE.
C***    Loop over Input-Parameters
   10   CONTINUE
          READ (hINPUT, '(A)', END=20) ACTION
   22     CONTINUE
          BFIRST = .FALSE.
          IF (ACTION /= 'INFO-D') THEN
            CNAME = ' '
            CFORMAT = ' '
          ELSE
            READ (hINPUT, '(A)', END=120) CNAME
            READ (hINPUT, '(A)', END=120) CFORMAT
          ENDIF

          CALL STORAGE (ICHANNEL, IADR, MAXADR, CNAME, CFORMAT, DUMMY, 1, 
     >                ACTION, MODE, IERR)
          WRITE (*,*) '*'
          WRITE (*,*) '*'

        GOTO 10

   20   CONTINUE
        IF (BFIRST) THEN
          ACTION = 'INFO'
          GOTO 22
        ENDIF      

        CLOSE(hINPUT)
      ENDIF
   29 WRITE (*,'(6A)') ' ### MSINFO from ', LINK_DATE, ' created by '
     >                 , LINK_USER(:IDX(LINK_USER))
     >     ,' at host ', LINK_HOST(:IDX(LINK_HOST))
      WRITE (*,*) '### MSINFO has worked correctly'

      RETURN

  100 WRITE (0,*) 'Error when Input-File msinfo.input was opened'
      STOP 'Error in MSINFO'

  110 WRITE (0,*) 'Input-File was too short'
      STOP 'Error in MSINFO'

  120 WRITE (0,*) 'Parameters for INFO-D where missing'
      STOP 'Error in MSINFO'

      END
