      SUBROUTINE AS2ST

      PARAMETER (IADRL = 128)
c      PARAMETER (MAXADR = 12800)
      PARAMETER (MAXADR = 25600)

      DIMENSION X(100000), J(100000)
      DIMENSION IADR(MAXADR)

      CHARACTER*8 NAME, LEER
      CHARACTER*104 MODHEAD
      CHARACTER LINE*132, ACTPAR*60, ACTPAR2*60, C*1

      LOGICAL BINT, OUTPUT

      OUTPUT = .FALSE.

      ICHANNEL = 20
      LEER = '        '

      OPEN (UNIT=4, FILE='ASCII_MODEL', STATUS='OLD')
      CALL STORAGE (ICHANNEL, IADR, MAXADR, DUMMY, DUMMY, X, 0, 
     >              'OPEN', 'CRAY', IERR)

   10 CONTINUE
        READ (4,'(A)') LINE
        CALL SARGV (LINE, 1, ACTPAR)
        IF (ACTPAR .EQ. 'MODHEAD') THEN
          READ (4,'(A)') LINE
          CALL STORAGE (ICHANNEL, IADR, MAXADR, 'MODHEAD ', DUMMY, 
     >                  LINE, 13, 'WRITE', 'CRAY', IERR)
          READ (4,'(A)') LINE
        ELSE
          DO I=1, 7
            IF (ACTPAR(I:I) .NE. ' ') GOTO 60
          ENDDO
   60     CONTINUE
          
          NAME(1: ) = ACTPAR(I:IDX(ACTPAR)) // LEER
          IF (NAME(1:1) .EQ. 'N' .OR. NAME(1:1) .EQ. 'M' .OR.
     >        NAME(1:1) .EQ. 'J') THEN
            BINT = .TRUE.
          ELSE
            BINT = .FALSE.
          ENDIF
          IF (NAME(1:3) .EQ. 'BLF') NAME = NAME(1:3) // ' ' // NAME(5:7)
          if (output) then
            write (0,*) 'name=',name
          endif

          NDAT = 0
   20     READ (4, '(A)', END=70) LINE
          CALL SARGV (LINE, 1, ACTPAR)
          IF (ACTPAR(1:1) .EQ. '*') GOTO 30
          CALL SARGC(LINE,N)
          DO I=1, N
            CALL SARGV (LINE, I, ACTPAR)
            NDAT = NDAT + 1
            IF (BINT) THEN
              IF (NAME .NE. 'MODHIST' .OR. NDAT .EQ. 1) THEN
                READ (ACTPAR, '(I25)') J(NDAT)
              ELSE
                DECODE (8,65,ACTPAR) J(NDAT)
   65           FORMAT(A8)
C!!!            READ (ACTPAR, '(I25)') J(NDAT)
              ENDIF
c       write (0,*) 'Einlesen Integer'
            ELSE
c       write (*,*) 'actpar=',actpar
              READ (ACTPAR, '(F25.0)', ERR=50) X(NDAT)
   50         CONTINUE
c       write (0,*) 'Einlesen Real'
            ENDIF
          ENDDO
          GOTO 20
   30     CONTINUE
          IF (BINT) THEN
            IF (NAME .EQ. 'MODHIST') THEN
c     >write (0,*) 'modhist:ndat=',ndat
c!!!              WRITE (0,*) 'DIMENSION OF MODHIST SET TO 4000'
              NDAT = 4000
            ENDIF
            CALL STORAGE (ICHANNEL, IADR, MAXADR, NAME, DUMMY, J, NDAT, 
     >                    'WRITE', 'CRAY', IERR)
       if (output) write (0,*) 'Speichern Integer'
          ELSE
            CALL STORAGE (ICHANNEL, IADR, MAXADR, NAME, DUMMY, X, NDAT, 
     >                    'WRITE', 'CRAY', IERR)
       if (output) write (0,*) 'Speichern real'
          ENDIF
        ENDIF

      GOTO 10

   70 CONTINUE


      WRITE (0,'(A)') 'AS2ST has worked correctly'

      CLOSE (4)
      CALL STORAGE (ICHANNEL, IADR, MAXADR, DUMMY, DUMMY, X, 0, 
     >              'CLOSE', 'CRAY', IERR)

      RETURN

      END
