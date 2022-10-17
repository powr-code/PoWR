C***  MAIN PROGRAM NJN *********************************************
      SUBROUTINE NJN
C***********************************************************************
C***  THIS PROGRAM RESETS THE JOBNUMBER TO 100
C***   required to continue calculations after JOBMAX has been reached
C***   starting with 100 ensures that no start approximations or early
C***   calculation methods are used
CCCC  abarnisk: Remove IADR and MAXADR, these Parameters are 
CCCC  no longer needed.
C***********************************************************************
 
      IMPLICIT NONE

      INTEGER, PARAMETER :: MAXHIST =  4000 
      CHARACTER, DIMENSION(8*MAXHIST) :: MODHIST

      INTEGER :: IDUMMY, IERR, JOBNUM, JOB2, LAST,
     >           IADR, MAXADR, LASTUPD
      CHARACTER(8) :: BUFFER8
      CHARACTER(32) :: BUFFER32
      CHARACTER(40) :: BUFFER40

      !File and channel handles (=KANAL)
      INTEGER, PARAMETER :: hCPR = 0        !write to wruniqX.cpr (stderr)
      INTEGER, PARAMETER :: hMODEL = 3      !MODEL file (fort.3)
      INTEGER, PARAMETER :: hEDDI = 17      !EDDI file (fort.17)
      INTEGER, PARAMETER :: hOUT = 6        !write to wruniqX.out (stdout)

C***  READING OF THE MODEL FILE ----------------------------------------
      CALL OPENMS (hMODEL, IDUMMY, IDUMMY, 1, IERR)

      CALL READMS(hMODEL,LAST,1,'MODHIST ' , IERR)
      WRITE (hCPR,*) 'MODHIST LESEN', LAST
      IF (LAST+5 >= MAXHIST) THEN
         WRITE(hCPR,*)'MODHIST DIMENSION INSUFFICIENT'
         STOP 'ERROR IN SUBR. NJN'
      ENDIF
      CALL READMS (hMODEL,MODHIST,LAST,'MODHIST ', IERR)

      CALL READMS (hMODEL, JOBNUM ,  1, 'JOBNUM  ' , IERR)
      JOB2 = JOBNUM + 1
      WRITE (hOUT,'(A,I7)') 'Old JOBNUM was   : ', JOBNUM
      JOBNUM = 100

      WRITE (UNIT=BUFFER40,FMT=100) JOB2, JOBNUM
  100 FORMAT  ('/', I7, '. Jobnumber set to ', I7, '  ')
      CALL ADDHISTENTRY(MODHIST, -1, MAXHIST, 40, BUFFER40)
      LAST = LAST + 5
      WRITE (hCPR,*) 'MODHIST SCHREIBEN', LAST

      CALL WRITMS (hMODEL,MODHIST,MAXHIST,'MODHIST ',-1, IDUMMY, IERR)
      WRITE (hOUT,'(A,I7)')
     >  'New JOBNUM is now : ', JOBNUM

      CALL WRITMS (hMODEL, JOBNUM ,  1,  'JOBNUM   ', -1, IDUMMY, IERR)
      CALL CLOSMS (hMODEL, IERR)

      CALL OPENMS (hEDDI, IADR, MAXADR, 1, IERR)
      LASTUPD = JOBNUM - 1
      WRITE (hCPR,*) 'LASTUPD SCHREIBEN'
      CALL WRITMS (hEDDI, LASTUPD ,  1,  'LASTUPD ', -1, IDUMMY, IERR)
      CALL CLOSMS (hEDDI, IERR)
      CALL JSYMSET ('G0','0')

      STOP 'O.K.'
      END
