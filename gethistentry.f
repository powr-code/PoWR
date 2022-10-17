      SUBROUTINE GETHISTENTRY(HISTENTRY,JOBNR,MODHIST,MAXHIST)
C***********************************************************************
C***  RETURNS AN ENTRY FOR A GIVEN JOB NUMER FROM THE MODEL HISTORY 
C       CHARACTER ARRAY
C***********************************************************************

      IMPLICIT NONE

      INTEGER, INTENT(IN) :: JOBNR, MAXHIST
      CHARACTER(8*MAXHIST), INTENT(IN) :: MODHIST
      CHARACTER(*), INTENT(OUT) :: HISTENTRY

      INTEGER :: LASTCHAR, LAST, BUFFERINT, ISTART, IEND, ILEN,
     >           NRFOUND, IFROM, ITO, I
      CHARACTER(16) :: BUFFER8, JOBNRSTRING, JOBCANDSTR

      BUFFER8 = MODHIST(1:8)
      READ(UNIT=BUFFER8, FMT='(A8)') BUFFERINT
      LAST = BUFFERINT
      LASTCHAR = LAST * 8

      NRFOUND = -1
      ISTART = 0            !ISTART and IEND are used for reading out job numbers
      IEND = 0
      IFROM = 0             !IFROM and ITO are used to read out the history entry

      WRITE(UNIT=JOBNRSTRING, FMT='(I16)') JOBNR
      JOBNRSTRING = ADJUSTL(JOBNRSTRING)

      findjobnr: DO I=9, LASTCHAR
        SELECTCASE (MODHIST(I:I))
          CASE (" ")
            !Blank => Nothing happens
          CASE ("0", "1":"9", "A")
            !Nothing happens
            IF (NRFOUND >= 0) THEN
              NRFOUND = NRFOUND + 1
            ENDIF
          CASE (".")
            IF (NRFOUND > 0) THEN
              !success: number found and this is really a job number
              IF (IFROM == 0) THEN
                IEND = I - 1
                JOBCANDSTR = ADJUSTL(MODHIST(ISTART:IEND))
                IF (TRIM(JOBCANDSTR) == TRIM(JOBNRSTRING)) THEN
                  !This is the job we are searching for
                  IFROM = ISTART - 1
                ENDIF
                NRFOUND = -1
                ISTART = 0
                IEND = 0
              ELSE
                !This is already the jobnumber after the searched job
                ITO = ISTART - 2
                EXIT findjobnr
              ENDIF
            ELSE
              !Dot found at wrong positon => This is not a job number
              NRFOUND = -1
            ENDIF
          CASE ("/")
            IF (NRFOUND < 0) THEN
              ISTART = I + 1
              NRFOUND = 0
            ELSEIF (NRFOUND > 0) THEN
              !Second Slash found => This is not a job number
              ISTART = 0
              IEND = 0
              NRFOUND = -1
            ENDIF
          CASE DEFAULT
            NRFOUND = -1
        ENDSELECT        
        IF ((I == LASTCHAR) .AND. (IFROM /= 0)) THEN
          !Last job was the one that was searched for, end is
          !  determined by end of history instead of next entry
          ITO = LASTCHAR 
        ENDIF
      ENDDO findjobnr

      IF ((IFROM > 0) .AND. (ITO > 0)) THEN
        HISTENTRY = MODHIST(IFROM:ITO)
      ELSE
        HISTENTRY = 'FAILED'
      ENDIF

      RETURN
      END
