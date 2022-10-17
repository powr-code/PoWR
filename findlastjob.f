      SUBROUTINE FINDLASTJOB (JOBNAME,JOBNUM,MODHIST,MAXHIST,
     >                        HISTINDEX, OFFSET)
C**********************************************************************
C***  SEARCH ROUTINE TO FIND THE LAST JOBNUMBER OF THE GIVEN JOB
C***   HISTINDEX returns the first index of the entry that
C***       has been found by this search
C***   OFFSET defines the start entry, use -1 for default (last char)
C**********************************************************************

      IMPLICIT NONE

      INTEGER, INTENT(IN) :: OFFSET, MAXHIST
      INTEGER, INTENT(OUT) :: JOBNUM, HISTINDEX
      CHARACTER(8*MAXHIST), INTENT(IN) :: MODHIST
      CHARACTER(8) :: JOBNAME

      INTEGER :: LAST, LASTCHAR, LASTWHATEVER, NRFOUND, STARTINDEX,
     >           I, J, ISTART, IEND, ILEN, JNLEN, BUFFERINT

      CHARACTER(1) :: NAM1
      CHARACTER(6) :: NAM2
      CHARACTER(8) :: NEXTNAME, BUFFER8
      CHARACTER(20) :: BUFFER20
      CHARACTER(40) :: JOBFMT

      IF (LAST < 0) THEN
        BUFFER8 = MODHIST(1:8)
        READ(UNIT=BUFFER8, FMT='(A8)') BUFFERINT
        LAST = BUFFERINT        
      ENDIF
      LASTCHAR = LAST * 8

      IF ((OFFSET < 8) .OR. (OFFSET > LASTCHAR)) THEN
        !for useless offsets search from last char
        STARTINDEX = LASTCHAR
      ELSE
        STARTINDEX = OFFSET
      ENDIF

      ISTART = -1
      IEND = -1

      HISTINDEX = -1

      JNLEN = LEN_TRIM(JOBNAME)

C***  FIND LASTWHATEVER = JOBNUMBER OF LAST <JOBNAME>
      LASTWHATEVER=0
      lastcandloop: DO I=STARTINDEX-5,10,-1
        !First step: search for JOBNAME string
        IF (MODHIST(I:I+JNLEN-1) == TRIM(ADJUSTL(JOBNAME))) THEN
          J = I
          NRFOUND = -1
          !Second step: Try to find Jobnumber
          findjobnr: DO !- - - - - - -
            J = J - 1
            SELECTCASE (MODHIST(J:J))
              CASE (" ")
                !Blank => Nothing happens
                CYCLE findjobnr
              CASE ("0", "1":"9")
                !Nothing happens
                NRFOUND = NRFOUND + 1
                CYCLE findjobnr
              CASE (".")
                IF (NRFOUND < 0) THEN
                  NRFOUND = 0
                ELSEIF (NRFOUND > 0) THEN
                  !Second Dot found => This is not a job number
                  NRFOUND = -1
                  EXIT findjobnr
                ENDIF
                IEND = J - 1
                CYCLE findjobnr
              CASE ("/")
                IF (NRFOUND > 0) THEN
                  !success: number found and this is really a job number
                  ISTART = J + 1
                ELSE
                  !Slash found at wrong positon => This is not a job number
                  NRFOUND = -1
                ENDIF
                EXIT findjobnr
              CASE DEFAULT
                NRFOUND = -1
                EXIT findjobnr
            ENDSELECT
        
          ENDDO findjobnr !- - - - - -

          IF ((NRFOUND > 0) .AND. (ISTART > 0) .AND. (IEND > 0)) THEN
            ILEN = IEND - ISTART + 1
            WRITE(UNIT=JOBFMT, FMT='(I20)') ILEN
            BUFFER20 = '(I' // TRIM(ADJUSTL(JOBFMT))
            JOBFMT = TRIM(BUFFER20) // ')'
            READ(UNIT=MODHIST(ISTART:IEND), FMT=JOBFMT) LASTWHATEVER
            IF (MODHIST(I:I+JNLEN-1) == TRIM(ADJUSTL(JOBNAME))) THEN
              HISTINDEX = ISTART - 1
              JOBNUM = LASTWHATEVER
              EXIT lastcandloop
            ELSE
              !Something went really wrong here
            ENDIF
          ENDIF

        ENDIF

      ENDDO lastcandloop

      IF (NRFOUND <= 0) THEN
        JOBNUM = -1
        HISTINDEX = -1
      ENDIF

      RETURN
      END

