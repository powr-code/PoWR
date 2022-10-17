      SUBROUTINE NEXTJOB (JOBNUM,JOBMAX,MODHIST,LAST,CORMAX,EPSILON,
     $                    NEWWRC,MOREJOBS,CONVERG,NOEXTRAP, NOCON, 
     >                    BPRICORR, COREX, BCOREX, NCOLIP, BAG, BGFIN, 
     >                    BAUTO_ABORT, FLUXEPS, ND, FLUXERR, IHSSTATUS)
C**********************************************************************
C***  CALLED FROM: STEAL
C***  THIS SUBROUTINE DECIDES UPON THE NEXT JOB TO BE EXECUTED.
C***  THE NAME OF THE NEXTJOB IS WRITTEN INTO "G1".
C***    - "G1" IS A JCL VARIABLE (COS) OR A FILE NAME (UNICOS) -
C***  POSSIBLE NEXTJOBS: WRCONT, EXTRAP, REPEAT, MODEL.
C***  INDEPENDENTLY, IT WRITES "G3"='MOREJOBS', IF THE
C***    MODEL IS NOT FINALLY CONVERGED, AND JOBMAX IS NOT EXCEEDED
C**********************************************************************

      IMPLICIT NONE

      INTEGER, INTENT(IN) :: JOBNUM, JOBMAX, LAST, NEWWRC, NCOLIP, 
     >                       IHSSTATUS, ND
      REAL, INTENT(IN) :: COREX, CORMAX, EPSILON, FLUXEPS
      REAL, DIMENSION(ND), INTENT(IN) :: FLUXERR
      LOGICAL, INTENT(INOUT) :: CONVERG, MOREJOBS
      LOGICAL, INTENT(IN) :: NOEXTRAP, NOCON, BPRICORR, BCOREX,
     >                       BAG, BGFIN, BAUTO_ABORT
      CHARACTER(8*LAST), INTENT(IN) :: MODHIST

      REAL, PARAMETER :: CORLIMIT = 1E-8 !used to avoid LOG if CORMAX = 0
      REAL :: LOGCORMAX

      INTEGER :: JOBDIFF, LASTEX, LASTCHAR, LASTWRC, LASTMOD, NGDIFF,
     >           MODDIFF, LASTWHATEVER, NRFOUND,
     >           I, J, ISTART, IEND, ILEN, IMAX

      CHARACTER(1) :: NAM1
      CHARACTER(6) :: NAM2
      CHARACTER(8) :: NEXTNAME
      CHARACTER(20) :: BUFFER20
      CHARACTER(40) :: JOBFMT
      CHARACTER(255) :: HISTENTRY

      LOGICAL :: bCONVCRIT

C***  File and channel handles (=KANAL)
      INTEGER, PARAMETER :: hOUT = 6        !write to wruniqX.out (stdout)
      INTEGER, PARAMETER :: hCPR = 0        !write to wruniqX.cpr (stderr)
      
      LASTCHAR = LAST * 8
      ISTART = -1
      IEND = -1

C***  FIND LASTWRC = JOBNUMBER OF LAST WRCONT JOB
C***  and FIND LASTEX = JOBNUMBER OF LAST EXTRAP JOB
C***  and FIND LASTMOD = JOBNUMBER OF LAST MODIFY JOB
      LASTWRC=0
      LASTEX=0
      LASTMOD=0
      lastcandloop: DO I=LASTCHAR-5,10,-1
        !First step: search for "WRCONT", "EXTRAP" or "MODIFY" string
        ! (do this only once per jobtype)
        IF (((MODHIST(I:I+5) == "WRCONT") .AND. (LASTWRC == 0)) .OR. 
     >      ((MODHIST(I:I+5) == "EXTRAP") .AND. (LASTEX == 0))  .OR. 
     >      ((MODHIST(I:I+5) == "MODIFY") .AND. (LASTMOD == 0))) THEN

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
            IF (MODHIST(I:I+5) == "WRCONT") THEN
              LASTWRC = LASTWHATEVER
            ELSEIF (MODHIST(I:I+5) == "EXTRAP") THEN
              LASTEX = LASTWHATEVER
            ELSEIF (MODHIST(I:I+5) == "MODIFY") THEN    
              LASTMOD = LASTWHATEVER
            ENDIF
            IF ((LASTWRC > 0) .AND. (LASTEX > 0) 
     >                        .AND. (LASTMOD > 0)) THEN
              !exit main loop if all jobs are found
              EXIT lastcandloop
            ENDIF
          ENDIF

        ENDIF

      ENDDO lastcandloop
   
      JOBDIFF=JOBNUM-LASTWRC
      IF (JOBDIFF < 0) JOBDIFF=JOBDIFF+1000
 
      NGDIFF=JOBNUM-LASTEX
      IF (NGDIFF < 0) NGDIFF = NGDIFF+1000

      MODDIFF=JOBNUM-LASTMOD
      IF (MODDIFF < 0) MODDIFF=MODDIFF+1000
 
C***  DEFINITION OF NEXTNAME (NAME OF THE NEXT JOB), AND MESSAGES *****
C----------------------------------------------------------------------

C***  GENERAL CASE: NEXTJOB = REPEAT-CYCLE
      NEXTNAME='REPEAT'

C***  FIRST JOB AFTER WRSTART: WRCONT
      IF (JOBNUM == 1) THEN
         CALL REMARK ('STEAL: EDDIES NOT YET EXISTING')
         PRINT *,'STEAL: EDDIES NOT YET EXISTING'
         NEXTNAME='WRCONT'
         GOTO 15
      ENDIF 

C***  EDDIES TOO OLD: WRCONT
      IF (.FALSE.) THEN
C***  Old Version
        IF (JOBDIFF >= 3*NEWWRC) THEN
          CALL REMARK ('STEAL: EDDIES TOO OLD')
          PRINT *, 'STEAL: EDDIES TOO OLD'
          NEXTNAME='WRCONT'
        ENDIF
      ELSE
C***  New Version for COLI
        IF (NCOLIP >= NEWWRC .OR. NCOLIP == 0 .OR. NCOLIP == -1) THEN
          CALL REMARK ('STEAL: EDDIES TOO OLD')
          PRINT *, 'STEAL: EDDIES TOO OLD'
          NEXTNAME='WRCONT'
        ENDIF
      ENDIF
       
C***  Print out convergence criteria to cpr file
      WRITE (hCPR,'(A,F10.6,A,F10.6,A)')
     >       'STEAL: current max. popnum corrections = ',
     >       CORMAX, ' (EPSILON=', EPSILON, ')'
      IF (FLUXEPS > 0.) THEN
        WRITE (hCPR,'(A,F10.6,A,F10.6,A)')
     >       'STEAL: current max. flux error = ',
     >       MAXVAL(FLUXERR), ' (FLUXEPS=', FLUXEPS, ')'
      ELSE
        WRITE (hCPR,'(A,F10.6,A,F10.6,A)')
     >       'STEAL: current max. flux error = ',
     >       MAXVAL(FLUXERR), ' (no FLUXEPS criterion defined)'
      ENDIF      
       
C***  CONVERGENCE TEST
      bCONVCRIT = CORMAX < EPSILON .AND. 
     >            BPRICORR .AND. 
     >            .NOT.  (BAG .AND. .NOT. BGFIN)

      IF (bCONVCRIT .AND. FLUXEPS .GT. .0) THEN
        IF (MAXVAL(FLUXERR) .GT. FLUXEPS) THEN
          bCONVCRIT = .FALSE.
          WRITE (0,'(A,F7.3,A,F7.3)')
     >       'STEAL: not converged: FLUX not conserved; FLUXERR=',
     >       MAXVAL(FLUXERR), ' > FLUXEPS=', FLUXEPS 
          WRITE (*,'(A,F7.3,A,F7.3)')
     >       'STEAL: not converged: FLUX not conserved; FLUXERR=',
     >       MAXVAL(FLUXERR), ' > FLUXEPS=', FLUXEPS 
        ENDIF
      ENDIF            

      IF (bCONVCRIT .AND. IHSSTATUS >= 0) THEN
        IF (IHSSTATUS == 0) THEN
          bCONVCRIT = .FALSE.
          WRITE (0,'(A)')
     >       'STEAL: not converged: HYDROSTATIC PART not consistent'
          WRITE (*,'(A)')
     >       'STEAL: not converged: HYDROSTATIC PART not consistent'
        ENDIF
      ENDIF            
            
      IF (bCONVCRIT) THEN
         CALL REMARK ('STEAL: REPEAT CYCLE IS CONVERGED')
         PRINT *,'STEAL: REPEAT CYCLE IS CONVERGED'
         NEXTNAME='WRCONT'
         DO J=LASTWRC+1, JOBNUM-1
           !new convergence check: only COLI and COMO since last WRCONT allowed
           CALL GETHISTENTRY(HISTENTRY,J,MODHIST,LAST)
           IMAX=LEN_TRIM(HISTENTRY)
           jnloop: DO I=2, IMAX 
             !jnloop needed to skip jobnumber of undefined length at the start
             IF (HISTENTRY(I:I) == ".") THEN
               !end of jobnumber is reached with first dot, jobname is next
               HISTENTRY = ADJUSTL(HISTENTRY(I+1:IMAX))
               IF ( (HISTENTRY(1:4) /= 'COLI') .AND. 
     >              (HISTENTRY(1:4) /= 'COMO') ) THEN
                 bCONVCRIT = .FALSE.
               ENDIF
               EXIT jnloop
             ENDIF
           ENDDO jnloop
         ENDDO

         IF (bCONVCRIT) THEN
            CALL REMARK ('STEAL: MODEL FINALLY CONVERGED!')
            PRINT *, 'STEAL: MODEL FINALLY CONVERGED!'
            NEXTNAME='MODEL'
         ENDIF
      ENDIF

      !prevent log(0) error
      IF (CORMAX <= CORLIMIT) THEN
        LOGCORMAX = ALOG10(CORLIMIT)
      ELSE 
        LOGCORMAX = ALOG10(CORMAX)
      ENDIF
      
C***  EXTRAPOLATION
C***  Old Version (NEWWRC=6)
      IF (NEXTNAME == 'REPEAT' .AND.
     $       .NOT. NOEXTRAP      .AND.
     $      NEWWRC >= 5        .AND.
     $     JOBDIFF >= 12       .AND.
     $      NGDIFF >= 12       .AND.
     >    (.NOT. BCOREX .OR. 
     >       (LOGCORMAX <= COREX))) THEN
        NEXTNAME='EXTRAP'
      ENDIF
 
C***  EXTRAPOLATION
C***  New Version (NEWWRC=1)
      IF (
     >    .NOT. NOEXTRAP      .AND.
     >    NEWWRC == 1       .AND. 
     >    NGDIFF >= 23      .AND.
     >    MODDIFF >= 23     .AND.
     >    JOBNUM > 20      .AND.
     >    (.NOT. BCOREX .OR. 
     >       (LOGCORMAX <= COREX))) THEN
        NEXTNAME='EXTRAP'
      ENDIF


   15 CONTINUE

C***  EXECUTION OF THE ABOVE DECISIONS  *******************************
C----------------------------------------------------------------------

C***  BRANCH FOR NEXTJOB = WRCONT 
      IF (NEXTNAME == 'WRCONT') THEN
         CALL REMARK ('STEAL: NEXTJOB=WRCONT')
         PRINT *,'STEAL: NEXTJOB=WRCONT'
         CALL JSYMSET ('G1','WRCONT')
      ENDIF
 
C***  BRANCH FOR NEXTJOB = REPEAT
      IF (NEXTNAME == 'REPEAT') THEN
         CALL REMARK ('STEAL: NEXTJOB=REPEAT')
         PRINT *,'STEAL: NEXTJOB=REPEAT'
         CALL JSYMSET ('G1','REPEAT')
      ENDIF

C***  BRANCH FOR NEXTJOB = EXTRAP
      IF (NEXTNAME == 'EXTRAP') THEN
         CALL REMARK ('STEAL: NEXTJOB=EXTRAP')
         PRINT *, 'STEAL: NEXTJOB=EXTRAP'
         CALL JSYMSET ('G1','EXTRAP')
      ENDIF

C***  BRANCH FOR NEXTJOB = MODEL (... FINALLY CONVERGED)
      IF (NEXTNAME == 'MODEL') THEN
         CALL REMARK ('MODEL FINALLY CONVERGED')
         PRINT *,' -------  MODEL FINALLY CONVERGED ]  ---------'
         CALL JSYMSET ('G1','MODEL')
         CONVERG=.TRUE.
      ENDIF
      
C***  MAX. NUMBER OF JOBS EXCEEDED? 
      IF (JOBNUM >= JOBMAX) THEN
         CALL REMARK ('STEAL: MAX. NUMBER OF JOBS EXCEEDED')
         PRINT *,'STEAL: MAX. NUMBER OF JOBS EXCEEDED'
         MOREJOBS=.FALSE.
         CALL JSYMSET ('G3','ENDJOB')
      ELSE
         MOREJOBS=.TRUE.
         CALL JSYMSET ('G3','MOREJOBS')
      ENDIF
            

C***  ONE ORE MORE DEPTH POINTS NOT CONVERGED?
      IF (NOCON) THEN
         CALL REMARK ('STEAL: ONE OR MORE DEPTH POINTS NOT CONVERGED')
         PRINT *,'STEAL: ONE OR MORE DEPTH POINTS NOT CONVERGED'
      ENDIF
      IF (NOCON .AND. BAUTO_ABORT) THEN
         MOREJOBS = .FALSE.
         CALL JSYMSET ('G3','ENDJOB')
      ENDIF
 
      RETURN
      END
