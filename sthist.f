      SUBROUTINE STHIST (MODHIST, LAST, MAXHIST, GAMMAC, GAMMAL, GAMMAR,
     >                   GAMMAD, MODHEAD, JOBNUM, CORMAX, RTCMAX,
     >                   REDUCE, NSCHAR, BUNLU, 
     >                   DUNLU_LOC, DUNLU_INT, DUNLU_RMAX, DUNLU_TB, 
     >                   BXJLAPPNEW, BXJCAPPNEW, bBLOCKINVERSION,
     >                   XLAM_FINE_START, XLAM_FINE_END, bHydroDone, 
     >                   bTauUpdated)
C***********************************************************************
C***  UPDATING THE MODEL HISTORY FOR MAIN PROGRAM STEAL
C***********************************************************************

      IMPLICIT NONE

C***  MAXENTRY = MAXIMUM SPACE NEEDED BY THE NEXT REPEAT CYCLE.
C***   (MAXENTRY SHOULD BE GREATER THAN THE NUMBER OF NON-RUDIMENTAL LINES)
      INTEGER, PARAMETER :: MAXENTRY = 1000
      INTEGER :: JOBNUM, NSCHAR, MAXHIST, LAST, LENGTH, IDX,
     >           LBACK, LASTWRC, LASTCOM, LASTNEW, LNEW,
     >           IOLD, INEW, NRFOUND,
     >           ISTART, IEND, LASTWHATEVER, LASTCHAR, I, J



      REAL :: GAMMAC, GAMMAL, GAMMAR, GAMMAD,
     >        DELTAC, CORMAX, RTCMAX, REDUCE,
     >        CORLOG, 
     >        DUNLU_LOC, DUNLU_INT, DUNLU_RMAX, DUNLU_TB, 
     >        XLAM_FINE_START, XLAM_FINE_END

      CHARACTER(8*MAXHIST) :: MODHIST
      CHARACTER(100) :: MODHEAD
      CHARACTER(88) :: BUFFER88
      CHARACTER(32) :: BUFFER32
      CHARACTER(24) :: BUFFER24
      CHARACTER(16) :: BUFFER16
      CHARACTER(8) :: MAXCOR, BUFFER8
      CHARACTER(2) :: C2
      LOGICAL :: BUNLU, BXJCAPPNEW, BXJLAPPNEW, bHydroDone, 
     >           bBLOCKINVERSION, bTauUpdated
 
      MAXCOR='UNDEF.  '
      IF (JOBNUM > 1) THEN
         IF (CORMAX > 1.E-100) THEN
            CORLOG=LOG10(CORMAX)
            WRITE (UNIT=MAXCOR, FMT='(F8.4)') CORLOG
         ENDIF
      ENDIF
      
      IF (GAMMAC /= .0 .OR. GAMMAL /= .0) THEN
            BUFFER88 = ' '
            WRITE(UNIT=BUFFER88, FMT=2) JOBNUM,
     >      NINT(GAMMAC), NINT(GAMMAL), NINT(GAMMAR), NINT(GAMMAD), 
     >      MAXCOR
    2       FORMAT ('/',I7,'. STEAL', 
     >            '  GAMMA(CLRD)=',4(I5,1X),'  COR.=',A8,2X )
            LENGTH = IDX(BUFFER88)
            LENGTH = (LENGTH/8 + 1) * 8
            CALL ADDHISTENTRY
     >           (MODHIST,LAST,MAXHIST,LENGTH,BUFFER88(:LENGTH))
      ELSE
            WRITE(UNIT=BUFFER32, FMT=12) JOBNUM,MAXCOR
   12       FORMAT ('/',I7,'. STEAL   COR.=',A8 )
            CALL ADDHISTENTRY(MODHIST,LAST,MAXHIST,32,BUFFER32)
      ENDIF
 
C***  Fine Integration of Sourcefunction
      IF (BXJLAPPNEW .OR. BXJCAPPNEW) THEN
            C2 = '  '
            IF (BXJCAPPNEW) C2(1:1) = 'C'
            IF (BXJLAPPNEW) C2(2:2) = 'L'
            WRITE(UNIT=BUFFER24, FMT=22) 
     >        C2, XLAM_FINE_START, XLAM_FINE_END
   22       FORMAT ('JAPPNEW:', A2, F6.1, 1X, F6.1, 1X)
            CALL ADDHISTENTRY(MODHIST,LAST,MAXHIST,24,BUFFER24)
      ENDIF
            
      
C***  TEMPERATURE CORRECTIONS AND WEIGHT FOR FLUX CONSERVATION
C***   OR REMARK 'NOTEMP'

      IF (BUNLU) THEN
        WRITE(UNIT=BUFFER88, FMT=17) RTCMAX, 
     >             DUNLU_LOC, DUNLU_INT, DUNLU_RMAX, DUNLU_TB
   17   FORMAT ('TC=',F7.4,' UNLU=',4(F4.2,1X))
        LENGTH = IDX(BUFFER88)
        LENGTH = (LENGTH/8 + 1) * 8
        CALL ADDHISTENTRY
     >           (MODHIST,LAST,MAXHIST,LENGTH,BUFFER88(:LENGTH))
      ELSE
        WRITE(UNIT=BUFFER8, FMT=8)
    8   FORMAT ('NOTEMP')
        CALL ADDHISTENTRY(MODHIST,LAST,MAXHIST,8,BUFFER8)
      ENDIF

      IF (NSCHAR /= -1) THEN
        WRITE(UNIT=BUFFER8, FMT=5) NSCHAR
    5   FORMAT ('NSC=',I2)
        CALL ADDHISTENTRY(MODHIST,LAST,MAXHIST,8,BUFFER8)
      ENDIF

      IF (REDUCE /= 1.) THEN
        WRITE(UNIT=BUFFER16, FMT=4) REDUCE
    4   FORMAT ('REDUCE=',F3.2)
        CALL ADDHISTENTRY(MODHIST,LAST,MAXHIST,16,BUFFER16)
      ENDIF

      IF (bBLOCKINVERSION) THEN
        WRITE(UNIT=BUFFER8, FMT='(A8)') ' SPLIT  '
        CALL ADDHISTENTRY(MODHIST,LAST,MAXHIST,8,BUFFER8)
      ENDIF

      IF (bTauUpdated) THEN
        WRITE(UNIT=BUFFER8, FMT='(A8)') ' TAUFIX '
        CALL ADDHISTENTRY(MODHIST,LAST,MAXHIST,8,BUFFER8)
      ENDIF

      IF (bHydroDone) THEN
        WRITE(UNIT=BUFFER8, FMT='(A8)') ' HYDRO  '
        CALL ADDHISTENTRY(MODHIST,LAST,MAXHIST,8,BUFFER8)
      ENDIF
      
C***  SAFETY BRANCH PREVENTING MODEL HISTORY OVERFLOW:
      IF (LAST >= MAXHIST-MAXENTRY) THEN
            CALL PRIHIST (MODHIST,LAST,MODHEAD,JOBNUM)
C***        RESTORE HISTORY FROM LAST ENTRY OF REPEAT CYCLE BEFORE WRCONT-JOB
C***        FOR USE BY SUBR. DECNOT

            !These two are just needed if you want to read the job number
            ISTART = -1
            IEND = -1

            LASTCHAR = LAST*8

            LASTWRC=0
            LASTCOM=0
            LASTWHATEVER = -1
            lastcandloop: DO I=LASTCHAR-5,10,-1 !--------------------------
            !1. FIND LAST ENTRY OF 'WRCONT':
            !2. FIND LAST ENTRY OF 'COMO' before this last WRCONT:
            !Method> search for string plus valid job number
            !        and save "/"-position before the job number
              IF (((MODHIST(I:I+5) == "WRCONT") .AND. (LASTWRC == 0)) 
     >          .OR. 
     >            ((MODHIST(I:I+3) == "COMO") 
     >             .AND. (LASTCOM == 0) .AND. (LASTWRC /= 0) )) THEN
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
                          LASTWHATEVER = J
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

                  IF ((NRFOUND > 0) .AND. (LASTWHATEVER > 0)
     >              .AND. (ISTART > 0) .AND. (IEND > 0)) THEN
C** reinclude this if you need the job number
C                    READ(UNIT=MODHIST(ISTART:IEND), FMT='(I)') 
C     >                  LASTWHATEVERJOBNR
                    IF (MODHIST(I:I+5) == "WRCONT") THEN
                      LASTWRC = LASTWHATEVER
                    ELSEIF (MODHIST(I:I+3) == "COMO") THEN
                      LASTCOM = LASTWHATEVER
                    ENDIF
                    LASTWHATEVER = -1
                    IF ((LASTWRC > 0) .AND. (LASTCOM > 0)) THEN
                      !exit main loop if all jobs are found
                      EXIT lastcandloop
                    ENDIF
                  ENDIF

              ENDIF

            ENDDO lastcandloop !-------------------------------------------                 

            IF ((LASTWRC == 0) .OR. (LASTCOM == 0)) THEN
              !Search failed => end subroutine
              RETURN
            ENDIF


C***        3. SHIFT ENTRY TO THE BEGINNING OF NEW HISTORY:
            LASTNEW=LASTCHAR-LASTCOM+9
            MODHIST(9:LASTNEW) = MODHIST(LASTCOM:LASTCHAR)
            LAST=LASTNEW / 8
            IF (MOD(LASTNEW,8) /= 0) THEN
              !Failsafe if LASTNEW is not a multiple of 8 
              !(Usually this should not happen because LASTCHAR is always a multiple of 8)
              LAST = LAST + 1
            ENDIF
            WRITE(UNIT=BUFFER8, FMT='(A8)') LAST
            MODHIST(1:8) = BUFFER8 
C            WRITE(*,*) "new LAST: ", LAST
      ENDIF

      RETURN
      END
