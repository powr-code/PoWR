      SUBROUTINE CLOPENE (KANAL, MODHEAD, JOBNUM, NCOLIP, 
     >                    NDEDDIA, BCLEERR)
C**********************************************************************
C***  OPEN THE EDDI-FILE 
C***  (used for the calculation of the radiation field by COLI with 
C***   Variable Eddington Factors)
C***  AND CHECK WHETHER THESE DATA EXIST AND BELONG TO THE PRESENT MODEL
C***  CALLED FROM: COLI
C**********************************************************************
      IMPLICIT NONE
      
      INTEGER, INTENT(IN) :: JOBNUM, NCOLIP, NDEDDIA, KANAL

      CHARACTER(100), INTENT(IN) :: MODHEAD
      CHARACTER(100) :: MODEL
      LOGICAL, INTENT(INOUT) :: BCLEERR
      
      INTEGER :: NDEDDIAtest, IERR, IERRND, IDUMMY, LASTUPD

      BCLEERR = .FALSE.

      CALL OPENMS (KANAL,IDUMMY, IDUMMY, 1, IERR)

      IERR=1
      CALL READMS (KANAL,MODEL,13,'MODEL   ',IERR)

C***  CASE: FILE IS PRESENT BUT DOES NOT CONTAIN THE "MODEL" RECORD 
C***        -> ERROR STOP
      IF (IERR .LT. 0 .AND. IERR. NE. -10) THEN
            WRITE (0,*) 'ERROR WHEN READING MODHEAD FROM FILE'
            STOP 'ERROR IN CLOPENE'
            ENDIF

      IF (IERR /= -10) THEN
C***      Read length for NDEDDIA storage array in current EDDI file
          CALL READMS (KANAL,NDEDDIAtest,1,'NDEDDIA ',IERRND)
          IF (IERRND == -10) THEN
            NDEDDIAtest = 0
          ENDIF
      ENDIF
            
C***  CASE: FILE IS NOT PRESENT OR FROM WRONG MODEL
C***        OR COLI++ IS FORCED VIA NCOLIP = -1 (e.g. after hydro step)
C***        -> SET FLAG, WRITE MODHEAD
      IF (IERR == -10 .OR. MODEL /= MODHEAD  
     >            .OR. NCOLIP == -1 .OR. NDEDDIA /= NDEDDIAtest) THEN
      
         IF (IERR .EQ. -10) WRITE (*,*) 
     >      ' COLI: CLOPENE> File EDDI not present -> COLI++'
         IF (MODEL .NE. MODHEAD) WRITE (*,*) 
     >      ' COLI: CLOPENE> File EDDI from wrong model -> COLI++'
         IF (NCOLIP .EQ. -1) WRITE (*,*) 
     >      ' COLI: CLOPENE> Renewed EDDI file enforced -> COLI++'
          IF (NDEDDIA /= NDEDDIAtest) THEN
C***  Length check for NDEDDIA storage array
C           if length does not match the EDDI file was created with a different
C           COLI version and cannot be read => create new one in this case  
            WRITE (*,*) 
     >      ' COLI: CLOPENE> EDDI created with different COLI version '
     >         // '-> COLI++'
            BCLEERR = .TRUE.
          ENDIF

         BCLEERR = .TRUE.
         !Replace with new, empty file and start with MODEL information
         CALL CLOSMS (KANAL, IERR)
         CALL OPENEXMS (KANAL, IDUMMY, IDUMMY, 1, 'REPLACE', IERR)
         CALL WRITMS (KANAL,MODHEAD,13,'MODEL   ',-1, IDUMMY, IERR)
         CALL WRITMS (KANAL,JOBNUM,1,  'LASTUPD ',-1, IDUMMY, IERR)
         CALL WRITMS (KANAL,NDEDDIA,1, 'NDEDDIA ',-1, IDUMMY, IERR)

C*** CASE: FILE IS FROM CORRECT MODEL -> CHECK AND UPDATE JOBNUM
      ELSE
          CALL READMS (KANAL,LASTUPD,1, 'LASTUPD ', IERR)

          IF (.NOT. BCLEERR .AND. LASTUPD .GE. JOBNUM) THEN 
            BCLEERR = .TRUE.
            PRINT *, 
     >      ' COLI:CLOPENE> File EDDI younger than present model -> '
     >        // 'COLI++'
          ENDIF

          IF (.NOT. BCLEERR .AND. LASTUPD .LT. JOBNUM-6) THEN
C         NOTE: It seems that these lines are never called because
C               LASTUPD is always written, even if no new EDDIs are stored (normal COLI)
            BCLEERR = .TRUE.
            PRINT *,
     >      ' COLI: CLOPENE> File EDDI older than 6 jobs -> '
     >        // 'COLI++'
          ENDIF

                     
           CALL WRITMS (KANAL,JOBNUM,1,'LASTUPD ',-1, IDUMMY, IERR)

      ENDIF



      RETURN
      END
