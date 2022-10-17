      SUBROUTINE LOAD_FF (ND, NDDIM, MODHEAD, JOBNUM, bFFASSET,
     >                    FF_INFO, IFF_DK, IFF_N_MS, IFF_MAX_MS)
C***  Read stuff for Fine-Integration in SETXJFINE

      IMPLICIT NONE
     
      INTEGER, PARAMETER :: TINYINT = SELECTED_INT_KIND(2)
      
      INTEGER, INTENT(IN) :: ND, NDDIM, IFF_MAX_MS, JOBNUM
      INTEGER, INTENT(OUT) :: IFF_N_MS
      LOGICAL, INTENT(OUT) :: bFFASSET

      INTEGER, DIMENSION(2*NDDIM+8) :: IADR16
      
      REAL, DIMENSION(10) :: FF_INFO
      INTEGER(KIND=TINYINT), DIMENSION(IFF_MAX_MS) :: IFF_DK   

      CHARACTER(8) :: CKIND, CNAME
      CHARACTER(100) :: MODHEAD, MODINFO
      
      INTEGER :: IERR, LASTUPD
      
      !File and channel handles (=KANAL)
      INTEGER, PARAMETER :: hOUT = 6        !write to wruniqX.out (stdout)
      INTEGER, PARAMETER :: hCPR = 0        !write to wruniqX.cpr (stderr)
      INTEGER, PARAMETER :: hFF = 27        !file handle for FFASSET
     
      IFF_N_MS = 0.     !means: better read nothing if nothing is found
     
      CALL OPENEXMS (hFF, IADR16, 2*ND+8, 'FFASSET', 'UNKNOWN', IERR)

      CALL READMS (hFF, MODINFO, 13, 'MODEL   ', IERR)
      IF (IERR < 0  .AND.  IERR /= -10) THEN
         CALL REMARK ('LOAD_FF: ERROR WHEN READING FFASSET')
         STOP 'ERROR'
      ENDIF

      IF (IERR == -10) THEN
        WRITE(hOUT,*) 'STEAL> FFASSET NOT YET EXISTING'
        bFFASSET = .FALSE.
        RETURN
      ENDIF

      IF (MODINFO /= MODHEAD) THEN
         PRINT *, 'STEAL> FFASSET BELONGS TO A DIFFERENT MODEL'
        bFFASSET = .FALSE.
        RETURN
      ENDIF

C***  CHECK JOBNUM OF LAST UPDATE
      CALL READMS (hFF, LASTUPD, 1, 'LASTUPD ', IERR)
      IF (LASTUPD >= JOBNUM) THEN
        WRITE(hOUT,*) 'STEAL> FFASSET YOUNGER THAN PRESENT MODEL'
        bFFASSET = .FALSE.
        RETURN
      ENDIF
      
      
      CALL READMS   (hFF,FF_INFO,     10, 'FF_INFO ', IERR)
      IF (IERR == -10) THEN
        WRITE(hOUT,*) 'STEAL> FFASSET META INFORMATION NOT FOUND'
        bFFASSET = .FALSE.
        RETURN
      ENDIF

      bFFASSET = .TRUE.
      
      CALL LENGTHMS(hFF, IFF_N_MS, 'IFF_DK  ', IERR)
      IF (IFF_N_MS > IFF_MAX_MS) THEN
        WRITE (hCPR,*) 
     >      'Number of fine frequencies lower than in MODEL file',
     >      'Try version with higher dim. (e.g. vd20 instead of vd50)'
        WRITE (hCPR,'(A,I8)') 'Needed IFF_N_MS = ', IFF_N_MS
        WRITE (hCPR,'(A,I8)') 'Present IFF_MAX_MS = ', IFF_MAX_MS
        STOP 'ERROR in Subroutine STEAL>LOAD_FF'          
      ENDIF

      CALL READMS (hFF,IFF_DK, IFF_N_MS,  'IFF_DK  ', IERR)
      
      RETURN
      END
