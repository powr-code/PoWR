      SUBROUTINE OPENEXMS(ICHANNEL, IADR, MAXADR, IFNAME, CSTATUS, IERR)
C***********************************************************************
C***  ROUTINE BY ANDREAS SANDER      6-Jun-2012 15:16:12
C***   extended version of OPENMS with additional STATUS flag
C***********************************************************************

      IMPLICIT NONE

      INTEGER :: IFNAME, ICHANNEL, IADR, MAXADR, IDUMMY, IERR
      REAL :: DUMMY
      CHARACTER(7) :: CSTATUS
      CHARACTER(8) :: FNAME     !IMPORTANT: only the first 7 characters can be used!

      INTEGER, EXTERNAL :: IDX

      IF (IDX(CSTATUS) == 0) THEN
        CSTATUS = 'UNKNOWN'
      ENDIF

      IF ((IFNAME /= 0) .AND. (IFNAME /= 1)) THEN
        WRITE(UNIT=FNAME, FMT='(A8)') IFNAME
      ELSE
        FNAME = '        '
      ENDIF

      CALL CMSSTORE (ICHANNEL, IADR, MAXADR, CSTATUS, FNAME, 
     >              DUMMY, IDUMMY, 'OPEN', IERR)

      RETURN
      END
