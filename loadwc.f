      SUBROUTINE LOADWC(IFF_WCHARM, IFF_MAX_MS, L)
C*******************************************************************
C***  Read the WCHARM factors from FFASSET file (Channel 27)
C***    at first Iteration of each depth point.
C***  Called by SETXJFINE
C*******************************************************************

      IMPLICIT NONE

      INTEGER, PARAMETER :: TINYINT = SELECTED_INT_KIND(2)
      INTEGER, INTENT(IN) :: IFF_MAX_MS

      INTEGER (KIND=TINYINT), DIMENSION(IFF_MAX_MS) :: IFF_WCHARM

      INTEGER :: L, IERR
      CHARACTER(8) :: CNAME

      INTEGER, PARAMETER :: hFF = 27        !file handle for FFASSET
      
      WRITE (UNIT=CNAME, FMT='(A5, I3)') 'IFFWC', L      
      CALL READMS (hFF, IFF_WCHARM, IFF_MAX_MS, CNAME, IERR)

      RETURN
      END
