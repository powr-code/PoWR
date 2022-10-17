      SUBROUTINE STORE_FF(MODHEAD, JOBNUM, ND, NDDIM,
     >                    FF_INFO, IFF_DK, IFF_MAX, IFF_WCHARM)
C***********************************************************************
C*** Creates the FFASSET file and stores the fine-frequency weights
C***  only called by COLI main program
C***********************************************************************
      IMPLICIT NONE

      INTEGER, PARAMETER :: TINYINT = SELECTED_INT_KIND(2)
      
      INTEGER, INTENT(IN) :: ND, NDDIM, IFF_MAX, JOBNUM
      INTEGER, DIMENSION(2*NDDIM+8) :: IADR16
      
      REAL, DIMENSION(10) :: FF_INFO
      INTEGER(KIND=TINYINT), DIMENSION(IFF_MAX) :: IFF_DK   
      INTEGER(KIND=TINYINT), DIMENSION(IFF_MAX,ND) :: IFF_WCHARM

      CHARACTER(8) :: CKIND, CNAME
      CHARACTER(100) :: MODHEAD
      
      INTEGER :: L, IERR, IDUMMY
      
      INTEGER, PARAMETER :: hFF = 27
      
      !Stuff is only used in next STEAL, so always create new FFASSET file
      CALL OPENEXMS (hFF, IADR16, 2*ND+8, 'FFASSET', 'REPLACE', IERR)

      CALL WRITMS (hFF, MODHEAD, 13, 'MODEL   ', -1, IDUMMY, IERR)
      CALL WRITMS (hFF, JOBNUM ,  1, 'LASTUPD ', -1, IDUMMY, IERR)

      CALL WRITMS (hFF, FF_INFO, 10, 'FF_INFO ', -1, IDUMMY, IERR)
      
      !Storage of the one-byte-integer arrays needs KIND-Parameter (instead default -1)
      WRITE(UNIT=CKIND, FMT='("I",I1)') TINYINT
      CALL WRITMS (hFF, IFF_DK, IFF_MAX, 'IFF_DK  ', 
     >                                        CKIND, IDUMMY, IERR)

      DO L=1, ND
        WRITE (UNIT=CNAME, FMT='(A5, I3)') 'IFFWC', L
        CALL WRITMS (hFF, IFF_WCHARM(1,L), IFF_MAX, CNAME,
     >                                        CKIND, IDUMMY, IERR)
      ENDDO

      CALL CLOSMS (hFF)

      RETURN
      END
