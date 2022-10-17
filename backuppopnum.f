      SUBROUTINE BACKUPPOPNUM(ND, N, POPNUM, SCRATCH)
C***********************************************************************
C***  copies the POPNUM array into the SCRATCH array
C***********************************************************************
      
      IMPLICIT NONE

      INTEGER, INTENT(IN) :: ND, N
      REAL, DIMENSION(ND, N), INTENT(IN) :: POPNUM
      REAL, DIMENSION(ND, N), INTENT(OUT) :: SCRATCH
      
      INTEGER :: L, J
      
      SCRATCH = 0.
      DO J=1, N
        DO L=1, ND
          SCRATCH(L,J) = POPNUM(L,J)
        ENDDO
      ENDDO

      RETURN
      END
      