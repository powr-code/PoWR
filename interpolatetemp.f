      SUBROUTINE INTERPOLATETEMP(Tnew, Told, Rnew, Rold, ND)
C***********************************************************************
C***  Interpolation of the electron temperature on a new grid
C***  
C***  called by ENSURETAUMAX
C***********************************************************************

      IMPLICIT NONE
      INCLUDE 'interfacebib.inc'

      INTEGER, INTENT(IN) :: ND
      REAL, DIMENSION(ND), INTENT(INOUT) :: Tnew
      REAL, DIMENSION(ND), INTENT(IN) :: Rnew, Rold, Told

      REAL :: RnewLOGL
      INTEGER :: L
      LOGICAL, PARAMETER :: bLOG = .TRUE.
      
C***  Heap variables      
      REAL, DIMENSION(ND) :: RoldLOG      


      IF (bLOG) THEN
C***    Interpolation over log radius  
        DO L=1, ND
          RoldLOG(L) = LOG10(Rold(L))
        ENDDO                 
        DO L=1, ND
          RnewLOGL = LOG10(Rnew(L))
          IF (RnewLOGL > RoldLOG(1)) THEN
            Tnew(L) = Told(1)
          ELSEIF (RnewLOGL < RoldLOG(ND)) THEN
            Tnew(L) = Told(ND) 
          ELSE
            CALL SPLINPOX(Tnew(L), RnewLOGL, Told, RoldLOG, ND)
          ENDIF
        ENDDO
      ELSE
C***    Interpolation over radius      
        DO L=1, ND
          IF (Rnew(L) > Rold(1)) THEN
            Tnew(L) = Told(1)
          ELSEIF (Rnew(L) < Rold(ND)) THEN
            Tnew(L) = Told(ND) 
          ELSE
            CALL SPLINPOX(Tnew(L), Rnew(L), Told, Rold, ND)
          ENDIF
        ENDDO
      ENDIF

      RETURN

      END
