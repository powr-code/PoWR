      SUBROUTINE INTERPOLATEXJC(XJCnew,  !old (and new) XJC array
     >                          XJCold,  !old XJC array
     >                          Rnew,    !new radius grid vector
     >                          Rold,    !old radius grid vector
     >                          NF,      !number of frequencies
     >                          ND)      !number of depth points
C**********************************************************************
C***
C***    Interpolation of XJC on new Radius-Grid 
C***    called from: HYDROSOLVE, ENSURETAUMAX
C***
C**********************************************************************

      IMPLICIT NONE
      INCLUDE 'interfacebib.inc'

      INTEGER, INTENT(IN) :: NF, ND
      REAL, DIMENSION(ND), INTENT(IN) :: Rnew, Rold
      REAL, DIMENSION(ND, NF), INTENT(IN) :: XJCold
      REAL, DIMENSION(ND, NF), INTENT(INOUT) :: XJCnew
        
      INTEGER :: L, K

      DO K=1, NF
        DO L=1, ND
          IF (Rnew(L) > Rold(1)) THEN
             XJCnew(L,K) = XJCold(1,K)  !If we are far more out than in the old grid, use old outer boundary value
          ELSEIF (Rnew(L) < Rold(ND)) THEN
             XJCnew(L,K) = XJCold(ND,K) !If we are far more inside than in the old grid, use old inner boundary
          ELSE
            CALL SPLINPOX(XJCnew(L,K),Rnew(L),XJCold(1,K),Rold,ND)
          ENDIF
        ENDDO
      ENDDO

      RETURN
      END
