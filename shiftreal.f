      SUBROUTINE SHIFTREAL (RealArray,L,ND)
C***********************************************************************
C***  SHIFTS ARRAY ELEMENTS R(L) TO R(ND) BY ONE INDEX
C***    uses INTERFACE to overload shift routine to handle
C***     arrays of types real, integer and string(character)
C***********************************************************************
      
      IMPLICIT NONE
      INTEGER, INTENT(IN) :: L, ND
      REAL, DIMENSION(ND), INTENT(INOUT) :: RealArray
      INTEGER :: II, I

      IF (L.LT.1 .OR. L.GT.ND) THEN
        WRITE (0,'(A,2I4)') 'ERROR IN SUBR. SHIFTREAL, L,ND=',L, ND
        STOP 'SHIFT'
      ENDIF
      DO II=L,ND
        I=ND+L-II
        RealArray(I+1)=RealArray(I)
      ENDDO

      RETURN
      END SUBROUTINE SHIFTREAL
