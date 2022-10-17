      SUBROUTINE SHIFTSTRING (StringArray,L,ND)
C***********************************************************************
C***  SHIFTS ARRAY ELEMENTS R(L) TO R(ND) BY ONE INDEX
C***    uses INTERFACE to overload shift routine to handle
C***     arrays of types real, integer and string(character)
C***********************************************************************

      IMPLICIT NONE
      INTEGER, INTENT(IN) :: L, ND
      CHARACTER(*), DIMENSION(ND), INTENT(INOUT) :: StringArray
      INTEGER :: II, I

      IF (L.LT.1 .OR. L.GT.ND) THEN
        WRITE (0,'(A,2I4)') 'ERROR IN SUBR. SHIFTSTRING, L,ND=',L, ND
        STOP 'SHIFT'
      ENDIF
      DO II=L,ND
        I=ND+L-II
        StringArray(I+1)=StringArray(I)
      ENDDO

      RETURN
      END SUBROUTINE SHIFTSTRING
