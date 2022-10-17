      SUBROUTINE FOLR (N, X, IX, A, IA)

      DIMENSION X(2), A(2)

      IF (IX .NE. 1 .OR. IA .NE. 1) THEN
        STOP 'ERROR IN FOLR: INCREMENTS NOT ALLOWED'
      ENDIF

      DO I=2, N
        A(I) = A(I) - X(I)*A(I-1)
      ENDDO

      RETURN
      END
