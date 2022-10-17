      SUBROUTINE BRNORM2(FN,F,NDIM)
C*****************************************************************
C***  FN := SQR (NORM(F))
C*****************************************************************
      DIMENSION F(NDIM)
      FN=0.
      DO I=1,NDIM
        IF (ABS(F(I)) .GT. 1.E100) THEN
          FN = 1.E300
          GOTO 2
        ENDIF
        FN=FN+F(I)*F(I)
      ENDDO

    2 CONTINUE

      RETURN
      END
