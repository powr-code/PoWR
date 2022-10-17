      LOGICAL FUNCTION CHRINSTR (TESTCHAR, STRING)
C********************************************************************
C***  Returns .TRUE. if the character TESTCHAR is member of STRING 
C********************************************************************
      CHARACTER TESTCHAR*1, STRING*(*)

      CHRINSTR = .FALSE.
      L = IDX(STRING)

      DO I=1, L
      CHRINSTR = CHRINSTR .OR. (TESTCHAR .EQ. STRING(I:I))
      ENDDO

      RETURN
      END
