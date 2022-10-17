      SUBROUTINE FLUSHLEFT (STRING)
C***  Macht einen String linksbuendig
      CHARACTER STRING*(*)
      DO I=1, LEN(STRING)
         IFIRST = I
	 IF (STRING(I:I) .NE. ' ') GOTO 1
      ENDDO
      
      RETURN
      
    1 STRING = STRING(IFIRST:)
      RETURN
      END
      
      
