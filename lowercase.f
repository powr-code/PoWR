
        SUBROUTINE LOWERCASE (TEXT)
C ***
C ******************************************************************************
C ***   converts string TEXT to lower case.
C ******************************************************************************

      CHARACTER*(*) TEXT, TAB(2)*26
      DATA TAB /'ABCDEFGHIJKLMNOPQRSTUVWXYZ', 
     >	        'abcdefghijklmnopqrstuvwxyz'/

      DO L = 1, LEN(TEXT)
         J = INDEX( TAB(1), TEXT(L:L) )
         IF( J .NE. 0 ) TEXT(L:L) = TAB(2)(J:J)
      ENDDO

      RETURN
      END
