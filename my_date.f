      SUBROUTINE MY_DATE (STRING)

      CHARACTER STRING*8, STRING2*10

      CALL DATE(STRING2)

C***  Create STRING as yy/mm/dd
      STRING = '  /  /  '

C***  Year
      STRING(1:2) = STRING2(8:9)

C***  Month
      IF      (STRING2(4:6) .EQ. 'Jan') THEN
        STRING(4:5) = '01'
      ELSE IF (STRING2(4:6) .EQ. 'Feb') THEN
        STRING(4:5) = '02'
      ELSE IF (STRING2(4:6) .EQ. 'Mar') THEN
        STRING(4:5) = '03'
      ELSE IF (STRING2(4:6) .EQ. 'Apr') THEN
        STRING(4:5) = '04'
      ELSE IF (STRING2(4:6) .EQ. 'May') THEN
        STRING(4:5) = '05'
      ELSE IF (STRING2(4:6) .EQ. 'Jun') THEN
        STRING(4:5) = '06'
      ELSE IF (STRING2(4:6) .EQ. 'Jul') THEN
        STRING(4:5) = '07'
      ELSE IF (STRING2(4:6) .EQ. 'Aug') THEN
        STRING(4:5) = '08'
      ELSE IF (STRING2(4:6) .EQ. 'Sep') THEN
        STRING(4:5) = '09'
      ELSE IF (STRING2(4:6) .EQ. 'Oct') THEN
        STRING(4:5) = '10'
      ELSE IF (STRING2(4:6) .EQ. 'Nov') THEN
        STRING(4:5) = '11'
      ELSE IF (STRING2(4:6) .EQ. 'Dec') THEN
        STRING(4:5) = '12'
      ELSE
        WRITE (0,*) 'MONTH NOT RECOGNIZED'
        WRITE (*,*) 'STRING2(4:6)=',STRING2(4:6)
        STOP 'ERROR IN MY_DATE'
      ENDIF

C***  day
      STRING(7:8) = STRING2(1:2) 


      RETURN
      END
