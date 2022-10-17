      SUBROUTINE SARGREST (TEXT, N, I, IFIRST, LAST)
C**********************************************************************
C***  Ermittelt den Beginn des i-ten Arguments und das nichtleere Ende 
C***  gesamten restlichen Strings. 
C***  Balancierte "..." am Anfang und Ende werden entfernt, sonst 
C***  bleiben " erhalten. Faengt des Argument mit einem Trennzeichen
C***  an (=,:/), so muss (!!!, sonst Fehlerabbruch!) der Reststring in 
C***  "..." eingeschlossen werden, ebenso natuerlich wenn der 
C***  Reststring mit einem Blank beginnen soll. An spaeterer Position sind
C***  Trennzeichen ohne Wirkung. 
C**********************************************************************

      CHARACTER TEXT*(*)

      CALL SARGP (TEXT, N, I, IFIRST, ILAST)
      LAST = IDX(TEXT)

C***  Leerer Reststring
      IF (IFIRST .EQ. -1) THEN
         IFIRST = 1
         LAST   = 1

C***  Der Reststring  beginnt mit einem Trennzeichen ohne " davor
      ELSE IF (IFIRST .EQ. 0) THEN
         WRITE (6, '(A)') 'Error when parsing the following string:'
         WRITE (6, '(A)') TEXT
         WRITE (6, '(A, I3, A)') 'Argument ', I, 
     >                        ' begins with delimiter -> use "..."'
         STOP '>>> ERROR IN SUBROUTINE SARGREST <<<'
      ENDIF

      IF (IFIRST .GT. 1) THEN
C***  Reststring begins with "
         IF (TEXT(IFIRST-1:IFIRST-1) .EQ. '"') THEN 
C***     Remove balanced closing quote if present, else restore leading quote 
            IF (TEXT(LAST:LAST) .EQ. '"' .AND. LAST .GT. IFIRST) THEN
               LAST = LAST - 1
               ELSE
               IFIRST = IFIRST - 1
            ENDIF
         ENDIF
      ENDIF

c      print *, 'Test: ', text
c      print *, 'Test: ', ifirst 
c      print *, 'Test: ', text(ifirst:ifirst)
c      print *, 'Test: ', last
c      print *, 'Test: ', text(last:last)

      RETURN
      END
