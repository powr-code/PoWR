      SUBROUTINE PLOTCON (KANAL,X,Y,N,ISYMBOL)
C***********************************************************************
C***  DIESE ROUTINE DIENST ZUM EINTRAGEN EINER 2. ODER 3. FUNKTION IN DEN
C***  BEGONNENEN PLOT ( PLOTCON ) , BZW. AB ENTRY PLOTTAB ZUM EINTRAGEN
C***  DER ERSTEN FUNKTION
C***********************************************************************
      DIMENSION X(N),Y(N)

C*** Backspace last Line (current ENDE)
      BACKSPACE KANAL

      ENTRY PLOTTAB (KANAL,X,Y,N,ISYMBOL)

C***  CHECK FOR SMALL X- AND Y-VALUES
      DO I=1, N
        IF (ABS(X(I)) .LT. 1.E-30) X(I) = 0.
        IF (ABS(Y(I)) .LT. 1.E-30) Y(I) = 0.
      ENDDO

      WRITE (KANAL,1) N,ISYMBOL
    1 FORMAT (' N=',I5,'   PLOTSYMBOL=',I3)

      DO 3 I=1,N,5
      J=MIN0( 5,N-I+1)
      WRITE (KANAL,2) (X(I+L-1),L=1,J)
      WRITE (KANAL,2) (Y(I+L-1),L=1,J)
    2 FORMAT (1X, 5(G15.8,1X))
    3 CONTINUE

      WRITE (KANAL,4)
    4 FORMAT (' ENDE')

      RETURN
      END
