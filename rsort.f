      SUBROUTINE RSORT(MAXIND,RFELD,IFELD)

C     **************************************************************
C     *  ROUTINE RSORT    ERSTELLT 7.MAR.89   VON    G. DUENNEBEIL *
C     *                   GEANDERT  .......   VON        .....     *
C     *                                                            *
C     *  Die Routine sortiert ein Feld von Realzahlen              *
C     *  Das Feld wird aufsteigend (RFELD(1) = kleinste Zahl)      *
C     *  sortiert. Dabei laufen die Feldgrenzen von 1..MAXIND.     *
C     *  Das Feld IFELD wird dabei mitsortiert, so das es als      *
C     *  auf weitere nicht mitsortierte Felder dienen kann.        *
C     *  IFELD muss vorbesetzt sein, wie es die aufrufende Routine *
C     *  benoetigt. RSORT stellt hier keinerlei Bedingungen        *
C     **************************************************************

C     **************************************************************
C     *  Beschreibung des Verfahrens:                              *
C     *  Das gesamte Feld zwischen ITOP und MAXIND wird als        *
C     *  unsortiert angenommen. Im unsortierten Feld wird das      *
C     *  kleinste Element gesucht und mit dem obersten vertauscht. *
C     *  Danach ist also das Feld von 1 bis ITOP sortiert.         *
C     *  Dementsprechend wird ITOP um eins erhoeht und das         *
C     *  Verfahren startet erneut.                                 *
C     *  Das Feld ist garantiert sortiert, wenn ITOP=MAXIND-1      *
C     **************************************************************

      DIMENSION RFELD(MAXIND)
      DIMENSION IFELD(MAXIND)

      DO 20 ITOP=1,MAXIND-1
C     Diese Schleife bringt das jeweils kleinste Element nach oben
C     und verschiebt 'oben' dabei nach unten

        XMIN=RFELD(ITOP)
        NMIN=ITOP

        DO 10 LINDEX=ITOP+1,MAXIND
C         Diese Schleife sucht das kleinste Element im Feld
          IF ( RFELD(LINDEX) .LT. XMIN ) THEN
             XMIN=RFELD(LINDEX)
             NMIN=LINDEX
          ENDIF
10      CONTINUE

C       Hier werden das Top-Element und das kleinste Element vertauscht
        XMIN=RFELD(NMIN)
        RFELD(NMIN)=RFELD(ITOP)
        RFELD(ITOP)=XMIN

        ITEMP=IFELD(NMIN)
        IFELD(NMIN)=IFELD(ITOP)
        IFELD(ITOP)=ITEMP
20    CONTINUE

      RETURN
      END
