      SUBROUTINE NEWDATOMION (ELEMENT,NION,IONDEGREE,PATH,
     >           DRTRANSIT,NLEVEL,LEVEL,NDIM, MAXION, LEVELCOUNT,
     >           FEHLERZEILE, KSHELL)
C**********************************************************************
C***
C***  Subroutine called by NEWDATOM for each element which is to be included, 
C***  reads data from separate DATOM files for each ion
C***  and writes it into the new DATOM output-file in the following order:
C***  LEVELS, LINES, CONTINUUM, DRTRANSIT, K-SHELL 
C***  
C**********************************************************************     
     
      CHARACTER ZEILE*80, NAMEION*5(20), FILENAME*255 
      CHARACTER ELEMENT*2, IONDEGREE(MAXION)*10
      DIMENSION NLEVEL(MAXION)
      LOGICAL, DIMENSION(MAXION) :: DRTRANSIT, KSHELL
      INTEGER LEVELCOUNT(MAXION)
      CHARACTER LEVEL(NDIM, MAXION)*(*)
      LOGICAL LEVELFOUND
      LOGICAL LEVELFOUNDA
      LOGICAL LEVELFOUNDZ
      CHARACTER ZEILE2*80
      CHARACTER ZEILE3*80
      CHARACTER ZEILE4*80
      CHARACTER PATH(MAXION)*255
      CHARACTER FEHLERZEILE(MAXION)*80
      CHARACTER ZIELNIVEAU(MAXION, NDIM)*10
      LOGICAL ZIEL
      INTEGER Z


C     Oeffnen aller Ionendateien fuer ein Element
      
   1  DO ION=1, NION    
         FILENAME=PATH(ION)(:IDX(PATH(ION))) // '/' // 'DATOM.' // 
     >   ELEMENT(:IDX(ELEMENT)) // '_' // IONDEGREE(ION) 
         OPEN (ION,FILE=FILENAME, ACTION = 'READ', ERR=901)
         WRITE (0,'(2A)') 'Opened: ', FILENAME(:IDX(FILENAME))
      ENDDO

      
C     Um als Zielniveaus benoetigte Level vorher zu finden: 
C     Vorbelegen des Arrays mit Leerstrings
      DO ION=1, NION
       DO I=1,NDIM
        ZIELNIVEAU(ION-1,I) = ''
       END DO
      END DO      
      
C     Um vorher zu wissen, welche Level ich zusaetzlich brauche, muss
C     ich schon vorher die Continuumskarten lesen und alle benoetigten
C     Zielniveaus fuer jedes Ion und jedes Level zwischenspeichern.     

      DO ION=1, NION
       Z=0
       DO
        READ (ION,'(A)', END=8 ) ZEILE2
        IF (ZEILE2(:10) .EQ. 'CONTINUUM ') THEN
         Z=Z+1  !Nummer der aktuellen Levelzeile innerhalb eines Ions  
         ZIELNIVEAU(ION,Z)=ZEILE2(70:79) !Zielniveau abspeichern
        END IF       
       END DO
    8  REWIND(ION) !Pointer wieder zurueck auf den Anfang setzen
      END DO
      
C     Hier ist der eigentliche Anfang      
     
C     Level
      WRITE (99,'(A)') !Header fuer Levelkarten ausgeben
     >      '*KEYWORD--  ---NAME--- CH WEIG--ENERGY-- EION----- QN'
      DO ION=1, NION !Schleife ueber alle Ionen
C        WRITE(0,*) 'Levelmax =', NLEVEL(ION)    !Testausgabe   
        LEVELCOUNT(ION) = 0 !Zaehler fuer die Anzahl der Levels
       DO  !Schleife fuer jedes einzelne Ion       
         READ (ION,'(A)', END=2 ) ZEILE2 !Fuer jedes Ion seine
C                                         Datendatei lesen
         IF (ZEILE2(:10) .EQ. 'LEVEL     ') THEN !Levelzeilen raussuchen
          IF (LEVELCOUNT(ION)+1 .GT. NDIM) GOTO 900 !Fehler durch zu
C                                                   viele Levels abfangen
           
C          Pruefen, ob noch Levels als Ziel benoetigt werden oder ob das
C          Array mit den Zielniveaus schon leer geraeumt wurde
           ZIEL = .TRUE.
           DO I=1,NDIM
            ZIEL = ZIEL .AND. ZIELNIVEAU(ION-1,I) .EQ. ''
           END DO 
   
C         Fuer das erste Ion gibt es keine niedrigeren, deshalb prueft 
C         man nur, ob die gewuenschte Zahl an Levels abgearbeitet wurde:
          IF (ION .EQ. 1) THEN 
           IF (LEVELCOUNT(ION)+1 .GT. NLEVEL(ION) .AND. 
     >      NLEVEL(ION) .NE. -1) GOTO 2

C         Fuer alle anderen muss man pruefen, ob die Levelzahl erreicht
C         ist und alle benoetigen Levels aus dem niedrigeren Ion
C         geschrieben wurden
          ELSE    
            IF (LEVELCOUNT(ION)+1 .GT. NLEVEL(ION) .AND. 
     >       NLEVEL(ION) .NE. -1 .AND. ZIEL) GOTO 2
          END IF
   
C          Mitzaehlen, wie viele Levels geschrieben werden:
           LEVELCOUNT(ION) = LEVELCOUNT(ION) + 1 
C          Merken, welches Level gleich geschrieben wird
C           (wird fuer Lines und Continuum benoetigt):
           LEVEL (LEVELCOUNT(ION),ION) = ZEILE2(13:22) 
           WRITE (99,'(A)') ZEILE2(:IDX(ZEILE2)) !Levelzeile

C          Wenn ein Zielniveau geschrieben wurde, wird es aus dem Array
C          geloescht und durch einen Leerstring ersetzt
           DO J=1, NDIM
            IF (ZEILE2(13:22) .EQ. ZIELNIVEAU (ION-1,J))THEN
             ZIELNIVEAU (ION-1,J) = '' 
            END IF   
           END DO

         ELSE IF (ZEILE2(:10) .EQ. 'LINE      ') THEN 
C          Jetzt sind die Levels zuende und es kommen Lines
           IF (LEVELCOUNT(ION) .LT. NLEVEL(ION)) THEN
             WRITE(0,*) 'WARNING: MORE LEVELS REQUESTED THAN AVAILABLE'
             WRITE(0,'(5A,I3,A,I3)') 
     >            '  Ion: ', ELEMENT,' ', IONDEGREE(ION),  
     >            '  Requested: ', NLEVEL(ION), 
     >            '  Available: ', LEVELCOUNT(ION)
           ENDIF
 
           BACKSPACE (ION) !Pointer einen zurueck um auch die erste
C                           Line-Zeile zu erwischen
           EXIT !Wenn die Linien anfangen, soll er hier die
C                Level-Schleife verlassen
         
         ELSE IF (ZEILE2(:10) .EQ. 'CONTINUUM ') THEN !Falls keine
                              !Linien da sind, kommt jetzt das Continuum
           BACKSPACE (ION) !Pointer einen zurueck um auch die erste
C                           Continuum-Zeile zu erwischen
           EXIT !Wenn es keine Linien mehr gibt und nach den Levels das
C           Continuum anfaengt, soll er hier die Line-Schleife verlassen
  
         ENDIF !Ende der Unterscheidung zwischen Levels, Lines,
C               Continuum
          
         CYCLE !Weitere Zeilen fuer ein Ion einlesen
    2    EXIT !Ende der Datei erreicht         
       END DO !Ende der Schleife fuer jedes einzelne Ion

      END DO !Ende der Schleife ueber alle Ionen
      
C     Line 
      WRITE (99,'(A)') !Header fuer Line-Karten ausgeben
     >'*KEYWORD--UPPERLEVEL  LOWERLEVEL--EINSTEIN  RUD-CEY ''
     >--COLLISIONAL COEFFICIENTS--'
      DO ION=1, NION !Schleife ueber alle Ionen
       DO !Einlesen der Datei fuer jedes einzelne Ion
        READ (ION,'(A)',END=3) ZEILE2 
        IF (ZEILE2(:10) .EQ. 'LINE      ') THEN  !Line-Zeilen

C     Ausgangs- und Zielniveau muessen existieren
         LEVELFOUNDA = .FALSE.
         LEVELFOUNDZ = .FALSE.
      
         DO I=1, LEVELCOUNT(ION) !Nach Ausgangs- und Zielniveau suchen
          LEVELFOUNDA = LEVELFOUNDA .OR. LEVEL(I,ION) .EQ. ZEILE2(23:32)
          LEVELFOUNDZ = LEVELFOUNDZ .OR. LEVEL(I,ION) .EQ. ZEILE2(11:20)  
         END DO
      
         IF (LEVELFOUNDA .AND. LEVELFOUNDZ) THEN !Nur wenn beide
C                                                 gefunden wurden, schreibe die Line-Zeile
         WRITE (99,'(A)') ZEILE2(:IDX(ZEILE2))
         END IF
        ELSE IF (ZEILE2(:10) .EQ. 'CONTINUUM ') THEN !Lines sind vorbei,
C                                                     Continuum faengt an
          BACKSPACE (ION)
          EXIT ! Wenn er eine Continuumszeile erreicht, soll er die
C                Lines-Schleife verlassen
        ENDIF !Ende der Fallunterscheidung zw. Lines und Continuum
        CYCLE !gleiche Datei, naechste Zeile
    3  EXIT !Ende der Datei erreicht
       END DO !Ende der Schleife fuer jedes einzelne Ion
      END DO !Ende der Schleife ueber alle Ionen
      
C     Continuum           
      WRITE (99,'(A)') '*KEYWORD  LOWERLEVEL ----SIGMA ----ALPHA
     >----SEXPO -IGAUNT- -KEYCBF- --IONLEV--' !Header schreiben
      DO ION=1, NION-1 !Schleife ueber Ionen, aber beim letzten Ion 
C      duerfen keine Kontinuumskarten mehr geschrieben werden, 
C      weil dann die Zielniveaus fehlen.
       Z=0
       DO !Schleife ueber die Zeilen jeder einzelnen Iondatei
        READ (ION,'(A)',END=4) ZEILE2 !Datei fuer das Ion einlesen
        IF (ZEILE2(:10) .EQ. 'CONTINUUM ') THEN !Continuumszeilen

C     Folgezeilen mitnehmen   
         IF (ZEILE2(52:57) .EQ. 'PIKB12' .OR. ZEILE2(52:59) .EQ. 'OPAPROIX') THEN
          READ (ION,'(A)',END=4) ZEILE3 !naechste Zeile zwischenspeichern
          BACKSPACE(ION) !Cursor wieder zurueck
         END IF 

C     Ausgangslevel muss existieren und gewollt sein  
        LEVELFOUND = .FALSE.
         DO I=1, LEVELCOUNT(ION) !Sucht nach dem Ausgangslevel
          LEVELFOUND = LEVELFOUND .OR. LEVEL(I,ION) .EQ. ZEILE2(11:20)  
         END DO
 
         IF (LEVELFOUND) THEN 
          WRITE (99,'(A)') ZEILE2(:IDX(ZEILE2)) !Schreibe die Zeile
  
          IF (ZEILE2(52:57) .EQ. 'PIKB12' .OR. 
     >     ZEILE2(52:59) .EQ. 'OPAPROIX') THEN
             WRITE (99,'(A)') ZEILE3(:IDX(ZEILE3)) !Schreibe die Folgezeile
          END IF
     
         END IF 
        ELSE IF (ZEILE2(:10) .EQ. 'DRTRANSIT ') THEN !Continuum vorbei,
C                                                     DRTRANSIT beginnt
          BACKSPACE (ION) !Cursor zuruecksetzen um keine Zeile
                          !auszulassen
          EXIT !Continuumsschleife verlassen
        ELSE IF (ZEILE2(:10) .EQ. 'K-SHELL   ') THEN !Kein DRTRANSIT
C                                                     vorhanden, also gleich zu K-SHELL
          BACKSPACE (ION)
          EXIT !Continuumsschleife verlassen
        ENDIF !Ende der Fallunterscheidung zw. Continuum, DRTRANSIT und
C              K-SHELL
        CYCLE !naechste Zeile
    4   EXIT !Ende der Datei erreicht
       END DO !Ende der Schleife ueber Zeilen fuer einzelnes Ion
      END DO !Ende der Schleife ueber alle Ionen

C     DRTRANSIT
      J=1 !Hilfsvariable, nur fuer die Kopfzeile von Bedeutung
      DO ION=1, NION !Schleife ueber alle Ionen
        IF (DRTRANSIT(ION) .EQ. .TRUE.) THEN !Nur wenn DRTRANSIT an ist
         DO !Schleife ueber alle Zeilen eines Ions
          READ (ION,'(A)',END=5) ZEILE2 !Einlesen
          IF (ZEILE2(:10) .EQ. 'DRTRANSIT ') THEN !nur DRTRANSIT-Zeilen
C                                                  nehmen  
C     Ausgangslevel muss existieren und gewollt sein
           LEVELFOUND = .FALSE.
           DO I=1, LEVELCOUNT(ION) !Suche nach dem Ausgangslevel
            LEVELFOUND = LEVELFOUND .OR. LEVEL(I,ION) .EQ. ZEILE2(11:20)  
           END DO
           IF (LEVELFOUND) THEN  !Wenn Ausgangslevel gefunden
C           Schreibe beim ersten Mal den Header   
            IF (J .EQ. 1) WRITE (99,'(A)') '*KEYWORD  LOWERLEVEL' 
     >       // '  UPPERLEVEL --G- --ENERGY-- -EINSTEIN-'  
            WRITE (99,'(A)') ZEILE2(:IDX(ZEILE2)) !Schreibe die Zeile
            J = J+1 !Hochzaehlen, damit der Header nur einmal
C                    geschrieben wird
            END IF
    
           ELSE IF (ZEILE2(:10) .EQ. 'K-SHELL   ') THEN !DRTRANSIT
C                                           vorbei, K-SHELL beginnt
            BACKSPACE (ION)
            EXIT !DRTRANSIT-Schleife verlassen
           END IF
          CYCLE !naechste Zeile
    5     EXIT !Ende der Datei erreicht
         END DO !Ende der Schleife fuer einzelnes Ion
       ELSE IF (ZEILE2(:10) .EQ. 'K-SHELL   ') THEN !Auch wenn es kein
C            DRTRANSIT gab, koennte irgendwann ein K-SHELL-Block kommen
           BACKSPACE (ION)
           EXIT !DRTRANSIT-Schleife verlassen
       END IF  !Ende der Fallunterscheidung zwischen DRTRANSIT und
C               K-SHELL 
      END DO !Ende der Schleife ueber alle Ionen
       
C     K-SHELL 
    7 DO ION=1, NION !Fuer jedes Ion seine Datendatei einlesen
        IF (KSHELL(ION) == .FALSE.) CYCLE
        DO
         READ (ION,'(A)',END=6) ZEILE2 
         IF (ZEILE2(:10) .EQ. 'K-SHELL ') THEN 
C         Falls es eine K-Shell-Zeile gibt, soll einmal der Header 
C         geschrieben werden
          IF (ION .EQ. 1) WRITE (99,'(A)') 
     >      '*KEYWORD--*****SY*I*<-K-SIGMA><-K-SEXPO>***<-K-EION->'
          WRITE (99,'(A)') ZEILE2(:IDX(ZEILE2)) !Zeile schreiben
         ENDIF 
         CYCLE !naechste Zeile
    6    EXIT !Ende der Datei erreicht
        END DO !Ende der Schleife fuer ein Ion
        CYCLE !naechstes Ion
      END DO !Ende der Schleife ueber alle Ionen
                
      
C     Alle Ionendateien wieder schliessen         
   9  DO ION=1, NION
       CLOSE (ION)
      END DO

   99 RETURN
   
   
C********  ERROR EXITS **************************************
  900 WRITE(0,*) 'ERROR: MORE LEVELS THAN DIMENSIONED'
      GOTO 990
      
  901 WRITE(0,*) 'ERROR: Cannot open file:'
      WRITE(0,*) FILENAME(:IDX(FILENAME))
      GOTO 990
      
  903 WRITE(0,*) 'ERROR: cannot read levelnumber'
      GOTO 990        
          
            
  990 WRITE(0,*) '*** The error occured in the following line:'
      WRITE(0,*) FEHLERZEILE(ION)(:IDX(FEHLERZEILE(ION)))
      WRITE(0,*) 'Element was '//ELEMENT
      STOP '*** ERROR DETECTED IN SUBROUTINE NEWDATOMION'
      
         
      END 
