      SUBROUTINE NEWDATOM
C**********************************************************************
C***
C***  Program to create a DATOM (= atomic data) file for the PoWR code
C***      The program needs a database of separate DATOM files 
C***      for each ion to be included. 
C***
C***  This program requires an input file NEWDATOM_INPUT
C***  to specify which elements and ions and how many levels 
C***  are to be included.
C***
C**********************************************************************
      PARAMETER (NDIM=1000) !Maximale Zahl der Levels
      PARAMETER (MAXION=20)
      CHARACTER LEVEL(NDIM, MAXION)*10
      CHARACTER ZEILE*80, NAMEION*5(MAXION), FILENAME*255 
      CHARACTER ELEMENT*2 !Elementsymbol
      CHARACTER IONDEGREE(MAXION)*10 !vorkommende Ionen
      CHARACTER HIGHESTION*10, ROMAN(20)*5
      CHARACTER ACTPAR*20, ACTPARION*20
      CHARACTER ATOM*80
      CHARACTER DEGREE*20
      LOGICAL DRTRANSITGLOBAL, KSHELLGLOBAL, DRTRANSITELEM, KSHELLELEM
      LOGICAL, DIMENSION(MAXION) :: DRTRANSIT, KSHELL
      CHARACTER*10 NAME !Elementname
      REAL STAGE
      DIMENSION NLEVEL(MAXION)
      CHARACTER ZEILE1*80
      CHARACTER*10 CDATE, CTIME
      REAL MASSE
      DIMENSION LEVELCOUNT(MAXION) 
      CHARACTER PATHSAVE*255
      CHARACTER PATH(MAXION)*255
      INTEGER LOW
      INTEGER NUP
      CHARACTER FEHLERZEILE(MAXION)*80

C***  Link data to identify program version
      CHARACTER LINK_DATE*30, LINK_USER*10, LINK_HOST*60
      COMMON / COM_LINKINFO / LINK_DATE, LINK_USER, LINK_HOST

C***  Converting roman numerals to arabic numerals
      DATA ROMAN / 'I', 'II', 'III', 'IV', 'V', 'VI', 'VII', 'VIII',
     >             'IX', 'X', 'XI', 'XII', 'XIII', 'XIV', 'XV', 
     >             'XVI', 'XVII', 'XVIII', 'XIX', 'XX' /  

C***  Write Link Data (Program Version) tp CPR file
      WRITE(0,'(2A)') '>>> NEWDATOM: Program Version from ', 
     >                LINK_DATE
      WRITE(0,'(4A)') '>>> created by ', LINK_USER(:IDX(LINK_USER)),
     >      ' at host ', LINK_HOST(:IDX(LINK_HOST))
    

C***  Error default for PATH to the data
      PATHSAVE = '..UNDEFINED..'


C     Ausgabe in einer neuen Datei 
      OPEN (99,FILE='DATOM', ERR = 904)  
    
C     Einlesen des Steuerfiles
      OPEN (98,FILE = 'NEWDATOM_INPUT', ACTION = 'READ', ERR = 905) 
      
C     Gesamtes Steuerfile als Kommentar oben drueber schreiben, d.h.
C     einlesen und ausgeben
      CALL DATE(CDATE)
      CALL TIME(CTIME)
      WRITE (99,'(A)') 
     > '* This DATOM-FILE has been created at ' // CDATE 
     >                  // '' // CTIME
      WRITE (99,'(A)') '* with the program NEWDATOM' 
      WRITE (99,'(A)') '*'
      WRITE (99,'(A)') '* NEWDATOM_INPUT was'       
    7 DO
       READ (98,'(A)', END=1) ZEILE1
       IF (ZEILE1(:1) .EQ. '*' .OR. ZEILE1(:1) .EQ. '-') GOTO 7
       WRITE (99,'(A)') '* ' // ZEILE1(:IDX(ZEILE1))
      END DO
    1 CONTINUE
      WRITE (99,'(A)') '*'      

      REWIND(98)   ! Am Schluss den Pointer wieder an den Anfang setzen
      
C     Defaults 
      DRTRANSITGLOBAL = .FALSE.
      KSHELLGLOBAL = .TRUE.
      

      
C     Schleife fuer die Elemente    
      DO 
    3  READ (98,'(A)', END = 5) ZEILE !Datei zeilenweise durchlesen
       IF (ZEILE .EQ. '') GOTO 3 !Leerzeilen abfangen 
       CALL SARGV(ZEILE,1,ACTPAR)   !Erstes Wort der Zeile (Keyword)
C      einlesen und anschliessend untersuchen (Fallunterscheidung)

C      Path-Variable einlesen in PATHSAVE
       IF (ACTPAR .EQ. 'PATH') THEN
         CALL SARGV(ZEILE,2,PATHSAVE) 

C      Hier wird DRTRANSIT global fuer alle Ionen der folgenden 
C       Elemente angeschaltet  
       ELSEIF (ACTPAR .EQ. 'DRTRANSIT') THEN 
        DRTRANSITGLOBAL = .TRUE. 

C      Hier wird DRTRANSIT global fuer alle Ionen der folgenden 
C       Elemente abgeschaltet  
       ELSEIF (ACTPAR .EQ. 'NO-DRTRANSIT') THEN 
        DRTRANSITGLOBAL = .FALSE. 

       ELSE IF (ACTPAR .EQ. 'NO-K-SHELL') THEN
        KSHELLGLOBAL = .FALSE.       

       ELSE IF (ACTPAR .EQ. 'K-SHELL') THEN
        KSHELLGLOBAL = .TRUE.       
        
C      Einzelnes Element wird hier drin abgearbeitet
       ELSE IF (ACTPAR .EQ. 'ELEMENT') THEN
          NION = 0 !Zaehlt, wie viele Ionen angegeben sind
          CALL SARGC (ZEILE, NPAR) !Anzahl der Parameter ermitteln
          IF (NPAR .LT. 2) GOTO 900 !Fehler abfangen
          CALL SARGV(ZEILE,2,ELEMENT) !Element = zweites Wort der Zeile
  
C         Ausnahme fuer Generic: Hier wird nur Ionlow und Iontop aus dem
C         Steuerfile eingelesen und hingeschrieben, alles andere ist fest
  
          IF (ELEMENT .EQ. 'G') THEN
           CALL SARGC (ZEILE, NPAR) !Anzahl der Parameter ermitteln
           IF (NPAR .LT. 4) GOTO 903 !Fehler abfangen
           CALL SARGV(ZEILE,3,ACTPAR) !Ionlow steht an dritter Stelle in
C                                 der entsprechenden Zeile im Steuerfile
           READ (ACTPAR,'(I20)') LOW  !Typumwandlung
           CALL SARGV(ZEILE,4,ACTPAR) !Iontop steht an vierter Stelle
           READ (ACTPAR,'(I20)') NUP  !Typumwandlung
C          Ausgabe der Daten fuer Generic
           WRITE (99,'(A)')
     >        '*KEYWORD--  ---NAME--- SYMB   IONLOW  IONTOP'
           WRITE (99,'(A)')
     >        '**************************************************'
           WRITE (99,'(A27,7X,I2,6X,I2)')
     >       'ELEMENT     GENERIC    (G )', LOW, NUP
           WRITE (99,'(A)')
     >       '*           =======                              *'
           WRITE (99,'(A)')
     >       '**************************************************'
   
C         Andere Elemente
          ELSE

C          We have to find the highest ionization stage
           NBACKSPACE = 0
           DO
             READ (98,'(A)',END=41) ZEILE
             IF (ZEILE .EQ. '') GOTO 42
             CALL SARGV (ZEILE,1,ACTPARION)
             IF (ACTPARION .EQ. 'ION') THEN
                CALL SARGC (ZEILE,NPAR) 
                IF (NPAR .LT. 2) GOTO 902 ! Fehler
                CALL SARGV(ZEILE,2,HIGHESTION) 
                NBACKSPACE = NBACKSPACE + 1
             ELSEIF (ACTPARION .EQ. 'ELEMENT') THEN
 41             NBACKSPACE = NBACKSPACE + 1
                EXIT
             ELSE 
 42             NBACKSPACE = NBACKSPACE + 1
             ENDIF             
           ENDDO
C          Move REC-Pointer back to first line of ELEMENT
           DO IBACKSPACE = 1, NBACKSPACE
              BACKSPACE 98
           ENDDO
           DO I=1, 20
              IF (HIGHESTION .EQ. ROMAN(I)) THEN
                 MAINSTAGE = MAX((I-1),1)
                 EXIT
              ENDIF
           ENDDO
              
C          Aufruf einer Subroutine, die aus dem Elementsymbol auf
C          weitere Eigenschaften schliesst  
           CALL FINDELEMENT(ELEMENT, NAME, STAGE, ATOMICMASS)
C          Ausgabe des Headers fuer das aktuelle Element
           WRITE (99,'(A)')
     >       '*======================================================='
           WRITE (99,'(A)')
     >       '*KEYWORD--  ---NAME--- SYMB   ATMASS   STAGE' 
           WRITE (99,'(A)')
     >       '************************************************'
C              USE MAINSTAGE INSTEAD OF STAGE
           WRITE (99, '(A7,5X,A10,1X,A1,A2,A1,3X,F6.2,3X,I3)') 
     >       'ELEMENT', NAME, '(', ELEMENT, ')', ATOMICMASS, MAINSTAGE
  
           WRITE (99,'(A)')
     >       '*           ======                             *'
           WRITE (99,'(A)')
     >       '************************************************'

C          default: global settings for DR-TRANSIT and K-SHELL in element
           DRTRANSITELEM = DRTRANSITGLOBAL
           KSHELLELEM = KSHELLGLOBAL
           
C***       Read ION options for each ion
            
    4      READ (98,'(A)', END = 2) ZEILE !Ionenzeilen nacheinander
C                                          einlesen
           IF (ZEILE .EQ. '') GOTO 4 !Leerzeilen abfangen
           CALL SARGV(ZEILE,1,ACTPARION) !erstes Wort ist Schluesselwort
   
C          Pfadvariable aendert sich
           IF (ACTPARION .EQ. 'PATH') THEN
C                              ====
            CALL SARGV(ZEILE,2,PATHSAVE) !PATH neu setzen
    
C          Kommandos zwischen ELEMENT und vor einem ION-Kommando 
C           koennen DR-TRANSIT und K-SHELL fuer alle folgenden Ionen
C           desselben Elements ein- oder ausschalten
           ELSEIF (ACTPARION == 'DRTRANSIT') THEN 
             DRTRANSITELEM = .TRUE. 

           ELSEIF (ACTPARION == 'NO-DRTRANSIT') THEN 
             DRTRANSITELEM = .FALSE. 

           ELSEIF (ACTPARION == 'NO-K-SHELL') THEN
             KSHELLELEM = .FALSE.       

           ELSEIF (ACTPARION == 'K-SHELL') THEN
             KSHELLELEM = .TRUE.       
    
C          Ionenzeile, hier werden Ionname und Parameter fuer einzelne 
C          Ionen angegeben
           ELSE IF (ACTPARION .EQ. 'ION') THEN
C                                   ===
             CALL SARGC (ZEILE, NPAR) !Zaehle die Parameter der Zeile
             IF (NPAR .LT. 2) GOTO 902 !Fehler abfangen
             NION = NION + 1 !Ionenzaehler hochsetzen fuer Schleifen in
C                             newdatomion
             FEHLERZEILE(NION) = ZEILE

             IF (PATHSAVE .EQ. '..UNDEFINED..') GOTO 908

             PATH(NION) = PATHSAVE !aktueller Pfad           

             CALL SARGV(ZEILE,2,IONDEGREE(NION)) !IONDEGREE einlesen     

C            entscheiden, ob DRTRANSIT mitgenommen wird (default = element setting)
             DRTRANSIT(NION) = DRTRANSITELEM

C            entscheiden, ob K-SHELL mitgenommen wird (default = element setting)
             KSHELL(NION) = KSHELLELEM

C            Zahl der Levels einschraenken oder (wenn keine
C            Einschraenkung gegeben ist) alle nehmen
             NLEVEL(NION) = -1 !Default: Wenn kein Level angegeben
C                               ist, nimmt er alle     
             DO ICOUNT = 1, NPAR !Gehe alle Parameter durch und pruefe,
C                                 ob irgendwo DRTRANSIT, NO-DRTRANSIT
C                                 oder NLEVEL steht      
              CALL SARGV(ZEILE,ICOUNT,ACTPAR) !Zeile wird in Parameter
C                                              zerlegt
              IF (ACTPAR .EQ. 'NO-DRTRANSIT') THEN
               DRTRANSIT(NION) = .FALSE. !DRTRANSIT fuer das ION
C                                         ausschalten
              ELSE IF (ACTPAR .EQ. 'DRTRANSIT') THEN
               DRTRANSIT(NION) = .TRUE. !DRTRANSIT fuer dieses ION
C                                        einschalten       
              ELSE IF (ACTPAR .EQ. 'NO-K-SHELL') THEN
               KSHELL(NION) = .FALSE.   !K-SHELL fuer das ION
C                                         ausschalten
              ELSE IF (ACTPAR .EQ. 'K-SHELL') THEN
               KSHELL(NION) = .TRUE.    !K-SHELL fuer dieses ION
C                                        einschalten       
              ELSE IF (ACTPAR .EQ. 'NLEVEL') THEN !Zahl der Level
C                                                  beschraenken
               IF (ICOUNT+1 .GT. NPAR) GOTO 906
               CALL SARGV(ZEILE,ICOUNT+1,ACTPAR) !Levelzahl wird
C                           fuer jedes Ion als Character abgespeichert
               READ (ACTPAR, '(I20)', ERR = 907) NLEVEL(NION)

              END IF 
             END DO      
           
           ELSEIF (ACTPARION .EQ. 'ELEMENT') THEN 
C             Reading for one Element is finished: 
C             Reset the pointer for continuing later with next element
              BACKSPACE (98)     
              GOTO 2
     
           ENDIF  

C***       Read next line within ELEMENT block
           GOTO 4 

C      Reading input options for last ELEMENT completed
C      (next ELEMENT or E-O-F encountered) 
C      Subroutine NEWDATOMION wird fuer jedes Element ein Mal
C      aufgerufen, Informationen fuer einzelne Ionen stehen in Arrays 
C      und werden hier uebergeben   

    2      CALL NEWDATOMION (ELEMENT,NION,IONDEGREE,
     >          PATH,DRTRANSIT,NLEVEL,LEVEL,NDIM,MAXION, 
     >          LEVELCOUNT, FEHLERZEILE, KSHELL)
         
          END IF !Element finished
                 
       ENDIF  
               
      END DO !End of loop reading input lines
    5 CONTINUE 
      
      CLOSE (99) !Schliessen der output-DATOM-Datei
      CLOSE (98) !Schliessen des Steuerfiles

C**   Regular program end - remove error code 
   99 CONTINUE
      CALL JSYMSET ('G0', '0')
      STOP 'O.K.'

C********  ERROR EXITS **************************************
  900 WRITE(0,*) 'ERROR: ELEMENT option without parameter'
      GOTO 990
      
  901 WRITE(0,*) 'ERROR: ELEMENT GENERIC NOT FOUND'
      GOTO 990    
      
  902 WRITE(0,*) 'ERROR: ION option without parameter'
      GOTO 990     
  
  903 WRITE(0,*) 'ERROR: ELEMENT GENERIC without parameter'
      GOTO 990         
      
  904 WRITE(0,*) 'ERROR: CANNOT WRITE IN FILE DATOM'
      GOTO 990     
      
  905 WRITE(0,*) 'ERROR: NEWDATOM_INPUT NOT FOUND'
      GOTO 990       

  906 WRITE(0,*) 'ERROR: LEVEL option without parameter'
      GOTO 990   
  
  907 WRITE(0,*) 'ERROR: invalid parameter for NLEVEL'
      GOTO 990      
     
  908 WRITE(0,*) 'ERROR: PATH to atomic data not (yet) defined'
      GOTO 990
          
            
  990 WRITE(0,*) '*** The error occured in the following line:'
      WRITE(0,*) ZEILE(:IDX(ZEILE))
      STOP '*** ERROR DETECTED IN PROGRAM NEWDATOM'
      
      
      END
