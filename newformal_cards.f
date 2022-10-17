      SUBROUTINE NEWFORMAL_CARDS
C*******************************************************************************
C***
C***  Program to create a FORMAL_CARDS file for the PoWR code
C***      The program needs a database of separate FORMAL_CARDS files 
C***      for each ion to be included. 
C***
C***  This program uses the existing DATOM file and requires an
C***  inputfile NEWFORMAL_CARDS_INPUT to specify the ranges.
C***
C*******************************************************************************

      IMPLICIT NONE
      INCLUDE 'interfacebib.inc'

C***  SET ARRAY DIMENSION PARAMETERS
      INTEGER, PARAMETER :: MAXATOM =          26 
      INTEGER, PARAMETER :: NDIM    =        1560 
      INTEGER, PARAMETER :: MAXAUTO =        2850 
      INTEGER, PARAMETER :: MAXIND  =       20000 
      INTEGER, PARAMETER :: MAXKONT =        NDIM 
      INTEGER, PARAMETER :: MAXKODR =        NDIM 
      INTEGER, PARAMETER :: MAXMULTI=       10000 
      INTEGER, PARAMETER :: MAXDRTRANSIT =   1000 
      INTEGER, PARAMETER :: MAXIONRESTRICT =  200
      INTEGER, PARAMETER :: LENIONNAME =       10
      
      CHARACTER(80) :: IONZEILE, HILFSZEILE, IONZEILE2,
     >                 MULTIZEILE, RANGENAME
      CHARACTER(30) :: CDATETIME
      CHARACTER(10) :: CDATE, CTIME, CELEM, CION,
     >                 UPPERLEVEL,  !oberes Level eines Multiplets
     >                 LOWERLEVEL   !unteres Level eines Multiplets
      CHARACTER(20) :: INPUTKEYWORD, IONKEYWORD, UNTERGRENZE, 
     >                 OBERGRENZE, LAMBDA, LAMBDA1,
     >                 UPNAME, LOWNAME, NAME,       !Levelbezeichnungen
     >                 UPENERGIE, LOWENERGIE, WORD, RESTRICTNAME
      CHARACTER(20), DIMENSION(10) :: RANGENAMES    !Array erlaubt Aliasnamen fuer Ranges
      CHARACTER(200) :: TESTPATH, STANDARDPATH,  INPUTZEILE
      CHARACTER(200), DIMENSION(100) :: FILENAME
      REAL :: MIN, MAX                              !Entspricht Ober- und Untergrenze
      INTEGER :: K,                                 !Zaehler fuer Ionen, ua auch Kanalnummer
     >           L, Z,                              !Indizes fuer Eintraege im Multipletarray
     >           ILEVEL                             !Laufindex fuer Level aus DATOM (1 bis N)
      CHARACTER(80), DIMENSION(MAXMULTI) :: MULTI   !Array fuer Multiplet-Bloecke
      CHARACTER(80), DIMENSION(MAXDRTRANSIT) :: DRTRANSIT !Array fuer DRTRANSIT-Bloecke
      REAL WL, UPE, LOWE     !Wellenlaenge, obere Energie, untere Energie
      REAL DELTAE           !Energiedifferenz
      INTEGER :: LMAX            !Anzahl der Zeilen im Multiplet-Array
      INTEGER :: NMAX            !Gesamtanzahl der Ionen
      INTEGER :: IION            !Schleifenindex
      INTEGER :: INDEXL, INDEXU  !Indizes fuer oberes und unteres Level
      LOGICAL :: PASST,    !true, wenn Wellenlaenge im Band
     >           OPENED,   !true, wenn Ionendateien geoeffnet sind
     >           ERROR,    !true, wenn Fehler beim Oeffnen mit Testpfad
     >           HEADER,   !true, wenn Header fuer Ion geschrieben werden soll
     >           LEVELUP,  !true, wenn oberes Level im DATOM vorhanden
     >           LEVELLOW, !true, wenn unteres Level im DATOM vorhanden
     >           ANFANG, BEXIST,
     >           bSKIPION, bRESTRICT, bUseDRTRANSIT


      REAL :: DXFE, XLAM0FE, VDOPFE
      
      INTEGER :: I, IATOM, IFIRST, ILAST, LASTSELF, 
     >           N, NPAR, NION, NATOM, NAUTO, NWORDS, IWORD,
     >           LASTFE, LASTKDR, LASTKON, LASTIND, IALIAS, NALIAS,
     >           IRESTRICT, LASTRESTRICT, LENRES

      INTEGER, EXTERNAL :: IDX

C***  HANDLING OF DIELECTRONIC RECOMBINATION / AUTOIONIZATION (SUBR. DATOM)
      INTEGER, DIMENSION(MAXAUTO) :: LOWAUTO, IONAUTO, KRUDAUT
      INTEGER, DIMENSION(MAXKODR) :: KODRNUP, KODRLOW
      REAL, DIMENSION(MAXAUTO) :: WAUTO, EAUTO, AAUTO
      
      INTEGER, DIMENSION(NDIM) :: NCHARG, MAINQN, NOM, IONGRND
      INTEGER, DIMENSION(MAXATOM) :: KODAT, NFIRST, NLAST
      REAL, DIMENSION(NDIM) :: WEIGHT, ELEVEL, EION
      REAL, DIMENSION(NDIM,NDIM) :: EINST
      REAL, DIMENSION(MAXATOM) :: ATMASS, STAGE

      REAL, DIMENSION(MAXATOM,MAXATOM) :: SIGMATHK, SEXPOK, EDGEK
      REAL, DIMENSION(4,NDIM) :: ALTESUM
      
      REAL, DIMENSION(MAXKONT) :: ALPHA, SEXPO, 
     >                            ADDCON1, ADDCON2, ADDCON3
      INTEGER, DIMENSION(MAXKONT) :: KONTNUP, KONTLOW
      CHARACTER*8 IGAUNT(MAXKONT), KEYCBF(MAXKONT)
      CHARACTER(10), DIMENSION(NDIM) ::    LEVEL
      CHARACTER(10), DIMENSION(MAXATOM) :: ELEMENT
      CHARACTER(4), DIMENSION(MAXIND) ::   KEYCBB
      CHARACTER(2), DIMENSION(MAXATOM) ::  SYMBOL
      INTEGER, DIMENSION(MAXIND) :: INDNUP, INDLOW
      REAL, DIMENSION(4,MAXIND) :: COCO

C***  IRON: 
C***  include "dimblock"
      INTEGER, PARAMETER :: INDEXMAX = 1E7, NFEREADMAX = 3E5    !std
C      INTEGER, PARAMETER :: INDEXMAX = 4E7, NFEREADMAX = 5E5     !vd20
C      INTEGER, PARAMETER :: INDEXMAX = 1E8, NFEREADMAX = 6E5     !xxl

      INTEGER, PARAMETER :: MAXFEIND  =       1500 
 
      INTEGER, DIMENSION(MAXFEIND) :: INDRB, INDRF, IFRBSTA, IFRBEND,
     >                                IFENUP, IFELOW
      REAL, DIMENSION(MAXFEIND) :: SIGMAINT

      REAL, DIMENSION(NFEREADMAX) :: FEDUMMY
      REAL, DIMENSION(INDEXMAX) :: SIGMAFE
      LOGICAL :: BFEMODEL
      
      INTEGER, DIMENSION(MAXIONRESTRICT) :: LISTRESTRICT 
      
      CHARACTER(10) :: TIM1

C***  Operating system:
      CHARACTER(8) :: OPSYS
      COMMON / COMOS / OPSYS

C***  Dimensions for the GETION programm
      LOGICAL :: BNEWION, BIRON, BIONSTRICT
      INTEGER, PARAMETER :: MAXIONNAMES = 200
      CHARACTER(LENIONNAME), DIMENSION(MAXIONNAMES) :: IONNAME
      CHARACTER(10), DIMENSION(MAXIONNAMES) :: ELEMENTNAME
      CHARACTER(7), DIMENSION(MAXIONNAMES) :: NUMMER
      CHARACTER(2), DIMENSION(MAXIONNAMES) :: SYM

      INTEGER, PARAMETER :: MAXION = 27
      CHARACTER(4), DIMENSION(0:MAXION) :: ROMION
      
      !File and channel handles (=KANAL)
      INTEGER, PARAMETER :: hOUT   = 6    !write to wruniqX.out (stdout)
      INTEGER, PARAMETER :: hCPR   = 0    !write to wruniqX.cpr (stderr)
      INTEGER, PARAMETER :: hFC    = 999  !write to FORMAL_CARDS
      INTEGER, PARAMETER :: hNFCIN = 998  !read from NEWFORMAL_CARDS_INPUT
      
            
C***  Ende der Variablen fuer getion ***********************************

C***  Link data to identify program version
      CHARACTER(30) :: LINK_DATE
      CHARACTER(10) :: LINK_USER
      CHARACTER(60) :: LINK_HOST
      COMMON / COM_LINKINFO / LINK_DATE, LINK_USER, LINK_HOST

C***  Write Link Data (Program Version) tp CPR file
      WRITE(hCPR,'(2A)') '>>> NEWFORMAL_CARDS: Program Version from ', 
     >                LINK_DATE
      WRITE(hCPR,'(4A)') '>>> created by ', LINK_USER(:IDX(LINK_USER)),
     >      ' at host ', LINK_HOST(:IDX(LINK_HOST))
    
      WRITE(*,*) 'BFEMODEL', BFEMODEL
C***  UNTERPROGRAMM GETION  *****************************************************
C     SUBROUTINE GETION (MAXIONNAMES, NION, IONNAME)
C*******************************************************************************
C***  Test program for finding which Ions are in a DATOM file 
C*******************************************************************************

      DATA (ROMION(I),I=0,MAXION) /'I   ','II  ','III ','IV  ','V   ',
     >                             'VI  ','VII ','VIII','IX  ','X   ',
     >                             'XI  ','XII ','XIII','XIV ','XV  ',
     >                             'XVI ','XVII','18  ','XIX ','XX  ',
     >                             'XXI ','XXII','23  ','XXIV','XXV ',
     >                             'XXVI','27  ','28  '/

      CALL INSTALL
                                 
      IF (OPSYS .EQ. 'CRAY' .OR. OPSYS .EQ. 'SGI') THEN
        CALL CLOCK(TIM1)
      ELSE
        CALL TIME(TIM1)
      ENDIF
      
* Eisenabfrage      
      BIRON = .FALSE.
      WRITE(*,*) 'BIRON=', BIRON

C***  READING THE ATOMIC DATA FROM FILE DATOM
      CALL       DATOM (NDIM,N,LEVEL,NCHARG,WEIGHT,ELEVEL,EION,MAINQN,
     $                  EINST,ALPHA,SEXPO,
     $                  ADDCON1, ADDCON2, ADDCON3, 
     $                  IGAUNT,COCO,KEYCBB,ALTESUM,
     $                  INDNUP,INDLOW,LASTIND,MAXIND,MAXATOM,NATOM,
     $                  ELEMENT,SYMBOL,NOM,KODAT,ATMASS,STAGE,
     $                  SIGMATHK,SEXPOK,EDGEK,NFIRST,
     $                  NLAST,NAUTO,MAXAUTO,LOWAUTO,WAUTO,EAUTO,AAUTO,
     $                  IONAUTO,KRUDAUT,KONTNUP,KONTLOW,LASTKON,MAXKONT,
     $                  IONGRND,KODRNUP,KODRLOW,LASTKDR,MAXKODR,KEYCBF,
C***  IRON: ADDITIONAL PARAMETERS FOR IRON-GROUP LINE BLANKETING
     >                  'NOIRON', INDEXMAX, NFEREADMAX, MAXFEIND,
     >                  LASTFE, SIGMAFE, INDRB, INDRF,
     >                  IFENUP, IFELOW, IFRBSTA, IFRBEND, FEDUMMY, 
     >                  VDOPFE, DXFE, XLAM0FE, SIGMAINT, BFEMODEL)
 
      WRITE(*,*) 'BFEMODEL', BFEMODEL
      IF (BFEMODEL) THEN
        BIRON = .TRUE.
        WRITE(*,*) 'BIRON=', BIRON
      ENDIF

C***  Make a list of all ions (exept iron) that are encountered in DATOM
      NION = 0 
      DO IATOM=1, NATOM
C        WRITE(*,*) 'BIRON=', BIRON
        WRITE(*,*) 'IATOM=', IATOM
        WRITE(*,*) 'NATOM=', NATOM
        WRITE(*,*) 'KODAT(26)=', KODAT(26)      
C*      Skip iron 
        IF (IATOM .EQ. KODAT(26)) THEN
          BIRON = .TRUE.
          WRITE(*,*) 'BIRON=', BIRON
          CYCLE
        ENDIF

        DO I = NFIRST(IATOM), NLAST(IATOM)
          BNEWION = (I == NFIRST(IATOM))
          IF (.NOT. BNEWION) BNEWION = (NCHARG(I-1) /= NCHARG(I))
              
          IF (BNEWION) THEN
            NION = NION + 1
            IF (NION > MAXIONNAMES) GOTO 907
            WRITE(*,*) 'ELEMENT=', ELEMENT(IATOM)
    
            IF (ELEMENT(IATOM) /= 'GENERIC') THEN
              ELEMENTNAME(NION) = ELEMENT(IATOM)
              SYM(NION) = SYMBOL(IATOM)
              IONNAME(NION) = SYMBOL(IATOM) 
              IONNAME(NION) = IONNAME(NION)(:IDX(IONNAME(NION))) // 
     >                       '_' // ROMION(NCHARG(I))
              NUMMER(NION) = ROMION(NCHARG(I))
      
              WRITE(*,*) 'BIRON=', BIRON
     
            ELSE
              BIRON = .TRUE.
              WRITE(*,*) 'BIRON=', BIRON
            ENDIF
          ENDIF
        ENDDO !NFIRST-NLAST loop
      ENDDO !NATOM loop


C*****************Ende von getion*****************************************

C***  Error default for PATH to the data
      STANDARDPATH = '..UNDEFINED..'

C     Ausgabe in einer neuen Datei 
      OPEN (hFC,FILE='FORMAL_CARDS', ERR = 904)  
      
C     Oeffnen des Steuerfiles NEWFORMAL_CARDS_INPUT 
      OPEN (hNFCIN,FILE='NEWFORMAL_CARDS_INPUT',ACTION = 'READ',ERR=903)
      
C     Header erstellen
      CALL DATE(CDATE)
      CALL TIME(CTIME)
      WRITE (hFC,'(A)') 
     > '* This FORMAL_CARDS-FILE has been created at ' 
     >         // CDATE(:IDX(CDATE)) // ' ' // CTIME
      WRITE (hFC,'(A)') '* with the program NEWFORMAL_CARDS' 
      WRITE (hFC,'(A)') '* and the following input from file'
      WRITE (hFC,'(A)') '* NEWFORMAL_CARDS_INPUT:' 
      WRITE (hFC,'(A)') '*' 
   11   READ (hNFCIN,'(A)', END=12) INPUTZEILE 
        WRITE (hFC,'(A)') '* ' // INPUTZEILE(:IDX(INPUTZEILE)) 
        GOTO 11
   12 WRITE (hFC,'(A)') '********** End of NEWFORMAL_CARDS_INPUT ***'  
      WRITE (hFC,'(A)') ' '
      REWIND (hNFCIN)
   
C     Defaults
      K=1  
      PASST = .FALSE.   
      UPENERGIE = ''
      LOWENERGIE = '' 
      UPE = 0
      LOWE = 0
      DO I=1, MAXMULTI
        MULTI(I) = ''
      ENDDO
      DO I=1, MAXDRTRANSIT
        DRTRANSIT(I) = ''
      ENDDO
      bUseDRTRANSIT = .FALSE.       !default: NO-DRTRANSIT

      TESTPATH = ''
      OPENED = .FALSE.
      NMAX = 0
      ANFANG = .TRUE.
      LAMBDA = ''
      INDEXL = 0
      INDEXU = 0
      HEADER = .FALSE.
      IONZEILE = 'not opened --> no Line read'
      IONZEILE2 = 'not opened --> no Line read'
      LEVELUP = .FALSE.
      LEVELLOW = .FALSE.
      UPPERLEVEL = ''
      LOWERLEVEL = ''
      bRESTRICT = .FALSE.
      LASTRESTRICT = 0
      DO IRESTRICT=1, MAXIONRESTRICT
        LISTRESTRICT(IRESTRICT) = 0
      ENDDO

C    Output von getion, welche Ionen gebraucht werden
      DO IION = 1, NION      
C       WRITE (0,*) IONNAME(IION)
C       WRITE (0,*) ELEMENTNAME(IION), NUMMER(IION)
        NMAX = NMAX+1
      ENDDO 
C      WRITE (0,*) NMAX
       
      CONTINUE

      
C     Steuerfile einlesen    
      DO !Schleife ueber Zeilen des Steuerfiles
    1   READ (hNFCIN,'(A)', END = 8) INPUTZEILE !Datei zeilenweise durchlesen
        IF (INPUTZEILE == '') GOTO 1 !Leerzeilen abfangen 
        CALL SARGV(INPUTZEILE,1,INPUTKEYWORD)   !Erstes Wort der Zeile (Keyword)
C      einlesen und anschliessend untersuchen (Fallunterscheidung)
C       WRITE (0,*) INPUTKEYWORD
C      Kommentarzeilen ignorieren
        IF (INPUTKEYWORD(:1) == '*' .OR. INPUTKEYWORD(:1) == '-') THEN
C          GOTO 1
          CYCLE
        ELSEIF (INPUTKEYWORD == 'STANDARDPATH') THEN
          CALL SARGV(INPUTZEILE,2,STANDARDPATH) !Standardpfad neu setzen
                                                !(sonst vordefiniert) 
          WRITE (hFC,'(A)') '* Standardpfad: '
          WRITE (hFC,'(2A)') '* ',STANDARDPATH(:IDX(STANDARDPATH))
        ELSEIF (INPUTKEYWORD == 'TESTPATH') THEN
          TESTPATH = ''                     !make sure that testpath will be reset
          CALL SARGV(INPUTZEILE,2,TESTPATH) !Testpfad einlesen
          WRITE (hFC,'(A)') '* Testpfad: '
          WRITE (hFC,'(2A)') '* ',TESTPATH(:IDX(TESTPATH))
        ELSEIF (INPUTKEYWORD == 'NOTESTPATH' .OR. 
     >          INPUTKEYWORD == 'NO_TESTPATH') THEN
          TESTPATH = ''                     !make sure that testpath will be reset
        ELSEIF (INPUTKEYWORD == 'DRTRANSIT') THEN
          bUseDRTRANSIT = .TRUE.
        ELSEIF (INPUTKEYWORD == 'NO-DRTRANSIT') THEN
          bUseDRTRANSIT = .FALSE.
        ELSEIF (INPUTKEYWORD == 'RESTRICT') THEN
C***      Restrict FORMAL_CARDS to specified Element or Ion (or remove restrictions)        
          CALL SARGV(INPUTZEILE,2,WORD) 
          IF (WORD == 'NONE') THEN
            bRESTRICT = .FALSE.
C***        Erase previous restriction list            
            DO IRESTRICT = 1, LASTRESTRICT
              LISTRESTRICT(IRESTRICT) = 0
            ENDDO
            LASTRESTRICT = 0
          ELSEIF (WORD == 'SED') THEN
            bRESTRICT = .TRUE.
C***        Erase previous restriction list            
            DO IRESTRICT = 1, LASTRESTRICT
              LISTRESTRICT(IRESTRICT) = 0
            ENDDO
            LASTRESTRICT = 0
          ELSE
            bRESTRICT = .TRUE.
C***        Read in restriction and add it (if possible)
            CELEM = ''
            IF (LEN_TRIM(WORD) > 2) THEN
C***          Long input => assumes that full element name is given
              DO IATOM=1, NATOM
                IF (WORD == ELEMENT(IATOM)) THEN
                  CELEM = SYMBOL(IATOM)
                  EXIT
                ENDIF
              ENDDO
            ELSE 
C***          Short input => symbol is given
              DO IATOM=1, NATOM
                IF (WORD == SYMBOL(IATOM)) THEN
                  CELEM = SYMBOL(IATOM)
                  EXIT
                ENDIF
              ENDDO
            ENDIF
C***        Stop if no valid/used element was given            
            IF (LEN_TRIM(CELEM) < 1) THEN
              WRITE (hCPR,'(A)')
     >           'ERROR: Invalid or unused element specified!'
              GOTO 990
            ENDIF

            RESTRICTNAME = ''
            RESTRICTNAME = TRIM(CELEM) // '_'
            BIONSTRICT = .FALSE. ! all ions of an element
C***        Is the restriction limited to a certain ion?            
            CALL SARGC(INPUTZEILE, NWORDS)
            IF (NWORDS > 2) THEN
              CALL SARGV(INPUTZEILE,3,CION)
              RESTRICTNAME = TRIM(RESTRICTNAME) // TRIM(CION)
              BIONSTRICT=.TRUE. ! only the given ion
            ENDIF
C***        Fill restriction array with "allowed" ions
C***        The list is generated via name comparison
C***        Name matching for RESTRICT ion names
            IF (BIONSTRICT) THEN
              LENRES=LENIONNAME ! exact ion name match
            ELSE
              LENRES = LEN_TRIM(RESTRICTNAME) ! only element match
            ENDIF
            DO K=1, NMAX
             IF (TRIM(IONNAME(K)(1:LENRES)) == TRIM(RESTRICTNAME)) THEN
                 LASTRESTRICT = LASTRESTRICT + 1
                 LISTRESTRICT(LASTRESTRICT) = K
             ENDIF
            ENDDO
         ENDIF

        ELSEIF (INPUTKEYWORD == 'RANGE') THEN
c          IF (ANFANG .EQ. .TRUE.) THEN
c            WRITE (hFC,'(A)') '* Standardpfad: '
c            WRITE (hFC,'(2A)') '* ',STANDARDPATH(:IDX(STANDARDPATH))
c          ENDIF 
          ANFANG = .FALSE. 

          !Definition der Baender einlesen        
          CALL SARGV(INPUTZEILE,2,UNTERGRENZE)
          READ (UNTERGRENZE,'(F10.0)') MIN               !Typumwandlung       
          CALL SARGV(INPUTZEILE,3,OBERGRENZE)
          READ (OBERGRENZE,'(F10.0)') MAX                !Typumwandlung
          CALL SARGREST(INPUTZEILE,NPAR,4,IFIRST,ILAST)  !Alles ab dem vierten Wort wird der Name
          RANGENAME = INPUTZEILE(IFIRST:ILAST)
          !Allow alias names:
          NALIAS = 0
          alsearch: DO I=IFIRST, ILAST
            IF (INPUTZEILE(I:I+5) == '#ALIAS') THEN
              RANGENAMES(1) = INPUTZEILE(IFIRST:I-1)
              RANGENAME = INPUTZEILE(I+6:ILAST)
              CALL SARGC(RANGENAME,NALIAS)
              DO IALIAS=1, NALIAS
                CALL SARGV(RANGENAME, IALIAS, RANGENAMES(IALIAS+1))
              ENDDO
              NALIAS = NALIAS + 1                        !include original name in counter
              EXIT alsearch
            ENDIF
          ENDDO alsearch
                    
          IF (NALIAS == 0) THEN
            RANGENAMES(1) = RANGENAME
            NALIAS = 1
          ENDIF
          
          WRITE (hCPR,'(A)') ''                     !Leerzeile trennt RANGEs in der CPR-Datei
          WRITE (hCPR,'(A,A)') 'RANGE: ', RANGENAMES(1)  
          WRITE (hCPR,'(A,F10.0)') ' from: ', MIN   
          WRITE (hCPR,'(A,F10.0)') ' to:   ', MAX     
          IF (NALIAS > 1) THEN
            WRITE (hCPR,'(A,$)') ' alias: '
            DO IALIAS=2, NALIAS
              WRITE (hCPR,'(A,$)') RANGENAMES(IALIAS)
            ENDDO
            WRITE (hCPR,'(A)') ''
          ENDIF

          !Ueberpruefen, ob Obergrenze groesser als Untergrenze ist
          IF (MAX < MIN) THEN
            WRITE (hCPR,*) 
     >          'WARNING: lower limit greater then upper limit'
            WRITE (hCPR,*) 'lower limit: ', MIN
            WRITE (hCPR,*) 'upper limit: ', MAX
          ENDIF

          !Oeffnen der Ionendateien        
C***      Achtung: K ist nicht nur Laufvariable, sondern auch file handle
C***               d.h. fort.1 ist Datei zu K = 1 etc.
          DO K = 1, NMAX !Schleife ueber Ionen 
            IF (TESTPATH /= '') THEN
              !Testpfad probieren falls gesetzt
              FILENAME(K) = TESTPATH(:IDX(TESTPATH)) // '/' //
     >                        'FORMAL_CARDS.' // IONNAME(K)
              INQUIRE (FILE=FILENAME(K), EXIST=BEXIST)
            ELSE
              BEXIST = .FALSE.
            ENDIF
            IF (BEXIST) THEN
              OPEN (K,FILE=FILENAME(K), ACTION = 'READ', ERR=901)
            ELSE
              FILENAME(K)=STANDARDPATH(:IDX(STANDARDPATH)) // '/' //
     >                        'FORMAL_CARDS.' // IONNAME(K)
              OPEN (K,FILE=FILENAME(K), ACTION = 'READ', ERR=901)
            ENDIF
            WRITE (hCPR,'(2A)') 'Opened: ', 
     >                   FILENAME(K)(:IDX(FILENAME(K)))  
          ENDDO
          OPENED = .TRUE. 


          !Header fuer das Band schreiben
          WRITE (hFC,'(A)') 
     >          '************************************************'
          DO IALIAS=1, NALIAS
            WRITE (hFC,'(A18,A20)')   
     >          'STRING COMMENT  * ', RANGENAMES(IALIAS)
          ENDDO
          WRITE (hFC,'(A6,F10.0,F10.0)') 'RANGE ', MIN, MAX
          WRITE (hFC,'(A)') 
     >          '************************************************'
          WRITE (hFC,'(A)') 'BLEND'
          WRITE (hFC,'(A)') ''

C      Fuer jedes Ion FORMAL_CARDS-Datei durcharbeiten
          DO K = 1, NMAX !Schleife ueber Ionen
            bSKIPION = .FALSE.
            IF (bRESTRICT) THEN
C***          Restricted output: Check if current ion is in allowed list
              bSKIPION = .TRUE.
              resloop: DO IRESTRICT=1, LASTRESTRICT
                IF (LISTRESTRICT(IRESTRICT) == K) THEN
                  bSKIPION = .FALSE.
                  EXIT resloop
                ENDIF
              ENDDO resloop
            ENDIF
            IF (bSKIPION) CYCLE
            HEADER = .TRUE.
            REWIND(K)
            singleionloop: DO !Schleife fuer ein einzelnes Ion
    4         READ (K,'(A)', END = 7) IONZEILE !Datei zeilenweise durchlesen
              IF (IONZEILE == '') GOTO 4 !Leerzeilen abfangen 
              IONKEYWORD = ''                     !ensure empty keyword variable here
              CALL SARGV(IONZEILE,1,IONKEYWORD)   !Erstes Wort der Zeile (Keyword)
C               einlesen und LINES und MULTIPLETS suchen

C ******* Einzelne Lines **********************************************
              IF ((IONKEYWORD == 'LINE') .OR. 
     >            (IONKEYWORD == '+LINE')) THEN
                !Wellenlaenge steht entweder im dritten Eintrag oder muss berechnet werden    
                CALL SARGV(IONZEILE,3,LAMBDA)             

                !DATOM nach Energien fragen und Wellenlaenge berechnen
                READ (K,'(A)', END = 7) HILFSZEILE !Levelbezeichnungen stehen in der naechsten Zeile
                UPNAME = HILFSZEILE(12:21)
                LOWNAME = HILFSZEILE(34:43)
                !zugehoerige Energien im Datom suchen   
                INDEXL = 0          
                DO I=1, NDIM
                  IF (LEVEL(I) .EQ. LOWNAME) THEN !unteres Level
                    INDEXL=I !Index des richtigen Levels merken
                    LEVELLOW = .TRUE.
                  ENDIF
                ENDDO
                LOWE = ELEVEL(INDEXL) !an diesem Index Energie auslesen
                IF (LEVELLOW == .FALSE.) THEN
                  !Level nicht gefunden
                  WRITE (hCPR,*) 'WARNING: Lowerlevel not found: ',
     >                              LOWNAME
                  GOTO 4
                ENDIF 
            
                INDEXU = 0 
                DO I=1, NDIM
                  IF (LEVEL(I) .EQ. UPNAME) THEN !oberes Level
                    INDEXU=I !Index des richtigen Levels merken
                    LEVELUP = .TRUE.
                  ENDIF
                ENDDO
                UPE = ELEVEL(INDEXU) !an diesem Index Energie auslesen
C              WRITE (0,*) 'Upperlevel: ', UPNAME, UPE
                IF (LEVELUP == .FALSE.) THEN
                  !Level nicht gefunden
                  WRITE (hCPR,*) 'WARNING: Upperlevel not found: ', 
     >                              UPNAME
                  GOTO 4    
                ENDIF

                IF (LAMBDA == '') THEN

                  IF (LEVELLOW .AND. LEVELUP) THEN
                    !Wellenlaenge berechnen              
                    DELTAE = UPE - LOWE 
                    IF (DELTAE == 0.) GOTO 908
                    WL = 1E8/DELTAE !1E8 wegen Kayser und Angstrom
                    !WRITE (0,*) 'Wellenlaenge: ', WL

C Pointer zuruecksetzen, damit dann die richtigen Zeilen geschrieben werden
                    BACKSPACE(K)
                    BACKSPACE(K)
                    READ (K,'(A)', END = 7) IONZEILE
                  ELSE
                    GOTO 4 
                  ENDIF       
              
                ELSE !(means: lambda is given)
                  BACKSPACE(K)
                  READ (K,'(A)', END = 7) IONZEILE2
                  BACKSPACE(K)
                  READ (LAMBDA, '(F12.0)', ERR = 912) WL    !Typumwandlung
                  LAMBDA = '' !wieder auf Default setzen
                ENDIF

C          Pruefen, ob die Wellenlaenge im Range liegt    
                IF (MIN <= WL .AND. MAX >= WL .AND.
     >                 LEVELUP .AND. LEVELLOW) THEN

C          Header fuer das Ion schreiben
                  IF (HEADER .EQ. .TRUE.) THEN
                    WRITE (hFC,'(A)') 
     >               '************************************************'
                    WRITE (hFC,'(A)',advance='no') '*'
                    WRITE (hFC,*)
     >                   ELEMENTNAME(K), NUMMER(K), '(', 
     >                   SYM(K)(:IDX(SYM(K))),' ',
     >                   NUMMER(K)(:IDX(NUMMER(K))), ')'
                    WRITE (hFC,'(A)') 
     >               '************************************************'
                    HEADER = .FALSE. 
                  ENDIF

                  WRITE (hFC,'(A)') IONZEILE(:IDX(IONZEILE))
                  READ (K,'(A)', END = 1) IONZEILE !naechste Zeile auch
                  WRITE (hFC,'(A)') IONZEILE(:IDX(IONZEILE))
                  WRITE (hFC,'(A)') '' !Fuer Uebersichtlickeit
                ENDIF

C               Default wiederherstellen  
                LEVELLOW = .FALSE.
                LEVELUP = .FALSE.

C ******* Multiplets ***********************************************    
              ELSEIF ((IONKEYWORD == 'MULTIPLET') .OR. 
     >            (IONKEYWORD == '+MULTIPLET')) THEN
              
                MULTI(1) = IONZEILE !schreibe die Zeile ins Array
                L = 1
                DO !Schleife ueber Zeilen bis Multiplet-Ende
                  L = L+1 !Zeilen zaehlen
   
    5             READ (K,'(A)', END = 1) IONZEILE
                  IF (IONZEILE == '') GOTO 5 !Leerzeilen abfangen
                  IF (L > MAXMULTI) GOTO 906   !Pruefen, dass die Groesse des Arrays ausreicht
          
                  MULTI(L) = IONZEILE          !schreibe die Zeile ins Array

                  CALL SARGV(IONZEILE,1,IONKEYWORD)

                  IF (IONKEYWORD == 'UPPERLEVEL') THEN
C        Ueberpruefen, ob Levels im DATOM stehen
                    UPPERLEVEL = IONZEILE(12:21)
                    LOWERLEVEL = IONZEILE(34:43)
C                    WRITE (hCPR,*) "Multizeile: ", IONZEILE
C                    WRITE (hCPR,*) "Up, Low: ", UPPERLEVEL, LOWERLEVEL
C   WRITE (0,*) UPPERLEVEL, ' ', LOWERLEVEL
                    DO I=1, NDIM
                      IF (LEVEL(I) == LOWERLEVEL) THEN
                        LEVELLOW = .TRUE.
                      ENDIF 
                      IF (LEVEL(I) == UPPERLEVEL) THEN
                        LEVELUP = .TRUE.
                      ENDIF
                    ENDDO
C                    WRITE (hCPR,*) "Found Up, Low: ", LEVELUP, LEVELLOW
   
                  ELSEIF (IONKEYWORD == '/SUBLINE') THEN      !sublines suchen   
                    LAMBDA = IONZEILE(43:55)                  !Wellenlaenge raussuchen

C         Berechnen der Wellenlaenge, wenn sie nicht da steht
                    IF (LAMBDA == '') THEN
                      UPNAME = IONZEILE(10:19)   
                      LOWNAME = IONZEILE(22:31)
C           WRITE (0,*) UPNAME, LOWNAME 
C         In den bisherigen Multipletzeilen nach Levelnamen suchen
                      DO Z = 1, L
                        HILFSZEILE = MULTI(Z)
                        NAME = HILFSZEILE(13:22)
                        IF (NAME(:IDX(NAME)) == 
     >                                 UPNAME(:IDX(UPNAME))) THEN
                          UPENERGIE = HILFSZEILE(29:48)           !Energie gefunden
                          !WRITE (0,*) 'UPENERGIE:', UPENERGIE
                        ELSEIF (NAME(:IDX(NAME)) == LOWNAME) THEN
                          LOWENERGIE = HILFSZEILE(29:48)          !Energie gefunden
                          !WRITE (0,*) 'LOWENERGIE:', LOWENERGIE
                        ENDIF  
                      ENDDO    

C          Wenn beide Levels aufgepalten sind und Energien da stehen     
                      IF (UPENERGIE /= '' .AND. LOWENERGIE /= '') THEN
C            selbst einlesen und Wellenlaenge berechnen
                        READ (UPENERGIE, '(F10.0)', ERR = 902) UPE
                        READ (LOWENERGIE, '(F10.0)', ERR = 902) LOWE
                        !WRITE (0,*) UPE, ' ', LOWE
                        DELTAE = UPE - LOWE
                        WL = 1E8/DELTAE !1E8 wegen Kayser und Angstrom
                        !WRITE (0,*) WL
                      ELSEIF (UPENERGIE /= '' .AND.
     >                        LOWENERGIE == '') THEN
C                 Nur Energie fuer oberes Level steht da
                        READ (UPENERGIE, '(F10.0)', ERR = 902) UPE
C              WRITE (0,*) 'Upperlevel (bekannt)', UPNAME, UPE
C            Fuer unteres Level DATOM fragen 
C            Index des unteren Levels im LEVEL-Array suchen
                        DO I=1, NDIM
                          IF (LEVEL(I) == LOWNAME) THEN
                            INDEXL=I
                          ENDIF
                        ENDDO
C             An diesem Index im ELEVEL-Array Energie auslesen    
                        LOWE = ELEVEL(INDEXL)
C              WRITE (0,*) 'Lowerlevel:', LOWNAME, LOWE
                        DELTAE = UPE - LOWE
                        WL = 1E8/DELTAE !1E8 wegen Kayser und Angstrom
C              WRITE (0,*) 'Wellenlaenge: ', WL
              
                      ELSEIF (LOWENERGIE /= '' .AND.
     >                        UPENERGIE == '') THEN
C            Nur Energie fuer unteres Level steht da
                        READ (LOWENERGIE, '(F10.0)', ERR = 902) LOWE
C              WRITE (0,*) 'Lowerlevel (bekannt)', LOWNAME, LOWE
C            Fuer oberes Level DATOM fragen 
C            Index des oberen Levels im LEVEL-Array suchen
                        DO I=1, NDIM
                          IF (LEVEL(I) == UPNAME) THEN
                            INDEXU=I
                          ENDIF
                        ENDDO
C             An diesem Index im ELEVEL-Array Energie auslesen
                        UPE = ELEVEL(INDEXU)
C                    WRITE (0,*) 'Upperlevel', UPNAME, UPE
                        DELTAE = UPE - LOWE
                        IF (DELTAE .EQ. 0.) GOTO 913
                        WL = 1E8/DELTAE !1E8 wegen Kayser und Angstrom
C                  WRITE (0,*) 'Wellenlaenge', WL
                      ELSEIF (LOWENERGIE == '' .AND.
     >                        UPENERGIE == '') THEN
                        WRITE (hCPR,*) 
     >                      'WARNING: Sublevels not found: ',
     >                      LOWNAME(:IDX(LOWNAME)), ' ', 
     >                      UPNAME(:IDX(UPNAME))
                      ENDIF 
     
                    ELSE  !zu LAMBDA .EQ. ''

C          Wellenlaenge schon da
C          Absturz verhindern, wenn nach der Wellenlaenge noch ein
C          Argument steht, z.B. VAC, AIR, VOIGT
                      CALL SARGV(LAMBDA,1,LAMBDA1)
                      READ (LAMBDA1, '(F12.0)', ERR = 902) WL !Typumwandlung
C             WRITE (0,*) WL

                    ENDIF !zu LAMBDA .EQ. ''
    
                    IF (MIN <= WL .AND. MAX >= WL .AND.
     >                        LEVELLOW .AND. LEVELUP) THEN       !Vergleich mit Baendergrenzen
                      PASST = .TRUE.

                      !Defaults fuer naechste Zeile
                      WL = 0 
                      UPENERGIE = ''
                      LOWENERGIE = ''   

                    ELSEIF (MIN <= WL .AND. MAX >= WL .AND. 
     >                  (.NOT. LEVELLOW) .AND. LEVELUP )  THEN
                      WRITE (hCPR,*) 
     >                 'WARNING: Lowerlevel not found in DATOM: ', 
     >                    LOWERLEVEL
                    ELSEIF (MIN <= WL .AND. MAX >= WL .AND. 
     >                  LEVELLOW  .AND. (.NOT. LEVELUP) )  THEN
                      WRITE (hCPR,*) 
     >                 'WARNING: Upperlevel not found in DATOM: ', 
     >                    UPPERLEVEL
                    ELSEIF (MIN <= WL .AND. MAX >= WL .AND. 
     >                  (.NOT. LEVELLOW) .AND. (.NOT. LEVELUP) )  THEN
                      WRITE (hCPR,*) 
     >                 'WARNING: Upperlevel not found in DATOM: ', 
     >                    UPPERLEVEL
                      WRITE (hCPR,*) 
     >                 'WARNING: Lowerlevel not found in DATOM: ', 
     >                    LOWERLEVEL
                    ENDIF !Ende vom Vergleich mit Baendergrenzen
   
                  ELSEIF (IONKEYWORD == '-MULTIPLET') THEN  !Abbruchbed.
                    GOTO 6 !Ende des Multiplets
                  ENDIF !Ende von IONKEYWORD .EQ. '/SUBLINE'
                ENDDO !Ende des Multiplets
  
    6           LMAX = L !Zahl der Eintrage im Multipletarray merken
                IF (PASST) THEN !eine Wellenlaenge im Band
                  PASST = .FALSE.
C             Header fuer das Ion schreiben
                  IF (HEADER .EQ. .TRUE.) THEN
                    WRITE (hFC,'(A)') 
     >                '************************************************'
                    WRITE (hFC,'(A)',advance='no') '*'
                    WRITE (hFC,*)
     >                  ELEMENTNAME(K), NUMMER(K), '(', 
     >                  SYM(K)(:IDX(SYM(K))),' ',
     >                  NUMMER(K)(:IDX(NUMMER(K))), ')'
                    WRITE (hFC,'(A)') 
     >                '************************************************'
                    HEADER = .FALSE. 
                  ENDIF  

                  DO L = 1, LMAX
                    WRITE (hFC,'(A)') MULTI(L)(:IDX(MULTI(L)))  !Alle Eintraege schreiben
                  ENDDO
                  WRITE (hFC,'(A)') ''      !Leerzeile zur Uebersichtlickeit
                ENDIF
          
C ******* Default wiederherstellen fuer naechstes Multiplet ********
                UPENERGIE = ''
                LOWENERGIE = ''
                LAMBDA = ''         
                DO I=1, MAXMULTI
                  MULTI(I) = ''
                ENDDO
                PASST = .FALSE.
                INDEXL = 0
                INDEXU = 0
                LEVELUP = .FALSE.
                LEVELLOW = .FALSE.
                UPPERLEVEL = ''
                LOWERLEVEL = ''
  
C ******* DRTRANSIT (Achtung! Erstmal nur kopiertes Multiplet!) ***********

              ELSEIF ((IONKEYWORD == 'DRTRANSIT') .OR. 
     >            (IONKEYWORD == '+DRTRANSIT')) THEN
              
                DRTRANSIT(1) = IONZEILE !schreibe die Zeile ins Array
                L = 1
                DO !Schleife ueber Zeilen bis DRTRANSIT-Ende
                  L = L+1 !Zeilen zaehlen
    9             READ (K,'(A)', END = 1) IONZEILE
                  IF (IONZEILE == '') GOTO 9      !Leerzeilen abfangen
                  IF (L > MAXDRTRANSIT) GOTO 909  ! Pruefen, dass die Groesse des Arrays ausreicht
                  DRTRANSIT(L) = IONZEILE         !schreibe die Zeile ins Array        

                  CALL SARGV(IONZEILE,1,IONKEYWORD)
                  IF (IONKEYWORD .EQ. '/ADDLINE') THEN  !sublines suchen   
                    LAMBDA = IONZEILE(43:55)            !Wellenlaenge raussuchen
 
                    !Wellenlaenge schon da
                    READ (LAMBDA, '(F12.0)', ERR = 902) WL !Typumwandlung
C              WRITE (0,*) WL            
 
                    IF (MIN <= WL .AND. MAX >= WL) THEN !Vergleich mit Baendergrenzen
                      PASST = .TRUE.

                      !Defaults fuer naechste Zeile 
                      WL = 0 
                      UPENERGIE = ''
                      LOWENERGIE = ''   
                  
                    ENDIF !Ende vom Vergleich mit Baendergrenzen

                  ELSEIF (IONKEYWORD == '-DRTRANSIT') THEN !Abbruchbed.
                    GOTO 10 !Ende des DRTRANSITS
                  ENDIF !Ende von IONKEYWORD .EQ. '/ADDLINE'
                ENDDO !Ende des DRTRANSITS
  
   10           LMAX = L !Zahl der Eintrage im DRTRANSITarray merken
                IF (PASST .AND. bUseDRTRANSIT) THEN !eine Wellenlaenge im Band
                  PASST = .FALSE.
  
C          Header fuer das Ion schreiben
                  IF (HEADER) THEN
                    WRITE (hFC,'(A)') 
     >                '************************************************'
                    WRITE (hFC,'(A)',advance='no') '*'
                    WRITE (hFC,*)
     >                   ELEMENTNAME(K), NUMMER(K), '(', 
     >                   SYM(K)(:IDX(SYM(K))),' ',
     >                   NUMMER(K)(:IDX(NUMMER(K))), ')'
                    WRITE (hFC,'(A)')
     >                '************************************************'
                    HEADER = .FALSE. 
                  ENDIF  
 
                  DO L = 1, LMAX
                    WRITE (hFC,'(A)') DRTRANSIT(L)(:IDX(DRTRANSIT(L))) !Alle Eintraege schreiben
                  ENDDO
                  WRITE (hFC,'(A)') '' !Fuer Uebersichtlickeit
                ENDIF

C ******* Default wiederherstellen fuer naechstes DRTRANSIT ********
                UPENERGIE = ''
                LOWENERGIE = ''
                LAMBDA = ''         
                DO I=1, MAXDRTRANSIT
                  DRTRANSIT(I) = ''
                ENDDO
                PASST = .FALSE.
                INDEXL = 0
                INDEXU = 0

              ENDIF !Ende der Unterscheidung zw. Lines, Multplets und DRTRANSIT

            ENDDO singleionloop !Ende der Schleife fuer ein Ion

   7        CONTINUE
            PASST = .FALSE.

          ENDDO !Ende der Schleife ueber Ionen (K = 1, NMAX)
          WRITE (hFC,'(A)') '-BLEND'
          WRITE (hFC,'(A)') '' !Fuer Uebersichtlichkeit

        ELSE !zusaetzliche Optionen mitkopieren
          WRITE (hFC,'(A)') INPUTZEILE(:IDX(INPUTZEILE))

        ENDIF !Ende der Fallunterscheidung fuer INPUTKEYWORD     

      ENDDO !Ende der Schleife fuer Steuerfile-Zeilen
 
C      WRITE (99,'(A)') '-BLEND'
    8 CLOSE (hFC) !Schliessen der output-Datei
      CLOSE (hNFCIN) !Schliessen des Input files
C     Schliessen der Ionendateien  
      DO K=1, NMAX
        CLOSE (K)
      END DO    

C**   Regular program end - remove error code 
   99 CONTINUE
      CALL JSYMSET ('G0', '0')
      STOP 'O.K.'

C********  ERROR EXITS **************************************
  900 WRITE(0,*) 'ERROR: PATH to atomic data not (yet) defined'
      GOTO 990
      
  901 WRITE(0,*) 'ERROR: Cannot open file:'
      WRITE(0,*) FILENAME(K)(:IDX(FILENAME(K)))
      GOTO 990    
      
  902 WRITE(0,*) 'ERROR: String kann nicht in real umgewandelt werden'
      GOTO 992     
  
  903 WRITE(0,*) 'ERROR: NEWFORMAL_CARDS_INPUT NOT FOUND'
      GOTO 990         
      
  904 WRITE(0,*) 'ERROR: CANNOT WRITE IN FILE FORMAL_CARDS'
      GOTO 990     
      
  905 WRITE(0,*) 'ERROR: DATOM NOT FOUND'
      GOTO 990       

  906 WRITE(0,*) 'ERROR: MORE SUBLINES AND -LEVELS THAN DIMENSIONED'
      WRITE(0,*) 'Multiplet-Block hat zu viele Zeilen (L=', L, ')'
      GOTO 992    
      
  907 WRITE (0,*) 'More ions found in DATOM than dimensioned'
      WRITE (0,*) 'MAXIONNAMES =', MAXIONNAMES
      GOTO 990          
      
  908 WRITE (0,*) 'Delta E =0'
      GOTO 991
      
  909 WRITE(0,*) 'ERROR: MORE SUBLINES AND -LEVELS THAN DIMENSIONED'
      WRITE(0,*) 'DRTRANSIT-Block hat zu viele Zeilen(L=', L, ')'
      GOTO 993             
        
  912 WRITE(0,*) 'ERROR: String kann nicht in real umgewandelt werden'
      GOTO 991  
      
  913 WRITE (0,*) 'Delta E =0'
      GOTO 992                     

*Allgemeiner Fehlerausgang            
  990 WRITE(0,*) '*** The error occured in the following line:'
      WRITE(0,*) INPUTZEILE(:IDX(INPUTZEILE))
      WRITE(0,*) '*** Line in FORMAL_CARDS-SUBFILE was:'
      WRITE(0,*) IONZEILE(:IDX(IONZEILE))
      WRITE(0,*) '*** Ion: ', IONNAME(K)
      STOP '*** ERROR DETECTED IN PROGRAM NEWFORMAL_CARDS'

*Fehlerausgang fuer einzelne Lines      
  991 WRITE(0,*) '*** The error occured in the following line:'
      WRITE(0,*) INPUTZEILE(:IDX(INPUTZEILE))
      WRITE(0,*) '*** Line in NEWFORMAL_CARDS-FILE was:'
      WRITE(0,*) IONZEILE(:IDX(IONZEILE))
      WRITE(0,*) '*** Line was:'
      WRITE(0,*) IONZEILE2(:IDX(IONZEILE2))
      WRITE(0,*) '*** Ion: ', IONNAME(K)
      STOP '*** ERROR DETECTED IN PROGRAM NEWFORMAL_CARDS' 

*Fehlerausgang fuer Multiplets            
  992 WRITE(0,*) '*** The error occured in the following line:'
      WRITE(0,*) INPUTZEILE(:IDX(INPUTZEILE))
      WRITE(0,*) '*** Line in NEWFORMAL_CARDS-FILE was:'
      WRITE(0,*) IONZEILE(:IDX(IONZEILE))
      WRITE(0,*) '*** Second Line of MULTIPLET was'
      WRITE(0,*) MULTI(2)
      WRITE(0,*) '*** Ion: ', IONNAME(K)
      STOP '*** ERROR DETECTED IN PROGRAM NEWFORMAL_CARDS'   
      
*Fehlerausgang fuer DRTRANSIT           
  993 WRITE(0,*) '*** The error occured in the following line:'
      WRITE(0,*) INPUTZEILE(:IDX(INPUTZEILE))
      WRITE(0,*) '*** Line in NEWFORMAL_CARDS-FILE was:'
      WRITE(0,*) IONZEILE(:IDX(IONZEILE))
      WRITE(0,*) '*** Second Line of DRTRANSIT was'
      WRITE(0,*) DRTRANSIT(2)
      WRITE(0,*) '*** Ion: ', IONNAME(K)
      STOP '*** ERROR DETECTED IN PROGRAM NEWFORMAL_CARDS'               
      
      
      END
