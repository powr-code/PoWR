      SUBROUTINE FINDELEMENT (SYMBOL, NAME, STAGE, ATOMICMASS)
C*****************************************************************
C***  Hier wird einem Elementsymbol der Elementname, Masse und Stage
C***  zugeordnet.
C*****************************************************************

      CHARACTER SYMBOL*(*)
      PARAMETER (MAXELEM = 26)
      CHARACTER*10 ELEMNAME(MAXELEM)
      CHARACTER*2 ELEMSYMBOL(MAXELEM)
      REAL ELEMMASSE(MAXELEM)
      REAL ELEMSTAGE(MAXELEM)
C      CHARACTER*5 ELEMMASSE(MAXELEM)
C      CHARACTER*3 ELEMSTAGE(MAXELEM)
      INTEGER Z
      CHARACTER*10 NAME

      
      DATA ELEMSYMBOL /'H ', 'HE', 'LI', 
     2                 'BE', 'B ', 'C ', 
     3                 'N ', 'O ', 'F ',
     4                 'NE', 'NA', 'MG', 
     5                 'AL', 'SI', 'P ',
     6                 'S ', 'CL', 'AR', 
     7                 'K ', 'CA', 'SC', 
     8                 'TI', 'V ', 'CR', 
     9                 'MN', 'G '/

      
      DATA ELEMNAME /'HYDROGEN  ', 'HELIUM    ', 'LITHIUM   ', 
     2               'BERYLLIUM ', 'BORON     ', 'CARBON    ', 
     3               'NITROGEN  ', 'OXYGEN    ', 'FLUORINE  ',
     4               'NEON      ', 'SODIUM    ', 'MAGNESIUM ', 
     5               'ALUMINIUM ', 'SILICON   ', 'PHOSPHORUS',
     6               'SULFUR    ', 'CHLORINE  ', 'ARGON     ', 
     7               'POTASSIUM ', 'CALCIUM   ', 'SCANDIUM  ', 
     8               'TITANIUM  ', 'VANADIUM  ', 'CHROMIUM  ', 
     9               'MANGANESE ', 'GENERIC   '/

C    Fuer H, HE, C, N, O, SI, P, NE sind hier Werte eingetragen, die ich aus
C    DATOM-Files zusammengesucht habe, der Rest ist aus dem
C    Periodensystem abgeschrieben.   
      DATA ELEMMASSE / 1.00, 4.00, 6.94, 
     2               9.01, 10.81, 12.00, 
     3               14.00, 16.00, 19.00,
     4               20.18, 22.99, 24.31, 
     5               26.98, 28.09, 30.97,
     6               32.07, 35.45, 39.95, 
     7               39.10, 40.08, 44.96, 
     8               47.87, 50.94, 52.00, 
     9               54.94, 0.00 /
     
     
C    Fuer H, HE, C, N, O, SI, P, NE sind hier Werte eingetragen, die ich aus
C    DATOM-Files zusammengesucht habe, der Rest ist willkuerlich auf 2 gesetzt     
      DATA ELEMSTAGE / 2. ,3. ,2. , 
     2                 2. ,2. ,5. , 
     3                 6. ,7. ,2. ,
     4                 2. ,2. ,2. , 
     5                 2. ,5. ,5. , 
     6                 2. ,2. ,2. ,  
     7                 2. ,2. ,2. , 
     8                 2. ,2. ,2. , 
     9                 2. ,2. /   
     

      
      DO Z=1, MAXELEM
         IF (SYMBOL .EQ. ELEMSYMBOL(Z)) THEN
           WRITE(NAME,'(A)') ELEMNAME(Z)
C	   READ (ELEMSTAGE(Z), FMT='(F5.0)', ERR=900) STAGE
C	   READ (ELEMMASSE(Z), FMT='(F6.2)', ERR=900) MASSE
	   STAGE = ELEMSTAGE(Z)
	   ATOMICMASS = ELEMMASSE(Z)
           EXIT
         ENDIF
      END DO
      
      RETURN 
      
  900  PRINT *, 'ERROR: ' 
       STOP '*** ERROR DETECTED' 

      END
