      SUBROUTINE PLOTANF (KANAL,NPLOT,NHEAD,NX,NY,
     $ B1,B2,B3,B4,B5,B6,C1,C2,C3,C4,C5,C6,
     $ X,Y,N,ISYMBOL )
C***********************************************************************
C***  DIESE ROUTINE BEREITET EIN NEUES PLOTFILE ZUR UEBERTRAGUNG VOR
C***  DIE ANGEGEBENEN WERTE BESCHREIBEN DEN PLOTKASTEN
C***********************************************************************
      CHARACTER*(*) NPLOT,NHEAD,NX,NY

      WRITE (KANAL,1) NPLOT
      WRITE (KANAL, '(A)') 'KASDEF FONT=HELVET'
      WRITE (KANAL,2) NHEAD
      WRITE(KANAL,3) NX
      WRITE (KANAL,4) NY
    1 FORMAT (' PLOT   :',A)
    2 FORMAT (' HEADER :',A)
    3 FORMAT (' X-ACHSE:',A)
    4 FORMAT (' Y-ACHSE:',A)

C***  Special branch activating Lars' AUTO axes: set XMIN = XMAX!
      IF (B2 .EQ. B3) THEN
         WRITE (KANAL, 6)
    6   FORMAT (5X,
     $ 'MASSTAB    MINIMUM    MAXIMUM    TEILUNGEN  BESCHRIFT. DARUNTER'
     $  / ,' X: AUTO',/,' Y: ')
      ELSE IF (C2 .EQ. C3) THEN
        WRITE (KANAL,7) B1,B2,B3,B4,B5,B6
    7   FORMAT (5X,
     $ 'MASSTAB    MINIMUM    MAXIMUM    TEILUNGEN  BESCHRIFT. DARUNTER'
     $  / ,' X: ',6(G12.6,1X),/,' Y: AUTO')
      ELSE 
        WRITE (KANAL,5) B1,B2,B3,B4,B5,B6,C1,C2,C3,C4,C5,C6
    5   FORMAT (5X,
     $ 'MASSTAB    MINIMUM    MAXIMUM    TEILUNGEN  BESCHRIFT. DARUNTER'
     $  / ,' X: ',6(G12.6,1X),/,' Y: ',6(G12.6,1X))
      ENDIF

C***  No Dataset opened for Zero Length
      IF (N .GT. 0) THEN
        CALL PLOTTAB (KANAL,X,Y,N,ISYMBOL)
      ELSE
        WRITE (KANAL, '(A)') 'END'
      ENDIF

      RETURN
      END
