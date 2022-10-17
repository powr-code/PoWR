      SUBROUTINE PRELINE (NUP, LOW, IND, N, XLAM, ND, XJLMEAN,
     $                    ELEVEL, NDIM, INDNUP, INDLOW)
C*******************************************************************************
C***  Preparing some quantities for the considered line transition
C***  Called from Main Programm COLI
C*******************************************************************************

C     By Lars Koesterke,  7-Aug-1997 19:08:34

C***  Parameter
C***  In:

C*    N          aktuelle Zahl der Level
C*    ND         aktuelle Zahl der Tiefenpunkte
C*    LINE       Liste aller zu berechnender Linien
C*    ELEVEL     Energien aller Niveaus
C*    NL         Index der aktuelle Linie
C*    NDIM       Maximale Zahl der Levels
cC*    EINST      Einsteinkoefizienten und Rudimental-kennung
C*    INDNUP     Zu jeder Linie das obere  Niveau
C*    INDLOW     Zu jeder Linie das untere Niveau
C*    LASTIND    Anzahl aller Linien (auch der nicht zu berechnenden)

C***  Out:
C*    NUP        Index des oberen  Niveaus der aktuellen Linie
C*    LOW        Index des unteren Niveaus der aktuellen Linie
C*    IND        Nummer der aktuellen Linie
C*    XJLMEAN    Strahlungsfeldintegral ueber die Linie an jedem Tiefenpunkt

 
      DIMENSION ELEVEL(N)
      DIMENSION INDNUP(NDIM),INDLOW(NDIM)
      DIMENSION XJLMEAN(ND)
      LOGICAL BLASERL

C***  FIND THE LEVEL INDICES NUP, LOW
      NUP=INDNUP(IND)
      LOW=INDLOW(IND)

C***  WAVELENGTH OF PRESENT LINE:
      XLAM=1.E8/(ELEVEL(NUP)-ELEVEL(LOW))
 
C***  INITIALIZE XJLMEAN FOR LATER INTEGRATION 
      DO L=1,ND
        XJLMEAN(L)=.0
      ENDDO
 
      RETURN
      END
