      SUBROUTINE TABREAD (ND,RADIUS,VELO,GRADI,T,ITAB)
C***********************************************************************
C***  READS TEMPERATURE AND/OR VELOCITY STRUCTURE FROM FILE TAPE8=TABLE ********
C***  ITAB = 0: NO TABLE INPUT (DEFAULT)
C***  ITAB = 1: INPUT OF TABULATED TEMPERATURE STRUCTURE T(R)
C***  ITAB = 2: INPUT OF TABULATED VELOCITY FIELD V(R)
C***  ITAB = 3: INPUT OF T(R) AND V(R)
C***********************************************************************
 
      REAL RADIUS(ND),VELO(ND),GRADI(ND),T(ND)
      REAL R(70),V(70),G(70),TH(70),IP(70)
      CHARACTER KARTE*80
 
      DO 10 L=1,70
   10 IP(L)=0
 
      OPEN (8,FILE='TABLE', STATUS='UNKNOWN')
    1 READ(8,2, END=99) KARTE
    2 FORMAT(A)

      IF (KARTE(:1) .EQ. '*' ) GOTO 1
      IF (KARTE(:5).EQ.'TABLE') THEN
         DECODE (80,11,KARTE) NRP,VTAB,TTAB
   11    FORMAT (15X,I3,7X,F7.1,7X,F8.0)
         IF (NRP.LE.0) THEN
            CALL REMARK ('NO TABLE INPUT')
            STOP 'ERROR'
            ENDIF
         IF (NRP.GT.70) THEN
            CALL REMARK ('NR TAB PTS .GT. 70')
            STOP 'ERROR'
            ENDIF
         GOTO 1
      ENDIF
      IF (KARTE(:7).EQ.'ENDGRID') GOTO 30
 
C***  DECODE INPUT
      DECODE (80,12,KARTE) L,R(L),V(L),G(L),TH(L)
   12 FORMAT (I5,4F10.5)
      IP(L)=1
      GOTO 1
   30 CONTINUE
      CLOSE (8)
C***  END INPUT
 
      PRINT 31, VTAB,TTAB
   31 FORMAT (1H1,//,20X,' TABULAR INPUT OF THE VELOCITY AND/OR TEMPERATURE
     $ URE STRUCTURE',/,20X,63(1H=),//,10X,
     $ 'VELOCITY UNIT =',F8.1,10X,'TEMPERATURE UNIT=',F9.0,//,10X,
     $                                    '  NO   RADIUS    VELOCITY   GRAD
     $ RADIENT   TEMPERATURE',//)
      DO 33 L=1,NRP
      PRINT 32, L,R(L),V(L),G(L),TH(L)
   32 FORMAT (10X,I4,4F11.5)
      IF (IP(L).NE.1) THEN
         CALL REMARK ('GRID POINT MISSING')
         STOP 'ERROR'
         ENDIF
   33 CONTINUE
 
C***  INTERPOLATION OF THE TABLE GRID
C***  TEMPERATURE STRUCTURE T(R)
      IF (TTAB.GT.0.0) THEN
      ITAB=1
      DO 51 L=1,ND
      REF=RADIUS(L)
      CALL LIPO (F,REF,TH,R,NRP)
   51 T(L)=F*TTAB
      END IF
C***  VELOCITY FIELD V(R)
      IF (VTAB.GT.0.0) THEN
      IF (ITAB .GT. 0) THEN
         ITAB=3
      ELSE
         ITAB=2
      ENDIF
      DO 41 L=1,ND
      REF=RADIUS(L)
      CALL LIPO (F,REF,V,R,NRP)
      VELO(L)=F*VTAB
      CALL LIPO (F,REF,G,R,NRP)
      GRADI(L)=F*VTAB
   41 CONTINUE
      END IF
 
      RETURN

C***  ERROR EXIT:
   99 CONTINUE
      CALL REMARK ('NO ENDGRID FOUND')
      STOP 'ERROR'

      END
