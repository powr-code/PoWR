      SUBROUTINE PLOTRTAU1 (NF, XLAMBDA, K, ND, RADIUS, OPA, MAINPRO, 
     >                     MAINLEV, MODHEAD)
C***********************************************************************
C***  PLOT of the radius where taucont=1 ("Nick-Knatterton-Diagram")
C***  CALLED FROM: COMO
C***  Option: PLOT RTAU1
C***  wrh, 27-Aug-2002 
C***********************************************************************

      DIMENSION OPA(ND), RADIUS(ND), XLAMBDA(NF) 
      CHARACTER*10 MAINPRO(ND), MAINLEV(ND), LASTMAINPRO, LASTMAINLEV
      CHARACTER HEAD1*60, HEAD2*60, MODHEAD*100
      PARAMETER (MAXIDENT = 100)
      CHARACTER*45 IDLINE(MAXIDENT)
 
C***  INITIALIZATION
      IF (K .EQ. 1) THEN
         RTAU1LAST = 1.0
         NIDENT=0
         KANAL=2
         OPEN (KANAL, FILE='PLOT', STATUS='UNKNOWN')
         CALL JSYMSET ('G2','TRANSFER')
         CALL REMARK ('RTAU1-PLOT TO BE ROUTED')

C         XMIN=ALOG10(XLAMBDA(1))
C         XMAX=ALOG10(XLAMBDA(NF))
C***     AUTO-SCALING by WRplot
         XMIN = .0
         XMAX = .0
 
         XTICK = 0.1
         XABST = 1.
         YTICK = 10.
         YABST = 50.
         CALL PLOTANF (KANAL,'RTAU1', '&E'//MODHEAD
     >        ,'\CENTER\log #l#/\A'
     >        ,'\CENTER\Radius where #t# = 1'
     >        ,0., XMIN, XMAX, XTICK, XABST,.0
     >        ,0., 1.,   RADIUS(1), YTICK, YABST, .0
     >        ,DUMMY, DUMMY, 0, 0)
         BACKSPACE(KANAL)
         WRITE (KANAL, '(A)') 'N=? XYTABLE COLOR=2 PEN=3' 
      ENDIF

      TAU=.0
      RTAU1 = 1.
      LTAU1 = ND 

      DO 6 L=1, ND-1
         TAUOLD = TAU
         TAU=TAU+0.5*(OPA(L)+OPA(L+1))*(RADIUS(L)-RADIUS(L+1))
         IF (TAU .GE. 1.) THEN
            Q = (1.-TAUOLD) / (TAU-TAUOLD)
            RTAU1 = (1.-Q) * RADIUS(L) + Q * RADIUS(L+1)
            WRITE(KANAL, '(F12.4, 1X, F12.4)') 
     >           ALOG10(XLAMBDA(K)), RTAU1
            LTAU1 = L + 1            
            EXIT
         ENDIF
    6 CONTINUE

      IF (K. GT. 1) THEN
        IF (LASTMAINPRO .NE. MAINPRO(LTAU1) 
     >       .OR. LASTMAINLEV .NE. MAINLEV(LTAU1) ) THEN 
ccc     >       .OR. RTAU1 .LT. RTAU1LAST*0.95) THEN 
          IF (NIDENT .LT. MAXIDENT) THEN
           NIDENT = NIDENT + 1
           WRITE (IDLINE(NIDENT), '(A, 1X, 1PG14.6, 1X, A)')
     >      '\IDENT', XLAMBDA(K-1), '&E'//LASTMAINPRO//' '//LASTMAINLEV 
          ENDIF
        ENDIF
      ENDIF
      LASTMAINPRO = MAINPRO(LTAU1)
      LASTMAINLEV = MAINLEV(LTAU1)
      RTAU1LAST = RTAU1

      IF (K .EQ. NF) THEN 
         WRITE (KANAL, '(A)') 'FINISH'

         WRITE (KANAL, '(A)') ' '
         WRITE (KANAL, '(A)') '* Identification of Edges'
         WRITE (KANAL, '(A)') '\IDY 10'
         WRITE (KANAL, '(A)') '\IDSIZE 0.2'
         WRITE (KANAL, '(A)') '\IDLOG'
         DO I=1, NIDENT
            WRITE (KANAL, '(A)') IDLINE(I)
         ENDDO
         WRITE (KANAL, '(A)') 'END'
      ENDIF



      RETURN
      END
