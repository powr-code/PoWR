      SUBROUTINE CLUMP_STRUCT (DENSCON, FILLFAC, ND, DENSCON_FIX, VELO, 
     >                         TAUROSS, DENSCON_LINE, RADIUS, T, XMU)

C************************************************************************
C***  This routine allows to specify a density stratification of 
C***  the clumping constrast, DENSCON(L). 
C***  The subroutine interpretes directly the DENSCON-line from the  
C***  CARDS file, which is handed over from DECSTAR via the character
C***  variable DENSCON_LINE. 
C***  Default is depth-independent clumping with DENSCON_FIX, 
C***  which is the second parameter on the DENSCON-line and decoded 
C***  already in Subr. DECSTAR
C************************************************************************

      IMPLICIT NONE
      INCLUDE 'interfacebib.inc'

      REAL, DIMENSION(ND) :: DENSCON, FILLFAC, VELO, DCSMOOTH,
     >                       TAUROSS, RADIUS, T, XMU  !NOTE: Tauross is Tauross_cont
      CHARACTER(10) :: CLUMP_CRIT, CLUMP_CRIT2 
      CHARACTER(20) :: ACTPAR
      CHARACTER*(*) :: DENSCON_LINE

C***  Local arrays      
      INTEGER, PARAMETER :: NDIPMAX = 100
      REAL, DIMENSION(NDIPMAX) :: AMACH, VHELP 
      
      INTEGER :: L, IPAR, ND, NPAR, NDSTART, Lcand, LENPAR
      REAL :: DENSCON_FIX, PAR1, PAR2, PARL, D2, F1, F8, X1, X2, 
     >        DX, X, Q, Rsonic, Vsonic, TAUsonic
      
      !File and channel handles (=KANAL)
      INTEGER, PARAMETER :: hOUT = 6        !write to wruniqX.out (stdout)
      INTEGER, PARAMETER :: hCPR = 0        !write to wruniqX.cpr (stderr)

      REAL, PARAMETER :: PI = 3.141592654
      REAL, PARAMETER :: RGAS = 8.3145E7            !GAS CONSTANT in CGS 


      IF (ND > NDIPMAX) THEN
         WRITE (hCPR,'(A)') 'CLUMP_STRUCT: FATAL ERROR ******'
         WRITE (hCPR,'(A)') 'CLUMP_STRUCT: NDIPMAX INSUFFICIENT'
         WRITE (hCPR,'(2(A,I4))') 'ND = ', ND, ', NDIPMAX = ', NDIPMAX
         STOP 'FATAL ERROR IN CLUMP_STRUCT'
      ENDIF
      
C***  Decoding the input line with DENSCON specifications
      CALL SARGC (DENSCON_LINE, NPAR)
C***  Clumping not depth-dependent?
      IF (NPAR .LE. 2) THEN
         DO L = 1, ND
         DENSCON(L) = DENSCON_FIX
         FILLFAC(L) = 1. / DENSCON_FIX
         ENDDO
         GOTO 100
      ENDIF

      CALL SARGV (DENSCON_LINE, 3, CLUMP_CRIT) 


      Lcand = 0
      Rsonic = 1.
      Vsonic = VELO(ND)
      DO L=1, ND
C***    Calculate speed of sound for all depth points)
        AMACH(L) = SQRT(RGAS * T(L) / XMU(L)) / 1.E5
        VHELP(L) = VELO(L) - AMACH(L)
        IF ((VHELP(L) < 0.) .AND. (Lcand == 0)) THEN
          Lcand = L
        ENDIF
      ENDDO                
C***  Find sonic point parameters
      IF (Lcand > 1) THEN
        CALL SPLINPOX(Rsonic,0.,RADIUS,VHELP,ND,.FALSE.,Lcand)
        CALL SPLINPOX(Vsonic,Rsonic,AMACH,RADIUS,ND)
        CALL SPLINPOX(TAUsonic,Rsonic,TAUROSS,RADIUS,ND)
      ENDIF

C***  Branch to simulate Hillier's formula, given in Martins et al. 
C***  (2004, A&A 420, 1087; SMC-Paper)       
      IF (CLUMP_CRIT == 'HILLIER') THEN
         IF (NPAR .LT. 4) GOTO 94
            CALL SARGV (DENSCON_LINE, 4, ACTPAR) 
            LENPAR = LEN_TRIM(ACTPAR)
            IF (ACTPAR(LENPAR:LENPAR) == 'S') THEN
C***          Parameter is interpreted as a fraction of the sonic speed            
              READ (ACTPAR(1:LENPAR-1), '(F20.0)', ERR=97) PAR1
              PAR1 = PAR1 * Vsonic
            ELSEIF (ACTPAR == 'SONIC') THEN
              PAR1 = Vsonic
            ELSE
              READ (ACTPAR, '(F20.0)', ERR=97) PAR1
            ENDIF
            IF (PAR1 <= .0) GOTO 95 
            F8 = 1. / DENSCON_FIX
            DO L = 1, ND
               FILLFAC(L) = F8 + (1.-F8)* EXP(-VELO(L)/PAR1)
               DENSCON(L) = 1. / FILLFAC(L)
            ENDDO
         GOTO 100
      ENDIF

C***  Branch to simulate Paco's (F. Najarro's) FASTWIND clumping strat.
C***  from Najarro et al. (2009, ApJ 691, 1816)
C***  SYNTAX:   DENSCON [D1] PACO V1 D2 V2
C***      e.g.  DENSCON 10 PACO 100. 4  300.
C***            DENSCON 10 PACO 2.5  1. 2.0  (more typical use)
      IF (CLUMP_CRIT == 'NAJARRO' .OR. CLUMP_CRIT == 'PACO') THEN
         IF (NPAR < 6) GOTO 93
         CALL SARGV (DENSCON_LINE, 4, ACTPAR) 
         LENPAR = LEN_TRIM(ACTPAR)
         IF (ACTPAR(LENPAR:LENPAR) == 'S') THEN
C***        Parameter is interpreted as a fraction of the sonic speed            
            READ (ACTPAR(1:LENPAR-1), '(F20.0)', ERR=97) PAR1
            PAR1 = PAR1 * Vsonic
         ELSEIF (ACTPAR == 'SONIC') THEN
            PAR1 = Vsonic
         ELSE
            READ (ACTPAR, '(F20.0)', ERR=93) PAR1
         ENDIF
         IF (PAR1 <= .0) GOTO 95 

C***     Read outermost clumping factor (reached only in the limit)
         CALL SARGV (DENSCON_LINE, 5, ACTPAR)
         READ (ACTPAR, '(F20.0)', ERR=93) D2

         CALL SARGV (DENSCON_LINE, 6, ACTPAR) 
         LENPAR = LEN_TRIM(ACTPAR)
         IF (ACTPAR(LENPAR:LENPAR) == 'S') THEN
C***        Parameter is interpreted as a fraction of the sonic speed            
            READ (ACTPAR(1:LENPAR-1), '(F20.0)', ERR=97) PAR2
            PAR2 = PAR2 * Vsonic
         ELSEIF (ACTPAR == 'SONIC') THEN
            PAR2 = Vsonic
         ELSE
            READ (ACTPAR, '(F20.0)', ERR=93) PAR2
         ENDIF
         IF (PAR2 <= .0) GOTO 95 

         F1 = 1. / DENSCON_FIX
         F8 = 1. / D2
         DO L = 1, ND
            FILLFAC(L) = F1 + (1.-F1)*EXP(-VELO(L)/PAR1) 
     >                      + (F8-F1)*EXP((VELO(L)-VELO(1))/PAR2)
            DENSCON(L) = 1. / FILLFAC(L)
         ENDDO

         GOTO 100
      ENDIF

C***  Exponential clumping onset (like in Hillier's formula)
C***  but using the radius scale
      IF (CLUMP_CRIT == 'EXPRADIUS') THEN
         IF (NPAR .LT. 4) GOTO 94
            CALL SARGV (DENSCON_LINE, 4, ACTPAR) 
            LENPAR = LEN_TRIM(ACTPAR)
            IF (ACTPAR(LENPAR:LENPAR) == 'S') THEN
C***          Parameter is interpreted as a fraction of the 
C***          distance of the sonic radius from RSTAR
              READ (ACTPAR(1:LENPAR-1), '(F20.0)', ERR=97) PAR1
              PAR1 = 1. + PAR1 * (Rsonic - 1.)
            ELSEIF (ACTPAR == 'SONIC') THEN
              PAR1 = Rsonic
            ELSE
              READ (ACTPAR, '(F20.0)', ERR=97) PAR1
            ENDIF
            IF (PAR1 <= .0) GOTO 95 
            F8 = 1. / DENSCON_FIX
            DO L = 1, ND
               FILLFAC(L) = F8 + (1.-F8)*EXP(-(RADIUS(L)-1.)/(PAR1-1.))
               DENSCON(L) = 1. / FILLFAC(L)
            ENDDO
         GOTO 100
      ENDIF

C***  Exponential clumping onset (like in Hillier's formula)
C***  but using the Tauross_cont scale
      IF (CLUMP_CRIT == 'EXPTAU') THEN
         IF (NPAR .LT. 4) GOTO 94
            CALL SARGV (DENSCON_LINE, 4, ACTPAR) 
            LENPAR = LEN_TRIM(ACTPAR)
            IF (ACTPAR(LENPAR:LENPAR) == 'S') THEN
C***          Parameter is interpreted as a fraction of the 
C***          sonic Tau value
              READ (ACTPAR(1:LENPAR-1), '(F20.0)', ERR=97) PAR1
              PAR1 = PAR1 * TAUsonic
            ELSEIF (ACTPAR == 'SONIC') THEN
              PAR1 = TAUsonic
            ELSE
              READ (ACTPAR, '(F20.0)', ERR=97) PAR1
            ENDIF
            IF (PAR1 <= .0) GOTO 95 
            F8 = 1. / DENSCON_FIX
            DO L = 1, ND
              IF (TAUROSS(L) <= 1.E-30) THEN
                FILLFAC(L) = F8
              ELSE
                FILLFAC(L) = F8 + (1.-F8)*EXP(-PAR1/TAUROSS(L))
              ENDIF
              DENSCON(L) = 1. / FILLFAC(L)
            ENDDO
         GOTO 100
      ENDIF
      
      
C***  Stratification needs enough parameters!
      IF (NPAR .LT. 5) GOTO 98
      CALL SARGV (DENSCON_LINE, 4, ACTPAR) 
      IF (CLUMP_CRIT == 'VELO' .AND. ACTPAR == 'SONIC') THEN
        !Clumping starts at sonic point        
        PAR1 = Vsonic
      ELSEIF (CLUMP_CRIT == 'TAU' .AND. ACTPAR == 'SONIC') THEN
        CALL SPLINPOX(PAR1,Rsonic,TAUROSS,RADIUS,ND)
      ELSEIF (CLUMP_CRIT == 'RADIUS' .AND. ACTPAR == 'SONIC') THEN
        PAR1 = Rsonic
      ELSEIF (CLUMP_CRIT == 'TAU' .OR. CLUMP_CRIT == 'VELO' 
     >                        .OR. CLUMP_CRIT == 'RADIUS') THEN
        READ (ACTPAR, '(F20.0)', ERR=97) PAR1
      ENDIF
      CALL SARGV (DENSCON_LINE, 5, ACTPAR) 
      IF (CLUMP_CRIT == 'VELO' .AND. ACTPAR == 'SONIC') THEN
C***    Clumping ends at sonic point
        PAR2 = Vsonic
      ELSEIF (CLUMP_CRIT == 'TAU' .AND. ACTPAR == 'SONIC') THEN
        CALL SPLINPOX(PAR2,Rsonic,TAUROSS,RADIUS,ND)
      ELSEIF (CLUMP_CRIT == 'RADIUS' .AND. ACTPAR == 'SONIC') THEN
        PAR2 = Rsonic
      ELSEIF (CLUMP_CRIT == 'TAU' .OR. CLUMP_CRIT == 'VELO' 
     >                           .OR. CLUMP_CRIT == 'RADIUS') THEN
        READ (ACTPAR, '(F20.0)', ERR=97) PAR2
      ENDIF

C***  Depth dependent Clumping Factor
         
C***  Clumping increases from VDENS1 to VDENS2


C***  Branch: Clumping scales with Rosseland optical depth
      IF (CLUMP_CRIT .EQ. 'TAU') THEN
         X1 = AMIN1 (PAR1, PAR2)
         X1 = AMAX1 (TAUROSS(1), X1)
         X2 = AMAX1 (PAR1, PAR2)
         X2 = AMIN1 (TAUROSS(ND), X2)
         DX = X2 - X1
         DO L=1, ND
            IF (TAUROSS(L) .LE. X1) THEN
               DENSCON(L) = DENSCON_FIX
            ELSE IF (TAUROSS(L) .GE. X2) THEN
               DENSCON(L) = 1.
            ELSE
               X = PI * (TAUROSS(L) - X1 ) / DX
               Q = 0.5 + 0.5 * COS(X)
               DENSCON(L) = (1.-Q) + Q * DENSCON_FIX
            ENDIF
            FILLFAC(L) = 1. / DENSCON(L)
         ENDDO

C***  Branch: Clumping scales with velocity
      ELSE IF (CLUMP_CRIT .EQ. 'VELO') THEN
         X1 = AMIN1 (PAR1, PAR2)
         X1 = AMAX1 (VELO(ND), X1)
         X2 = AMAX1 (PAR1, PAR2)
         X2 = AMIN1 (VELO(1), X2)
         DX = X2 - X1
         DO L=1, ND
            IF (VELO(L) .LE. X1) THEN
               DENSCON(L) = 1.
            ELSE IF  (VELO(L) .GE. X2) THEN
               DENSCON(L) = DENSCON_FIX
            ELSE
               X = PI * (VELO(L) - X1 ) / DX
               Q = 0.5 + 0.5 * COS(X)
               DENSCON(L) = Q + (1.-Q) * DENSCON_FIX
            ENDIF
            FILLFAC(L) = 1. / DENSCON(L)
         ENDDO

      ELSE IF (CLUMP_CRIT == 'RADIUS') THEN
         X1 = MIN(PAR1, PAR2)
         X1 = MAX(RADIUS(ND), X1)
         X2 = MAX(PAR1, PAR2)
         X2 = MIN(RADIUS(1), X2)
         DX = X2 - X1
         DO L=1, ND
            IF (RADIUS(L) <= X1) THEN
               DENSCON(L) = 1.
            ELSE IF  (RADIUS(L) >= X2) THEN
               DENSCON(L) = DENSCON_FIX
            ELSE
               X = PI * (RADIUS(L) - X1 ) / DX
               Q = 0.5 + 0.5 * COS(X)
               DENSCON(L) = Q + (1.-Q) * DENSCON_FIX
            ENDIF
            FILLFAC(L) = 1. / DENSCON(L)
         ENDDO
         
C***  The LIST criterion allows multiple clumping steps
C***  which can be defined by velocity, optical depth or radius
C***  The first value before the LIST keyword defines the outermost clumping factor
C***  The keyword after LIST defines the criterion for the following list which
C***  can be of arbitrary length, but must follow in (parval dval)-pairs.
C***   
C***  SYNTAX:  DENSCON [Dout] LIST TAU|VELO|RADIUS p1 D1 p2 D2 p3 D3 ...
C***  Example: DENSCON 10 LIST TAU 0.1 1. 0.005 20. 
C***       means that (seen from R_star outwarts)
C***           - clumping starts with D=1. until TAUROSS_cont drops below 0.1
C***           - increases to D=20. until TAUROSS_cont drops below 0.005
C***           - decreases to D=10 in the outermost wind (TAUROSS_cont < 0.005)
      ELSE IF (CLUMP_CRIT == 'LIST') THEN
        DO L=1, ND
          DENSCON(L) = DENSCON_FIX
        ENDDO
        CALL SARGV (DENSCON_LINE, 4, CLUMP_CRIT2) 
        IF (CLUMP_CRIT2 == 'VELO') THEN
          PARL = 0.
          NDSTART = ND
          DO IPAR=5, NPAR, 2
            CALL SARGV (DENSCON_LINE, IPAR, ACTPAR) 
            IF (ACTPAR == 'SONIC') THEN
              PAR1 = Vsonic
            ELSE
              READ (ACTPAR, '(F20.0)', ERR=97) PAR1           
            ENDIF
            IF (PAR1 < PARL) THEN
              WRITE (hCPR,*) 'ERROR: ' 
              WRITE (hCPR,*) 'Velocity steps not in monotonic order'
              GOTO 99
            ENDIF
            CALL SARGV (DENSCON_LINE, IPAR+1, ACTPAR) 
            READ (ACTPAR, '(F20.0)', ERR=97) PAR2
            ndloop: DO L=NDSTART, 1, -1
              IF (VELO(L) <= PAR1) THEN
                DENSCON(L) = PAR2
              ELSE
                PARL = PAR1
                NDSTART = L
                EXIT ndloop
              ENDIF
            ENDDO ndloop
          ENDDO
        ELSEIF (CLUMP_CRIT2 == 'TAU') THEN
          PARL = TAUROSS(ND) * 1.1
          NDSTART = ND
          DO IPAR=5, NPAR, 2
            CALL SARGV (DENSCON_LINE, IPAR, ACTPAR) 
            READ (ACTPAR, '(F20.0)', ERR=97) PAR1           
            IF (PAR1 > PARL) THEN
              WRITE (hCPR,'(A)') ' ERROR: ' 
              WRITE (hCPR,'(A)') ' Tau steps not in decreasing order'
              WRITE (hCPR,'(2(A,F8.3))') 'Last = ',PARL, ' Next = ',PAR1
              GOTO 99
            ENDIF
            CALL SARGV (DENSCON_LINE, IPAR+1, ACTPAR) 
            READ (ACTPAR, '(F20.0)', ERR=97) PAR2
            ndloop2: DO L=NDSTART, 1, -1
              IF (TAUROSS(L) >= PAR1) THEN
                DENSCON(L) = PAR2
              ELSE
                PARL = PAR1
                NDSTART = L
                EXIT ndloop2
              ENDIF
            ENDDO ndloop2
          ENDDO
        ELSEIF (CLUMP_CRIT2 == 'RADIUS') THEN
          PARL = 0.
          NDSTART = ND
          DO IPAR=5, NPAR, 2
            CALL SARGV (DENSCON_LINE, IPAR, ACTPAR) 
            IF (ACTPAR == 'SONIC') THEN
              PAR1 = Rsonic
            ELSE
              READ (ACTPAR, '(F20.0)', ERR=97) PAR1           
            ENDIF
            IF (PAR1 < PARL) THEN
              WRITE (hCPR,*) 'ERROR: ' 
              WRITE (hCPR,*) 'Radius steps not in monotonic order'
              GOTO 99
            ENDIF
            CALL SARGV (DENSCON_LINE, IPAR+1, ACTPAR) 
            READ (ACTPAR, '(F20.0)', ERR=97) PAR2
            ndloop3: DO L=NDSTART, 1, -1
              IF (RADIUS(L) <= PAR1) THEN
                DENSCON(L) = PAR2
              ELSE
                PARL = PAR1
                NDSTART = L
                EXIT ndloop3
              ENDIF
            ENDDO ndloop3
          ENDDO
        ELSE
          WRITE (hCPR,*) 'ERROR: ' 
          WRITE (hCPR,*)
     >      'Invalid choice of clumping-structure criterion: ', 
     >       CLUMP_CRIT2
          GOTO 99
        ENDIF
        !Smoothing of clumping grid
        DCSMOOTH(ND) = DENSCON(ND)
        DCSMOOTH(1) = DENSCON(1)
        DO L=ND-1, 2, -1
          DCSMOOTH(L) = DENSCON(L) * 0.6 +
     >          DENSCON(L+1) * 0.2 + DENSCON(L-1) * 0.2
        ENDDO
        DO L=1, ND
          DENSCON(L) = DCSMOOTH(L)
          FILLFAC(L) = 1. / DENSCON(L)
        ENDDO
      ELSE
C***  Error: invalid choice of Clumping Criterion
         GOTO 96
      ENDIF
      
  100 RETURN


C***  ERROR branches *********************************************

  93  WRITE (0,*) 'ERROR: '
      WRITE (0,*) 'Najarro formula needs three additional parameters'
      GOTO 99

  94  WRITE (0,*) 'ERROR: ' 
      WRITE (0,*) 'Hillier formula needs one parameter: V_cl ' 
      GOTO 99

  95  WRITE (0,*) 'ERROR: ' 
      WRITE (0,*) 'Hillier/Najarro velocity parameters must be > 0'
      GOTO 99

  96  WRITE (0,*) 'ERROR: ' 
      WRITE (0,*) 'Invalid choice of Clumping-structure criterion: ', 
     >            CLUMP_CRIT       
      GOTO 99

  97  WRITE (0,*) 'ERROR: ',  ACTPAR, 
     >      ' could not be decodes as floating-point number!'
      GOTO 99    

  98  WRITE (0,*) 'ERROR: Clumping-structure criterion needs two', 
     >            'parameters'
      GOTO 99

  99  WRITE (0,*) 'The error occured in the following input line:'
      WRITE (0,*) DENSCON_LINE
      STOP 'FATAL ERROR detected by SUBROUTINE CLUMP_STRUCT'

      END


