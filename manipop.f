      SUBROUTINE MANIPOP (ENLTE, WEIGHT, NCHARG, EION, ELEVEL, NOM, 
     >              ABXYZ, NFIRST, NLAST, NATOM, POPNUM, RNE, ENTOT, 
     >              N, ND, MANIPOP_OPTIONS, DENSCON, T, level)
C************************************************************************
C***  Manipulation of popnumbers, in order to simulate the emission of 
C***  a hot component. 
C************************************************************************

      DIMENSION RNE(ND), ENTOT(ND), POPNUM(ND,N), ENLTE(N), T(ND)
      PARAMETER ( MAXPAR = 6 )
      CHARACTER MANIPOP_OPTIONS*(*), PAR(MAXPAR)*20, level(N)*(*)
      LOGICAL BTCOOL
      DIMENSION DENSCON(ND)

      BTCOOL = .FALSE.

C***  THOT    = Temperature of the hot component (optional 'TCOOL') 
C***  DENSHOT = Density of the hot component, relative to the 
C***            density of the cool component including the clumping factor
C***  HOTMASS = Mass fraction of the hot material

      CALL SARGC (MANIPOP_OPTIONS, NPAR)      
      IF (NPAR .EQ. 0) RETURN

      DO I=1, NPAR
         CALL SARGV (MANIPOP_OPTIONS, I, PAR(I))
      ENDDO       

      THOT    = -1. 
      DENSHOT = -1.
      HOTMASS = -1.
      DO I=1, NPAR-1
        IF (PAR(I) .EQ. 'THOT'   )  THEN
           IF (PAR(I+1) .EQ. 'TCOOL') THEN
              BTCOOL = .TRUE.
           ELSE
              READ (PAR(I+1),'(G20.0)') THOT
           ENDIF
        ENDIF
        IF (PAR(I) .EQ. 'DENSHOT') READ (PAR(I+1),'(G20.0)') DENSHOT
        IF (PAR(I) .EQ. 'HOTMASS') READ (PAR(I+1),'(G20.0)') HOTMASS
      ENDDO
      IF (THOT .LE. 0. .AND. .NOT. BTCOOL) THEN 
         WRITE (0,'(A)') '*** INVALID OR MISSING PARAMETER THOT'
         STOP '*** FATAL ERROR IN SUBROUTINE MANIPOP'
      ENDIF 
      IF (DENSHOT .LE. 0.) THEN 
         WRITE (0,'(A)') '*** INVALID OR MISSING PARAMETER DENSHOT'
         STOP '*** FATAL ERROR IN SUBROUTINE MANIPOP'
      ENDIF 
      IF (HOTMASS .LE. 0.) THEN 
         WRITE (0,'(A)') '*** INVALID OR MISSING PARAMETER HOTMASS'
         STOP '*** FATAL ERROR IN SUBROUTINE MANIPOP'
      ENDIF 

      DO L=1, ND
        ENE = RNE(L) * ENTOT(L) * DENSCON(L) * DENSHOT 
        IF (BTCOOL) THOT = T(L)
        
C***    No manipulation if 'HOT' component not hotter than the cool one 
        IF (THOT .LT. T(L)) CYCLE

        CALL LTEPOP (N,ENLTE,THOT,ENE,WEIGHT,NCHARG,EION,ELEVEL,NOM,
     $                  ABXYZ,NFIRST,NLAST,NATOM)

        DO I=1, N
           POPNUM(L,I) = POPNUM(L,I) + HOTMASS * ENLTE(I)
        ENDDO

cccc    test output
c        if (l .eq. 30 ) then
c           write (0, '(A,I3)') '*** L =', L
c           do i=1,n
c            write (0, '(A,X,G14.4)') level(i), enlte(i)
c           enddo
c        endif

      ENDDO

      RETURN
      END
