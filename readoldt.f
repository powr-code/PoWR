      SUBROUTINE READOLDT (ICH, ND, NDDIM, T, RADIUS,
     $                     TOLD, ROLD, MODOLD, JOBNOLD, TEFF, TEFFOLD,
     >                     TAURCONT, TAURCONTOLD, BTAUR, DTDRIN_OLD)
C**********************************************************************
C***  CALLED FROM: WRSTART
C***  READS TEMPEREATURE STRUCTURE FROM OLD MODEL FILE, IF REQUESTED
C***  AS THE RADIUS GRIDS MAY DEVIATE, THE NEW TEMPERATURE STRUCTURE
C***    IS OBTAINED BY INTERPOLATION
C**********************************************************************

      DIMENSION T(ND), RADIUS(ND), TOLD(ND), ROLD(ND)
      DIMENSION TAURCONT(ND), TAURCONTOLD(NDDIM)
      LOGICAL BTAUR, BTAUR_INTERPO

      BTAUR_INTERPO = BTAUR

C***  READ FROM CHANNEL ICH = OLD MODEL FILE ****************************
      CALL OPENMS (ICH, IDUMMY, IDUMMY, 1, IERR)
      IERR=1
      CALL READMS (ICH,MODOLD,13,'MODHEAD ',IERR)
      IF (IERR .LT. 0) THEN
         CALL REMARK (' OLD MODEL FILE NOT AVAILABLE')
         PRINT *,     ' OLD MODEL FILE NOT AVAILABLE'
         STOP 'ERROR'
         ENDIF
      CALL READMS (ICH,JOBNOLD,     1, 'JOBNUM  ' , IERR)
      CALL READMS (ICH,NDOLD  ,     1, 'ND      ' , IERR)
C***  ARRAY BOUND CHECK
      IF (NDOLD .GT. NDDIM) THEN
         CALL REMARK (' OLD MODEL HAS TOO MANY DEPTH POINTS')
         PRINT *,      'OLD MODEL HAS TOO MANY DEPTH POINTS'
         STOP 'ERROR'
      ENDIF

      CALL READMS (ICH,TEFFOLD,     1, 'TEFF    ' , IERR)
      CALL READMS (ICH,TOLD   , NDOLD, 'T       ' , IERR)
      CALL READMS (ICH,ROLD   , NDOLD, 'R       ' , IERR)

      IF (BTAUR) THEN
         CALL READMS (ICH,TAURCONTOLD, NDOLD, 'TAURCONT' , IERR)
         IF (IERR .EQ. -10) THEN
            WRITE (0,*) 
     >       '*** WARNING: TAURCONT not on MODEL file, take TAUROSS'
            CALL READMS (ICH,TAURCONTOLD, NDOLD, 'TAUROSS ' , IERR)
            IF (IERR .EQ. -10) THEN
               WRITE (0,*) 
     >          '*** WARNING: TAUROSS not on MODEL file,' 
     >          // ' TAU option in OLD T disabled'
               BTAUR_INTERPO = .FALSE.
            ENDIF        
         ENDIF
      ENDIF

      CALL READMS (ICH,DTDRIN_OLD,1,   'DTDRIN  ' , IERR)
      ! DTDRIN_OLD only used if OLD T + OLD STRAT/V is set

      CALL CLOSMS (ICH, IERR)


C***  INTERPOLATION ***************************************************
      BTAUR_INTERPO = BTAUR_INTERPO .AND. (TAURCONT(ND) .NE. 0.)

      IF (BTAUR_INTERPO) THEN

         WRITE (0,*) 'Interpolation of Temperature on Tau-Grid'

cC***     scale TAUROLD
c         DO L=1, NDOLD
c            TAURCONTOLD(L) = 
c     >         TAURCONTOLD(L) * (TAURCONT(ND)+1.E-8)/TAURCONTOLD(NDOLD)
c         ENDDO
c Auskommentiert, und testweise durch das nachstehende IF ersetzt:
c     wrh  4-Mar-2019

         DO L=1, ND
            TAURL = TAURCONT(L)
            IF (TAURL .GT. TAURCONTOLD(NDOLD)) THEN
               T(L) = TOLD(NDOLD)
            ELSE             
               CALL LIPO (T(L),TAURL,TOLD,TAURCONTOLD,NDOLD)
            ENDIF
         ENDDO

      ELSE

         WRITE (0,*) 'Interpolation of Temperature on Radius-Grid'
         DO 1 L=1, ND
            IF (RADIUS(L) .GT. ROLD(1)) THEN
               T(L)=TOLD(1)
            ELSE
               CALL LIPO (T(L),RADIUS(L),TOLD,ROLD,NDOLD)
            ENDIF
    1    CONTINUE 

      ENDIF

C***  SCALE TEMPERATURE WITH THE RATIO TEFF ( NEW / OLD )
      IF (BTAUR_INTERPO) THEN
      IF (TEFF .NE. TEFFOLD) THEN
         Q = TEFF / TEFFOLD
         DO L=1, ND
         FTAU = EXP(-TAURCONT(L)*2.)
C         FTAU = 1.-MIN(TAURCONT(L), 1.)
         T(L) = T(L) * (FTAU + (1.-FTAU) * Q)
         ENDDO
      ENDIF

      ELSE
      IF (TEFF .NE. TEFFOLD) THEN
         Q = TEFF / TEFFOLD
         DO 2 L=1, ND
    2    T(L) = T(L) * Q
         ENDIF

      ENDIF

      RETURN
      END
