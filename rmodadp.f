      SUBROUTINE RMODADP (NDDIM, OLDHEAD, N, NOLD, NDIM, 
     >             NATOM, MODHEAD, ND, NDOLD, ABXYZ,LAST,MODHIST,
     >             RADIUS,ROLD,POPHELP, POPNUM, TAURCONT, 
     >             TAURCONTOLD, POPLTE, POPLTE_OLD, BTAUR)
C**********************************************************************
C***  ... READS OLD AND NEW MODEL DATA FOR PROGRAM 'ADAPTER'
C***  this routine also interpolates the population numbers 
C***    for the radius grid of the new model  
C**********************************************************************

      LOGICAL BTAUR


C***  READING OF THE OLD MODEL FILE  -----------------------------------
      IERR=0
      CALL OPENMS (9,IDUMMY, IDUMMY,1,IERR)
      IF ((IERR .NE. 0) .AND. (IERR .NE. -5)) THEN
         CALL REMARK ('ADAPTER: ERROR WHEN OPENING OLD MODEL FILE')
         STOP 'ERROR'
      ENDIF

      CALL READMS (9,NDOLD,1,        'ND      ', IERR)
      IF (NDOLD .GT. NDDIM) THEN
         CALL REMARK ('TOO MANY DEPTH POINTS IN OLD MODEL FILE')
         STOP 'ERROR'
         ENDIF

C***  Use MS-routine STORAGE to find out the number of levels
      CALL LENGTHMS (9, NDNOLD, 'POPNUM  ', IERR)
      NOLD = NDNOLD / NDOLD 
c      write (0,*) 'NUMBER OF LEVELS IN OLD MODEL ATOM:', NOLD
      IF (NOLD .GT. NDIM) THEN
         WRITE (0,*) 'OLD MODEL HAS MORE LEVELS THAN DIMENSIONED'
         WRITE (0,'(A,I5)') 'AVAILABLE NDIM:', NDIM 
         WRITE (0,'(A,I5)') 'REQUIRED  NDIM:', NOLD
         STOP 'ERROR in ADAPTER:RMODADP'
      ENDIF 
      CALL READMS (9,OLDHEAD,13,     'MODHEAD ', IERR)
      CALL READMS (9,POPHELP,     NDNOLD,  'POPNUM  ', IERR)
      CALL READMS (9,POPLTE_OLD , NDNOLD,  'POPLTE  ', IERR)
      CALL READMS (9,ROLD   ,     NDOLD,   'R       ', IERR)

      IF (BTAUR) THEN
         CALL READMS (9,TAURCONTOLD, NDOLD, 'TAURCONT' , IERR)
         IF (IERR .EQ. -10) THEN
            WRITE (0,*) 
     >       '*** WARNING: TAURCONT not on MODEL file, take TAUROSS'
            CALL READMS (9,TAURCONTOLD, NDOLD, 'TAUROSS ' , IERR)
            IF (IERR .EQ. -10) THEN
               WRITE (0,*) 
     >          '*** WARNING: TAUROSS not on MODEL file,' 
     >          // ' TAU option in OLDSTART disabled'
               BTAUR = .FALSE.
            ENDIF        
         ENDIF
      ENDIF

      CALL CLOSMS (9, IERR)

C***  READING OF THE NEW MODEL FILE  -----------------------------------
      CALL OPENMS (3, IDUMMY, IDUMMY, 1, IERR)
      CALL READMS (3,ND,1,'ND      ', IERR)
      IF (ND .GT. NDDIM) THEN
         CALL REMARK ('INSUFFICIENT DIMENSION: NDDIM')
         STOP 'ERROR'
         ENDIF
      NDN = ND * N
      CALL READMS (3, MODHEAD, 13   , 'MODHEAD ', IERR)
      CALL READMS (3, POPNUM , NDN  , 'POPNUM  ', IERR)
      CALL READMS (3, POPLTE , NDN  , 'POPLTE  ', IERR)
      CALL READMS (3, ABXYZ  , NATOM, 'ABXYZ   ', IERR)
      CALL READMS (3, LAST   , 1    , 'MODHIST ', IERR)
      CALL READMS (3, MODHIST, LAST , 'MODHIST ', IERR)
      CALL READMS (3, RADIUS , ND   , 'R       ', IERR)
      CALL READMS (3, TAURCONT,ND   , 'TAURCONT', IERR)

      RETURN
      END
