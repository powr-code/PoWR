      SUBROUTINE REMOCO (ND, NDDIM, RADIUS, ENTOT,T , RNE, NF, NFDIM,
     $                  XLAMBDA, FWEIGHT, KEY, XJC, 
     $                  N, POPNUM, RSTAR,
     $                  MODHEAD, LAST, MAXHIST, MODHIST, JOBNUM,
     $                  NOTEMP, TEFF, HEDDI, EDDI, MAXXDAT, XDATA, 
     >                  DENSCON, FILLFAC, OPARND, ABXYZ, NATOM, 
     >                  ZERO_RATES, NDIM)
C*******************************************************************************
C***  READING THE MODEL FILE FOR MAIN PROGRAM "COMO"
C*******************************************************************************
 
      DIMENSION XJC(2), KEY(2),HEDDI(2),EDDI(2)
      DIMENSION ABXYZ(NATOM)
      DIMENSION XDATA(MAXXDAT)
      INTEGER :: ND, NDDIM, NF, NFDIM
      REAL, DIMENSION(NDDIM) :: RADIUS, ENTOT, T
      LOGICAL NOTEMP
      LOGICAL ZERO_RATES (NDIM*NDDIM)
C*** tiefenabh. clumping nach goetz
      DIMENSION DENSCON(NDDIM),FILLFAC(NDDIM)

C *** NOTE: ASSURE 64-BIT TYPE FOR USE IN WRITMS
      CHARACTER NAME*8   

      CALL OPENMS (3, IDUMMY, IDUMMY, 1, IERR)
      CALL READMS (3,ND,1,      'ND      ', IERR)
      IF(ND.GT.NDDIM) THEN
            CALL REMARK ('TOO MANY DEPTH POINTS')
            STOP 'ERROR IN REMOCO'
            ENDIF
      CALL READMS (3,RADIUS,ND, 'R       ', IERR)
      CALL READMS (3,ENTOT,ND,  'ENTOT   ', IERR)
      DENSCON(1:ND) = 0.0
      CALL READMS(3,DENSCON,ND, 'DENSCON ', IERR)
      IF (IERR .EQ. -10) DENSCON(1:ND) = 1.
      DO L=1,ND
         IF (DENSCON(L) .LE. 0. ) THEN
            IF (DENSCON(1) .LE. 0. ) THEN
               CALL REMARK ('Zero or negative depth-dep. clumping!')
               WRITE(0,*)'Error in depth-dep. clumping during COMO!'
               STOP 'Error in depth-dep. clumping during COMO!'
            ENDIF
            DENSCON(L) = DENSCON(1)
            FILLFAC(L) = 1. / DENSCON(L)
         ELSE
            FILLFAC(L) = 1. / DENSCON(L)
         ENDIF
      ENDDO
      CALL READMS (3,T,ND,      'T       ', IERR)
      CALL READMS (3,RNE,ND,    'RNE     ', IERR)
      CALL READMS (3,NF,1,      'NF      ', IERR)
      IF(NF.GT.NFDIM) THEN
            CALL REMARK ('TOO MANY FREQUENCY POINTS')
            STOP 'ERROR'
            ENDIF
      CALL READMS (3,XLAMBDA,NF,'XLAMBDA ', IERR)
      CALL READMS (3,FWEIGHT,NF,'FWEIGHT ', IERR)
      CALL READMS (3,KEY,NF,    'KEY     ', IERR)
      CALL READMS (3,POPNUM,ND*N,'POPNUM  ', IERR)
      CALL READMS (3,RSTAR,1,   'RSTAR   ', IERR)


C***  READ 'XDATA' AND CHECK WHETHER THE RECORD EXISTS
      IERR=1
      CALL READMS (3,XDATA,MAXXDAT,'XDATA   ',IERR)
      IF (IERR .LT. 0) THEN
         CALL REMARK ('ERROR WHEN READING XDATA FROM MODEL FILE')
C***     XFILL EQ 0.
         XDATA(1) = 0.  
      ENDIF
 
C***  Flags for the POPMIN levels
      CALL READMS (3,ZERO_RATES,  N*ND, 'ZERO_RAT', IERR)
C*    Default if variable does not exist yet
      IF (IERR .EQ. -10) THEN
        DO I=1, N*ND
         ZERO_RATES(I) = .FALSE.
        ENDDO
      ENDIF

      DO 15 K=1,NF
      WRITE (NAME,'(A3,I4,A1)') 'XJC',K, ' '
   15 CALL READMS (3,XJC(1+ND*(K-1)),ND,NAME, IERR)
  
C***  IF TEMPERATURE-CORRECTIONS APPLY: READ TEFF FROM MODEL FILE 
      IF (.NOT. NOTEMP) THEN
C***        CAUTION: TEFF RECORD MAY NOT EXIST IN VERY OLD BERLIN MODELS
            IERR=1
            CALL READMS (3,TEFF,1,'TEFF    ',IERR)
            IF (IERR .LT. 0) THEN
                CALL REMARK (' TEFF NOT ON MODEL FILE: USE OLDSTART')
                PRINT *,     ' TEFF NOT ON MODEL FILE: USE OLDSTART'
                STOP 'ERROR'
                ENDIF
C***  READ HEDDI = EDDINGTON FACTOR AT INNER BOUNDARY
         DO 14 K=1,NF
         IF (K <= 999) THEN
           WRITE (NAME,'(A4,I3,A1)') 'EDDI',K,' '
         ELSE
           WRITE (NAME,'(A4,I4)') 'EDDI', K
         ENDIF
         CALL READMS (3,EDDI,3*ND,NAME, IERR)
   14    HEDDI(K)=EDDI(3*ND)

         ENDIF

      CALL READMS (3,MODHEAD,13,'MODHEAD ', IERR)
      CALL READMS (3,LAST,1,    'MODHIST ', IERR)
      IF (LAST.GT.MAXHIST) THEN
            CALL REMARK ('MODEL HISTORY DIMENSION TOO SMALL')
            STOP 'ERROR'
            ENDIF
      CALL READMS (3,MODHIST,LAST,'MODHIST ', IERR)
 
      OPARND = 0.
      CALL READMS (3, OPARND, 1, 'OPARND  ', IERR)

      CALL READMS (3,JOBNUM,1,    'JOBNUM  ', IERR)
      JOBNUM=JOBNUM+1
C      IF (JOBNUM .GE. 1000) JOBNUM=JOBNUM-1000
      CALL WRITMS (3,JOBNUM,1,    'JOBNUM  ',-1, IDUMMY, IERR)

C***  READ 'ABXYZ' (REL. ABUNDANCES OF ALL ELEMENTS) AND CHECK WHETHER
C***  THE RECORD EXISTS
      IERR=1
      CALL READMS (3,ABXYZ,NATOM,'ABXYZ   ',IERR)
      IF ((IERR .LT. 0) .AND. (IERR .NE. -10)) THEN
         CALL REMARK ('ERROR WHEN READING ABXYZ FROM MODEL FILE')
         STOP 'ERROR'
      ENDIF

      RETURN
      END
