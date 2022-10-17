      SUBROUTINE FORMOSA (ND, R, NP, P, Z, ENTOT, RNE, ABXYZ, NATOM, 
     >                 T,VELO,NF,XLAMBDA,GRADI,
     >        POPNUM, RSTAR, VDOP, JOBNUM, N, NDDIM, NPDIM, NFDIM,
     >        MODHEAD,TEFF,MAXXDAT,XDATA, XJC, IMOD, 
     >        DENSCON, FILLFAC, TAURCONT, ZERO_RATES, RCON, NDIM, XMDOT)
C***********************************************************************
C***  CALLED FROM FORMAL, READS THE MODEL FILE
C***********************************************************************
 
      DIMENSION XDATA(MAXXDAT), ABXYZ(NATOM)
      DIMENSION XJC(2),EDDI(2),xlambda(2)
      CHARACTER*8 NAME, MODHEAD*(*)
      DIMENSION DENSCON(ND), FILLFAC(ND), TAURCONT(ND)
      LOGICAL, DIMENSION(NDIM*NDDIM) :: ZERO_RATES
      REAL :: RCON

      M = IMOD - 1

      CALL OPENMS(3+M, IDUMMY, IDUMMY, 1, IERR)
      CALL READMS(3+M, ND,1,          'ND      ', IERR)
      IF (ND .GT. NDDIM) THEN
            CALL REMARK ('TOO MANY DEPTH POINTS')
            STOP 'ERROR'
            ENDIF
      CALL READMS (3+M, R,ND,         'R       ', IERR)
      CALL READMS (3+M, NP,1,         'NP      ', IERR)
      IF (NP.GT.NPDIM) THEN
            CALL REMARK ('TOO MANY IMPACT-PARAMETER POINTS')
            STOP 'ERROR'
            ENDIF
      CALL READMS (3+M, P,NP,         'P       ', IERR)
      CALL READMS (3+M, Z,ND*NP,      'Z       ', IERR)
      CALL READMS(3+M, ENTOT,ND,      'ENTOT   ', IERR)
      CALL READMS(3+M, TAURCONT,ND,   'TAURCONT', IERR)
      DENSCON(1:ND) = 0.0
      CALL READMS(3+M,DENSCON,ND, 'DENSCON ', IERR)
      IF (IERR .EQ. -10) DENSCON(1:ND) = 1.
      DO L=1,ND
         IF (DENSCON(L) .LE. 0. ) THEN
            IF (DENSCON(1) .LE. 0. ) THEN
               CALL REMARK ('Zero or negative depth-dep. clumping!')
               WRITE(0,*)'Error in depth-dep. clumping during FORMAL!'
               STOP 'Error in depth-dep. clumping during FORMALLARGECL!'
            ENDIF
            DENSCON(L) = DENSCON(1)
            FILLFAC(L) = 1. / DENSCON(L)
         ELSE
            FILLFAC(L) = 1. / DENSCON(L)
         ENDIF
      ENDDO
      CALL READMS (3+M, RNE,ND,       'RNE     ', IERR)
      CALL READMS (3+M, T,ND,         'T       ', IERR)
      CALL READMS (3+M, VELO,ND,      'VELO    ', IERR)
      CALL READMS (3+M, GRADI,ND,     'GRADI   ', IERR)
      CALL READMS (3+M, NF,1,         'NF      ', IERR)
      IF (NF.GT.NFDIM) THEN
         CALL REMARK ('TOO MANY FREQUENCY POINTS')
         STOP 'ERROR'
      ENDIF
      CALL READMS (3+M, XLAMBDA,   NF,'XLAMBDA ', IERR)
      CALL READMS (3+M, POPNUM,ND*N,  'POPNUM  ', IERR)
      CALL READMS (3+M, RSTAR,1,      'RSTAR   ', IERR)
      CALL READMS (3+M, VDOP,1,       'VDOP    ', IERR)
      CALL READMS (3+M, JOBNUM,1,     'JOBNUM  ', IERR)
      CALL READMS (3+M, MODHEAD,13,   'MODHEAD ', IERR)
      CALL READMS (3+M, TEFF,1,       'TEFF    ', IERR)
      CALL READMS (3+M, XMDOT,1,      'XMDOT   ', IERR)

C***  Flags for the POPMIN levels
      CALL READMS (3+M,ZERO_RATES,  N*ND, 'ZERO_RAT', IERR)
C*    Default if variable does not exist yet
      IF (IERR .EQ. -10) THEN
        DO I=1, N*ND
          ZERO_RATES(I) = .FALSE.
        ENDDO
      ENDIF

      CALL READMS(3+M, RCON,1,        'RCON    ', IERR)
      
C***  'XDATA' IN FORMAL NOT NECESSARY
C***  READ 'XDATA' AND CHECK WHETHER THE RECORD EXISTS
      IERR=1
      CALL READMS (3+M, XDATA,MAXXDAT,'XDATA   ',IERR)
      IF (IERR .LT. 0) THEN
         CALL REMARK ('ERROR WHEN READING XDATA FROM MODEL FILE')
C***     XFILL EQ 0.
         XDATA(1) = 0.  
      ENDIF

C***  READ ALL CONTINUUM INTENSITIES
      ND3=3*ND
      DO 15 K=1,NF
        WRITE (NAME, '(A3, I4, A1)') 'XJC', K, ' '
        CALL READMS(3+M, XJC(1+ND*(K-1)),ND,NAME, IERR)
   15 CONTINUE
 
C***  ABXYZ only needed for calling LTEPOP in MANIPOP, and for PRI_PAR
C***  READ 'ABXYZ' (REL. ABUNDANCES OF ALL ELEMENTS) AND CHECK WHETHER
C***  THE RECORD EXISTS
      IERR=1
      CALL READMS (3+M, ABXYZ,NATOM,'ABXYZ   ',IERR)
      IF ((IERR .LT. 0) .AND. (IERR .NE. -10)) THEN
         CALL REMARK ('ERROR WHEN READING ABXYZ FROM MODEL FILE')
         STOP 'ERROR'
      ENDIF

C***  NOT EXISTING RECORD 'ABXYZ': DEFAULT IS AN ATOMIC DATA FILE "DATOM"
C***  CONTAINING "HELIUM" AS THE ONLY ELEMENT
      IF (IERR .EQ. -10) THEN
         IF (NATOM .EQ. 1) THEN
            ABXYZ(1)=1.
         ELSE
            CALL REMARK ('NOT EXISTING RECORD ABXYZ')
            STOP 'ERROR'
         ENDIF
      ENDIF
 
      CALL CLOSMS(3+M,  IERR)

      RETURN
      END
 
 
