      SUBROUTINE ST2AS

      PARAMETER (IADRL = 128)
      PARAMETER (MAXADR = 6000)
      PARAMETER (MAXXJL = 10000)

      DIMENSION X(100000), J(100000)
      DIMENSION IADR(MAXADR), MODHIST(3000)

      CHARACTER*8 NAME
      CHARACTER*104 MODHEAD

      LOGICAL BSHORT

C***  SHORT VERSION: DO NOT GIVE OLD POPNUMS AND TEMPERATURES AND
C***                 DO NOT GIVE RADIATION FIELDS, BLOCKING FACTORS AND
C***                 EDDIS
      BSHORT = .FALSE.

C***  NUMBER OF LEVELS = N : MUST BE INSERTED BY HAND
C***  MOG1 : 132
C***  MOG2 : 150
C***  MOG3 : 175
C***  MOG4 : 209
C***  MOG5 : 223
C***  MOG6 : 241
      N = 175

C***  NUMBER OF ATOMS = NATOM : MUST BE INSERTED BY HAND
      NATOM = 3

      OPEN (UNIT=4, FILE='ASCII_MODEL', STATUS='NEW')
      CALL OPENMS(3, IADR, MAXADR, 1, IERR)

C***  CONVERT MODHEAD
      IERR = -1
      CALL READMS(3, MODHEAD, 13, 'MODHEAD ', IERR)
      IF (IERR .NE. 0) GOTO 100
      WRITE (4,50) 'MODHEAD'
   50 FORMAT ('"', A7, '"')
      WRITE (4,'(A100)') MODHEAD
      WRITE (4,'(A1)') '*'

C***  WRITE N
      WRITE (4,50) 'NLEVEL'
      WRITE (4,'(I3)') N
      WRITE (4,'(A1)') '*'

C***  WRITE NATOM
      WRITE (4,50) 'NATOM'
      WRITE (4,'(I3)') NATOM
      WRITE (4,'(A1)') '*'

C***  CONVERT ND
      IERR = -1
      CALL READMS(3, ND, 1, 'ND      ', IERR)
      IF (IERR .NE. 0) GOTO 100
      WRITE (4,50) 'ND'
      WRITE (4,'(I3)') ND
      WRITE (4,'(A1)') '*'

C***  CONVERT R
      IERR = -1
      CALL READMS(3, X, ND, 'R       ', IERR)
      IF (IERR .NE. 0) GOTO 100
      WRITE (4,50) 'R'
      WRITE (4,'(3(G25.19,1X))') (X(I),I=1,ND)
      WRITE (4,'(A1)') '*'

C***  CONVERT NP
      IERR = -1
      CALL READMS(3, NP, 1, 'NP      ', IERR)
      IF (IERR .NE. 0) GOTO 100
      WRITE (4,50) 'NP'
      WRITE (4,'(I3)') NP
      WRITE (4,'(A1)') '*'

C***  CONVERT P
      IERR = -1
      CALL READMS(3, X, NP, 'P       ', IERR)
      IF (IERR .NE. 0) GOTO 100
      WRITE (4,50) 'P'
      WRITE (4,'(3(G25.19,1X))') (X(I),I=1,NP)
      WRITE (4,'(A1)') '*'

C***  CONVERT Z
      IERR = -1
      CALL READMS(3, X, ND*NP, 'Z       ', IERR)
      IF (IERR .NE. 0) GOTO 100
      WRITE (4,50) 'Z'
      WRITE (4,'(3(G25.19,1X))') (X(I),I=1,ND*NP)
      WRITE (4,'(A1)') '*'

C***  CONVERT ENTOT
      IERR = -1
      CALL READMS(3, X, ND, 'ENTOT   ', IERR)
      IF (IERR .NE. 0) GOTO 100
      WRITE (4,50) 'ENTOT'
      WRITE (4,'(3(G25.19,1X))') (X(I),I=1,ND)
      WRITE (4,'(A1)') '*'

C***  CONVERT TEFF
      IERR = -1
      CALL READMS(3, X, 1, 'TEFF    ', IERR)
      IF (IERR .NE. 0) GOTO 100
      WRITE (4,50) 'TEFF'
      WRITE (4,'(3(G25.19,1X))') (X(I),I=1,1)
      WRITE (4,'(A1)') '*'

C***  CONVERT T
      IERR = -1
      CALL READMS(3, X, ND, 'T       ', IERR)
      IF (IERR .NE. 0) GOTO 100
      WRITE (4,50) 'T'
      WRITE (4,'(3(G25.19,1X))') (X(I),I=1,ND)
      WRITE (4,'(A1)') '*'

C***  CONVERT TOLD1
      IF (.NOT. BSHORT) THEN
        IERR = -1
        CALL READMS(3, X, ND, 'TOLD1   ', IERR)
        IF (IERR .EQ. 0) THEN
          WRITE (4,50) 'TOLD1'
          WRITE (4,'(3(G25.19,1X))') (X(I),I=1,ND)
          WRITE (4,'(A1)') '*'
        ENDIF
      ENDIF

C***  CONVERT TOLD2
      IF (.NOT. BSHORT) THEN
        IERR = -1
        CALL READMS(3, X, ND, 'TOLD2   ', IERR)
        IF (IERR .EQ. 0) THEN
          WRITE (4,50) 'TOLD2'
          WRITE (4,'(3(G25.19,1X))') (X(I),I=1,ND)
          WRITE (4,'(A1)') '*'
        ENDIF
      ENDIF

C***  CONVERT VELO
      IERR = -1
      CALL READMS(3, X, ND, 'VELO    ', IERR)
      IF (IERR .NE. 0) GOTO 100
      WRITE (4,50) 'VELO'
      WRITE (4,'(3(G25.19,1X))') (X(I),I=1,ND)
      WRITE (4,'(A1)') '*'

C***  CONVERT GRADI
      IERR = -1
      CALL READMS(3, X, ND, 'GRADI   ', IERR)
      IF (IERR .NE. 0) GOTO 100
      WRITE (4,50) 'GRADI'
      WRITE (4,'(3(G25.19,1X))') (X(I),I=1,ND)
      WRITE (4,'(A1)') '*'

C***  CONVERT RSTAR
      IERR = -1
      CALL READMS(3, X, 1, 'RSTAR   ', IERR)
      IF (IERR .NE. 0) GOTO 100
      WRITE (4,50) 'RSTAR'
      WRITE (4,'(3(G25.19,1X))') (X(I),I=1,1)
      WRITE (4,'(A1)') '*'

C***  CONVERT VDOP
      IERR = -1
      CALL READMS(3, X, 1, 'VDOP    ', IERR)
      IF (IERR .NE. 0) GOTO 100
      WRITE (4,50) 'VDOP'
      WRITE (4,'(3(G25.19,1X))') (X(I),I=1,1)
      WRITE (4,'(A1)') '*'

C***  CONVERT NF
      IERR = -1
      CALL READMS(3, NF, 1, 'NF      ', IERR)
      IF (IERR .NE. 0) GOTO 100
      WRITE (4,50) 'NF'
      WRITE (4,'(I5)') NF
      WRITE (4,'(A1)') '*'

C***  CONVERT XLAMBDA
      IERR = -1
      CALL READMS(3, X, NF, 'XLAMBDA ', IERR)
      IF (IERR .NE. 0) GOTO 100
      WRITE (4,50) 'XLAMBDA'
      WRITE (4,'(3(G25.19,1X))') (X(I),I=1,NF)
      WRITE (4,'(A1)') '*'

C***  CONVERT FWEIGHT
      IERR = -1
      CALL READMS(3, X, NF, 'FWEIGHT ', IERR)
      IF (IERR .NE. 0) GOTO 100
      WRITE (4,50) 'FWEIGHT'
      WRITE (4,'(3(G25.19,1X))') (X(I),I=1,NF)
      WRITE (4,'(A1)') '*'

C***  CONVERT KEY
      IERR = -1
      CALL READMS(3, J, NF, 'KEY     ', IERR)
      IF (IERR .NE. 0) GOTO 100
      WRITE (4,50) 'KEY'
      WRITE (4,'(3(I25,1X))') (J(I),I=1,NF)
      WRITE (4,'(A1)') '*'

C***  CONVERT POPNUM
      IERR = -1
      CALL READMS(3, X, N*ND, 'POPNUM  ', IERR)
      IF (IERR .NE. 0) GOTO 100
      WRITE (4,50) 'POPNUM'
      WRITE (4,'(3(G25.19,1X))') (X(I),I=1,N*ND)
      WRITE (4,'(A1)') '*'

C***  CONVERT POP1
      IF (.NOT. BSHORT) THEN
        IERR = -1
        CALL READMS(3, X, N*ND, 'POP1    ', IERR)
        IF (IERR .EQ. 0) THEN
          WRITE (4,50) 'POP1'
          WRITE (4,'(3(G25.19,1X))') (X(I),I=1,N*ND)
          WRITE (4,'(A1)') '*'
        ENDIF
      ENDIF

C***  CONVERT POP2
      IF (.NOT. BSHORT) THEN
        IERR = -1
        CALL READMS(3, X, N*ND, 'POP2    ', IERR)
        IF (IERR .EQ. 0) THEN
          WRITE (4,50) 'POP2'
          WRITE (4,'(3(G25.19,1X))') (X(I),I=1,N*ND)
          WRITE (4,'(A1)') '*'
        ENDIF
      ENDIF

C***  CONVERT RNE
      IERR = -1
      CALL READMS(3, X, ND, 'RNE     ', IERR)
      IF (IERR .NE. 0) GOTO 100
      WRITE (4,50) 'RNE'
      WRITE (4,'(3(G25.19,1X))') (X(I),I=1,ND)
      WRITE (4,'(A1)') '*'

C***  CONVERT ABXYZ
      IERR = -1
      CALL READMS(3, X, NATOM, 'ABXYZ   ', IERR)
      IF (IERR .NE. 0) GOTO 100
      WRITE (4,50) 'ABXYZ'
      WRITE (4,'(3(G25.19,1X))') (X(I),I=1,NATOM)
      WRITE (4,'(A1)') '*'

C***  CONVERT NEXTK
      IERR = -1
      CALL READMS(3, J, 1, 'NEXTK   ', IERR)
      IF (IERR .NE. 0) GOTO 100
      WRITE (4,50) 'NEXTK'
      WRITE (4,'(3(I5,1X))') (J(I),I=1,1)
      WRITE (4,'(A1)') '*'

C***  CONVERT JOBNUM
      IERR = -1
      CALL READMS(3, J, 1, 'JOBNUM  ', IERR)
      IF (IERR .NE. 0) GOTO 100
      WRITE (4,50) 'JOBNUM'
      WRITE (4,'(3(I5,1X))') (J(I),I=1,1)
      WRITE (4,'(A1)') '*'

C***  CONVERT TOTOUT
      IERR = -1
      CALL READMS(3, X, 1, 'TOTOUT  ', IERR)
      IF (IERR .NE. 0) GOTO 100
      WRITE (4,50) 'TOTOUT'
      WRITE (4,'(3(G25.19,1X))') (X(I),I=1,1)
      WRITE (4,'(A1)') '*'

C***  CONVERT XDATA
      IERR = -1
      CALL READMS(3, X, 10, 'XDATA   ', IERR)
      IF (IERR .NE. 0) GOTO 100
      WRITE (4,50) 'XDATA'
      WRITE (4,'(3(G25.19,1X))') (X(I),I=1,10)
      WRITE (4,'(A1)') '*'

C***  CONVERT MODHIST
      IERR = -1
      CALL READMS(3, LAST, 1, 'MODHIST ', IERR)
      IF (IERR .NE. 0) GOTO 100
      IERR = -1
      CALL READMS(3, J, LAST, 'MODHIST ', IERR)
      IF (IERR .NE. 0) GOTO 100
      WRITE (4,50) 'MODHIST'
      WRITE (4,'(I5)') LAST
      WRITE (4,55) (J(I),I=2,LAST)
   55 FORMAT (3('"',A8,'"',10X))
      WRITE (4,'(A1)') '*'

C***  CONVERT EMFLUX
      IERR = -1
      CALL READMS(3, X, NF, 'EMFLUX  ', IERR)
      IF (IERR .NE. 0) GOTO 100
      WRITE (4,50) 'EMFLUX'
      WRITE (4,'(3(G25.19,1X))') (X(I),I=1,NF)
      WRITE (4,'(A1)') '*'

C***  CONVERT TOTIN
      IERR = -1
      CALL READMS(3, X, 1, 'TOTIN   ', IERR)
      IF (IERR .NE. 0) GOTO 100
      WRITE (4,50) 'TOTIN'
      WRITE (4,'(3(G25.19,1X))') (X(I),I=1,1)
      WRITE (4,'(A1)') '*'

C***  CONVERT XJC, EDDI, BLF
      IF (.NOT. BSHORT) THEN
        DO K=1, NF
          WRITE (NAME,'(A3,I4,A1)') 'XJC',K, ' '
C!!!          ENCODE (7,10,NAME) 'XJC ', K
C!!!   10     FORMAT(A4,I3)
          IERR = -1
          CALL READMS(3, X, ND, NAME, IERR)
          IF (IERR .NE. 0) GOTO 100
          WRITE (4,50) NAME
          WRITE (4,'(3(G25.19,1X))') (X(I),I=1,ND)
          WRITE (4,'(A1)') '*'
C***
          WRITE (NAME,'(A4,I3,A1)') 'EDDI',K,' '
C!!!          ENCODE (7,10,NAME) 4HEDDI,K
          IERR = -1
          CALL READMS(3, X, 3*ND, NAME, IERR)
          IF (IERR .NE. 0) GOTO 100
          WRITE (4,50) NAME
          WRITE (4,'(3(G25.19,1X))') (X(I),I=1,3*ND)
          WRITE (4,'(A1)') '*'
C***
          WRITE (NAME,'(A4,I3,A1)') 'BLF ',K, ' '
C!!!          ENCODE (7,10,NAME) 3HBLF,K
          IERR = -1
          CALL READMS(3, X, ND, NAME, IERR)
          IF (IERR .NE. 0) GOTO 100
          WRITE (4,50) NAME
          WRITE (4,'(3(G25.19,1X))') (X(I),I=1,ND)
          WRITE (4,'(A1)') '*'
        ENDDO
      ENDIF

C***  CONVERT XJL
      IF (.NOT. BSHORT) THEN
        DO K=1, MAXXJL
          WRITE (NAME, '(A3,I4,A1)') 'XJL',K,' '
C!!!          ENCODE (7,11,NAME) 3HXJL,K
C!!!   11     FORMAT(A3,I4)
          IERR = -1
          CALL READMS(3, X, ND, NAME, IERR)
          IF (IERR .EQ. -10) GOTO 20
            IF (IERR .NE. 0) GOTO 100
            WRITE (4,50) NAME
            WRITE (4,'(3(G25.19,1X))') (X(I),I=1,ND)
            WRITE (4,'(A1)') '*'
   20     CONTINUE
        ENDDO
      ENDIF

      WRITE (4,'(2A,L1)') 'CR2AS HAS WORKED CORRECTLY ',
     >                   'SHORT VERSION WAS ',BSHORT

      CALL CLOSMS(3, IERR)
      CLOSE (4)

      RETURN

  100 WRITE (4,*) 'ERROR WHILE READING FROM MASS-STORAGE'
      WRITE (4,'(A5,I4)') 'IERR=',IERR

      END
