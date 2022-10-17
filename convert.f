      SUBROUTINE CONVERT

      PARAMETER (IADRL = 128)
      PARAMETER (MAXADR = 25600)
      PARAMETER (MAXXJL = 10000)
      PARAMETER (MAXARRAY = 100000)

      DIMENSION X(MAXARRAY), J(MAXARRAY)
      DIMENSION IADR(MAXADR), MODHIST(3000)

      CHARACTER*8 NAME
      CHARACTER*104 MODHEAD

      CHARACTER LINE*132, ACTPAR*40, ACTPAR2*40, ACTPAR4*40, C1*1
      CHARACTER OFILE*20, IFILE*20

      LOGICAL BSHORT, BTOAS

C***  SHORT VERSION: DO NOT GIVE OLD POPNUMS AND TEMPERATURES AND
C***                 DO NOT GIVE RADIATION FIELDS, BLOCKING FACTORS AND
C***                 EDDIS

      ICH = 0
      BTOAS = .FALSE.
      BSHORT = .FALSE.

C***  Read the Input-File convert.inp
      OPEN (UNIT=1, FILE='convert.inp')
   10 CONTINUE
        READ (1, '(A)') LINE
        CALL SARGV(LINE, 1, ACTPAR)
        IF (ACTPAR(1:7) .EQ. 'HEADEND') THEN
          GOTO 12
        ENDIF
        IF (ACTPAR(1:1) .NE. '*') THEN
          CALL SARGV(LINE, 2, ACTPAR2)
          IF (ACTPAR .EQ. 'SHORT') THEN
            IF (ACTPAR2 .EQ. 'TRUE') THEN
              BSHORT = .TRUE.
            ENDIF
          ELSE IF (ACTPAR .EQ. 'MODUS') THEN
            IF (ACTPAR2 .EQ. 'TO-ASCII') THEN
              BTOAS = .TRUE.
            ELSEIF (ACTPAR2 .EQ. 'TO-MS') THEN
              BTOAS = .FALSE.
            ELSE
              WRITE (0,*) 'Modus not known : MODUS=', ACTPAR2
              STOP 'ERROR in Subr. CONVERT'
            ENDIF
          ELSE IF (ACTPAR .EQ. 'ICH') THEN
            READ (ACTPAR2,'(I40)') ICH
          ELSE IF (ACTPAR .EQ. 'OFILE') THEN
            OFILE = ACTPAR2
          ELSE IF (ACTPAR .EQ. 'IFILE') THEN
            IFILE = ACTPAR2
          ENDIF
        ENDIF
      GOTO 10
   12 CONTINUE

      IF (BTOAS) THEN
C***  FROM MS TO ASCII
        CALL OPENMS(ICH, IADR, MAXADR, 1, IERR)
        OPEN (UNIT=4, FILE=OFILE, STATUS='NEW')

   60   CONTINUE
          READ (1, '(A)', END=62) LINE
          CALL SARGV(LINE, 1, ACTPAR)
          CALL SARGV(LINE, 2, ACTPAR2)
          ACTPAR = LINE(1:8)
          ACTPAR2 = LINE(10:14)
          READ (ACTPAR2,'(I40)') NVAR
          IF (NVAR .GT. MAXARRAY) THEN
            WRITE (0,*) 'CONVERT : ARRAY LENGTH TOO LARGE!!!'
            STOP 'ERROR IN CONVERT'
          ENDIF
          C1 = ACTPAR(1:1)
          IF (ACTPAR .EQ. 'MODHEAD' .OR. ACTPAR .EQ. 'MODEL') THEN
            CALL READMS(ICH, J, NVAR, ACTPAR(1:8), IERR)
              IF (IERR .NE. 0) GOTO 100
            WRITE (4,'(A)') '*'
            WRITE (4,'(A,A8,A6,I5)') 
     >        'NAME=', ACTPAR(1:8), ' NVAR=', NVAR
            WRITE (4,'(15(A8))') (J(I),I=1,NVAR)
          ELSEIF (ACTPAR .EQ. 'REDISMO' .OR. 
     >            ACTPAR .EQ. 'IONNAME' .OR. 
     >            ACTPAR .EQ. 'FGRIDS' .OR. 
     >            ACTPAR(5:7) .EQ. 'NAM' .OR. 
     >            ACTPAR .EQ. 'GENERIC') THEN
            CALL READMS(ICH, J, NVAR, ACTPAR(1:8), IERR)
            IF (IERR .NE. 0) GOTO 100
            WRITE (4,'(A)') '*'
            WRITE (4,'(A,A8,A6,I5)') 
     >        'NAME=', ACTPAR(1:8), ' NVAR=', NVAR
            WRITE (4,'(3(A8))') (J(I),I=1,NVAR)
          ELSEIF (ACTPAR .EQ. 'MODHIST') THEN
            CALL READMS(ICH, J, NVAR, ACTPAR(1:8), IERR)
            IF (IERR .NE. 0) GOTO 100
            WRITE (4,'(A)') '*'
            WRITE (4,'(A,A8,A6,I5)') 
     >        'NAME=', ACTPAR(1:8), ' NVAR=', NVAR
            WRITE (4,'(I8,15(A8))') J(1), (J(I),I=2,15)
            WRITE (4,'(15(A8))') (J(I),I=16,J(1))
          ELSEIF (C1 .EQ. 'I' .OR.
     >            C1 .EQ. 'J' .OR.
     >            C1 .EQ. 'K' .OR.
     >            C1 .EQ. 'L' .OR.
     >            C1 .EQ. 'M' .OR.
     >            C1 .EQ. 'N') THEN
            CALL READMS(ICH, J, NVAR, ACTPAR(1:8), IERR)
            IF (IERR .NE. 0) GOTO 100
            WRITE (4,'(A)') '*'
            WRITE (4,'(A,A8,A6,I5)') 
     >        'NAME=', ACTPAR(1:8), ' NVAR=', NVAR
            WRITE (4,'(3(I20,1X))') (J(I),I=1,NVAR)
          ELSE
            CALL READMS(ICH, X, NVAR, ACTPAR(1:8), IERR)
            IF (IERR .NE. 0) GOTO 100
            WRITE (4,'(A)') '*'
            WRITE (4,'(A,A8,A6,I5)') 
     >        'NAME=', ACTPAR(1:8), ' NVAR=', NVAR
            WRITE (4,'(3(G25.19,1X))') (X(I),I=1,NVAR)
          ENDIF
        GOTO 60
   62   CONTINUE
      ELSE
C***  FROM ASCII TO MS
        WRITE (0,*) 'FROM ASCII TO MS: File=', IFILE
        CALL OPENMS(ICH, IADR, MAXADR, 1, IERR)
        OPEN (UNIT=4, FILE=IFILE, STATUS='OLD')
  120   CONTINUE
          READ (4, '(A)', END=140) LINE
          ACTPAR=' '
          CALL SARGV(LINE, 1, ACTPAR)
          IF (ACTPAR .EQ. 'NAME') THEN
            ACTPAR2 = LINE(6:13)
            C1 = ACTPAR2(1:1)
            ACTPAR4 = LINE(20:24)
            READ (ACTPAR4,'(I40)') NVAR
            IF (NVAR .EQ. 0) THEN
              CALL WRITMS(ICH, X, NVAR, ACTPAR2(1:8), IERR)
            ELSEIF (ACTPAR2 .EQ. 'MODHEAD' .OR. 
     >              ACTPAR2 .EQ. 'MODEL') THEN
              READ (4, '(A)') LINE
              READ (LINE,'(15(A8))') (J(I),I=1,NVAR)
              CALL WRITMS(ICH, J, NVAR, ACTPAR2(1:8), IERR)
            ELSEIF (ACTPAR2 .EQ. 'REDISMO' .OR. 
     >              ACTPAR2 .EQ. 'IONNAME' .OR. 
     >              ACTPAR2 .EQ. 'FGRIDS' .OR. 
     >              ACTPAR2(5:7) .EQ. 'NAM' .OR. 
     >              ACTPAR2 .EQ. 'GENERIC') THEN
              DO I=1, NVAR, 3
                K = MIN (I+2, NVAR)
                READ (4, '(A)') LINE
                READ (LINE,'(3(A8))') (J(L),L=I,K)
              ENDDO
C              READ (4, '(A)') LINE
C              READ (LINE,'(3(A8))') J(1), (J(L),L=1,NVAR)
              CALL WRITMS(ICH, J, NVAR, ACTPAR2(1:8), IERR)
            ELSEIF (ACTPAR2 .EQ. 'MODHIST') THEN
              READ (4, '(A)') LINE
              READ (LINE,'(I8,15(A8))') J(1), (J(L),L=2,15)
              DO I=16, J(1), 15
                K = MIN (I+14, J(1))
                READ (4, '(A)') LINE
                READ (LINE,'(15(A8))') (J(L),L=I,K)
              ENDDO
              CALL WRITMS(ICH, J, NVAR, ACTPAR2(1:8), IERR)
            ELSEIF (C1 .EQ. 'I' .OR.
     >              C1 .EQ. 'J' .OR.
     >              C1 .EQ. 'K' .OR.
     >              C1 .EQ. 'L' .OR.
     >              C1 .EQ. 'M' .OR.
     >              C1 .EQ. 'N') THEN
              DO I=1, NVAR, 3
                READ (4, '(A)') LINE
                K = MIN (I+2, NVAR)
                READ (LINE,'(3(I20,1X))') (J(L),L=I,K)
              ENDDO
              CALL WRITMS(ICH, J, NVAR, ACTPAR2(1:8), IERR)
            ELSE
              DO I=1, NVAR, 3
                READ (4, '(A)') LINE
                K = MIN (I+2, NVAR)
                READ (LINE,'(3(F25.0,1X))') (X(L),L=I,K)
              ENDDO
              CALL WRITMS(ICH, X, NVAR, ACTPAR2(1:8), IERR)
            ENDIF
          ENDIF
        GOTO 120
  140   CONTINUE
      ENDIF

      CALL CLOSMS(ICH, IERR)
      WRITE (0,*) 'CONVERT has worked correctly'
      STOP 

  100 write (0,*) 'ERROR IN READMS'
      STOP 'ERROR in Subr. CONVERT'

      END
