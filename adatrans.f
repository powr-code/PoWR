      SUBROUTINE ADATRANS (NTRANS, ADPWEIGHT, N, NOLD, NEWATOM, POPMIN,
     >                LEVELCARD, NLEVELCARD, MAXLEVELCARD)

C*****************************************************************
C***  Interprets the LEVEL-Cards for the ADAPTER and writes 
C***  the vector NTRANS 
C***     index: level index in new model atom
C***     value: -1 --> level is not overwritten 
C***             0 --> level is set to POPMIN
C***            >0 --> corresponding index in old model atom 
C***
C***  LEVEL options modernized by wrh 29-Jul-2007 
C***  This subroutine was separated from DECADP, because only 
C***  after reading both MODEL files, NOLD becomes available
C***      wrh, 27-Nov-2010 
C*****************************************************************

      DIMENSION NTRANS(N), ADPWEIGHT(N)
      LOGICAL NEWATOM, CHRINSTR
      CHARACTER KARTE*80, REST*80, ACTPAR1*20, ACTPAR2*20
      CHARACTER LEVELCARD*80(MAXLEVELCARD)

C***  DEFAULT VALUES: 
      NEWATOM=.FALSE.
      DO I=1,N
        ADPWEIGHT(I) = 1.
        NTRANS(I)    = -1
      ENDDO

      NEWATOM = NLEVELCARD .GT. 0

C***  Simple OLDSTART  (NO "LEVEL" CARDS GIVEN)
      IF (.NOT. NEWATOM) THEN
         DO J=1, N
           NTRANS(J) = J
         ENDDO
C***    ...  but model atom became bigger? --> Set new levels ZERO
         IF (N .GT. NOLD) THEN
            NEWATOM = .TRUE.
            DO J=NOLD+1, N
               NTRANS(J) = 0
            ENDDO
            WRITE (0,'(A,1PG8.1)') 
     >            'NOTE (ADAPTER): NEW MODEL HAS MORE ' //
     >            'LEVELS THAN OLD --> SET to POPMIN:',POPMIN  
         ENDIF
      ENDIF


C***  Oldstart with LEVEL cards

C***  Loop over all LEVEL cards
      DO ILC =1, NLEVELCARD
         KARTE = LEVELCARD(ILC)
         REST = KARTE(6:)       
         CALL SARGC (REST, NPAR)
         IF (NPAR .LE. 0) GOTO 90        
         CALL SARGV (REST,1,ACTPAR1)
         LAST1 = IDX (ACTPAR1)

C***     If there are no more parameters, skip the next part!
         IF (NPAR .GE.2) THEN       
            CALL SARGREST (REST, NPAR, 2, IA, IE)
            REST = REST(IA:)
            CALL SARGV (REST,1,ACTPAR2)

C***        If last character of ACTPAR1 is "-" 
C***         -->  directly append next parameter 
            IF (ACTPAR1(LAST1:LAST1) .EQ. '-') THEN
               ACTPAR1(LAST1+1:) = ACTPAR2
               CALL SARGREST (REST, NPAR, 2, IA, IE)
               REST = REST(IA:)
            ENDIF

C***        If first character of next parameter is "-"
            IF (ACTPAR2(1:1) .EQ. '-') THEN
C***        ... and this is the whole parameter --> glue next parameter
               IF (IDX(ACTPAR2) .EQ. 1) THEN
                  CALL SARGV (REST,2,ACTPAR2)
                  ACTPAR1 = ACTPAR1(:LAST1) // '-' // ACTPAR2
                  CALL SARGC (REST, NPAR)
                  IF (NPAR .GE. 3) THEN 
                     CALL SARGREST (REST, NPAR, 3, IA, IE)
                     REST = REST(IA:)
                  ELSE
                     REST = ' '
                  ENDIF
               ELSE
                  ACTPAR1 = ACTPAR1(:LAST1) // ACTPAR2
                  IF (NPAR .GE. 2) THEN 
                     CALL SARGREST (REST, NPAR, 2, IA, IE)
                     REST = REST(IA:)
                  ELSE
                     REST = ' '
                  ENDIF
               ENDIF
            ENDIF
         ELSE
            REST = ' '
         ENDIF

C***     Now ACTPAR1 should contain the compact string nnn-mmm
C***        --> decode nnn and mmm

         IF (CHRINSTR ('-', ACTPAR1)) THEN
            LAST1 = IDX(ACTPAR1)
            DO I=1, LAST1
               IF (ACTPAR1(I:I) .EQ. '-') THEN
                  MPOS = I
                  EXIT
               ENDIF
            ENDDO
            READ (ACTPAR1(:MPOS-1), '(I10)', ERR=91) NSTART
            READ (ACTPAR1(MPOS+1:), '(I10)', ERR=91) NSTOP
         ELSE
            READ (ACTPAR1, '(I10)', ERR=91) NSTART
            NSTOP = NSTART
         ENDIF

         NSTOPOLD = MIN(NSTOP, NOLD)  !this value is needed for SHIFT cards
         NSTOP = MIN0 (NSTOP, N)

         IF (NSTART .LE. 0) GOTO 96   ! Error exit

C***     Decode the further parameters

C***     No further parameters -- identify the levels
         IF (REST .EQ. ' ') THEN
            DO I = NSTART, NSTOP 
               IF (I .GT. NOLD) THEN
                  WRITE (0,'(A)') 'NOTE FROM ADAPTER:ADATRANS:'
                  WRITE (0,*) KARTE(:IDX(KARTE))
                  WRITE (0,'(A,I5)') ' LOOP TRUNCATED AT', I-1 
                  EXIT
               ENDIF
               NTRANS(I)=I
            ENDDO

         ELSE IF (REST(1:5) .EQ. 'SHIFT') THEN
C                                 =====
            REST = REST(6:)
            IF (REST .EQ. ' ') GOTO 92
            CALL SARGV (REST, 1, ACTPAR1)
            READ (ACTPAR1, '(I10)', ERR=92) ISHIFT
            NSTOPMAX = MIN0(NOLD, N-ISHIFT)
            NSTOP = MIN0(NSTOPOLD,NSTOPMAX)  ! corrected 1-Mar-2017 ansander
            DO I = NSTART, NSTOP 
               IF (I+ISHIFT .GT. N .OR. I .GT. NOLD) THEN
                  WRITE (0,'(A)') 'NOTE FROM ADAPTER:ADATRANS:'
                  WRITE (0,*) KARTE(:IDX(KARTE))
                  WRITE (0,'(A,I5)') ' LOOP TRUNCATED AT', I-1 
                  EXIT
               ENDIF
               NTRANS(I+ISHIFT)=I
            ENDDO


         ELSE IF (REST(1:4) .EQ. 'NULL') THEN
C                                 =====
C***       SET NEW LEVEL TO ZERO
            DO I = NSTART, NSTOP 
               NTRANS(I) = 0
            ENDDO

         ELSE IF (REST(1:6) .EQ. 'WEIGHT') THEN
            CALL SARGC (REST, NPAR)
            IF (NPAR .LT. 2) GOTO 94
            CALL SARGV (REST, 2, ACTPAR1)
            READ (ACTPAR1, '(F10.0)', ERR=94) ADPW
            DO I = NSTART, NSTOP 
               ADPWEIGHT(I) = ADPW
            ENDDO

         ELSE
C***        the only remaining possibility: next parameter is old level
            IF (NSTART .EQ. NSTOP) THEN
               CALL SARGV (REST, 1, ACTPAR1)
               READ (ACTPAR1, '(I10)', ERR=91) NOLDLEVEL
               NTRANS(NSTART) = NOLDLEVEL
            ELSE
C***           Unidentified further parameter
               GOTO 93
            ENDIF


         ENDIF
      ENDDO !--------- end of loop over all LEVEL cards


C***  Check assigned oldlevel for existence
      DO J=1, N
         IF (NTRANS(J) .GT. NOLD) THEN
            WRITE (0, '(A,/,A, I4, A, I4, A)')
     >          'ADAPTER: WARNING OF INTERNAL INCONSISTENCY',
     >          'non-existing old level:', NTRANS(J),
     >          ' asigned to new level:', J, ' - replaced by NULL'
            NTRANS(J) = 0
         ENDIF
      ENDDO

      RETURN

C***  ERROR EXITS *******************************

   90 WRITE (0,*) 'LEVEL option needs parameters!'
      GOTO 95

   91 WRITE (0,*) 'Error when decoding LEVEL index'
      GOTO 95 

   92 WRITE (0,*) 'Error when decoding SHIFT index'
      GOTO 95

   93 WRITE (0,*) 'Unidentified parameter:', REST(:IDX(REST))
      GOTO 95

   94 WRITE (0,*) 'A number must be given for the WEIGHT'
      GOTO 95

   96 WRITE (0,*) 'Invalid level index (.LE. 0) encountered'
      GOTO 95

   95 WRITE (0,*) 'THE ERROR OCCURED IN THE FOLLOWING LINE:'
      WRITE (0,*) KARTE 
      STOP 'FATAL ERROR IN ADAPTER:ADPTRANS'

      END 
