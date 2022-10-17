      SUBROUTINE DECMULTI(KARTE, NMOD, VERSION, PMERGE)
C***********************************************************************
C***  DECODING INPUT CARD "MULTI" FOR PROGRAM "FORMAL" 
C***********************************************************************
 
      CHARACTER KARTE*80, ACTPAR*10
      CHARACTER*(*) VERSION
 
      NMOD = 1

    1 READ (2, 2, END=94) KARTE
    2 FORMAT (A)
      CALL SARGV (KARTE, 1, ACTPAR)

C*** END OF INPUT DATA REACHED
      IF (IOSTAT .LT. 0 .OR. KARTE(:3) .EQ. 'END') THEN
        GOTO 94

      ELSE IF ( ACTPAR .EQ. 'MULTI' ) THEN
C                            =====
        CALL SARGV (KARTE, 2, ACTPAR)
        READ (ACTPAR,'(I10)', ERR=96) NMOD
        CALL SARGC(KARTE, NPAR)
        IF (NPAR .GT. 2) THEN
          CALL SARGV (KARTE, 3, ACTPAR)
          VERSION = ACTPAR
          IF (VERSION .EQ. 'MERGE') THEN
            CALL SARGV (KARTE, 4, ACTPAR)
            PMERGE=0.5
            READ (ACTPAR,'(F10.0)', ERR=80) PMERGE
   80       CONTINUE
          ENDIF
          IF (VERSION .NE. 'CONE' .AND. VERSION .NE. 'MERGE' .AND.
     >        VERSION .NE. 'SPHERE' .AND.
     >        VERSION .NE. 'TEST') THEN
            WRITE (0,*) 'VERSION NOT KNOWN, VERSION=', VERSION
            STOP 'ERROR IN SUBR. DECMULTI'
          ENDIF
        ELSE
          VERSION = 'CONE'
        ENDIF
        IF (VERSION .EQ. 'MERGE') THEN
          WRITE (0,'(A,A,A,F8.4)') 'MULTI-CARD DETECTED: VERSION=', 
     >                 VERSION( :IDX(VERSION)), '  PMERGE=', PMERGE
        ELSE
          WRITE (0,*) 'MULTI-CARD DETECTED: VERSION=', 
     >                 VERSION( :IDX(VERSION))
        ENDIF

      ENDIF
 
      GOTO 1
 
   94 RETURN

   96 WRITE (0,*) 'ERROR WHEN READING MULTI'
      STOP 'ERROR IN SUBR. DECMULTI'

      END
