      SUBROUTINE PRIHIST (MODHIST,LAST,MODHEAD,JOBNUM)
C***********************************************************************
C***  PRINTOUT OF THE MODEL HISTORY
C***  called by STEAL, ...?
C***********************************************************************

      IMPLICIT NONE

      !DIMENSION MODHIST(LAST)
      INTEGER, INTENT(IN) :: LAST, JOBNUM
      CHARACTER(LAST*8), INTENT(IN) :: MODHIST
      CHARACTER(100), INTENT(IN) :: MODHEAD

      INTEGER :: LASTCHAR, ICHAR1, ICHAR2,
     >           IA, IB, IFROM, ITO, INOW, J,
     >           ISTART, IEND,
     >           LINELEN, MAXLINELEN
      CHARACTER(1) :: SLASH, DOT

      LOGICAL :: NEWJOBLINE

      CHARACTER(40) :: FMTFIRSTLINE, FMTOTHERLINES

      !Output format for model history
      MAXLINELEN = 120
      FMTFIRSTLINE = '(11X,(A))'
      FMTOTHERLINES = '(27X,(A))'


      PRINT 1,MODHEAD,JOBNUM
    1 FORMAT ('1',//, 10X, A,10X, 'JOB NO.', I7, //, 11X,
     $       'M O D E L   H I S T O R Y', /, 11X, 25('='), //)

      LASTCHAR = LAST * 8
      IFROM = 9
      ITO = 8   !Just to have an intial value

      IF (LAST <= 1) THEN
        !no history saved yet => exit subroutine
        PRINT FMTFIRSTLINE, 'no history has been saved yet'
        RETURN
      ENDIF

      histmainloop: DO INOW=9, LASTCHAR, 1 !-----------------------

        !Find the beginning of the next new line (= next job)
        SLASH = MODHIST(INOW:INOW)
        NEWJOBLINE = .FALSE.
        IF (SLASH == '/') THEN
          !Slash gefunden => pruefe ob dies der Anfang eines neuen Jobeintrags ist
          !Nach dem Slash duerfen nur Leerzeichen oder 0-9 plus folgen. Die Folge
          !wird mit einem Punkt abgeschlossen
          J=INOW
          newlineloop: DO !- - - - - - - - - - - - - - - -
            J = J+1
            IF (J >= LASTCHAR) THEN
              !prevent running behind the last character
              ITO = LASTCHAR
              NEWJOBLINE = .TRUE.
              EXIT newlineloop
            ENDIF

            SELECTCASE(MODHIST(J:J))
              CASE (' ', '0', '1':'9')
                CYCLE newlineloop
              CASE ('.')
                ITO = INOW - 1
                NEWJOBLINE = .TRUE.
                EXIT newlineloop
              CASE DEFAULT
                EXIT newlineloop
            ENDSELECT
          
          ENDDO newlineloop !- - - - - - - - - - - - - - -

        ENDIF

        !Failsafe: If end of MODhist is reached, printout the rest
        IF ((.NOT. NEWJOBLINE) .AND. (INOW == LASTCHAR)) THEN
            NEWJOBLINE = .TRUE.
            ITO = LASTCHAR
        ENDIF

        IF (NEWJOBLINE .AND. (ITO >= IFROM)) THEN    
            !Ausgabe der aktuellen Zeile
            LINELEN = ITO - IFROM + 1
            IF (LINELEN > MAXLINELEN) THEN
               !Entry is too long for one line, output in several lines
               ISTART = IFROM
               IEND = ISTART + MAXLINELEN - 1
               PRINT FMTFIRSTLINE, MODHIST(ISTART:IEND)                 
               ISTART = IEND + 1
               DO WHILE(ISTART <= ITO)
                 IEND = ISTART + MAXLINELEN - 1
                 IF (IEND > ITO) THEN
                   IEND = ITO
                 ENDIF
                 PRINT FMTOTHERLINES, MODHIST(ISTART:IEND)
                 ISTART = IEND + 1                 
               ENDDO
            ELSE
               PRINT FMTFIRSTLINE, MODHIST(IFROM:ITO)
            ENDIF

            IFROM = INOW
        ENDIF

      ENDDO histmainloop !----------------------------------------
 
      RETURN
      END
