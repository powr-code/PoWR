      SUBROUTINE LIMBDARK_OUTPUT (KANAL, LIMB_LINE, NFOBS, EMINT_P, 
     >                            NP_MERGED, PGRID_MERGED) 
C******************************************************************************
C***  PLOT OF THE LIMB DARKENING 
C******************************************************************************
      IMPLICIT NONE
      CHARACTER :: LIMB_LINE*(*), LIMBFILENAME*100
      INTEGER, INTENT(IN) :: NP_MERGED, NFOBS
      REAL, DIMENSION(NP_MERGED) :: EMINT_P, PGRID_MERGED, TAURCONT

C***  Local variables
      REAL, DIMENSION(NP_MERGED) :: XPLOT, YPLOT
      CHARACTER(100) :: HEAD1, HEAD2, STYLE, ACTPAR,  
     >                 XAXISSTR, YAXISSTR
      INTEGER :: KANAL, JP, NPAR, IPAR, IDX, NP_CUT, ISRCHFGT, ISRCHFLT
      INTEGER :: LIMB_UNIT
      REAL PCUT, EMINT_MIN

C***  Defaults: 
      LIMBFILENAME = ''
      XAXISSTR = '\CENTER\p [R&T*&M]'
      YAXISSTR = 
     >   '\CENTER\I&T#n#&M [erg cm&H-2&M s&H-1&M Hz&H-1&M sterad&H-1&M]'

      DO JP = 1, NP_MERGED
         XPLOT(JP) = PGRID_MERGED(JP)
         YPLOT(JP) = EMINT_P(JP)
      ENDDO

      NP_CUT = NP_MERGED


C***  Search for plot options 
      CALL SARGC(LIMB_LINE, NPAR)

      DO IPAR=2, NPAR
         CALL SARGV(LIMB_LINE, IPAR, ACTPAR)

         IF (ACTPAR(:4) .EQ. 'NORM') THEN
C                             ====
             DO JP=1, NP_MERGED
                YPLOT(JP) = YPLOT(JP) / EMINT_P(1)
             ENDDO
                YAXISSTR = '\CENTER\I&T#n#&M / I&T#n#&M(p=0)'

         ELSEIF (ACTPAR .EQ. 'MU') THEN
C                             ==
            NP_CUT = MIN (NP_CUT, 
     >               ISRCHFGT(NP_MERGED,PGRID_MERGED,1, 1.) - 1 )
            DO JP=1, NP_CUT
               XPLOT(JP) = SQRT(1. - PGRID_MERGED(JP)**2)
            ENDDO
            XAXISSTR = '\CENTER\#m#'

C***     Plot ends before impact parameter PCUT
         ELSEIF (ACTPAR .EQ. 'PCUT') THEN
C                             ====
            IF (IPAR+1 .GT. NPAR) GOTO 100
            CALL SARGV(LIMB_LINE, IPAR+1, ACTPAR)
            READ (ACTPAR, '(F12.0)',ERR=100) PCUT
            NP_CUT = MIN (NP_CUT, 
     >               ISRCHFGT(NP_MERGED,PGRID_MERGED,1,PCUT) - 1)

         ELSEIF (ACTPAR .EQ. 'MIN') THEN
C                             ====
            IF (IPAR+1 .GT. NPAR) GOTO 100
            CALL SARGV(LIMB_LINE, IPAR+1, ACTPAR)
            READ (ACTPAR, '(F12.0)',ERR=100) EMINT_MIN
            EMINT_MIN = EMINT_MIN * EMINT_P(1)
            NP_CUT = MIN (NP_CUT, 
     >               ISRCHFLT(NP_MERGED, EMINT_P, 1, EMINT_MIN) -1)

         ELSEIF (ACTPAR .EQ. 'FILE') THEN
C                             ====
            IF (IPAR+1 .GT. NPAR) GOTO 100
            CALL SARGV(LIMB_LINE, IPAR+1, LIMBFILENAME)

         ENDIF
      ENDDO
     

  110 CONTINUE
C***  HEADER  ------------------------------------------------------
      HEAD1 = 'LIMB DARKENING'
      HEAD2 = LIMB_LINE(:IDX(LIMB_LINE))

      WRITE (KANAL, '(A,A)') 'PLOT: ', LIMB_LINE
      WRITE (KANAL, '(A)') '\INBOX'

      CALL PLOTANFS (KANAL,HEAD1,HEAD2
     >        ,XAXISSTR
     >        ,YAXISSTR
     >        , 0., 0., 0., 0., 0.,.0
     >        , 0., 0., 0., 0., 0.,.0
     >        , XPLOT ,YPLOT, NP_CUT, 'COLOR=2')
      
      WRITE (0,'(A)') 'LIMBDARKENING plot written'

C***  ASCII file output (optional)
      IF (LIMBFILENAME .NE. '') THEN
         LIMB_UNIT = 52
         OPEN (UNIT=LIMB_UNIT, FILE=LIMBFILENAME, ACTION="write", 
     >            status="unknown")
         WRITE (LIMB_UNIT, '(A)') '# (1) ' // XAXISSTR(9:IDX(XAXISSTR)) 
         WRITE (LIMB_UNIT, '(A)') '# (2) ' // YAXISSTR(9:IDX(YAXISSTR)) 
         DO JP=1, NP_CUT
            WRITE (LIMB_UNIT, '(1P,G14.4,1X,G14.4)') 
     >             XPLOT(JP), YPLOT(JP)
         ENDDO
         CLOSE (UNIT=LIMB_UNIT)
      ENDIF

      RETURN

C***  ERROR branch ************************************************
  100 WRITE (0, 101) LIMB_LINE(:IDX(LIMB_LINE)) 
  101 FORMAT  ("*** LIMBDARK: SYNTAX ERROR with plot options" /
     >         "The error occured in the following line:" /, A, /,
     >         "Possible plot options:" /
     >         "LIMBDARK ... [PCUT=x.x] [MIN=x.x] [MU] [NORM]" /
     >         'Non-fatal error detected in subr. LIMBDARK_OUTPUT' /
     >         ' --> LIMBDARK plot option(s) ignored' )
      GOTO 110

      END
