      SUBROUTINE SECONDMODEL_DEFINE (SECONDMODEL_LINE, NMOD, 
     >             SECONDMODEL_PATH, SECONDMODEL_CHANGED,
     >             RMAX, SECMOD_RRANGE, IGNORE_DIFF)
C******************************************************************
C***  Reads the pathname from the SECONDMODEL_LINE, and checks whether
C***  this name has changed or the option is swiched OFF
C***  Moreover, the radius-range in which the radius-grid needs to be 
C***  refined is estimated here:
C***  For CONE: 
C***      SECMOD_RRANGE = [RMAX, 1]
C***  For SPHERE:
C***      SECMOD_RRANGE = [DSPHERE+RSPHERE, DSPHERE-RSPHERE]
C******************************************************************

      CHARACTER SECONDMODEL_LINE*(*), ACTPAR*100, ACTPAR2*100
      CHARACTER SECONDMODEL_PATH*(*), NEW_SECONDMODEL_PATH*400
      CHARACTER SHAPE*10
      LOGICAL SECONDMODEL_CHANGED, IGNORE_DIFF
      DIMENSION SECMOD_RRANGE(2)

      DATA NEW_SECONDMODEL_PATH / 'UNDEFINED' /

C***  Secondmodel active
      NMOD = 2

C***  Decode parameters from input line
      CALL SARGC (SECONDMODEL_LINE, NPAR)
      DO IPAR=1, NPAR
         CALL SARGV (SECONDMODEL_LINE, IPAR, ACTPAR)

         IF (ACTPAR .EQ. 'SHAPE') THEN
C*                        =====
            IF (NPAR .LT. IPAR+1) THEN
               GOTO 90
            ELSE
               CALL SARGV (SECONDMODEL_LINE, IPAR+1,SHAPE)
            ENDIF

         ELSEIF (ACTPAR .EQ. 'RSPHERE') THEN
C*                            =====
            IF (NPAR .LT. IPAR+1) THEN
               GOTO 90
            ELSE
               CALL SARGV (SECONDMODEL_LINE, IPAR+1,ACTPAR2)
               READ (ACTPAR2, '(F10.0)', ERR=98) RSPHERE
            ENDIF

         ELSEIF (ACTPAR .EQ. 'DSPHERE') THEN
C*                            =====
            IF (NPAR .LT. IPAR+1) THEN
               GOTO 90
            ELSE
               CALL SARGV (SECONDMODEL_LINE, IPAR+1,ACTPAR2)
               READ (ACTPAR2, '(F10.0)', ERR=98) DSPHERE
            ENDIF

         ELSEIF (ACTPAR .EQ. 'PATH') THEN
C*                            ====
            IF (NPAR .LT. IPAR+1) THEN
               GOTO 90
            ELSE
               CALL SARGV (SECONDMODEL_LINE, IPAR+1, 
     >                     NEW_SECONDMODEL_PATH)
            ENDIF

         ELSEIF (ACTPAR(:11) .EQ. 'IGNORE_DIFF') THEN
C*                                 ===========
            IGNORE_DIFF = .TRUE.

         ELSEIF (ACTPAR .EQ. 'OFF') THEN
C*                            ===
            NMOD = 1
         ENDIF
      ENDDO         

C***  Check for mandatory parameters
      IF (NEW_SECONDMODEL_PATH .EQ. 'UNDEFINED') THEN
         GOTO 901
      ELSEIF (SHAPE .EQ. 'UNDEFINED') THEN
         GOTO 92
      ELSEIF (SHAPE .EQ. 'SPHERE') THEN
         IF (RSPHERE .EQ. -999.) GOTO 97
         IF (DSPHERE .EQ. -999.) GOTO 971
         IF (RSPHERE .GE. DSPHERE+RMAX) GOTO 972
      ENDIF

      IF (SHAPE .EQ. 'CONE') THEN
         SECMOD_RRANGE (1) = RMAX 
         SECMOD_RRANGE (2) = 1. 
      ELSEIF (SHAPE .EQ. 'SPHERE') THEN
         SECMOD_RRANGE (1) = MIN(RMAX, DSPHERE+RSPHERE) 
         SECMOD_RRANGE (2) = MAX(1.  , DSPHERE-RSPHERE) 
      ENDIF

C***  Check if PATH has changed
      IF (NEW_SECONDMODEL_PATH .NE. SECONDMODEL_PATH) THEN
         SECONDMODEL_CHANGED = .TRUE.
         SECONDMODEL_PATH = NEW_SECONDMODEL_PATH
      ELSE
         SECONDMODEL_CHANGED = .FALSE.
      ENDIF

      RETURN

C*********************************************************************
C***  ERROR branches  ************************************************
C*********************************************************************
   90 WRITE (0,*) '*** ERROR: Option ', ACTPAR(:IDX(ACTPAR)), 
     >            ' needs a value (keyword)'
      GOTO 99

  901 WRITE (0,*) '*** ERROR: PATH to the SECONDMODEL is not defined'
      GOTO 99

   92 WRITE (0,*) '*** ERROR: mandatory parameter SHAPE is missing'
      GOTO 99

   97 WRITE (0,*) '*** ERROR: mandatory parameter RSPHERE is missing'
      GOTO 99

  971 WRITE (0,*) '*** ERROR: mandatory parameter DSPHERE is missing'
      GOTO 99

  972 WRITE (0,*) '*** ERROR: SPHERE covers the whole atmosphere'
      GOTO 99

   98 WRITE (0,*) '*** ERROR: Parameter cannot be decoded as number:'
      WRITE (0,*) '*** ERROR: ', ACTPAR(:IDX(ACTPAR)), '=', 
     >                           ACTPAR2(:IDX(ACTPAR2))
      GOTO 99


   99 STOP '*** FATAL ERROR in subroutine SECONDMODEL_DEFINE'

      END
