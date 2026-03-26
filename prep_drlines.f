      SUBROUTINE PREP_DRLINES (DRLINES_CARD, NAUTO, KRUDAUT, EAUTO)
C******************************************************************
C***  Prepares which of the DRTRANSITs (i.e. stabilizing transitions
C***  from auto-ionisation levels) should be treated as: 
C***    - "normal" lines:     KRUDAUT = 0 
C***    - "rudimental" lines: KRUDAUT = 1 
C***  This depends on the CARDS line DRLINES
C***  called from: WRSTART, COLI, STEAL
C******************************************************************
      CHARACTER*80 DRLINES_CARD
      CHARACTER*4 ACTPAR
      DIMENSION KRUDAUT(NAUTO), EAUTO(NAUTO)

C***  Default: no DRLINES_CARD was encountered
      IF (DRLINES_CARD .EQ. '') THEN
C***     Only lines that are truely auto-ionizing are set rudimantel         
         DO IND=1, NAUTO
            IF (EAUTO(IND) < .0) THEN
               KRUDAUT(IND) = 0
            ELSE
               KRUDAUT(IND) = 1
            ENDIF
         ENDDO
      ELSE
         CALL SARGC (DRLINES_CARD, NPAR)
         IF (NPAR .LT. 2) GOTO 97
         CALL SARGV (DRLINES_CARD, 2, ACTPAR)

         IF (ACTPAR .EQ. 'ALL') THEN
C***                      ===
C***        All lines set non-rudimental:
            DO IND=1, NAUTO
               KRUDAUT(IND) = 0
            ENDDO
          
         ELSEIF (ACTPAR .EQ. 'NONE') THEN
C***                          ====
C***        All lines set rudimental:
         write (0,*) 'in PREP_DRLINES: parameter NONE detected'
            DO IND=1, NAUTO
               KRUDAUT(IND) = 1
            ENDDO

         ELSE
            GOTO 98
         ENDIF 
      ENDIF   

      RETURN 

C***  ERROR exits
   97 WRITE (0,*) '*** ERROR: DRLEVELS needs one parameter:' 
      WRITE (0,*) '*** allowed are: ALL or NONE'  
      GOTO 99

   98 WRITE (0,*) '*** ERROR: DRLEVELS has invalid parameter:'  
      WRITE (0,*) '*** allowed are: ALL or NONE'  
      GOTO 99

   99 WRITE (0,*) '*** the error occured on the following CARDS line:'
      WRITE (0,*) DRLINES_CARD(:IDX(DRLINES_CARD)) 
      STOP 'Fatal ERROR detected by subroutine PREP_DRLINES'

      END
