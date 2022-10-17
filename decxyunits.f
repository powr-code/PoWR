      SUBROUTINE DECXYUNITS (LINE, IARG, B, C, VALUE)
C**********************************************************************
C***  Decodes a floting-point type argument from LINE
C***  ALLOWS THE USE OF THE "SAME" KEYWORD, "VARIABLES"
C***  XMIN, XMAX, YMIN, YMAX INSTEAD OF NUMERICAL VALUES
C***  and the prefixes LOG and DEX and is therefore good for X, Y in units
C**********************************************************************
      CHARACTER LINE*(*), ACTPAR*40
      DIMENSION B(6), C(6)

      COMMON / COMCMPT / CMPT
      COMMON / COMCLIP / ICLIP, CXMIN(2), CXMAX(2), CYMIN(2), CYMAX(2)
      COMMON / COMPLT / KPS, XOFFSET, YOFFSET, FORMATFA

      CALL SARGV (LINE, IARG, ACTPAR)
      LAST = IDX(ACTPAR)
      
      IF (ACTPAR(:4) .EQ. 'SAME') RETURN
      
      IF (ACTPAR .EQ. 'XMIN') THEN
         VALUE = B(2)
      ELSE IF (ACTPAR .EQ. 'XMAX') THEN
         VALUE = B(3)
      ELSE IF (ACTPAR .EQ. 'YMIN') THEN
         VALUE = C(2)
      ELSE IF (ACTPAR .EQ. 'YMAX') THEN
         VALUE = C(3)
      ELSE IF (ACTPAR .EQ. 'PXMIN') THEN
         XPT = 1.
         VALUE = (XPT/CMPT - XOFFSET) / (FORMATFA * B(1)) + B(2)
      ELSE IF (ACTPAR .EQ. 'PXMAX') THEN
         XPT = CXMAX(1) - 1.
         VALUE = (XPT/CMPT - XOFFSET) / (FORMATFA * B(1)) + B(2)
      ELSE IF (ACTPAR .EQ. 'PYMIN') THEN
         YPT = 1.
         VALUE = (YPT/CMPT - YOFFSET) / (FORMATFA * C(1)) + C(2)
      ELSE IF (ACTPAR .EQ. 'PYMAX') THEN
         YPT = CYMAX(1) - 1.
         VALUE = (YPT/CMPT - YOFFSET) / (FORMATFA * C(1)) + C(2)
      ELSE IF (ACTPAR(1:3) .EQ. 'LOG') THEN
         READ (ACTPAR(4: ), '(F10.0)') VALUE
         VALUE = ALOG10 (VALUE)
      ELSE IF (ACTPAR(1:3) .EQ. 'DEX') THEN
         READ (ACTPAR(4: ), '(F10.0)') VALUE
         VALUE = 10. ** (VALUE)
      ELSE IF (ACTPAR(1:1) .EQ. '(' .AND. 
     >        ACTPAR(LAST:LAST) .EQ. ')' ) THEN 
C***     Arithmetic expression
         CALL ARIPAR (ACTPAR, VALUE, IERR)
         IF (IERR .NE. 0) GOTO 99
      ELSE
         READ (ACTPAR, '(F20.0)', ERR=99) VALUE
      ENDIF

      RETURN

   99 PRINT *, '*** ERROR IN KASDEF OPTION: ***'
      PRINT *, '>>> ', ACTPAR(:IDX(ACTPAR)), 
     >         ' <<<< CANNOT BE DECODED AS FLOATING POINT NUMBER!'
      PRINT *, 'THE ERROR OCCURED IN THE FOLLOWING LINE:'
      PRINT *, LINE(:IDX(LINE))
      STOP '*** ERROR IN SUBROUTINE DECXYUNITS ***'

      END

