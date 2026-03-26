      SUBROUTINE INTEGRATE (SUM, X, Y, NDAT, XMIN, XMAX)
C**********************************************************************
C***  Integrates the tabulated function Y(X) 
C***  in the interval (XMIN, XMAX) 
C***  the X-values must be monotonically ascending 
C**********************************************************************

      DIMENSION X(NDAT), Y(NDAT)
      LOGICAL ASCENDING

C***  make sure that XMIN, XMAX are sorted in same sense as the X-values
      ASCENDING = X(1) .LT. X(NDAT) 

      IF (.NOT. ASCENDING) GOTO 90

      IF (XMIN .GE. XMAX) GOTO 91

      XMIN_LOC = MAX (XMIN, X(1))
      XMAX_LOC = MIN (XMAX, X(NDAT))

C***  FIND INDEX RANGE INSIDE (XMIN,XMAX)
      IMIN = ISRCHFGT (NDAT, X, 1, XMIN)
      IF (IMIN .EQ. 0) GOTO 92
      IMAX = ISRCHFGT (NDAT, X, 1, XMAX) - 1
      IF (IMAX .EQ. 0) GOTO 93

      SUM = .0

C***  First (fractional) interval (XMIN, X(IMIN))
      CALL SPLINPO (YMIN, XMIN, Y, X, NDAT)
      SUM = SUM + 0.5 * (X(IMIN)-XMIN) * (Y(IMIN)+YMIN)

C***  Last (fractional) interval (X(IMAX), XMAX)
      CALL SPLINPO (YMAX, XMAX, Y, X, NDAT)
      SUM = SUM + 0.5 * (XMAX-X(IMAX)) * (YMAX+Y(IMAX))

C***  No inner intervals?
      IF (IMIN .GE. IMAX) RETURN

C***  Loop 
      DO I= IMIN, IMAX
         IF (I .EQ. IMIN) THEN
            W =  0.5 * (X(IMIN+1)-X(IMIN)) 
         ELSEIF (I .EQ. IMAX) THEN
            W = 0.5 * (X(IMAX)-X(IMAX-1))
         ELSE
            W = 0.5 * (X(I+1)-X(I-1)) 
         ENDIF

         SUM = SUM + W * Y(I)
      ENDDO

      RETURN

C***  ERROR messages **************************************************

   90 WRITE (0,*) 'BRANCH FOR DESCENDING ORDER NOT YET CODED'
      GOTO 99

   91 WRITE (0,*) 'Integration boundaries XMIN .GE. XMAX'
      GOTO 99

   92 WRITE (0,*) 'Index search failed for XMIN'
      GOTO 99

   93 WRITE (0,*) 'Index search failed for XMAX'
      GOTO 99

   99 WRITE (0,*) 'FATAL ERROR in subr. INTEGRATE:'
      CALL TRBK

      END
