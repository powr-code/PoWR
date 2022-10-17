      SUBROUTINE PRIMAT (A, N, NDIM, NAME)
C**********************************************************************
C***  TEST PRINTOUT OF MATRIX A
C**********************************************************************

      DIMENSION A(NDIM, 1)
      PARAMETER (MAXCOLUMNS = 126)
      CHARACTER*1 C(MAXCOLUMNS)
      CHARACTER*(*) NAME

      PRINT 10, NAME, N, N
   10 FORMAT (//, 1X, 40('-'), '   MATRIX: ',A,'(',I3,',',I3,')   ',
     $        40('-'),//)

      jStart = 1
      jEnd = MIN(N, MAXCOLUMNS)
      NC = jEnd

    7 CONTINUE  
      PRINT 1, (J/100, J=jStart, jEnd)
      PRINT 1, ((J-J/100*100)/10, J=jStart, jEnd)
      PRINT 1, (J-J/10*10, J=jStart, jEnd)
    1 FORMAT (5X,126I1)
      PRINT 2
    2 FORMAT (1X)

      DO 3 I=1, N

      DO 5 J=jStart, jEnd
      IF (A(I,J) .EQ. .0) THEN
         C(J + 1 - jStart)=' '
      ELSE IF (A(I,J) .EQ. 1.) THEN
         C(J + 1 - jStart)='1'
      ELSE IF (A(I,J) .LT. 0) THEN
         C(J + 1 - jStart)='-'
      ELSE
         C(J + 1 - jStart)='+'
      ENDIF
    5 CONTINUE

      PRINT 4, I, (C(J), J=1, NC)
    4 FORMAT (I4, 1X, 126 A1)

    3 CONTINUE

      IF (N .GT. jEnd) THEN 
         jStart = jStart + MAXCOLUMNS
         jEnd = MIN(jEnd + MAXCOLUMNS, N)
         NC = jEnd - jStart + 1
         GOTO 7
      ENDIF

      RETURN
      END
