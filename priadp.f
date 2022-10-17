      SUBROUTINE PRIADP (MODHEAD, OLDHEAD, NEWATOM, LEVEL, NTRANS,
     >                   ADPWEIGHT, N, NOLD)
C**********************************************************************
C***  PRINTOUT FOR PROGRAM 'ADAPTER'
C**********************************************************************

      CHARACTER*100 MODHEAD, OLDHEAD
      CHARACTER*10 LEVEL(N), MESSAGE(N)*4
      DIMENSION NTRANS(N), ADPWEIGHT(N)
      LOGICAL NEWATOM

C**********************************************************************
C***  PRINT HEADER                                                  ***
C**********************************************************************
      PRINT 11, MODHEAD, OLDHEAD
   11 FORMAT (1H1,/,1X,A,20X,'JOB NO. 1A',///,50X,
     $  'A D A P T E R',/,50X,13('='),//,
     $  10X,'POPULATION NUMBERS TAKEN FROM OLD MODEL:',/,10X,A,//)

C**********************************************************************
C***  PRINTOUT OF LEVEL CORRESPONDENCE (IN 5 COLUMNS), IF CHANGED   ***
C**********************************************************************
      IF (NEWATOM) THEN
         PRINT 13
   13    FORMAT (5(6X,'OLD  NEW      LEVEL'),/)
         DO 30 K=1, N
           IF (ADPWEIGHT(K) .NE. 1.) THEN
              WRITE (MESSAGE(K), '(1X,F3.1)') ADPWEIGHT(K)
           ELSE
              MESSAGE(K) = '    '
           ENDIF
   30    CONTINUE
C***     NUMBER OF ROWS REQUIRED:
         NROW= (N-1)/5 + 1
         DO 14 IROW=1,NROW
            ILAST=MIN0(N,IROW+NROW*4)
            PRINT 15, (NTRANS(I), I, LEVEL(I),MESSAGE(I),
     >              I=IROW,ILAST,NROW)
   15       FORMAT (4X,5(1X,I4,I5,1X,A10,A4))
   14    CONTINUE
      ELSE IF (N .LT. NOLD) THEN
         WRITE (*,*) '*************************************************'
         WRITE (0,'(A,I4,A,I4)')  'WARNING: N0', N, '   NOLD= ', NOLD
         WRITE (*,'(A,I4,A,I4)')  'WARNING: N0', N, '   NOLD= ', NOLD
         WRITE (*,*) 'NEW MODEL HAS LESS LEVELS THAN OLDSTART MODEL'
         WRITE (*,*) '*************************************************'
      ENDIF

      PRINT 97
   97 FORMAT (//)


      RETURN
      END
