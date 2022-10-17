      SUBROUTINE LINSOL (X, A, B, N, NDIM, SCRATCH, ATEST, BTEST, DTEST, 
     >                   VERSION)
C**********************************************************************
C***  SOLVES THE LINEAR SYSTEM:  >>  X * A = B
C***  NDIM  = ROW DIMENSION  OF A
C***  N     = RANK OF THE SYSTEM
C***  A     = COEFFICIENT MATRIX 
C***  B     = RIGHT=HAND SIDE VECTOR
C***  X     = UNKNOWN VECTOR  --> SOLUTION
C***  Call from Subroutine LINPOP
C**********************************************************************

      DIMENSION A(NDIM,NDIM), X(NDIM), B(NDIM) 
      DIMENSION SCRATCH(2*NDIM)
      DIMENSION ATEST(NDIM,NDIM), BTEST(NDIM), DTEST(NDIM)
      CHARACTER*4 VERSION

C***  Operating system:
      COMMON / COMOS / OPSYS
      CHARACTER*8 OPSYS

C***  DEFINE THE INVERSION SUBROUTINE: 'FREE' OR 'OWN' OR 'OWNL'
C***  The difference between OWN and OWNL is, that in case of a singularity
C***    the Subroutine OWNINV stops for VERSION = OWN while it returns
C***    SING in the other case. LINSOL can then skip this depth point and 
C***    continues with a MODIFY.
      VERSION = 'OWNL'

C***  NUMBER OF IMPROVEMENT-ITERATIONS
      IMPMAX = 1

      IF (N .GT. NDIM) THEN
         PRINT *,'ERROR IN SUBROUTINE LINSOL: N .GT. NDIM'
         STOP 'ERROR'
         ENDIF
       
C***  SAVE MATRIX A BEFORE INVERSION
      DO 1 I=1, N
      DO 1 J=1, N
    1 ATEST(I,J) = A(I,J)

C***  MATRIX INVERSION -----------------------------
      CALL INV (N, NDIM, A, VERSION)
      IF (VERSION .EQ. 'SING') THEN
        RETURN
      ENDIF

C***  SOLUTION BY:   X = B * A-INV
      CALL VMF (X, B, A, N, NDIM)

C***  ITERATIVE IMPROVEMENT:  ----------------------------------
      DO 30 IMPRO=1, IMPMAX
C***  ACTUAL RIGHT-HAND SIDE VECTOR:  BTEST := X * A
      CALL VMF (BTEST,X,ATEST,N,NDIM)
C***  RIGHT-HAND SIDE ERROR:  BTEST := BTEST - B
      CALL VSUB (BTEST, B, N)
C***  SOLUTION OF:   DTEST * A = BTEST    BY:  DTEST = BTEST * A-INV
      CALL VMF (DTEST, BTEST, A, N, NDIM)
C***  SUBTRACTION OF THE ERROR TERM: X(NEW) := X(OLD) - DTEST
      CALL VSUB (X, DTEST, N)
   30 CONTINUE

C***  TEST PRINTOUT: RIGHT-HAND SIDE VECTOR (OPTINAL!!!)  .............
      IF (.FALSE.) THEN
      CALL VMF (BTEST, X, ATEST, N, NDIM)

      PRINT 11, VERSION, IMPMAX
   11 FORMAT (//, 10X, 'VERSION: ', A, 10X, 'ITERATIONS:',I2,
     $  //,10X,'INDEX', 11X, 'B(I)', 6X,'BTEST(I)',7X,'DEVIATION')

      DO 2 I=1,N
      DIFF = B(I) - BTEST(I)
    2 PRINT 3, I, B(I), BTEST(I), DIFF
    3 FORMAT (10X, I5, 3(3X,G12.3))
      ENDIF
C......................................................................

      RETURN
      END
