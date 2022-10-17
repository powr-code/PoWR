      SUBROUTINE LINSOL_SPLIT (X, A, B, N, NDIM,
     >                         NFIRST, NLAST, NATOM, AOFF,
     >                  SCRATCH, ATEST, BTEST, DTEST, VERSION, ELEMENT)
C**********************************************************************
C***  SOLVES THE LINEAR SYSTEM:  >>  X * A = B
C***    by using the atomic block structure of A
C***    (this means of course that A must be block-diagonal)
C***  NDIM  = ROW DIMENSION  OF A
C***  N     = RANK OF THE SYSTEM
C***  A     = COEFFICIENT MATRIX 
C***  B     = RIGHT=HAND SIDE VECTOR
C***  X     = UNKNOWN VECTOR  --> SOLUTION
C***  Call from Subroutine LINPOP
C**********************************************************************

C     This version is based on the ideas from the subroutine 
C      linsol_split in the Goetz branch. However, there are certain
C      things in Goetz code that are not in the wrh version. These are:

C     BFAST - Logical to switch between two methods, today always
C        set to BFAST = .TRUE. because BFAST = .FALSE. simply means
C        that no block inversion is used
C
C     B_SKIP_POP(N) - Array of Logicals
C     NSTART(100), NSTOP(100)
C        Allows to skip certain levels, therefore NSTART and NSTOP
C        had to be introduced to store the "reduced" block start
C        This version works like all entries of B_SKIP_POP are .FALSE.
C        which is equal to 
C      NSTART(NA)= NFIRST(NA)
C      NSTOP(NA) = NLAST(NA)
C        for all entries in the Goetz branch. However, when it comes
C        to the matrix inversion part, NSTART and NSTOP are not used,
C        so this stuff might have been discarded

C**********************************************************************

      IMPLICIT NONE

      INTEGER, INTENT(IN) :: NDIM, N, NATOM
      REAL, DIMENSION(NDIM,NDIM), INTENT(INOUT) :: A, AOFF
      REAL, DIMENSION(NDIM), INTENT(IN) :: B
      REAL, DIMENSION(NDIM), INTENT(OUT) :: X
      
C*** As we do not know the size of the blocks we init SCRATCH as a vector
C     (Note: SCRATCH is used for one block in the matrix, this can be
C            up to the dimension of the original matrix)
      REAL, DIMENSION(NDIM*NDIM), INTENT(INOUT) :: SCRATCH

      REAL, DIMENSION(NDIM,NDIM), INTENT(OUT) :: ATEST
      REAL, DIMENSION(NDIM), INTENT(OUT) :: BTEST, DTEST
      INTEGER, DIMENSION(NATOM), INTENT(IN) :: NFIRST, NLAST
      CHARACTER(4), INTENT(INOUT) :: VERSION
      CHARACTER*10 ELEMENT(NATOM)

      INTEGER :: IMPMAX, na, nn, i, j, nBlockDim, 
     >           indBlockI, indBlockJ, nFirstNA, nLastNA,
     >           nVectorIndex, IMPRO, iPrimat
      REAL :: DIFF

      LOGICAL :: bForceBlock, bPrintRightVector

C***  Operating system:
      COMMON / COMOS / OPSYS
      CHARACTER(8) :: OPSYS

C***  set iPrimat to zero to print the inverted matrix structure
      iPrimat = 1

      bForceBlock = .TRUE.          !Force block matrix structure 
      bPrintRightVector = .FALSE.   !Print B, BTEST and difference (debug option)


C***  DEFINE THE INVERSION SUBROUTINE: 'FREE' OR 'OWN' OR 'OWNL'
C***  The difference between OWN and OWNL is, that in case of a singularity
C***    the Subroutine OWNINV stops for VERSION = OWN while it returns
C***    SING in the other case. LINSOL can then skip this depth point and 
C***    continues with a MODIFY.
      VERSION = 'OWNL'

C***  NUMBER OF IMPROVEMENT-ITERATIONS
      IMPMAX = 1

      IF (N > NDIM) THEN
         PRINT *,'ERROR IN SUBROUTINE LINSOL: N .GT. NDIM'
         STOP 'ERROR'
      ENDIF
      
C***  SAVE MATRIX A BEFORE INVERSION 
C***  (OPTIONAL: ENSURE BLOCK STRUCTURE BY forceblock = .TRUE.)
C***  08.06.2010 forceblock MUST be used to avoid SUM-Error at some depth points
C     Side note: Only if GAMMA=0 has been set in the CARDS file
C      the matrix DM is fully block-diagonal. In all other cases
C      there are (usually only few) non-zero elements outside of
C      the blocks. There forceblock = .TRUE. should always be used
C      as otherwise these elements would be simply copied (without
C      any inversion) to the inverted matrix and could cause
C      strange errors
      IF (bForceBlock) THEN
C***     init ATEST as identity matrix
         ATEST = 0.
C         DO nn=1, NLAST(NATOM)
         DO nn=1, NDIM
            ATEST(nn, nn) = 1.
         ENDDO
C***     copy atom blocks
         DO na=1, NATOM
            nFirstNA = NFIRST(na)
            nLastNA = NLAST(na)
            nBlockDim = nLastNA - nFirstNA + 1
            DO j=1, nBlockDim
               indBlockJ = nFirstNA + j - 1
               DO i=1, nBlockDim
                  indBlockI = nFirstNA + i - 1
                  ATEST(indBlockI, indBlockJ) = A(indBlockI, indBlockJ)
               ENDDO
            ENDDO
         ENDDO
C***     copy additional diagonal elements (if there are any)
         IF (NDIM .GT. NLAST(NATOM)) THEN
            DO nn=NLAST(NATOM)+1, NDIM
               ATEST(nn, nn) = A(nn, nn)
            ENDDO
         ENDIF

C***  PRINT Matrix debug option
         IF (iPrimat .eq. 0) THEN
            write (*,*) 'LINSOL_SPLIT> PRIMAT'
            CALL PRIMAT(A, N, NDIM, 'Matrix A')
            write (*,*) 'LINSOL_SPLIT> end of PRIMAT'
            iPrimat = 1
         ENDIF
C***  end of block structure test

         DO i=1, NDIM
            DO j=1, NDIM
               A(i,j) = ATEST(i,j)
            ENDDO
         ENDDO
      ELSE      
         DO i=1, N
            DO j=1, N
               ATEST(i,j) = A(i,j)
            ENDDO
         ENDDO
      ENDIF

C***  MATRIX INVERSION -------------------------------
      DO na=1, NATOM         
         SCRATCH = 0.
C***     extract blocks
         nFirstNA = NFIRST(na)
         nLastNA = NLAST(na)
         nBlockDim = nLastNA - nFirstNA + 1
         DO j=1, nBlockDim
            indBlockJ = nFirstNA + j - 1
            DO i=1, nBlockDim
               indBlockI = nFirstNA + i - 1
C***           store column vectors of current block in one giant vector
C***             because the dimension of the blocks is unknown
C***             when SCRATCH has to be initialized
C***           subroutine INV will treat SCRATCH as a Matrix
               nVectorIndex = i + (j - 1) * nBlockDim
               SCRATCH(nVectorIndex) = A(indBlockI, indBlockJ)
            ENDDO
         ENDDO

C***     call INV and cancel everything if singularity encounters
C***       => returns inverted block in SCRATCH
         CALL INV (nBlockDim, nBlockDim, SCRATCH, VERSION)

C***     replace original block with inverted block in the matrix
         DO j=1, nBlockDim
            indBlockJ = nFirstNA + j - 1
            DO i=1, nBlockDim
               indBlockI = nFirstNA + i - 1
               nVectorIndex = i + (j - 1) * nBlockDim
               A(indBlockI, indBlockJ) = SCRATCH(nVectorIndex)
            ENDDO
         ENDDO

         IF (VERSION .EQ. 'SING') THEN
            write (0,*) 'SING in block for element ' // ELEMENT(NA)
            RETURN
         ENDIF
      ENDDO

C***  invert additional diagonal elements
      IF (NDIM .GT. NLAST(NATOM)) THEN
         DO nn=NLAST(NATOM)+1, NDIM
            A(nn, nn) = 1. / A(nn, nn)
         ENDDO
      ENDIF
C***  End of MATRIX INVERSION ---------------------------


C***  -------/ the following is unchanged from LINSOL /---------

C***  SOLUTION BY:   X = B * A-INV
      X = MATMUL(B, A)

C***  ITERATIVE IMPROVEMENT: ("Nachbrenner") ------------------------
      DO IMPRO=1, IMPMAX
        BTEST = MATMUL(X, ATEST)  !ACTUAL RIGHT-HAND SIDE VECTOR:  BTEST := X * A
        BTEST = BTEST - B         ! RIGHT-HAND SIDE ERROR:  BTEST := BTEST - B
        DTEST = MATMUL(BTEST, A)  !SOLUTION OF:   DTEST * A = BTEST    BY:  DTEST = BTEST * A-INV
        X = X - DTEST             !SUBTRACTION OF THE ERROR TERM: X(NEW) := X(OLD) - DTEST
      ENDDO
      
C***  TEST PRINTOUT: RIGHT-HAND SIDE VECTOR (OPTINAL!!!)  .............
      IF (bPrintRightVector) THEN
        !CALL VMF (BTEST, X, ATEST, N, NDIM)
        BTEST = MATMUL(X, ATEST)

        PRINT 11, VERSION, IMPMAX
   11   FORMAT (//, 10X, 'VERSION: ', A, 10X, 'ITERATIONS:',I2,
     $    //,10X,'INDEX', 11X, 'B(I)', 6X,'BTEST(I)',7X,'DEVIATION')


    3   FORMAT (10X, I5, 3(3X,G12.3))
        DO i=1,N
          DIFF = B(i) - BTEST(i)
          PRINT 3, i, B(i), BTEST(i), DIFF
        ENDDO
      ENDIF
C......................................................................

      RETURN
      END
