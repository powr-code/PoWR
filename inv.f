      SUBROUTINE INV (N,NDIM,A,CKEY)
C**********************************************************************
C***  MATRIX INVERSION BY CALLING THE CFT-LIBRARY SUBR. MINV
C***  A      = MATRIX
C***  N      = RANK OF SUBMATRIX (UPPER-LEFT CORNER) TO BE INVERTED
C***  NDIM   = ROW DIMENSION OF TOTAL MATRIX
C***  NMAX   = MAXIMUM VALUE OF NDIM
C**********************************************************************

      PARAMETER (NMAX=2000)
      DIMENSION A(NDIM,NDIM)
      CHARACTER*4 CKEY

C***  Dimensum must be greater-equal 2*NMAX
C***  If DEC-DXML Routines are used, it should be 64*NPDIM
C***  Here NPDIM=94 is assumed ==> 94 * 64 = 6016
C***  Or, for solving the statistical equations NDIM=300 ==> 300 * 64 = 19200
      PARAMETER (NDIMSC = 19200)
      DIMENSION SCRATCH(NDIMSC)

C***  Array for special for the Use in DEC-DXML-Routines 
C***  DGETRF and DGETRI
      DIMENSION IPIV(NMAX)

C***  Operating system:
      COMMON / COMOS / OPSYS
      CHARACTER*8 OPSYS

C***  OUTPUT of matrix dimension parameters (testing only)
ccc      write (0,*) '**** test INV> N=', N, 'NDIM=', NDIM

      IF (N .GT. NDIM) THEN
         PRINT *,'ERROR IN SUBROUTINE INV: N .GT. NDIM'
         CALL REMARK ('ERROR IN SUBROUTINE INV: N .GT. NDIM')
         STOP 'ERROR'
         ENDIF
      IF (NDIM .GT. NMAX) THEN
         WRITE (0,'(A)') '*** FATAL ERROR: NDIM .GT. NMAX'
         WRITE (0,'(A, I4, A, I4)') 'NDIM=', NDIM, ' .GT. NMAX=', NMAX
         STOP 'ERROR IN SUBROUTINE INV'
         ENDIF

C***  disable owninv for testing
      IF (CKEY(1:3) .EQ. 'OWN' .OR. CKEY(1:4) .EQ. 'OWNL') THEN
        CALL OWNINV(N, NDIM, A, CKEY) 
      ELSE
          if (ndim .gt. 300) then
            write (0,*) '*** obsolete program branch'
            write (0,*) '*** check DIMENSION SCRATCH'
            stop '*** FATAL ERROR in SUBR. INV' 
          endif 
          CALL DGETRF(N, N, A, NDIM, IPIV, INFO)
          CALL DGETRI(N, A, NDIM, IPIV, SCRATCH, NDIMSC, INFO)
      ENDIF

C***  Nur zum Loesen ginge auch folgendes
C***     Inverse Matrix stuende dann aber nicht zur Verfuegung
C***     Wert von Parameter 8 (LDB) ist unklar
C!!!          CALL DGETRS('N', N, N, A, NDIM, IPIV, B, NDIM, INFO)
C!!!          WRITE (0,*) 'INFO, SCRATCH(1)=', INFO, SCRATCH(1)

      RETURN
      END
