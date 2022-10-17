      SUBROUTINE OWNINV (N, NDIM, A, CKEY)
C*******************************************************************************
C***  MATRIX INVERSION (CRAY FORTRAN)
C***  MATRIX A
C***  N = RANK OF SUBMATRIX (LEFT UPPER BLOCK) TO BE INVERTED
C***  NDIM = ROW DIMENSION OF TOTAL MATRIX
C***  N .LE. NDIM .LE. NMAX
C***  Most time consuming routine in STEAL 
C***  in WRSTART-Job : (69 percent), (Loop 80: 66 percent)
C***  Most time consuming routine in WRCONT
C***                   (80 percent), (Loop 80: 76 percent)
C***  Tested on 30-Jan-1997 20:54:25, Lars
C***
C***  Now a Normalization of the rows and columns is possible (BNORM = .TRUE.)
C***  Row-Norm is Quadratic Norm
C***  Column-Norm is Maximum Norm
C***
C***
C*******************************************************************************

      LOGICAL BNORM

      PARAMETER (BNORM = .TRUE.)
C      PARAMETER (NMAX = 560)  old version before split
      PARAMETER (NMAX = 2002)

      DIMENSION A(NDIM,NDIM),IK(NMAX),JK(NMAX)
      DIMENSION CNORM(NMAX), RNORM(NMAX)
      CHARACTER*4 CKEY

C***  OUTPUT of matrix dimension parameters (testing only)
C      WRITE (0,*) 'OWNINV> N=', N, 'NDIM=', NDIM

      IF (N .GT. NDIM) THEN
            CALL REMARK ('N .GT. NDIM')
            STOP 'ERROR IN OWNINV'
            ENDIF
      IF (NDIM .GT. NMAX) THEN
            CALL REMARK ('NDIM .GT. NMAX')
            WRITE (0,'(A,i4,1x,a,i4)') 'NDIM=',NDIM, 'NMAX=',NMAX
            STOP 'ERROR IN OWNINV'
            ENDIF

C*** New Branch to normalize Columns and Rows
C*** Calculation of the Row-Norm (Quadratic Norm)
      IF (BNORM) THEN
        DO K=1, N
          R = 0.
          DO L=1, N
            XM = A(L,K)
            IF (ABS(XM) .GT. 1.E100) THEN
              R = 1.
              EXIT
            ENDIF
            R = R + XM*XM
          ENDDO
          RNORM(K) = 1. / SQRT(R)
        ENDDO

C*** Calculation of the Column-Norm (Maximum Norm)
        DO L=1, N
          C = A(L,1)
          DO K=2, N
            XM = ABS(A(L,K))
            IF (XM .GT. C) C = XM
          ENDDO
          CNORM(L) = 1. / C
        ENDDO

C*** Applikation of CNORM and RNORM
        DO K=1, N
          DO L=1, N
            A(L,K) = A(L,K) * CNORM(L) * RNORM(K)
          ENDDO
        ENDDO
      ENDIF

      DO 100  K=1,N
C
C     SUCHE MAXIMALES MATRIXELEMENT
      AMAX=0.
      ABSAMAX=0.
      DO 30 J=K,N
        L=N+1-K
        IMAX=MY_ISAMAX(L,A(K,J),1)+K-1
        IF(ABSAMAX.GT.ABS(A(IMAX,J))) GOTO 30
        AMAX=A(IMAX,J)
        ABSAMAX=ABS(AMAX)
        IK(K)=IMAX
        JK(K)=J
   30 CONTINUE
C
C***  WARNING IN CASE OF SINGULARITY (PIVOT ELEMENT = 0 )
C***  Exit the Routine with CKEY = 'SING'
      IF (AMAX .EQ. .0) THEN
        CALL REMARK ('Subr. OWNINV: Singularity discovered')
        IF (CKEY .EQ. 'OWNL') THEN
          CKEY = 'SING'
          RETURN
        ELSE
          STOP 'ERROR in Subr. OWNINV'
        ENDIF
      ENDIF
C
C     VERTAUSCHEN DER ZEILEN I UND K
      I=IK(K)
      IF(I.EQ.K) GOTO 51
      DO 50 J=1,N
      SAVE =A(K,J)
      A(K,J)=A(I,J)
   50 A(I,J)=-SAVE
C
C     VERTAUSCHEN DER SPALTEN J UND K
   51 J=JK(K)
      IF(J.EQ.K) GOTO 61
      DO 60 I=1,N
      SAVE = A(I,K)
      A(I,K)=A(I,J)
   60 A(I,J)=-SAVE
C
C     DIVISION DER SPALTE K DURCH AMAX
   61 DO 70  I=1,N
      A(I,K)=-A(I,K)/AMAX
   70 CONTINUE
      A(K,K)=0.
C
C     UMFORMEN DER ZEILEN UND SPALTEN
C***  ELEMENTARE UMFORMUNG:  ZEILE I = ZEILE I + A(I,K) * ZEILE K
      DO 81 I=1,N
        IF (I .EQ. K) GOTO 81
        AIK=A(I,K)
        DO 80 J=1,N
          A(I,J)=A(I,J)+AIK*A(K,J)
   80 CONTINUE
C!!!  Tried to replace Loop 80
C!!!  It works, but it is slower!
c      CALL DGEMA ('N', 'N', 1, N, 1., A(I,1), NDIM, AIK, A(K,1), NDIM, 
c     >            A(I,1), NDIM)
   81 CONTINUE
C!!!  Tried to replace Loops 80 and 81 by
C!!!  It works, but it is also slower!
c      FORALL (I=1:N, J=1:N, I/=K) A(I,J) = A(I,J) + A(I,K)*A(K,J)


C
C***  SPALTE K: DIVISION DUCH AMAX
      DO 90 J=1,N
      A(K,J)=A(K,J)/AMAX
  90  CONTINUE
C
C***  DIAGONALELEMENT
      A(K,K)=1./AMAX
  100 CONTINUE
C
C
C     ZEILEN UND SPALTEN RUECKTAUSCHOPERATIONEN
      DO 130 L=1,N
      K=N-L+1
      J=IK(K)
      IF(J.LE.K) GOTO 111
C***  VERTAUSCHEN DER SPALTEN J UND K
      DO 110 I=1,N
       SAVE=A(I,K)
      A(I,K)=-A(I,J)
  110 A(I,J)=SAVE
  111 I=JK(K)
      IF(I .LE. K) GOTO 130
C***  VERTAUSCHEN DER ZEILEN I UND K
      DO 120 J=1,N
      SAVE=A(K,J)
      A(K,J)=-A(I,J)
  120 A(I,J)=SAVE
  130 CONTINUE
C

C*** Applikation of RNORM and CNORM
      IF (BNORM) THEN
        DO K=1, N
          DO L=1, N
            A(L,K) = A(L,K) * RNORM(L) * CNORM(K)
          ENDDO
        ENDDO
      ENDIF

      RETURN
      END
