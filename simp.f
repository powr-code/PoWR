      SUBROUTINE SIMP(N,X,Y,F)
C************************************************************************
C***  SIMPSON-INTEGRATION OF THE FUNCTION ( N POINTS X(I),Y(I) )
C***  THE INTERVALLS X(I)-X(I-1) CAN BE DIFFERENT
C***  THE X(I)'S MUST BE SORTED (UP OR DOWN)
C************************************************************************
 
      DIMENSION X(N),Y(N)
      INTEGER ON

      IF (N .LT. 3) STOP 'SIMP'

      F=0.
      ON=INT((N-1)/2)

C***  INTERVALL I-1,I+1

      DO 1 J=1,ON
      Y13=Y(2*J-1)-Y(2*J+1)
      Y23=Y(2*J)-Y(2*J+1)
      X12=X(2*J-1)-X(2*J)
      X13=X(2*J-1)-X(2*J+1)
      X23=X(2*J)-X(2*J+1)

C***  COEFFICIENTS OF THE POLYNOM A2*X*X+A1*X+A0

      A2=Y13/(X12*X13)-Y23/(X12*X23)
      A1=Y13/X13-A2*(X(2*J-1)+X(2*J+1))
      A0=Y(2*J+1)-(A2*X(2*J+1)+A1)*X(2*J+1)

C***  INTEGRAL OF THE INTERVAL I-1,I+1

      XQ1=X(2*J-1)*X(2*J-1)
      XC1=XQ1*X(2*J-1)
      XQ3=X(2*J+1)*X(2*J+1)
      XC3=XQ3*X(2*J+1)
      F=F+A2*(XC3-XC1)/3+A1*(XQ3-XQ1)/2+A0*(X(2*J+1)-X(2*J-1))

    1 CONTINUE

C***  IF N EVEN, ADD PART OF THE INTERVAL N-1,N 

      IF ((2*ON+1).EQ.N) GOTO 3

      Y13=Y(N-2)-Y(N)
      Y23=Y(N-1)-Y(N)
      X12=X(N-2)-X(N-1)
      X13=X(N-2)-X(N)
      X23=X(N-1)-X(N)

C***  COEFFICIENTS OF THE POLYNOM A2*X*X+A1*X+A0

      A2=Y13/(X12*X13)-Y23/(X12*X23)
      A1=Y13/X13-A2*(X(N-2)+X(N))
      A0=Y(N)-(A2*X(N)+A1)*X(N)

C***  INTEGRAL OF THE INTERVAL N-1,N

      XQ2=X(N-1)*X(N-1)
      XC2=XQ2*X(N-1)
      XQ3=X(N)*X(N)
      XC3=XQ3*X(N)
      F=F+A2*(XC3-XC2)/3+A1*(XQ3-XQ2)/2+A0*(X(N)-X(N-1))

    3 CONTINUE

      IF (X(1) .GT. X(N)) F=-F

      RETURN

      END    
