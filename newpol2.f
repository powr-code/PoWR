      SUBROUTINE NEWPOL2(Y,X,YY,XX,ND)

C **** INTERPOLATION DURCH GEMITTELTE NEWTONPOLYNOME 2. GRADES.
C ***
C ***   
C ***
      IMPLICIT NONE
      INTEGER ND
      REAL Y(ND), X(ND), YY, XX

      INTEGER I, IS
      REAL Z, C21, C22, C23, C31, C32

C *** find XX
      DO I = 1, ND
          IF( XX .lt. X(I) ) GO TO 1
      end do
1     CONTINUE

      IS = I
      IF( IS .EQ. 1 ) THEN		! no extrapolation
	  YY = Y(1)
      ELSE IF( IS .GT. ND ) THEN	! still no extrapolation
	  YY = Y(ND)
      ELSE IF( IS .EQ. 2 ) THEN		! second point has problems
          C22 = (Y(2)-Y(1))/ (X(2)-X(1))
          C23 = (Y(3)-Y(2))/ (X(3)-X(2))
          C32 = (C23-C22)/ (X(3)-X(1))
          Z = XX - X(1)
          YY = Y(1) + C22*Z + C32*Z* (XX-X(2))
      ELSE IF( IS .LT. ND ) THEN	! normal interpolation
          C21 = (Y(IS-1)-Y(IS-2))/ (X(IS-1)-X(IS-2))
          C22 = (Y(IS)-Y(IS-1))/ (X(IS)-X(IS-1))
          C23 = (Y(IS+1)-Y(IS))/ (X(IS+1)-X(IS))
          C31 = (C22-C21)/ (X(IS)-X(IS-2))
          C32 = (C23-C22)/ (X(IS+1)-X(IS-1))
          Z = XX - X(IS-1)
          YY = Y(IS-1) + C22*Z + 0.5* (C31+C32)*Z* (XX-X(IS))
      ELSE ! IS .EQ. ND 		! some trouble again
          C21 = (Y(ND-1)-Y(ND-2))/ (X(ND-1)-X(ND-2))
          C22 = (Y(ND)-Y(ND-1))/ (X(ND)-X(ND-1))
          C31 = (C22-C21)/ (X(ND)-X(ND-2))
          Z = XX - X(ND-1)
          YY = Y(ND-1) + C22*Z + C31*Z* (XX-X(ND))
      END IF

      END
