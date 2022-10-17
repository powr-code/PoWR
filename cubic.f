      SUBROUTINE CUBIC (L,LTOT,ZRAY,XCMF,P1,P2,P3,P4)
C***********************************************************************
C***  CUBIC INTERPOLATION WITH FIXED DERIVATIVES AT THE MESH POINTS
C***  THESE DERIVATIVES ARE CALCULATED BY LINEAR I5TERPOLATION BETWEEN
C***  NEIGHBOURING POINTS
C***  XCMF(L) = GRID OF MESH POINTS
C***  ZRAY(L) = FUNCTION TO BE INTERPOLATED
C***  P1,P2,P3,P4 = RESULTING COEFFICIENTS
C***  EVALUATION OF THE INTERPOLATED VALUE Z(X) : 
C***  Z(X)=P1*(X-XCMF(L-1))**3 + P2*(X-XCMF(L-1))
C***                + P3*(XCMF(L)-X)**3 + P4*(XCMF(L)-X)
C***  THE COEFFICIENT MATRIX H HAS BEEN INVERTED ANALYTICALLY
C***********************************************************************
 
      DIMENSION ZRAY(LTOT),XCMF(LTOT)
 
      IF (L .LT. 2 .OR. L .GT. LTOT) THEN
        WRITE (0,'(A,2(I4,1X))') 'L, LTOT=',L, LTOT
        STOP 'ERROR in subroutine CUBIC'
      ENDIF
 
C***  SET UP THE COEFFICIENT MATRIX
      D1=1./(XCMF(L)-XCMF(L-1))
      D2=D1*D1
      D3=D1*D2
      D23=D2/3.
      H11=D3
      H12=-D3
      H13=D23
      H14=2.*D23
      H21=-D1
      H22=2.*D1
      H23=-0.333333333333333
      H24=-0.666666666666666
      H31=-D3
      H32=D3
      H33=-2.*D23
      H34=-D23
      H41=2.*D1
      H42=-D1
      H43=0.666666666666666
      H44=0.333333333333333
C***  FOR THE BOUNDARY INTERVALS THE DERIVATIVE CANNOT EXTEND OVER THE BOUNDARY
      LA=MAX0(L-2,1)
      LB=MIN0(L+1,LTOT)
C***  FUNCTION TO BE INTERPOLATED: ZRAY
      F1=ZRAY(L-1)
      F2=ZRAY(L)
      F3=(ZRAY(L)-ZRAY(LA))/(XCMF(L)-XCMF(LA))
      F4=(ZRAY(LB)-ZRAY(L-1))/(XCMF(LB)-XCMF(L-1))
 
C***  CALCULATE POLYNOMIAL COEFFICIENTS: P(VECTOR) = H(MATRIX) * F(VECTOR)
      P1=H11*F1+H12*F2+H13*F3+H14*F4
      P2=H21*F1+H22*F2+H23*F3+H24*F4
      P3=H31*F1+H32*F2+H33*F3+H34*F4
      P4=H41*F1+H42*F2+H43*F3+H44*F4
 
      RETURN
      END
