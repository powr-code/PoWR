      SUBROUTINE MOMO (OPA, ETA, THOMSON, EDDI, R, XJC, A, B, C, W, ND,
     $                XLAM, T, DTDR, TEFF, NOTEMP, HNUE) 
C***********************************************************************
C***  SOLUTION OF THE MOMENT EQUATION **********************************
C***  CALLED FROM: (MAIN PROGRAM) COMO
C***********************************************************************

      DIMENSION OPA(ND),ETA(ND),THOMSON(ND),EDDI(3,ND),R(ND),XJC(ND)
      DIMENSION A(ND),B(ND),C(ND),W(ND),HNUE(ND)
      LOGICAL NOTEMP
 
C***  OUTER BOUNDARY
      FL=EDDI(1,1)
      QL=EDDI(2,1)
      H =EDDI(3,1)
      FP=EDDI(1,2)
      QP=EDDI(2,2)
      RL=R(1)
      X=OPA(1)
      DX=(QL+QP)*(X+OPA(2))*(RL-R(2))/4.
      B(1)=2.*QL*QL*FL*X/DX/DX+2.*QL*X*H/DX+X*(1.-THOMSON(1))
      C(1)=2.*QL*QP*FP*X/DX/DX
      W(1)=RL*RL*ETA(1)
 
C***  NON-BOUNDARY POINTS
      NDM=ND-1
      DO 1 L=2,NDM
      FL=EDDI(1,L)
      QL=EDDI(2,L)
      FP=EDDI(1,L+1)
      QP=EDDI(2,L+1)
      FM=EDDI(1,L-1)
      QM=EDDI(2,L-1)
      X=OPA(L)
      RL=R(L)
      DR=(R(L-1)-R(L+1))/2.
      DA=DR*(R(L-1)-RL)*(QM+QL)*(OPA(L-1)+X)/4.
      DC=DR*(RL-R(L+1))*(QP+QL)*(OPA(L+1)+X)/4.
      A(L)=QM*FM/DA
      C(L)=QP*FP/DC
      B(L)=QL*FL*(1./DA+1./DC)+X*(1.-THOMSON(L))
    1 W(L)=RL*RL*ETA(L)
 
C***  INNER BOUNDARY
      FL=EDDI(1,ND)
      QL=EDDI(2,ND)
      H =EDDI(3,ND)
C***  NEW CALCULATION OF CURRENT "HPLUS" (IN CASE OF TEMP. CORRECTIONS):
C***  -- NOTE: "DTDR" IS PROVIDES BY SUBR. DIFDTDR (CALLED FROM COMO)
      IF (.NOT. NOTEMP) THEN
         CALL DIFFUS (XLAM,T,R,ND,BCORE,DBDR,DTDR,TEFF,NOTEMP)
         EDDI(3,ND-1)=BCORE/2. + DBDR/3./OPA(ND)
         ENDIF
      HPLUS=EDDI(3,ND-1)
      FM=EDDI(1,ND-1)
      QM=EDDI(2,ND-1)
      RL=R(ND)
      X=OPA(ND)
      DX=(QM+QL)*(X+OPA(ND-1))*(R(ND-1)-R(ND))/4.
      A(ND)=2.*QL*QM*FM*X/DX/DX
      B(ND)=2.*QL*QL*FL*X/DX/DX+2.*QL*X*H/DX+X*(1.-THOMSON(ND))
      W(ND)=RL*RL*(ETA(ND)+2.*QL*X*HPLUS/DX)
 
      CALL INVTRI (A,B,C,W,ND)
      DO 2 L=1,ND
      RL=R(L)
    2 XJC(L)=W(L)/RL/RL

C***  CALCULATION OF HNUE (USED IN STEAL FOR TEMPERATURE EQUATION)
C***  Note that HNUE(L) refers to the interstice L,L-1
      DO 3 L=2,ND
      FL=EDDI(1,L)
      QL=EDDI(2,L)
      FM=EDDI(1,L-1)
      QM=EDDI(2,L-1)
      X=0.5*(OPA(L-1)+OPA(L))
      Q=0.5*(QM+QL)
      DR=R(L-1)-R(L)
      HNUE(L)=(FL*QL*W(L)-FM*QM*W(L-1))/(X*Q*DR)
    3 CONTINUE

      RETURN
      END
