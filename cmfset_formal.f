      SUBROUTINE CMFSET_FORMAL(Z,
     $          ND,LMAX,TA,TB,TC,UB,VA,VB,GA,H,S,OPA,ETA,
     $          PP,BCORE,DBDR,RADIUS,XIMINUS,DXI,OPAK,ETAK,
     >          OPAFE, ETAFE, BWITHLINES)
C***********************************************************************
C***  THIS SUBROUTINE IS TO SET UP THE ARRAY ELEMENTS FOR THE CMF FORMALISM
C***********************************************************************

      DIMENSION S(ND),OPA(ND),ETA(ND),VA(ND),VB(ND)
      DIMENSION TA(ND),TB(ND),TC(ND),UB(ND),GA(ND),H(ND)
      DIMENSION PP(ND),Z(ND),OPAK(ND),ETAK(ND)
      DIMENSION OPAFE(ND), ETAFE(ND)
      LOGICAL BWITHLINES
 
      LZ=LMAX-1
  
C***  OUTER BOUNDARY CONDITION  -  FIRST ORDER
      IF (BWITHLINES) THEN
         AK=OPA(1)+OPAK(1)+OPAFE(1) ! PHIK*OPAL(1)
         AZ=0.5*(AK+OPA(2)+OPAK(2)+OPAFE(2)) ! PHIK*OPAL(2))
      ELSE
         AK=OPA(1)
         AZ=0.5*(AK+OPA(2))
      ENDIF  
      TAUZ=Z(1)-Z(2)
      DX=PP(1)/AK
      S(1)=XIMINUS+DX*DXI
      TC(1)=1./(AK*TAUZ)
      TB(1)=TC(1)+DX+1.
      UB(1)=DX
      VB(1)=0.0
C***  FOR G AND H, THE MATRIX ELEMENT S ARE NOT DIFFERENT FROM INNER POINTS
      DTZM=1./(AZ*TAUZ)
      DXZM=(PP(1)+PP(2))/AZ/2.
      DAZM=DTZM/(1.+DXZM)
      DBZM=DXZM/(1.+DXZM)
      GA(1)=-DAZM
      H(1)=DBZM
      IF(LZ.LT.2) GOTO 2
 
C***  NON-BOUNDARY POINTS
      DO 1 L=2,LZ
         IF (BWITHLINES) THEN
            AK=OPA(L)+OPAK(L)+OPAFE(L) ! PHIK*OPAL(L)
            EK=ETA(L)+ETAK(L)+ETAFE(L) ! PHIK*ETAL(L)
            AZ=0.5*(AK+OPA(L+1)+OPAK(L+1)+OPAFE(L+1)) ! PHIK*OPAL(L+1))
            AZM=0.5*(AK+OPA(L-1)+OPAK(L-1)+OPAFE(L-1)) ! PHIK*OPAL(L-1))
         ELSE
            AK=OPA(L)
            EK=ETA(L) 
            AZ=0.5*(AK+OPA(L+1)) 
            AZM=0.5*(AK+OPA(L-1))
         ENDIF
         S(L)=EK/AK
         TAU=0.5*(Z(L-1)-Z(L+1))
         TAUZ=Z(L)-Z(L+1)
         TAUZM=Z(L-1)-Z(L)
         DT=1./(AK*TAU)
         DTZ=1./(AZ*TAUZ)
         DTZM=1./(AZM*TAUZM)
         DX=PP(L)/AK
         DXZ=(PP(L)+PP(L+1))*0.5/AZ
         DXZM=(PP(L)+PP(L-1))*0.5/AZM
         DAZ=DTZ/(1.+DXZ)
         DAZM=DTZM/(1.+DXZM)
         DBZ=DXZ/(1.+DXZ)
         DBZM=DXZM/(1.+DXZM)
         TA(L)=DT*DAZM
         TC(L)=DT*DAZ
         TB(L)=TA(L)+TC(L)+DX+1.
         UB(L)=DX
         VA(L)=-DT*DBZM
         VB(L)=DT*DBZ
         GA(L)=-DAZ
         H(L)=DBZ
C     DTZM=DTZ
C     DAZM=DAZ
C     DBZM=DBZ
    1 CONTINUE
 
    2 L=LMAX
      IF(LMAX.LT.ND) GOTO 4
      
C***  INNER BOUNDARY CONDITION (CORE RAYS)  -  ONLY TO FIRST ORDER
      IF (BWITHLINES) THEN
         AK=OPA(ND)+OPAK(ND)+OPAFE(ND) ! PHIK*OPAL(ND)
      ELSE
         AK=OPA(ND)
      ENDIF  
      S(ND)=BCORE+DBDR*Z(ND)/AK
      TAUZ=Z(ND-1)-Z(ND)
      DT=1./(TAUZ*AK)
      DX=PP(L)/AK
      TA(L)=DT
      TB(L)=DT+DX+1.
      UB(L)=DX
      VA(L)=0.0
      RETURN
 
C***  INNER BOUNDARY CONDITION (NON-CORE RAYS)  -  SECOND ORDER
    4 CONTINUE
      IF (BWITHLINES) THEN      
         AK=OPA(L)+OPAK(L)+OPAFE(L) ! PHIK*OPAL(LMAX)
         EK=ETA(L)+ETAK(L)+ETAFE(L) ! PHIK*ETAL(L)
      ELSE
         AK=OPA(L)
         EK=ETA(L)
      ENDIF        
      S(L)=EK/AK
      TAUZ=Z(LZ)
      DT=1./(AK*TAUZ)
      DX=PP(L)/AK
      DA=DT/(1.+DX)
      DB=DX/(1.+DX)
      TA(L)=2*DT*DA
      TB(L)=TA(L)+DX+1.
      UB(L)=DX
      VA(L)=-2*DT*DB
      RETURN
      END



