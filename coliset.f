      SUBROUTINE COLISET (Z,
     >             ND,LMAX,TA,TB,TC,UB,VA,VB,GA,H,S,OPAK,ETAK,
     >             PP,BCORE,DBDR,RADIUS,XIMINUS,XIPLUS,JP)
C***********************************************************************
C***  THIS SUBROUTINE IS TO SET UP THE ARRAY ELEMENTS FOR THE CMF FORMALISM
C***********************************************************************

      DIMENSION S(ND),OPA(ND),ETA(ND),VA(ND),VB(ND)
      DIMENSION TA(ND),TB(ND),TC(ND),UB(ND),GA(ND),H(ND)
      DIMENSION PP(ND),Z(ND),OPAK(ND),ETAK(ND)

      LZ=LMAX-1
      NC=JP-ND
 
C***  OUTER BOUNDARY CONDITION  -  SECOND ORDER

      EK=ETAK(1)
      AK=OPAK(1)
      AZ=0.5*(AK+OPAK(2))
      TAU=Z(1)-Z(2)
      TAUZ=TAU
      DT=1./(AK*TAU)
      DTZ=1./(AZ*TAUZ)
      DX=PP(1)/AK
      DXZ=(PP(1)+PP(2))*0.5/AZ
      DAZ=DTZ/(1.+DXZ)
      DBZ=DXZ/(1.+DXZ)

      S(1) = EK/AK + 2.*DT*XIMINUS
      TC(1)= 2.*DT*DAZ
      TB(1)= TC(1) + 2.*DT + DX + 1.
      UB(1)= DX
      VB(1)= 2.*DT*DBZ

C***  FOR G AND H, THE MATRIX ELEMENT S ARE NOT DIFFERENT FROM INNER POINTS
C      DTZM=1./(AZ*TAUZ)
C      DXZM=(PP(1)+PP(2))/AZ/2.
C      DAZM=DTZM/(1.+DXZM)
C      DBZM=DXZM/(1.+DXZM)
      GA(1)=-DAZ
      H(1) = DBZ 
      
      IF(LZ.LT.2) GOTO 2
 
C***  NON-BOUNDARY POINTS
      DO 1 L=2,LZ
      AK=OPAK(L)                              ! PHIK*OPAL(L)
      EK=ETAK(L)                              ! PHIK*ETAL(L)
      S(L)=EK/AK
      AZ=0.5*(AK+OPAK(L+1))                 ! PHIK*OPAL(L+1))
      AZM=0.5*(AK+OPAK(L-1))                ! PHIK*OPAL(L-1))
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
C      IF(LMAX.LT.ND) GOTO 4
      IF(JP.GT.NC) GOTO 4
 
C***  INNER BOUNDARY CONDITION (CORE RAYS)  -  SECOND ORDER
      EK=ETAK(ND)
      AK=OPAK(ND)
      AZM=0.5*(OPAK(LZ)+AK)
      TAU=Z(LZ)-Z(ND)
      TAUZM=TAU
      DT=1./(AK*TAU)
      DTZM=1./(AZM*TAUZM)
      DX=PP(ND)/AK
      DXZM=(PP(LZ)+PP(ND))*0.5/AZM
      DAZM=DTZM/(1.+DXZM)
      DBZM=DXZM/(1.+DXZM)

C***  INNER BOUNDARY CONDITION (CORE RAYS)
C***  SECOND ORDER
      S(ND) = EK/AK + 2.*DT*XIPLUS
      TA(ND)= 2.*DT*DAZM
      TB(ND)= TA(ND) + 2.*DT + DX + 1.
      UB(ND)= DX
      VA(ND)= -2.*DT*DBZM
C***  FIRST ORDER
C      S(ND) = BCORE+DBDR*Z(ND)/AK
C      TA(ND)= DT
C      TB(ND)= DT+DX+1.
C      UB(ND)= DX
C      VA(ND)= 0.

      RETURN

C***  INNER BOUNDARY CONDITION (NON-CORE RAYS)  -  SECOND ORDER
    4 AK=OPAK(LMAX)
      EK=ETAK(LMAX)
      AZM=0.5*(OPAK(LZ)+AK)
      TAU=Z(LZ)
      TAUZM=TAU
      DT=1./(AK*TAU)
      DTZM=1./(AZM*TAUZM)
      DX=PP(LMAX)/AK
      DXZM=(PP(LZ)+PP(LMAX))*0.5/AZM
      DAZM=DTZM/(1.+DXZM)
      DBZM=DXZM/(1.+DXZM)

      S(LMAX) = EK/AK
      TA(LMAX)= 2.*DT*DAZM
      TB(LMAX)= TA(LMAX) + DX + 1.
      UB(LMAX)= DX
      VA(LMAX)= -2.*DT*DBZM

      RETURN
      END




