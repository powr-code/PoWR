      SUBROUTINE GRADIFF  (ND,VELO,GRADI,RADIUS)
C***********************************************************************
C***  FOR THE VELOCITY FIELD GIVEN BY VECTOR VELO(L), THE GRADIENTS GRADI(L)
C***  ARE COMPUTED BY LINEAR INTERPOLATION BETWEEN THE NEIGHBORING POINTS
C***********************************************************************
      DIMENSION VELO(ND),GRADI(ND),RADIUS(ND)
      NDM=ND-1
      DO 1 L=2,NDM
      GRADI(L)=(VELO(L+1)-VELO(L-1))/(RADIUS(L+1)-RADIUS(L-1))
      IF (L.EQ.2) GRADI(1)=GRADI(2)
      IF (L.EQ.NDM) GRADI(ND)=GRADI(NDM)
    1 CONTINUE

      RETURN
      END
