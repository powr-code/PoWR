      SUBROUTINE PREP_PPP(ND, NP, NDDIM, NPDIM, 
     >              PPP, Z, RADIUS, GRADI, VELO,
     >              K, BSTATIC)
C****************************************************************
C***  The Array PPP is either set to Zero (Static branch) 
C***    or calculated. 
C***    Called by COLI
C****************************************************************

      DIMENSION PPP(NDDIM, NPDIM), Z(NDDIM, NPDIM)
      DIMENSION RADIUS(ND), GRADI(ND), VELO(ND)

      LOGICAL BSTATIC

      
      IF (BSTATIC .OR. K .EQ. 0) THEN
        DO JP = 1, NP-1
          LMAX=MIN0(NP+1-JP,ND)
          DO L=1, LMAX
            IRAY= L + ND * (JP-1)
            PPP(IRAY,1) =  0.
          ENDDO
        ENDDO
      ELSE
        DO JP = 1, NP-1
          LMAX=MIN0(NP+1-JP,ND)
          DO L=1, LMAX
            RL=RADIUS(L)
            IRAY= L + ND * (JP-1)
            Y=Z(IRAY,1)/RL
            YY=Y*Y
            PPP(IRAY,1) =  YY*GRADI(L) + (1.-YY)*VELO(L)/RL
          ENDDO
        ENDDO
      ENDIF

      RETURN
      END
