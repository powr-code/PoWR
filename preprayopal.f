      SUBROUTINE PREPRAYOPAL(NP, ND, JP, NLOBLFI, NLOBLLA, 
     >                       NBLINE, NDDIM, NDADDIM, CORE, 
     >                       OPAL, ETAL, OPALRAY, ETALRAY)
C***********************************************************************
C***  Copies the necessary opacity and emissivity from the standard 
C***  arrays (1..ND) into the ray arrays (1...2*LMAX-1) which also 
C***  cover the negative z part and will be used in ZONEINT
C***
C***  called by OBSFRAM
C***********************************************************************
      IMPLICIT NONE
      
      INTEGER, INTENT(IN) :: NP, ND, JP, NLOBLFI, NLOBLLA, 
     >                       NDDIM, NDADDIM, NBLINE
      LOGICAL, INTENT(IN) :: CORE
      
      REAL, DIMENSION(NDADDIM,NBLINE) :: OPALRAY, ETALRAY
      REAL, DIMENSION(NDDIM,NBLINE) :: OPAL, ETAL
      
      INTEGER :: LMAX, L, LL, NBL
      
      LMAX=MIN0(NP+1-JP,ND)

c      DO NBL=1,NBLINE   !original version in PREPRAY
      DO NBL=NLOBLFI, NLOBLLA
        DO L=1, LMAX
          OPALRAY(L,NBL)=OPAL(L,NBL)
          ETALRAY(L,NBL)=ETAL(L,NBL)
          IF (.NOT. CORE) THEN
C***        For non-core rays we also need to fill the negative z part
            LL=2*LMAX-L
            OPALRAY(LL,NBL)=OPAL(L,NBL)
            ETALRAY(LL,NBL)=ETAL(L,NBL)
          ENDIF
        ENDDO
      ENDDO            
            
      RETURN
      
      END
      