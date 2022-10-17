      SUBROUTINE REGRIDOLDPOP(POPNUM, Told, SCRATCH, ND, N, T,
     >                        POPMIN, TAUTHICK, bNoRGrid,
     >                        RADIUS, RADIUSorg, ENTOT, ENTOTorg,
     >                        RNE, RNEorg, ABXYZ, NFIRST, NLAST, 
     >                        NATOM, WEIGHT, NCHARG, EION, ELEVEL, NOM, 
     >                        TAUROSScont, TAUSCALorg,
     >                        ENLTE, DEPARTNDorg, bUseENTOT)
C***********************************************************************      
C***  Adjusts old population numbers (i.e. POP1, etc. from older 
C***  STEAL jobs) to a new stratification (R, V)
C***
C***  Interpolation is performed over radius (if grid is changed) 
C***  in the outer part and always over TAUROSScont
C***  in the optical thick part (TAU > TAUTHICK)
C***
C***  Inner extrapolation is done using LTE population numbers
C***  multiplied with the old innermost departure coefficients
C***
C***  called from ENSURETAUMAX
C***********************************************************************      

      IMPLICIT NONE
      INCLUDE 'interfacebib.inc'
      
      INTEGER, INTENT(IN) :: ND, N, NATOM
      REAL, INTENT(IN) :: POPMIN, TAUTHICK
      LOGICAL, INTENT(IN) :: bNoRGrid, bUseENTOT
      
      REAL, DIMENSION(ND), INTENT(IN) :: RNEorg, RNE, ENTOT, ENTOTorg, 
     >                                   RADIUSorg, RADIUS, T, Told,
     >                                   TAUROSScont, TAUSCALorg
      REAL, DIMENSION(N), INTENT(IN) :: WEIGHT, ELEVEL, EION
      REAL, DIMENSION(NATOM) :: ABXYZ

      INTEGER, DIMENSION(N), INTENT(IN) :: NCHARG, NOM
      INTEGER, DIMENSION(NATOM), INTENT(IN) :: NFIRST, NLAST
                
      REAL, DIMENSION(ND, N), INTENT(INOUT) :: POPNUM, SCRATCH      
      REAL, DIMENSION(N), INTENT(INOUT) :: ENLTE, DEPARTNDorg
                  
      REAL :: ENEL
      INTEGER :: L, J
      
      CALL BACKUPPOPNUM(ND, N, POPNUM, SCRATCH)      
      IF (.NOT. bNoRGrid) THEN
        CALL INTERPOLATEPOPNUM(POPNUM, SCRATCH, POPMIN,
     >                         RADIUS, RADIUSorg, ENTOT, ENTOTorg,
     >                         N, ND, ABXYZ, NFIRST, NLAST, NATOM,
     >                         bUseENTOT)
      ENDIF     
C***  Calculate departure coefficients for the innermost
C***  depth point of the old stratification
      ENEL = RNEorg(ND) * ENTOTorg(ND)
      CALL LTEPOP (N, ENLTE, Told(ND), ENEL, 
     >             WEIGHT, NCHARG, EION, ELEVEL, NOM,
     >             ABXYZ, NFIRST, NLAST, NATOM) 
      DO J=1, N
        DEPARTNDorg(J) = SCRATCH(ND,J)/ENLTE(J)
      ENDDO      
      DO L=1, ND        
        IF (TAUROSScont(L) < TAUTHICK) CYCLE
        IF (TAUROSScont(L) <= TAUSCALorg(ND)) THEN
          DO J=1, ND
            CALL SPLINPOX(POPNUM(L,J), TAUROSScont(L), 
     >                                SCRATCH(1,J), TAUSCALorg, ND)
          ENDDO
        ELSE
C***      Inner extrapolation necessary
          ENEL = RNE(L) * ENTOT(L)
          CALL LTEPOP (N, ENLTE, T(L), ENEL, 
     >                 WEIGHT, NCHARG, EION, ELEVEL, NOM,
     >                 ABXYZ, NFIRST, NLAST, NATOM) 
          DO J=1, N
            POPNUM(L,J) = ENLTE(J) * DEPARTNDorg(J)
          ENDDO
        ENDIF        
      ENDDO      
      
      RETURN
      END
