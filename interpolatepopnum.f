      SUBROUTINE INTERPOLATEPOPNUM(POPNUM,  !old (and new) popnumber array
     >                             POPOLD,  !old POPNUM array
     >                             POPMIN,  !POPMIN value (replaces zeroes)
     >                             Rnew,    !new radius grid vector
     >                             Rold,    !old radius grid vector
     >                             ENTOTnew, !old total particle number
     >                             ENTOTold, !new total particle number
     >                             N,       !number of levels in DATOM
     >                             ND,      !number of depth points
     >                             ABXYZ,
     >                             NFIRST,  !Array with first level number of element blocks in level list
     >                             NLAST,   !similar to NFIRST, but with last level number
     >                             NATOM,   !Number of different elements in the model
     >                             bUseENTOT)     
C**********************************************************************
C***
C***    Interpolation of Popnumbers on new Radius-Grid
C***     crucial if radius grid has been updated in-between iterations
C***    called from: ENSURETAUMAX, HYDROSOLVE, HYDRO_REGRID
C***
C***    The interpolation is performed in  (log n_i) over (log n_tot)
C***    The older version (interpolation over radius) 
C***    can be activated by setting bUseENTOT to .FALSE.
C***
C**********************************************************************

      IMPLICIT NONE
      INCLUDE 'interfacebib.inc'

      INTEGER, INTENT(IN) :: N, ND, NATOM
      INTEGER, DIMENSION(NATOM) :: NFIRST, NLAST
      REAL, DIMENSION(NATOM) :: ABXYZ
      REAL, DIMENSION(ND), INTENT(IN) :: Rnew, Rold, ENTOTold, ENTOTnew
      REAL, DIMENSION(ND, N), INTENT(IN) :: POPOLD
      REAL, DIMENSION(ND, N), INTENT(INOUT) :: POPNUM
      REAL, INTENT(IN) :: POPMIN
      
      REAL, DIMENSION(ND) :: ENTOToldLOG, POPJLOG, RoldLOG
        
      INTEGER :: L, J, NA, NFIRNA, NLANA
      REAL :: SUMME, POPJLOGnewL, ENTOTnewLOGL, RnewLOGL

      LOGICAL, INTENT(IN) :: bUseENTOT
      
      
C***  Skip this routine if called for a not yet defined popnumber array
      IF (MAXVAL(POPOLD) <= 0.) THEN
        POPNUM = POPOLD
        RETURN
      ENDIF
      

C***  Recommended branch: interpolation of (log n_i) over (log n_tot)
      IF (bUseENTOT) THEN
        DO J=1, N
C         prepare logarithmic  vectors
          DO L=1, ND
            ENTOToldLOG(L) = LOG10(ENTOTold(L))
C           Note: minimum popnumber ist POPMIN
            POPJLOG(L) = LOG10(MAX(POPOLD(L,J),POPMIN))
          ENDDO          
C         Perform interpolation
          dploop: DO L=1, ND
            ENTOTnewLOGL = LOG10(ENTOTnew(L))
            IF (ENTOTnewLOGL > ENTOToldLOG(ND)) THEN
              !more dense than old innermost value => take old inner boundary value
              POPJLOGnewL = POPJLOG(ND)
            ELSEIF (ENTOTnewLOGL < ENTOToldLOG(1)) THEN
              !less dense than old outermost value => take old outer boundary value
              POPJLOGnewL = POPJLOG(1)
            ELSE
              CALL SPLINPOX(POPJLOGnewL, ENTOTnewLOGL,
     >                     POPJLOG, ENTOToldLOG, ND)
            ENDIF
            POPNUM(L,J) = 10**(POPJLOGnewL)
          ENDDO dploop
        ENDDO

      ELSE
C***   double-logarithmic interpolation over radius
        DO J=1, N
          DO L=1, ND
            RoldLOG(L) = LOG10(Rold(L))
C           Note: minimum popnumber ist POPMIN
            POPJLOG(L) = LOG10(MAX(POPOLD(L,J),POPMIN))
          ENDDO          
          DO L=1, ND            
            RnewLOGL = LOG10(Rnew(L))
            IF (RnewLOGL > RoldLOG(1)) THEN
C             If R outside the old grid, use old outermost value
              POPJLOGnewL = POPJLOG(1)
            ELSEIF (RnewLOGL < RoldLOG(ND)) THEN
C             If R < 1 use innermost value (should never happen)
              POPJLOGnewL = POPJLOG(ND)
            ELSE
              CALL SPLINPOX(POPJLOGnewL,RnewLOGL,POPJLOG,RoldLOG,ND)
            ENDIF
            POPNUM(L,J) = 10**(POPJLOGnewL)
          ENDDO
        ENDDO
      ENDIF

C***  Renormalization      
      DO L=1, ND
        DO NA=1, NATOM
          SUMME=0.0
          NFIRNA = NFIRST(NA)
          NLANA = NLAST(NA)
          DO J = NFIRNA, NLANA
            SUMME = SUMME + POPNUM(L,J)
          ENDDO
          SUMME = SUMME / ABXYZ(NA)
          IF (SUMME /= 0.) THEN
            DO J=NFIRNA,NLANA
              POPNUM(L,J) = POPNUM(L,J) / SUMME
            ENDDO
          ENDIF
        ENDDO

      ENDDO

      RETURN
      END
