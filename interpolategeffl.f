      SUBROUTINE INTERPOLATEGEFFL(GEFFLnew, GEFFL, ENTOTnew, ENTOTold,
     >                            RADIUS, RADIUSorg, ND, bUseENTOT)
C***********************************************************************
C***  Interpolation of GEFFL after grid changes
C***
C***  called by ENSURETAUMAX
C***********************************************************************
      IMPLICIT NONE
      INCLUDE 'interfacebib.inc'

      INTEGER, INTENT(IN) :: ND

      REAL, DIMENSION(ND), INTENT(IN) :: ENTOTnew, ENTOTold, 
     >                                   RADIUS, RADIUSorg, GEFFL
      REAL, DIMENSION(ND), INTENT(INOUT) :: GEFFLnew     
      REAL, DIMENSION(ND) :: ENTOTnewLOG, ENTOToldLOG

      LOGICAL, INTENT(IN) :: bUseENTOT
      
      INTEGER :: L

      IF (bUseENTOT) THEN
C***    Density interpolation        

        DO L=1, ND
          ENTOToldLOG(L) = LOG10(ENTOTold(L))
          ENTOTnewLOG(L) = LOG10(ENTOTnew(L))
        ENDDO          
        
        DO L=1, ND
          IF (ENTOTnewLOG(L) > ENTOToldLOG(ND)) THEN
            GEFFLnew(L) = GEFFL(ND)
          ELSEIF (ENTOTnewLOG(L) < ENTOToldLOG(1)) THEN            
            GEFFLnew(L) = GEFFL(1)
          ELSE
            CALL SPLINPOX(GEFFLnew(L), ENTOTnewLOG(L),
     >                    GEFFL, ENTOToldLOG, ND)
          ENDIF
        ENDDO
        
      ELSE
C***    Radius interpolation        
      
        DO L=1, ND
          IF (RADIUS(L) > RADIUSorg(1)) THEN
            GEFFLnew(L) = GEFFL(1)
          ELSE                 
            CALL SPLINPOX(GEFFLnew(L),RADIUS(L),GEFFL,RADIUSorg,ND)
          ENDIF
        ENDDO
        
      ENDIF
       
      RETURN

      END
       