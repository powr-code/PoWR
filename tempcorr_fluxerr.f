      SUBROUTINE TEMPCORR_FLUXERR(FLUXERR, ND, HTOTL, HTOTCMF0, 
     >                            FTEPS, TAUROSS, AUTOFLUXTAU)
C***********************************************************************
C***  Calculation of the current deviation from flux conservation
C***  and determination of an optical depth that can be used
C***  as a damping parameter for flux correction terms that should
C***  favor inner flux corrections.
C***
C***  called from TEMPCORR
C***********************************************************************
      IMPLICIT NONE
      
      INTEGER, INTENT(IN) :: ND
      REAL, INTENT(IN) :: FTEPS
      REAL, INTENT(OUT) :: AUTOFLUXTAU
     
      REAL, DIMENSION(ND-1) :: HTOTL, HTOTCMF0
      REAL, DIMENSION(ND) :: FLUXERR, TAUROSS
           
      
      INTEGER :: L
      
C***  Calculate current flux error      
      FLUXERR(1) = 0.      
      DO L=2, ND
          FLUXERR(L) = 2.*ABS(HTOTL(L-1) - HTOTCMF0(L-1))
     >                          /(HTOTL(L-1)+HTOTCMF0(L-1))     
      ENDDO
      
C***  Estimate automatic tau via innermost point where the
C***  flux deviates more than FTEPS
C***  => This point needs full correction, so put tau 
C***     at 0.1 times the current optical depth
      AUTOFLUXTAU = TAUROSS(ND) * 0.1
      DO L=ND, 2, -1
        IF (FLUXERR(L) > FTEPS) EXIT
        AUTOFLUXTAU = TAUROSS(L) * 0.1
      ENDDO
      
      RETURN
      
      END      
