      SUBROUTINE TEMPCORR_EXPDAMP(TAUDAMP, METHOD, TAUROSS, bSKIPCORR, 
     >                            L, ND, DE_LOC, DE_INT, DE_RMAX, DE_TB)
C***********************************************************************
C***  Optional exponential tau-based damping for temperature corrections
C***  suggested by Hubeny & Mihalas (book, 2015)
C***
C***  called from TEMPCORR
C***********************************************************************
      IMPLICIT NONE
      
      INTEGER, INTENT(IN) :: L, ND
      REAL, INTENT(IN) :: TAUDAMP
      LOGICAL, INTENT(IN) :: bSKIPCORR
      CHARACTER(4), INTENT(IN) :: METHOD

      REAL, INTENT(INOUT) :: DE_LOC, DE_INT, DE_RMAX, DE_TB
      
      REAL, DIMENSION(ND), INTENT(IN) :: TAUROSS
            
      REAL :: Q, DTAUL, DTAUD     
            
      REAL, PARAMETER :: PI = 3.141592654          
            
C***  File and channel handles (=KANAL)
      INTEGER, PARAMETER :: hOUT = 6        !write to wruniqX.out (stdout)
      INTEGER, PARAMETER :: hCPR = 0        !write to wruniqX.cpr (stderr)
      
C***  Default: No additional damping      
      DE_LOC = 1.
      DE_INT = 1.
      DE_RMAX = 1.
      DE_TB = 1.
      
C***  TAUDAMP can be given in the CARDS file as EXPTAU or COSTAU
C***  For a non-zero value, additional depth-dependent damping factors
C***  are calculated and applied
C***  FLUX CONSISTENCY terms are damped outwards, all other terms are
C***  damped inwards
      IF (TAUDAMP /= 0.) THEN
        IF (METHOD == 'EXP') THEN
          IF (TAUROSS(L) <= 0.) THEN
            DE_LOC = 1.
            DE_INT = 0.
            DE_RMAX = 0.
            DE_TB = 1.
          ELSE
            DE_LOC = EXP(-TAUROSS(L)/TAUDAMP)
            DE_INT = 1. - EXP(-TAUROSS(L)/TAUDAMP)
            DE_RMAX = 1. - EXP(-TAUROSS(L)/TAUDAMP)
            DE_TB = EXP(-TAUROSS(L)/TAUDAMP)
          ENDIF
        ELSEIF (METHOD == 'COS') THEN
C***      COS damping between TAUDAMP and 10*TAUDAMP    
          IF (TAUROSS(L) < TAUDAMP) THEN
            DE_LOC = 1.
            DE_INT = 0.
            DE_RMAX = 0.
            DE_TB = 1.
          ELSEIF (TAUROSS(L) > 10. * TAUDAMP) THEN
            DE_LOC = 0.
            DE_INT = 1.
            DE_RMAX = 1.
            DE_TB = 0.
          ELSE
            DTAUL = TAUROSS(L) - TAUDAMP
            DTAUD = 9. * TAUDAMP 
            Q = 0.5 + 0.5 * COS(PI*DTAUL/DTAUD)
            DE_LOC = Q          
            DE_INT = 1.-Q          
            DE_RMAX = 1.-Q          
            DE_TB = Q          
          ENDIF
        ELSE
          WRITE (hCPR,*) '*** ERROR: UNKNOWN EXPDAMP METHOD'
          STOP 'FATAL ERROR in STEAL->TEMPCORR_EXPDAMP'
        ENDIF  
      ENDIF
                  

C***    Option: In case of negative flux or arad, skip TCORR completely
        IF (bSKIPCORR) THEN
          WRITE (hCPR,*) 'Temperature corrections skipped to due to ' 
     >       // 'negative Flux or acceleration'
          DE_LOC = 0.
          DE_INT = 0.
          DE_RMAX = 0.
          DE_TB = 0.
        ENDIF

      RETURN
      END
        
