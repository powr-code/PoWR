      SUBROUTINE CONVOLOPAFE(YORIG, YSCRATCH, NDDIM, ND, NFL, DX, 
     >                       VDOP, VDOPFE, DD_VDOPDUFE)
C**********************************************************************
C***  CONVOLUTION OF PROFILE WITH GAUSS-Function exp(-x*x/(sigma^2))
C***
C***  called from FORMCMF
C**********************************************************************
      IMPLICIT NONE
      INCLUDE 'interfacebib.inc'

      INTEGER, INTENT(IN) :: NDDIM, NFL, ND

      REAL, DIMENSION(ND), INTENT(IN) :: DD_VDOPDUFE
      REAL, DIMENSION(NDDIM, NFL), INTENT(INOUT) :: YORIG
      REAL, DIMENSION(NFL), INTENT(INOUT) :: YSCRATCH

      INTEGER, PARAMETER :: KRELMAX = 250000
      INTEGER, PARAMETER :: NDCONVMAX =  200

      REAL, DIMENSION(KRELMAX) :: DOPWEIGHT
      INTEGER, DIMENSION(NDCONVMAX) :: KRELINDEX
      
      REAL, INTENT(IN) :: VDOP, VDOPFE, DX
      
      REAL :: SIGMA_SQRD, WEIGHT, WEIGHTSUM, X
      INTEGER :: L, KL, K, KREL, KR
      
C***  File and channel handles (=KANAL)
      INTEGER, PARAMETER :: hOUT = 6        !write to wruniqX.out (stdout)
      INTEGER, PARAMETER :: hCPR = 0        !write to wruniqX.cpr (stderr)
      
      IF (ND > NDCONVMAX) THEN
        WRITE (hCPR,'(A)') 'CONVOLOPAFE: FATAL ERROR ******'
        WRITE (hCPR,'(A)') 'CONVOLOPAFE: NDHDMAX INSUFFICIENT'
        WRITE (hCPR,'(2(A,I4))') 'ND = ', ND, ', NDCONVMAX = ', NDCONVMAX
        STOP 'FATAL ERROR IN FORMAL->CONVOLOPAFE'
      ENDIF
            
C***  LOOP OVER ALL DEPTH POINTS      
      DO L=1, ND
C***    Do not perform any convolution if current intended FE VDOP
C***      is below FEDAT VDOPFE
        IF (DD_VDOPDUFE(L) <= VDOPFE/VDOP) CYCLE
        
        SIGMA_SQRD = (DD_VDOPDUFE(L))**2 - (VDOPFE/VDOP)**2
        
C***    Prepare doppler weight vector        
        X = DX
        KREL = 0
        DO WHILE (X < 4.5*DD_VDOPDUFE(L))
          KREL = KREL + 1
          IF (KREL > KRELMAX) THEN
            WRITE (hCPR,'(A)') 'CONVOLOPAFE: FATAL ERROR ******'
            WRITE (hCPR,'(A)') 'CONVOLOPAFE: KRELMAX INSUFFICIENT'
            STOP 'FATAL ERROR IN FORMAL->CONVOLOPAFE'          
          ENDIF
          
          DOPWEIGHT(KREL) = EXP(-X*X/ SIGMA_SQRD)
          X = X + DX
        ENDDO
C***    Store maximum relative index
        KRELINDEX(L) = KREL
                  
C****   LOOP OVER ALL DATA POINTS **************************************
        DO KL=1, NFL

C***    CENTRAL POINT
          WEIGHTSUM = 1.
          YSCRATCH(KL) = YORIG(L,KL) 

C***      Loop over relative steps to the left          
          DO KR = 1, KRELINDEX(L)
            K = KL - KR
            IF (K < 1) EXIT
            YSCRATCH(KL) = YSCRATCH(KL) + YORIG(L,K) * DOPWEIGHT(KR)
            WEIGHTSUM = WEIGHTSUM + DOPWEIGHT(KR)
          ENDDO

C***      Loop over relative steps to the right          
          DO KR = 1, KRELINDEX(L)          
            K = KL + KR
            IF (K > NFL) EXIT
            YSCRATCH(KL) = YSCRATCH(KL) + YORIG(L,K) * DOPWEIGHT(KR)
            WEIGHTSUM = WEIGHTSUM + DOPWEIGHT(KR)
          ENDDO
          
C***      Normalization of Doppler profile
C***      (This cannot be done with DOPWEIGHT vector as they would
C***       not cover the boundary situations)
          YSCRATCH(KL) = YSCRATCH(KL) / WEIGHTSUM

        ENDDO
C*****  END OF DATAPOINT LOOP ******************************************

C***  The original data is overwritten with the convolved data
        DO KL=1, NFL
          YORIG(L,KL) = YSCRATCH(KL)
        ENDDO

      ENDDO

      RETURN
      END
