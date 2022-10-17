      SUBROUTINE CONVOLGAUSS (YORIG, YCONV, NDAT, DX)
C**********************************************************************
C***  CONVOLUTION OF PROFILE WITH GAUSS-Function exp(-x*x)
C**********************************************************************
      DIMENSION YORIG(NDAT), YCONV(NDAT)

C**** LOOP OVER ALL DATA POINTS *******************************
      DO L=1, NDAT

C***    CENTRAL POINT
        WEIGHTSUM = 1.
        YCONV(L) = YORIG(L) 

C***    LOOP OVER STEPS TO THE LEFT
        K = L
510     K = K - 1
        X = (L-K) * DX
        IF (K .LT. 1 .OR. X .GT. 4.5) GOTO 610
        WEIGHT = EXP(-X*X)
        WEIGHTSUM = WEIGHTSUM + WEIGHT
        YCONV(L) = YCONV(L) + YORIG(K) * WEIGHT
        GOTO 510

610     K = L

C***    LOOP OVER STEPS TO THE RIGHT
620     K = K + 1
        X = (K-L) * DX
        IF (K .GT. NDAT .OR. X .GT. 4.5) GOTO 720
        WEIGHT = EXP(-X*X)
        WEIGHTSUM = WEIGHTSUM + WEIGHT
        YCONV(L) = YCONV(L) + YORIG(K) * WEIGHT
        GOTO 620

C*** NACHNORMIERUNG
720     YCONV(L) = YCONV(L) / WEIGHTSUM

      ENDDO
C*****END OF LOOP ******************************************************

C***  The original date are finally overwritten by the convolved data
      DO L=1, NDAT
         YORIG(L) = YCONV(L)
      ENDDO



      RETURN
      END
