      SUBROUTINE SETUP_FF(FF_INFO, XLAM0, ALN,
     >                    XLAM_FINE_START, XLAM_FINE_END, IFF_N, KSPACE, 
     >                    IFF_DK, IFF_MAX)
C*******************************************************************
C***  Setup of the Array FF_INFO
C*******************************************************************

      IMPLICIT NONE

      INTEGER, PARAMETER :: TINYINT = SELECTED_INT_KIND(2)

      REAL, DIMENSION(10) :: FF_INFO
      INTEGER (KIND=TINYINT), DIMENSION(IFF_MAX) :: IFF_DK

      REAL :: XLAM0, ALN, XLAM_FINE_START, XLAM_FINE_END
      INTEGER :: IFF_MAX, IFF_N, KSPACE

      FF_INFO(1)  = XLAM0
      FF_INFO(2)  = ALN
      FF_INFO(3)  = XLAM_FINE_START
      FF_INFO(4)  = XLAM_FINE_END
      FF_INFO(7)  = FLOAT(IFF_N)
      FF_INFO(10) = FLOAT(KSPACE)

      IFF_DK(1) = 0 - 100

      RETURN
      END
