      SUBROUTINE PRI_PAR (TEFF, RSTAR, VFIN, DENSCON, XMDOT)
C     ************************************************************
C     * Short model parameter printout                           *
C     *  This routine is only called from FORMAL                 *
C     ************************************************************

      IMPLICIT NONE

      REAL, PARAMETER :: RSUN = 6.96E10 ! solar radius in cm
      REAL DENSCON, RT, RSTAR, VFIN, TEFF, RSTAR_S, XMDOT

      RSTAR_S = RSTAR / RSUN
      RT = (VFIN/2500. * 1.E-4 / 10.**XMDOT /SQRT(DENSCON))**(2./3.) 
     >      * RSTAR_S
      
      WRITE (*,*)
      WRITE (*,*) '--------------------------'
      WRITE (*,*) 'Model Parameters'
      WRITE (*,*) 
      WRITE (*,'(A, F7.0)') ' T* = ', TEFF
      WRITE (*,'(A, 1PG12.3,0P,A,F7.3,A)') 
     >      ' Rt = ', RT, ' = ', ALOG10(RT), 'dex'
      WRITE (*,'(A, F7.4)') ' R* = ', RSTAR_S
      WRITE (*,'(A, F7.3,A)') ' Md = ', XMDOT, 'dex'
      WRITE (*,'(A, F7.1)') ' V8 = ', VFIN
      WRITE (*,'(A, F7.2)') ' D(L=1)  = ', DENSCON
      WRITE (*,*) '--------------------------'
      WRITE (*,*) 

      RETURN
      END
