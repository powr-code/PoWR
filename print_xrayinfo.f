      SUBROUTINE PRINT_XRAYINFO (XDATA, MAXXDAT, XLOGLXLBOL)
C**************************************************************
C***  Prints a small block with XRAY data if specified
C***  Called from STEAL (only if model finally converged)
C**************************************************************

      DIMENSION XDATA(MAXXDAT)

      WRITE (*,'(/,A)') 'XRAY emission specifications:' 
      WRITE (*,'(A, G12.3)') 'XFILL = ', XDATA(1)
      WRITE (*,'(A, G12.3)') 'XRAYT = ', XDATA(2)
      WRITE (*,'(A, G12.3)') 'XRMIN = ', XDATA(3)

      IF (XDATA(5) .GT. .0) THEN
         WRITE (*,'(A)') 'Second plasma component:' 
         WRITE (*,'(A, G12.3)') 'XFILL2 = ', XDATA(5)
         WRITE (*,'(A, G12.3)') 'XRAYT2 = ', XDATA(6)
         WRITE (*,'(A, G12.3)') 'XRMIN2 = ', XDATA(7)
      ENDIF
 
      IF (XDATA(4) .GT. .0) THEN
         WRITE (*,'(A, G12.3)') 'Differential Emission Measure (DEM) '
     >      // 'with exponent = ', XDATA(4)
      ENDIF
 
      IF (XDATA(8) .LT. .0) THEN
         WRITE (*,'(A, G12.3)') 'XFILL was automatically adjusted '
     >      // 'in order to meet log (L_X / L_Bol) = ', XDATA(8)
         WRITE (*,'(A, G12.3)') 'The converged model has '
     >      // 'log (L_X / L_Bol) = ', XLOGLXLBOL
         WRITE (*,'(A,F7.3,A,F7.3,A)') 
     >       'L_X refers to a band from ', XDATA(9), 
     >       ' to ', XDATA(10), ' Angstroem'
         XMIN_KEV = 12.398 / XDATA(10)
         XMAX_KEV = 12.398 / XDATA(9)
         WRITE (*,'(A,F7.3,A,F7.3,A)') 
     >       'corresponding energies:   ', XMIN_KEV, 
     >       ' to ', XMAX_KEV, ' keV'
      ENDIF
 
C***  empty line
      WRITE (*,*)
      

      RETURN
      END
