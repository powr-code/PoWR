      SUBROUTINE PLOT_WINDROT_GRID (P, NPDIM, JPFIRST, JPLAST, 
     >      LPHISTA_ORIG, LPHIEND_ORIG, NPHI, 
     >      PHIARR, NPHIMAX, XPLOT, YPLOT)
C******************************************************************
C**   Plot of the ray positions, seen as facing the stellar disk 
C******************************************************************
      DIMENSION P(NPDIM), NPHI(NPDIM) 
      DIMENSION PHIARR(NPHIMAX,NPDIM)
      DIMENSION XPLOT(NPHIMAX), YPLOT(NPHIMAX)

      WRITE (65,*) 'PLOT: Grid for wind rotation'
      WRITE (65,*) '\NOBOX'
      NSET = 2 + JPLAST - JPFIRST
      IF (NSET .GT. 99) WRITE (65,'(A,I4)') '\SET_NSETMAX ', NSET  
    
      RMAX = 8.
      SCALE = 12./P(JPLAST)
      DO JP=JPFIRST, JPLAST
         WRITE (65,'(A,F10.4,A)') '\ARC 0 0 0 0 ', SCALE*P(JP), ' 0 180'
      ENDDO

      XPLOT(1)=.0
      YPLOT(1)=.0
      CALL PLOTANFS (65, '', '', '', '', 
     > SCALE, -P(JPLAST), P(JPLAST), 1., 10., 0.,
     > SCALE,     .0, P(JPLAST), 1., 10., 0.,
     > XPLOT, YPLOT, 1, 'SYMBOL=8 SIZE=-0.05 COLOR=1')

      DO JP=JPFIRST, JPLAST
C***     If specified, the plot is restricted to the points 
C***     between LPHISTA and LPHIEND
         LPHISTA = MAX(1,       LPHISTA_ORIG)
         LPHISTA = MIN(LPHISTA, NPHI(JP))
         LPHIEND = MIN(NPHI(JP),LPHIEND_ORIG)
         LPHIEND = MAX(LPHIEND, LPHISTA)

         DO LPHI=LPHISTA, LPHIEND
            XPLOT(LPHI) = P(JP) * COS(PHIARR(LPHI,JP))
            YPLOT(LPHI) = P(JP) * SIN(PHIARR(LPHI,JP))
         ENDDO

         CALL PLOTCONS (65, XPLOT(LPHISTA), YPLOT(LPHISTA), 
     >            1+LPHIEND-LPHISTA, 'SYMBOL=8 SIZE=-0.1 COLOR=2')
      ENDDO

      RETURN
      END
