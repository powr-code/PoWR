      SUBROUTINE PLOT_SECONDMODEL_GRID (P, NP, NPDIM, NPHI, PHIARR,
     >           NPHIMAX, JPFIRST, JPLAST, LPHISTA_ORIG, LPHIEND_ORIG, 
     >           ZINTER)
C******************************************************************
C**   Plot of the ray positions, seen as facing the stellar disk 
C******************************************************************
      DIMENSION P(NPDIM), NPHI(NPDIM) 
      DIMENSION PHIARR(NPHIMAX,NPDIM)
      PARAMETER ( MAXPLOT = 1000 )
      DIMENSION XPLOT(MAXPLOT), YPLOT(MAXPLOT)
      DIMENSION ZINTER (2,NPDIM, NPHIMAX)

      WRITE (66,*) 'PLOT: Grid for second-model geometry'
      WRITE (66,*) '\NOBOX'
      WRITE (66,*) '\OFS 4 2'
      IF (2*NP .GT. 99) WRITE (66,'(A,I4)') '\SET_NSETMAX ',  2*NP+1  
    
      SCALE = 9./ P(JPLAST)
      DO JP=1, JPLAST
         WRITE (66,'(A,F10.4,A)') '\ARC 0 0 0 0 ', SCALE*P(JP), ' 0 360'
         WRITE (66,'(A,F10.4,A, I3)') 
     >             '\LUN 0 0 ', SCALE*P(JP), ' 1 .2 ', JP
      ENDDO

      XPLOT(1)=.0
      YPLOT(1)=.0

      CALL PLOTANFS (66, '', '', '', '', 
     > SCALE, -P(JPLAST), P(JPLAST), 1., 10., 0.,
     > SCALE, -P(JPLAST), P(JPLAST), 1., 10., 0.,
     > XPLOT, YPLOT, 1, 'SYMBOL=8 SIZE=-0.05 COLOR=1')

      DO JP=JPFIRST, JPLAST
         IF (NPHI(JP) .GT. MAXPLOT) THEN
            WRITE (*,*) 'ERROR: INSUFFICIENT DIMENSION' 
            STOP '*** ERROR in subr. PLOT_SECONDMODEL_GRID'
         ENDIF
         LPHISTA = MAX(1,       LPHISTA_ORIG)
         LPHISTA = MIN(LPHISTA, NPHI(JP))
         LPHIEND = MIN(NPHI(JP),LPHIEND_ORIG)
         LPHIEND = MAX(LPHIEND, LPHISTA)

C***     Red dots: intersecting with SECONDMODEL domain
         NCOUNT = 0
         DO LPHI= LPHISTA, LPHIEND
            IF (ZINTER(1,JP,LPHI) .NE. ZINTER(2,JP,LPHI)) THEN
               NCOUNT = NCOUNT + 1 
               XPLOT(NCOUNT) = P(JP) * COS(PHIARR(LPHI,JP))
               YPLOT(NCOUNT) = P(JP) * SIN(PHIARR(LPHI,JP))
            ENDIF
         ENDDO
         IF (NCOUNT .GT. 0) CALL PLOTCONS (66, XPLOT, YPLOT, 
     >            NCOUNT, 'SYMBOL=8 SIZE=-0.1 COLOR=2')

C***     Blue dots: not intersecting with SECONDMODEL domain
         NCOUNT = 0
         DO LPHI= LPHISTA, LPHIEND
            IF (ZINTER(1,JP,LPHI) .EQ. ZINTER(2,JP,LPHI)) THEN
               NCOUNT = NCOUNT + 1 
               XPLOT(NCOUNT) = P(JP) * COS(PHIARR(LPHI,JP))
               YPLOT(NCOUNT) = P(JP) * SIN(PHIARR(LPHI,JP))
            ENDIF
         ENDDO

         IF (NCOUNT .GT. 0) CALL PLOTCONS (66, XPLOT, YPLOT, 
     >            NCOUNT, 'SYMBOL=8 SIZE=-0.1 COLOR=4')

      ENDDO

      RETURN
      END
