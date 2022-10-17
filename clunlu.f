      SUBROUTINE CLUNLU(ND, SUM_OPASMEAN, SUM_SMEAN, SUM_QFJMEAN, 
     >                  SUM_JMEAN, SUM_OPAJMEAN, SUM_QOPAHMEAN, 
     >                  SUM_HMEAN, SUM_HJMEAN, HTOTL, 
     >                  RADIUS, RSTAR, VDOP, TEFF, FTCOLI, T, TNEW, 
     >                  DUNLU, DUNLUR, RTCMAX, LRTCMAX, BPLOTUNLU)

      DIMENSION SUM_OPASMEAN(ND), SUM_SMEAN(ND), SUM_QFJMEAN(ND)
      DIMENSION SUM_JMEAN(ND), SUM_OPAJMEAN(ND), SUM_QOPAHMEAN(ND)
      DIMENSION SUM_HMEAN(ND), HTOTL(ND)
      DIMENSION FTCOLI(ND), RADIUS(ND), T(ND), TNEW(ND)
      LOGICAL BPLOTUNLU
      CHARACTER MODHEAD*100, HEADLINE*100, CENTER*8

C***  C1 = H * C / K        (CM * KELVIN)
      DATA C1 / 1.4388 /
C***  C2 = 2 * H * C     ( CGS UNITS )
      DATA C2 / 3.9724E-16 /
C***  DATA: STEBOL = STEFAN-BOLTZMANN CONSTANT (CGS-UNITS) / PI
      DATA STEBOL / 1.8046E-5 /

C***  PI8 = 8*PI
      DATA PI8 /25.1327412288 /

      DATA TMIN /3000./

      write (0,*) 'CLUNLU Start', bplotunlu

      IF (BPLOTUNLU) THEN
        OPEN (UNIT=100, FILE='coli_unlu.dat')
      ENDIF

        WRITE (0,'(7(A16))')
     >    'SUM_OPASMEAN', 'SUM_SMEAN', 'SUM_QFJMEAN',
     >    'SUM_JMEAN', 'SUM_OPAJMEAN', 'SUM_QOPAHMEAN',
     >    'SUM_HMEAN'
      DO L=1, ND
        WRITE (0,'(7(E15.7,1X))')
     >    SUM_OPASMEAN(L), SUM_SMEAN(L), SUM_QFJMEAN(L),
     >    SUM_JMEAN(L), SUM_OPAJMEAN(L), SUM_QOPAHMEAN(L),
     >    SUM_HMEAN(L)
      ENDDO
      WRITE (0,'(A,E15.7)') 'SUM_HJMEAN=', SUM_HJMEAN
      WRITE (0,*) '======================='
      WRITE (0,*) 

      VDOPCM = VDOP * 1.E5
      HNULL = 0.25 * STEBOL * TEFF*TEFF*TEFF*TEFF
      RTCMAX = 0.
      LRTCMAX = 0
      WRITE (0,*) 'DUNLU, DUNLUR=', DUNLU, DUNLUR
      WRITE (0,'(A8,A3,7A16)') 
     >     '        ', 'L', 'DTLOCAL', 'DTINT', 'DTRMAX', 'DTLO', 
     >     'DTL', 'T(L)', 'TNEW(L)'
      IF (BPLOTUNLU) THEN
      WRITE (100,'(A8,A3,7A16)') 
     >     '*       ', 'L', 'DTLOCAL', 'DTINT', 'DTRMAX', 'DTLO', 
     >     'DTL', 'T(L)', 'TNEW(L)'
      WRITE (100,*) '*DUNLU, DUNLUR=', DUNLU, DUNLUR
      ENDIF

      open (unit=199, file='unlu.doc')
      write (199,'(a,4(1x,e12.5))') 
     > 'htotl(1), teff, hnull, sum_jmean(1)=', 
     > htotl(1), teff, hnull, sum_jmean(1)

      DO L=1, ND
        RL2      = RADIUS(L) * RADIUS(L)
        TL3      = T(L) * T(L) * T(L)
        OPASMEAN = SUM_OPASMEAN(L) / (SUM_SMEAN(L) * RSTAR)
        QFJMEAN  = SUM_QFJMEAN(L) / SUM_JMEAN(L)
        OPAJMEAN = SUM_OPAJMEAN(L) / (SUM_JMEAN(L) * RSTAR)
C***  Local
        DTB = 1. / (4. * TL3 * STEBOL * OPASMEAN)
        DTLOCAL = -FTCOLI(L)/RSTAR * DTB
C***  Flux
        IF (L .GT. 1) THEN
          QOPAHMEAN = SUM_QOPAHMEAN(L-1) / HTOTL(L-1)
          DH = HTOTL(L-1) - HNULL
          WRADIUS = RADIUS(L) - RADIUS(L-1)
          DHINT = DHINT + DH * WRADIUS * QOPAHMEAN
        ELSE
          DHINT = 0.
          HJMEAN = SUM_HJMEAN / SUM_JMEAN(1)
          DHRMAX = (HNULL - HTOTL(1)) / HJMEAN - SUM_JMEAN(1) * RL2
        ENDIF

        DTINT = OPAJMEAN * DHINT / 
     >          (4. * TL3 * STEBOL * OPASMEAN * RL2 * QFJMEAN)
        DTRMAX = OPAJMEAN * DHRMAX / 
     >           (4. * TL3 * STEBOL * OPASMEAN * RL2)

        DTL = DTLOCAL + DUNLU*DTINT + DUNLUR*DTRMAX
        DTLO = DTL
C***  Maximum Correction 10 percent
        TRESH = 0.10
        IF (ABS(DTL) .GT. TRESH * T(L)) THEN
          IF (DTL .LT. 0.) THEN
            DTL = -TRESH * T(L)
          ELSE
            DTL = TRESH * T(L)
          ENDIF
        ENDIF
        TNEW(L) = T(L) + DTL
C***  REPLACE TEMPERATURES BELOW TMIN BY TMIN
        IF (TNEW(L) .LT. TMIN) THEN
          TNEW(L) = TMIN
        ENDIF

C***  Calculated Maximum of relative Corrections
        RTC = (TNEW(L)-T(L))/T(L)
        IF (ABS(RTC) .GE. ABS(RTCMAX)) THEN
          RTCMAX=RTC
          LRTCMAX=L
        ENDIF

        WRITE (0,'(A8,I3,7(E15.5,1X))') 
     >   '        ', L, DTLOCAL, DTINT, DTRMAX, DTLO, DTL, T(L), TNEW(L)

        IF (BPLOTUNLU) THEN
          WRITE (100,'(A8,I3,7(E15.5,1X))') 
     >     '        ', L, DTLOCAL, DTINT, DTRMAX, DTLO, DTL, 
     >     T(L), TNEW(L)
        ENDIF

        write (199,'(i3,9(1x,e12.5))')
     >    l, opasmean, qfjmean, opajmean, 
     >    t(l), dhint, dhrmax, ftcoli(l), hjmean, rl2

      ENDDO

      WRITE (0,35) RTCMAX, LRTCMAX, DUNLU, DUNLUR
      WRITE (*,35) RTCMAX, LRTCMAX, DUNLU, DUNLUR
      IF (BPLOTUNLU) THEN
        WRITE (100,36) RTCMAX, LRTCMAX, DUNLU, DUNLUR
      ENDIF
   35 FORMAT (15X,'COLI: MAX.REL.TEMP.CORR.=',F8.4,'  AT L=',I3,
     >        '    (DUNLU=', F5.2,1X,F5.2, ')'/)
   36 FORMAT ('*COLI: MAX.REL.TEMP.CORR.=',F8.4,'  AT L=',I3,
     >        '    (DUNLU=', F5.2,1X,F5.2, ')')

      IF (BPLOTUNLU) THEN
        CLOSE(100)
      ENDIF


      write (0,*) 'CLUNLU End'
      RETURN
      END
