      SUBROUTINE PRI_EPSG(BEMIX, BEMIXFIX, ITMAX, EPSG, NDDIM, ND, NIT, 
     >                    JOBNUM)
C****************************************************************
C***  An Overview on the EPSG Factors, which controls the EDDIMIX
C***    Scheme, is printed
C***    Called by COLI
C****************************************************************

      IMPLICIT NONE

      INTEGER, INTENT(IN) :: ITMAX, NDDIM, ND, NIT, JOBNUM
      REAL, DIMENSION(NDDIM,NIT), INTENT(IN) :: EPSG
      LOGICAL, INTENT(IN) :: BEMIX, BEMIXFIX

      INTEGER :: L, ITLOCAL
      REAL :: EMAX

        WRITE (0,*)
        WRITE (0,'(A,L1,A,L1,A)') 
     >    'Overview on current EPSG:    BEMIX=', BEMIX, 
     >    '   BEMIXFIX=', BEMIXFIX, 
     >    '      Note: Lines containing only zero are omitted'
        IF (ITMAX .EQ. 3) THEN
          WRITE (0,'(8X,5A8)') 
     >      'DEPTH', 'JOBNUM', 'IT=1', 'IT=2', 'IT=3'
        ELSE
          WRITE (0,'(8X,5A8)') 
     >      'DEPTH', 'JOBNUM', 'IT=1', 'IT=2'
        ENDIF
          DO L=1, ND
            EMAX = AMAX1(EPSG(L,1), EPSG(L,2), EPSG(L,3))
          IF (EMAX .EQ. 0.) CYCLE
            WRITE (0,'(A6,4X,I2,6X,I7,2X,3(1X,F7.2))')
     >        'EMIX: ', 
     >        L, JOBNUM, (EPSG(L,ITLOCAL), ITLOCAL=1, ITMAX)
          ENDDO
        WRITE (0,*)

      RETURN
      END
