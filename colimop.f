      SUBROUTINE COLIMOP(ND, R, GRADI, VELO, 
     >                   DLF, DLH, GLF2, GLH2, VLF2, VLH2, BPLOT, 
     >                   RADIUS2, RADIUSH, RADIUSH2)
C*************************************************************
C***  Precalculation of the the geometry for the Moment-equation
C***    solved in COLIMO
C*************************************************************

      DIMENSION R(ND), GRADI(ND), VELO(ND)
      DIMENSION RADIUS2(ND), RADIUSH(ND), RADIUSH2(ND)
      DIMENSION DLF(ND), GLF2(ND), VLF2(ND)
      DIMENSION DLH(ND-1), GLH2(ND-1), VLH2(ND-1)
      LOGICAL BPLOT

C***  Prepare R^2, R at Interstices ...
      DO L=1, ND
        RADIUS2(L)  = R(L) * R(L)
        IF (L .EQ. ND) CYCLE
        RADIUSH(L)  = 0.5 * (R(L) + R(L+1))
        RADIUSH2(L) = RADIUSH(L) * RADIUSH(L)
      ENDDO

C***  DLF, GLF2 and VLF2 at L=2,ND-1
      DO L=2, ND-1
        DLF(L) = 0.5 * (R(L-1)-R(L+1))
        VLF2(L) = VELO(L) / R(L)
        GLF2(L) = GRADI(L) - VLF2(L)
      ENDDO

C***  DLF, GLF2 and VLF2 at Boundaries
      DLF(1) = R(1) - R(2)
      VLF2(1) = VELO(1) / R(1)
      GLF2(1) = GRADI(1) - VLF2(1)
      DLF(ND) = R(ND-1) - R(ND)
      VLF2(ND) = VELO(ND) / R(ND)
      GLF2(ND) = GRADI(ND) - VLF2(ND)

C***  DLH, GLH2 and VLH2 at Interstices
      DO L=1, ND-1
        DLH(L) = R(L) - R(L+1)
C!!!        VLH2(L) = (VELO(L+1)+VELO(L)) / (R(L+1)+R(L))
C!!!        GLH2(L) = -(VELO(L+1)-VELO(L)) / DLH(L) - VLH2(L)
C!!!  Neu von Goetz
        VLH2(L) = 0.5*(VLF2(L)+VLF2(L+1))
        GLH2(L) = 0.5*(GLF2(L)+GLF2(L+1))
      ENDDO

      IF (BPLOT) THEN
        OPEN (UNIT=99, FILE='colimop.dat')
        WRITE (99,'(A6,6A16)')
     >    '* L', 'DL', 'VL2', 'GL2', 'VELO', 'GRADI', 'RADIUS'
        WRITE (99,*) '*FULL'
        WRITE (99,*) 'N=?'
        DO L=1, ND
          WRITE (99,'(F6.2,1X,6(E15.7,1X))')
     >      FLOAT(L), DLF(L), VLF2(L), GLF2(L), VELO(L), GRADI(L), R(L)
        ENDDO
        WRITE (99,*) '*HALF'
        WRITE (99,*) 'N=?'
        DO L=1, ND-1
          WRITE (99,'(F6.2,1X,3(E15.7,1X))')
     >      FLOAT(L)+0.5, DLH(L), VLH2(L), GLH2(L)
        ENDDO
      ENDIF

      RETURN
      END
