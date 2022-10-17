      SUBROUTINE GEOMESH (RADIUS, INCRIT,P,Z,ND,NDDIM,NP,NPDIM,RMAX, 
     >                    RadiusGridParameters, bNoRGrid, NC)
C***********************************************************************
C***  THIS SUBROUTINE GENERATES THE GEOMETRICAL POINT MESH IN R, P AND Z
C****   called only by WRSTART (wrh, goetz) and STEAL->HYDROSOLVE (wrh)
C***********************************************************************
  
      IMPLICIT NONE

      LOGICAL bNoRGrid
      INTEGER :: ND, NDDIM, NP, NPDIM, NC
      REAL, DIMENSION(NDDIM) :: R
      CHARACTER(8), DIMENSION(NDDIM) :: INCRIT      
      REAL, DIMENSION(ND) :: RADIUS, P
      REAL, DIMENSION(2) :: Z
      REAL :: RMAX, RR, PJ, PJPJ
      CHARACTER(80), DIMENSION(3) :: RadiusGridParameters
      
      INTEGER I, J, L, JMAX
      
      IF (.NOT. bNoRGrid) THEN
        CALL RGRID (NDDIM,ND,RADIUS,INCRIT,RMAX, 
     >              RadiusGridParameters)
      ENDIF

      CALL PGRID (NPDIM,NP,ND,RADIUS,P, NC)
      DO L=1,ND
        RR=RADIUS(L)*RADIUS(L)
        JMAX=NP+1-L
        DO J=1,JMAX
          PJ=P(J)
          PJPJ=PJ*PJ
          I=(J-1)*ND+L
          Z(I)=SQRT(RR-PJPJ)
        ENDDO
      ENDDO

      RETURN
      END
