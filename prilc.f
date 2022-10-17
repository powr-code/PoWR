      SUBROUTINE PRILC (IPRILC,LASTIND,XRED,XBLUE,TAUMIN,L,ND,ERXMIN,
     $                  MODHEAD,JOBNUM,GAMPRI)
C***********************************************************************
C***  PRINTOUT OF SCHARMER LINE CORES  *********************************
C***********************************************************************

      IMPLICIT NONE

      INTEGER, INTENT(IN) :: L, ND, JOBNUM, LASTIND
      INTEGER, INTENT(INOUT) :: IPRILC
      REAL, INTENT(IN) :: ERXMIN

      REAL, DIMENSION(2) :: XRED, XBLUE
      INTEGER, DIMENSION(2) :: NCHARG
      REAL, DIMENSION(0:2) :: TAUMIN, AMP, GAMPRI
      CHARACTER(100) :: MODHEAD

      INTEGER :: I
      REAL :: FC

      REAL, EXTERNAL:: ERF
 
      IPRILC=MIN0(IPRILC,LASTIND-2)
C***  PRINT THE HEADER OF THE TABLE
      IF (L .EQ. 1) THEN
      PRINT 1,MODHEAD,JOBNUM
    1 FORMAT (1X,  A  ,20X,'JOB NO.',I7,//)
      PRINT 10,(IPRILC+I,GAMPRI(I),I=0,2)
   10 FORMAT (       10X,'SCHARMER-CORES',//,5X,
     $   3('   -----  LINE',I4,'  --  GAMMA=',F5.1,'  -----'),//,
     $ '  L   ',3('  XRED   XBLUE  TAUMIN      AMP           '),/)
      ENDIF
 
      DO 2 I=0,2
      FC=(ERF(XBLUE(IPRILC+I))-ERF(XRED(IPRILC+I)))/(1.-2.*ERXMIN)
      AMP(I)=1./(1.-FC)
    2 CONTINUE
      PRINT 12, L,
     $        (XRED(IPRILC+I),XBLUE(IPRILC+I),TAUMIN(I),AMP(I),I=0,2)
   12 FORMAT (I3,2X,3(2F7.2,1P,G11.2,G11.2,0P,6X))

      RETURN
      END
