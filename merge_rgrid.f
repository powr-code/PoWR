         SUBROUTINE MERGE_RGRID (RADIUS, NDDIM, MAXMOD, ND, 
     >      RADIUS_MERGED, ND_MERGED, PGRID, NP, NPDIM, PGRID_MERGED,
     >                     NP_MERGED, ZGRID_MERGED, SECMOD_RRANGE)
C******************************************************
C***  The RADIUS-Grids of each model are optimized; 
C***  In order to secure the numerical accuracy of the integration
C***  in case of the combination of two models, these grids must
C***  be combined. 
C******************************************************

      DIMENSION RADIUS(NDDIM,MAXMOD), RADIUS_MERGED(NDDIM)
      DIMENSION ND(MAXMOD)
      DIMENSION PGRID(NPDIM), PGRID_MERGED(NPDIM) 
      DIMENSION ZGRID_MERGED(NDDIM*NPDIM)
      REAL, DIMENSION(2) :: SECMOD_RRANGE

      RMAX = MIN (RADIUS(1,1), RADIUS(1,2))
      RADIUS_MERGED(1) = MIN (RADIUS(1,1), RADIUS(1,2))

      IF (RADIUS(1,1) .NE. RADIUS(1,2)) THEN
         WRITE (0,'(A)')      'Note:  RMAX of the two models differ:'
         WRITE (0,'(A,F6.1)') '         Main model: RMAX=', RADIUS(1,1) 
         WRITE (0,'(A,F6.1)') '       Second model: RMAX=', RADIUS(1,2) 
         WRITE (0,'(A,F6.1)') '--> Minimum adopted: RMAX=', RMAX 
      ENDIF

C***  The "density" of RADIUS-Points is compared between the two models;
C***  the merged-grid points are taken from that model where the points
C***  are denser. 
C***  At radii outside the SECMOD_RRANGE, only model 1 (the main model) 
C***  is taken; this can happen when the second-model region has the
C***  SHAPE od a SHERE

      L = 2
C***  RADIUS grids are restricted to values below the new RMAX 
      L1= ISRCHFLT(ND(1), RADIUS(1,1), 1, RADIUS_MERGED(1))
      L2= ISRCHFLT(ND(2), RADIUS(1,2), 1, RADIUS_MERGED(1))

    1 CONTINUE
      IF (RADIUS(L1-1,1) - RADIUS(L1+1,1) .LT.
     >    RADIUS(L2-1,2) - RADIUS(L2+1,2) .OR. 
     >    (RADIUS(L1,1)-SECMOD_RRANGE(1)) * 
     >    (RADIUS(L1,1)-SECMOD_RRANGE(2)) .GE. .0) THEN
          RADIUS_MERGED(L) = RADIUS(L1,1)
ccc          write (0,*) 'from 1: L, RADIUS_MERGED =', L, RADIUS_MERGED(L)
      ELSE
          RADIUS_MERGED(L) = RADIUS(L2,2)
ccc          write (0,*) 'from 2: L, RADIUS_MERGED =', L, RADIUS_MERGED(L)
      ENDIF


      L1 = ISRCHFLT(ND(1), RADIUS(1,1), 1, RADIUS_MERGED(L))
      L2 = ISRCHFLT(ND(2), RADIUS(1,2), 1, RADIUS_MERGED(L))
      L  = L  + 1
      IF (L .GT. NDDIM) THEN
         WRITE (0,*) 'NDDIM insufficient for merging RADIUS-Grids'
         STOP 'FATAL ERROR in Subr. MERGE_RGRID'
      ENDIF
      IF (L1 .LT. ND(1) .AND. L2 .LT. ND(2)) GOTO 1

      RADIUS_MERGED(L) = 1.
      ND_MERGED = L

      WRITE (0,'(A)') 'MERGED RADIUS GRID CONSTRUCTED'
      WRITE (0,'(A,I4)') 'MODEL 1: ND =', ND(1)
      WRITE (0,'(A,I4)') 'MODEL 2: ND =', ND(2)
      WRITE (0,'(A,I4)') 'Merged : ND =', ND_MERGED

C***  Construct PGRID_MERGED on basis of the new RADIUS_MERGED
C***  Note: PGRID has actually two vectors for each model, 
C***        but only the first one is used here to maintain the core-rays 
C***        which might have been changed by ROTATION_PREP

      NCORE = NP - ND(1)
      DO JP = 1, NCORE
         PGRID_MERGED(JP) = PGRID(JP)
      ENDDO

      NP_MERGED = NCORE + ND_MERGED
      IF (NP_MERGED .GT. NPDIM) THEN
         WRITE (0,'(A)') '*** DIMENSION NPDIM insufficient!'
         WRITE (0,'(A,I4)') '*** presently: NPDIM =', NPDIM
         WRITE (0,'(A,I4)') '*** needed   : NPDIM =', NP_MERGED
         STOP '*** FATAL ERROR IN SUBR. MERGE_RGRID'
      ENDIF

      DO L=1,ND_MERGED
         JP = NP_MERGED + 1 - L
         PGRID_MERGED(JP) = RADIUS_MERGED(L)
      ENDDO


C***  Z-array must be updated with new RADIUS_MERGED and the new
C***        PGRID_MERGED
C***  Note: the Z array is filled in the upper-left corner 
C***  The impact-parameters P are taken from the first (main) model,
C***  where the core-rys might have been increased in ROTATION_PREP
      DO L = 1, ND_MERGED
        RR = RADIUS_MERGED(L) * RADIUS_MERGED(L)
        JMAX = NP_MERGED + 1 - L
        DO JP = 1, JMAX
          PJ = PGRID_MERGED(JP)
          PJPJ = PJ * PJ
          I=(JP-1)*ND_MERGED+L
          IF ( (RR-PJPJ) .GT. .0) THEN
             ZGRID_MERGED(I) = SQRT(RR-PJPJ)
          ELSE
             ZGRID_MERGED(I) = .0
          ENDIF
        ENDDO
      ENDDO

C***  test output
c      jp = 1
c      LMAX=MIN0(NP_MERGED+1-JP,ND_MERGED)
c      do l=1, lmax
c         I=(JP-1)*ND_MERGED+L
c         write (0,*) 'L, ZGRID_MERGED =', L, ZGRID_MERGED(i)
c      enddo 

      RETURN
      END
