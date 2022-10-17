C**********************************************************************
C***  If a second model is employed for a part of the atmosphere: 
C***    the radius points of both models were combined to a new 
C***    RADIUS_MERGED grid.
C***  This subroutine now interpolates all relevent arrays to the  
C***     new RADIUS_MERGED grid 
C***  Note: this routine is also abused after NOWIND for restricting 
C***        the radius grid, irrespective whether NMOD=1 or 2  
C***
C***  Called from: FORMAL
C**********************************************************************

      SUBROUTINE RESCALE_SECMOD (NDDIM, NFLDIM, MAXLAP, MAXMOD, MAXATOM, 
     >   ND, NATOM, NFL, NBLINE, RADIUS, RADIUS_MERGED, ND_MERGED,
     >   VDU, OPAL, ETAL, ETACK, ETACCK, OPACK, OPAFE, ETAFE, 
     >   DD_VDOPDU, DD_VDOP, DD_VMICDU, T, ENTOT, RNE, VEC_SECMOD, 
     >   RSTAR, NMOD)

      DIMENSION ND(MAXMOD), RSTAR(MAXMOD), VEC_SECMOD(NDDIM)  
      REAL, DIMENSION (NDDIM,MAXATOM,MAXMOD) :: DD_VDOPDU, DD_VDOP 
      REAL, DIMENSION (NDDIM, MAXMOD) :: 
     >          RADIUS, VDU, T, ENTOT, RNE, DD_VMICDU      
      REAL, DIMENSION (NDDIM,NFLDIM,MAXMOD) :: ETACK, ETACCK, OPACK,
     >                                         OPAFE, ETAFE
      REAL, DIMENSION (NDDIM,MAXLAP,MAXMOD) :: OPAL, ETAL
      REAL, DIMENSION (NDDIM) :: RADIUS_MERGED
      DIMENSION QRSTAR(2)
      REAL, PARAMETER :: RSUN = 6.96E10 ! solar radius in cm

ccc   temporary for test plots:
      DIMENSION XPLOT(100), YPLOT(100)

C***  In each of the models, opacities and emissivities are per  
C***  RSTAR of the respctive model, which might differ.
C***  The flux output (Subr. TRAPLO) multiplies the flux with 
C***  the square of RSTAR(1), i.e. adopting the radius of the first
C***  (main) model. Hence, the opacities and emisivities of the 
C***  second model must be re-scaled by the factor QRSTAR(2)

      QRSTAR(1) = 1.

      IF (NMOD .GT. 1) THEN
       QRSTAR(2) = RSTAR(1) / RSTAR(2)

       IF (ABS(QRSTAR(2)-1.0) .GT. 0.001 ) THEN
        WRITE (0,'(A)')      'WARNING: second model has different RSTAR'  
        WRITE (0,'(A,F8.3)') '         RSTAR(1) =', RSTAR(1)/RSUN   
        WRITE (0,'(A,F8.3)') '         RSTAR(2) =', RSTAR(2)/RSUN   
        WRITE (0,'(A,F8.3)') 'second model scaled by factor ', QRSTAR(2)  
        WRITE (0,'(A)')      'Note that this also affects log g and Mdot'  
       ENDIF
      ENDIF
  
      DO IMOD= 1, NMOD

      VEC_SECMOD = VDU(1:ND(IMOD),IMOD) 
      CALL TRANSFORM_RGRID (VDU(1,IMOD), ND_MERGED, VEC_SECMOD, ND(IMOD),
     >                            RADIUS_MERGED, RADIUS(1,IMOD))

C***  Note: There is only one VMIC option in FORMAL_CARDS which
C***        holds for both models, but since there are parameters 
C***        which refer tp RADIUS or VELO, vmic(r) might differ
C***        between the first and the second model
      VEC_SECMOD = DD_VMICDU(1:ND(IMOD),IMOD) 
      CALL TRANSFORM_RGRID (DD_VMICDU(1,IMOD), ND_MERGED, VEC_SECMOD,
     >                 ND(IMOD), RADIUS_MERGED, RADIUS(1,IMOD))

      VEC_SECMOD = T(1:ND(IMOD),IMOD) 
      CALL TRANSFORM_RGRID (T(1,IMOD), ND_MERGED, VEC_SECMOD, ND(IMOD),
     >                            RADIUS_MERGED, RADIUS(1,IMOD))

      VEC_SECMOD = ENTOT(1:ND(IMOD),IMOD) 
      CALL TRANSFORM_RGRID (ENTOT(1,IMOD), ND_MERGED, VEC_SECMOD, ND(IMOD),
     >                            RADIUS_MERGED, RADIUS(1,IMOD))

      VEC_SECMOD = RNE(1:ND(IMOD),IMOD) 
      CALL TRANSFORM_RGRID (RNE(1,IMOD), ND_MERGED, VEC_SECMOD, ND(IMOD),
     >                            RADIUS_MERGED, RADIUS(1,IMOD))

      DO NA=1, NATOM

         VEC_SECMOD = DD_VDOPDU(1:ND(IMOD),NA,IMOD) 
         CALL TRANSFORM_RGRID (DD_VDOPDU(1,NA,IMOD), ND_MERGED, VEC_SECMOD,
     >                 ND(IMOD), RADIUS_MERGED, RADIUS(1,IMOD))

         VEC_SECMOD = DD_VDOP(1:ND(IMOD),NA,IMOD) 
         CALL TRANSFORM_RGRID (DD_VDOP(1,NA,IMOD), ND_MERGED, VEC_SECMOD,
     >                 ND(IMOD), RADIUS_MERGED, RADIUS(1,IMOD))
      ENDDO

      DO NBL=1, NBLINE
     
         VEC_SECMOD = OPAL(1:ND(IMOD),NBL,IMOD) * QRSTAR(IMOD) 
         CALL TRANSFORM_RGRID (OPAL(1,NBL,IMOD), ND_MERGED, VEC_SECMOD,
     >                 ND(IMOD), RADIUS_MERGED, RADIUS(1,IMOD))

C******* Test plot facility
        IF (nbl .eq. 0) then
C          note: after rescale_secmod, all vectors are over radius of imod=1
c        do l=1, ND(2)
        do l=1, ND(1)
           xplot(L) = alog10(RADIUS(L,1)-.999)
c           yplot(L) = alog10(opal(L,nbl,1))
           yplot(L) = opal(L,nbl,1)
c           yplot(L) = VEC_SECMOD(L)
        enddo

        CALL PLOTANFS (77, '', '', '', '',
     >        0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.,
c     >        XPLOT, YPLOT, ND(2), 'SYMBOL=8 SIZE=-0.05 COLOR=2')
     >        XPLOT, YPLOT, ND(1), 'SYMBOL=8 SIZE=-0.05 COLOR=2')

        do l=1, ND(1)
           xplot(L) = alog10(RADIUS(L,1)-.999)
c           yplot(L) = alog10(opal(L,nbl,2))
           yplot(L) = opal(L,nbl,2)
c           yplot(L) = vdu(L,2)
        enddo

        CALL PLOTCONS (77, XPLOT, YPLOT, 
     >                ND(1), 'SYMBOL=5 COLOR=4')

      endif
***************************************************************


         VEC_SECMOD = ETAL(1:ND(IMOD),NBL,IMOD) * QRSTAR(IMOD)
         CALL TRANSFORM_RGRID (ETAL(1,NBL,IMOD), ND_MERGED, VEC_SECMOD,
     >                 ND(IMOD), RADIUS_MERGED, RADIUS(1,IMOD))

      ENDDO

C***  Loop over frequencies
      DO K=1, NFL

         VEC_SECMOD = ETACK(1:ND(IMOD),K,IMOD) * QRSTAR(IMOD)
         CALL TRANSFORM_RGRID (ETACK(1,K,IMOD), ND_MERGED, VEC_SECMOD,
     >                 ND(IMOD), RADIUS_MERGED, RADIUS(1,IMOD))

         VEC_SECMOD = ETACCK(1:ND(IMOD),K,IMOD) * QRSTAR(IMOD)
         CALL TRANSFORM_RGRID (ETACCK(1,K,IMOD), ND_MERGED, VEC_SECMOD,
     >                 ND(IMOD), RADIUS_MERGED, RADIUS(1,IMOD))

         VEC_SECMOD = OPACK(1:ND(IMOD),K,IMOD) * QRSTAR(IMOD)
         CALL TRANSFORM_RGRID (OPACK(1,K,IMOD), ND_MERGED, VEC_SECMOD,
     >                 ND(IMOD), RADIUS_MERGED, RADIUS(1,IMOD))

         VEC_SECMOD = OPAFE(1:ND(IMOD),K,IMOD) * QRSTAR(IMOD)
         CALL TRANSFORM_RGRID (OPAFE(1,K,IMOD), ND_MERGED, VEC_SECMOD,
     >                 ND(IMOD), RADIUS_MERGED, RADIUS(1,IMOD))

         VEC_SECMOD = ETAFE(1:ND(IMOD),K,IMOD) * QRSTAR(IMOD)
         CALL TRANSFORM_RGRID (ETAFE(1,K,IMOD), ND_MERGED, VEC_SECMOD,
     >                 ND(IMOD), RADIUS_MERGED, RADIUS(1,IMOD))

      ENDDO
C***  End of loop over frequencies

      ENDDO
C***  End of loop over IMOD (model 1, 2)

C******* Test plot facility
      k = 1
      IF (k .eq. 1) then
        do l=1, ND_MERGED
           xplot(L) = alog10(RADIUS_MERGED(L)-.999)
           yplot(L) = alog10(etacck(L,k,1))
        enddo

        CALL PLOTANFS (77, '', '&2mod1 rescaled  &4mod2 rescaled' , 
     >        'log ETACCK', 'log (radius - 0.999)',
     >        0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.,
     >        XPLOT, YPLOT, ND_MERGED, 'SYMBOL=8 SIZE=-0.05 COLOR=2')

        do l=1, ND_MERGED
            xplot(L) = alog10(RADIUS_MERGED(L)-.999)
            if (etacck(L,k,2) .gt. .0) then
               yplot(L) = alog10(etacck(L,k,2))
            else
               yplot(L) = .0
            endif
        enddo

        CALL PLOTCONS (77, XPLOT, YPLOT, 
     >                ND_MERGED, 'SYMBOL=5 COLOR=4')

      endif
***************************************************************


      RETURN

C**********************************************************************
C***  ERROR BRANCHES 
C**********************************************************************


      END
