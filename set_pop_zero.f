      SUBROUTINE SET_POP_ZERO(LEVEL, NDIM, N, SPZ1, SPZ2, POPNUM, 
     >                        ND, NMOD, MAXSPZ)
C************************************************************************
C***  Some selected Popnums are set to Zero 
C***  Called by FORMALCL
C************************************************************************


      DIMENSION POPNUM(ND,N,NMOD)
      CHARACTER*(*) SPZ1(MAXSPZ), SPZ2(MAXSPZ)
      CHARACTER*10 LEVEL(NDIM)

      !File and channel handles (=KANAL)
      INTEGER, PARAMETER :: hOUT   = 6    !write to wruniqX.out (stdout)
      INTEGER, PARAMETER :: hCPR   = 0    !write to wruniqX.cpr (stderr)

      INTEGER, EXTERNAL :: IDX

      DO ISPZ = 1, MAXSPZ

      ID_SPZ1 = IDX(SPZ1(ISPZ))
      ID_SPZ2 = IDX(SPZ2(ISPZ))
      IF (ID_SPZ1 <= 0) RETURN

c      write (0,*) 'SPZ1=', SPZ1(ISPZ)
c      write (0,*) 'SPZ2=', SPZ2(ISPZ)


      DO I=1, N
        IF ((LEVEL(I)(1:ID_SPZ1) /= SPZ1(ISPZ)(1:ID_SPZ1)) .OR.
     >      (ID_SPZ2 > 0 .AND.  
     >       LEVEL(I)(1:ID_SPZ2) == SPZ2(ISPZ)(1:ID_SPZ2))) THEN
          CYCLE
        ENDIF
        WRITE (0,'(A,A,A)') 'SET_POP_ZERO: Popnums for Level ',
     >                              LEVEL(I), ' have been nulled.'
        DO L=1, ND
          DO M=1, NMOD
            POPNUM(L,I,M) = 0.
          ENDDO
        ENDDO
      ENDDO

      ENDDO ! loop over all input cards

      RETURN
      END
