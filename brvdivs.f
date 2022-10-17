      SUBROUTINE BRVDIVS (V, S, N, BRDIVFAILED)
C********************************************************************
C***  VEKTOR V WIRD DURCHDEN SKALAR S GETEILT   --   V = V / S
C********************************************************************
      DIMENSION V(N)
      LOGICAL BRDIVFAILED

      IF (S .EQ. .0) THEN
         CALL REMARK ('BRVDIVS: DIVISION BY ZERO')
         WRITE (0,*) 'VECTOR V IS NOT DIIVIDED BY S'
C         STOP 'ERROR'
         BRDIVFAILED = .TRUE.
         GOTO 99
         ENDIF

      DO 1 I=1,N
        V(I) = V(I) / S
    1 CONTINUE

   99 CONTINUE

      RETURN
      END
