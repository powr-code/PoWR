      SUBROUTINE SHIFT (R,L,ND)
C***********************************************************************
C***  SHIFTS ARRAY ELEMENTS R(L) TO R(ND) BY ONE INDEX
C***********************************************************************
      DIMENSION R(ND)

      IF (L.LT.1 .OR. L.GT.ND) THEN
        WRITE (0,'(A,2I4)') 'ERROR IN SUBR. SHIFT, L,ND=',L, ND
        STOP 'SHIFT'
      ENDIF
      DO 1 II=L,ND
      I=ND+L-II
    1 R(I+1)=R(I)

      RETURN
      END
