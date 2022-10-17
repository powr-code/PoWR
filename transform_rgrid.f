      SUBROUTINE TRANSFORM_RGRID (VEC_NEW, N_NEW, VEC_OLD, N_OLD, 
     >                            R_NEW, R_OLD)
C********************************************************************
C***  Transformation of vector VEC_OLD (N_OLD) over radius-grid R_OLD
C***               into vector VEC_NEW (N_NEW) over radius-grid R_NEW
C***  If the new grid exceeds the old one: extrapolation as constant
C*********************************************************************

      DIMENSION VEC_NEW(N_NEW), R_NEW(N_NEW)
      DIMENSION VEC_OLD(N_OLD), R_OLD(N_OLD)

      DATA IWARN / 0 /

      IF (IWARN .EQ. 0) THEN 
         IF ( R_NEW(1) .GT. R_OLD(1) ) THEN 
         WRITE (0,'(A,/, A, F8.2,A,F8.2)') 
     >        '*** WARNING issued by subr. TRANSFORM_RGRID: ',  
     >        ' Extrapolation needed from ', R_NEW(1), ' TO ', R_OLD(1)
            IWARN = 1
         ENDIF
      ENDIF

      DO L=1, N_NEW
         IF (R_NEW(L) .GT. R_OLD(1)) THEN
            VEC_NEW(L) = VEC_OLD(1)
         ELSE
            CALL LIPO (VEC_NEW(L), R_NEW(L), VEC_OLD, R_OLD, N_OLD)
         ENDIF
      ENDDO

      RETURN
      END
