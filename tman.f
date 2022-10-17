      SUBROUTINE TMAN
C*** One-purpose program for manipulating the temperature structure on 
C*** the MODEL file (fort.3)

      DIMENSION T(200)

      CALL OPENMS(3,IADR,MAXADR,1, IERR)
      CALL READMS(3,ND,1,       'ND      ', IERR)

      CALL READMS(3,T,ND,       'T       ', IERR)

ccc      DO L=1, ND
      DO L=1, 25
        T(L) = T(L) * 0.8
cc        T(L) = 25000.
      ENDDO

      CALL WRITMS (3,T,ND,'T       ',-1, IDUMMY, IERR)
      CALL CLOSMS (3, IERR)

      RETURN
      END
