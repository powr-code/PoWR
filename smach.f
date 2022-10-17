      FUNCTION SMACH (I)
C*********************************************************
C***     int       result
C***
C***     1         0.7105427357601002E-14
C***               The machine epsilon (the smallest positive machine number
C***               epsilon for which 1.0 +- epsilon is not 1.0).
C*********************************************************

      IF (I .EQ. 1) THEN
        SMACH = 0.7105427357601002E-14
      ELSE
        STOP 'ERROR IN SMACH'
      ENDIF

      RETURN
      END
