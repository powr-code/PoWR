      SUBROUTINE INSTALL
C***********************************************************************
C***  Called from Main Programs
C***********************************************************************

C***  Operating system:
      COMMON / COMOS / OPSYS
      CHARACTER*8 OPSYS

C*********************************************
C***  Hier Hauptschalter:
C***  INST = 1 : Cray
C***  INST = 2 : Potsdam DEC/UNIX
C***  INST = 3 : Potsdam SGI Origin 2000
C*********************************************

C     vvvvvvvv
      INST = 2
C     ^^^^^^^^

C******  Cray  **************************************
      IF (INST .EQ. 1) THEN

      OPSYS = 'CRAY'

C******  DEC  **************************************
      ELSE IF (INST .EQ. 2) THEN

      OPSYS = 'DEC/UNIX'

C******  SGI  **************************************
      ELSE IF (INST .EQ. 3) THEN

      OPSYS = 'SGI'

C****** ERROR  **************************************
      ELSE
      STOP 'ERROR IN INSTALL'
      ENDIF

      RETURN
      END
