C***  MAIN PROGRAM GF *********************************************
      SUBROUTINE GF
C***********************************************************************
C***  Writes Variable GF in Model-file. The default is 100.
C***  This forces the AUTO-GAMMA procedure to set the GAMMAs to GF
C***    in the next ALI-Iteration
C***********************************************************************
 
      PARAMETER ( MAXADR  =       25600 )
      DIMENSION IADR(MAXADR)

C***  READING OF THE MODEL FILE ----------------------------------------
      CALL OPENMS (3, IADR, MAXADR, 1, IERR)

      GF_DEF = 100.
      CALL WRITMS (3,GF_DEF,1,'GF      ',-1, IDUMMY, IERR)

      CALL CLOSMS (3, IERR)

      CALL JSYMSET ('G0','0')

      STOP 'O.K.'
      END
