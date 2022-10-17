      SUBROUTINE KSIGMA(SIGMAK,SIGMATHK,EDGEK,WAVENUM,SEXPOK)
C******************************************************************
C***  CLCULATE KSIGMA, THE FREQUENZ-DEPENDENT K-SHELL CROSS SECTION
C***  CALLED FROM SOBROUTINE COOP
C******************************************************************

     
      X=EDGEK/WAVENUM
      
      SIGMAK = SIGMATHK * 1.E-18 * X ** SEXPOK
      
      RETURN
      END
