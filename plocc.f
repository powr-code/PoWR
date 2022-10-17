      SUBROUTINE PLOCC (LPLOCC, KPLOCC, 
     >                  ND, NF, WCHARM, GAMMAC, DELTAC, MODHEAD, JOBNUM, 
     >                  OPC, KANAL)
C***********************************************************************
C***  Plot of Scharmer Continuum Cores (weight function)
C***    at Depth LPLOCC at frequency KPLOCC
C***********************************************************************

      DIMENSION WCHARM(ND,NF)
      CHARACTER MODHEAD*100, OPC*8
      CHARACTER*80 NHEAD, NX, NY
 
      write (0,*) 'PLOCC: Kanal=', kanal

C***  L-PLOT
      WRITE (NHEAD,'(A,I3)') '\CENTER\CCORE-Plot at Depth ', LPLOCC
      NX = '\CENTER\Frequency'
      NY = '\CENTER\1-e&H-#t#&M'
      CALL PLOTANF(KANAL, NHEAD, NHEAD, NX, NY, 
     >             0., 1., FLOAT(NF), 5., 10., 0., 
     >             0., 0., 1., 0.1, 0.2, 0., 
     >             X, Y, 0, 5)

      WRITE (KANAL, '(A)') 'N=? XYTABLE'
      DO K=1, NF
        WRITE (KANAL,'(I7,E20.10)') K, WCHARM(LPLOCC,K)
      ENDDO        
      WRITE (KANAL,'(A)') 'FINISH'
      WRITE (KANAL,'(A)') 'END'

C***  K-PLOT
      WRITE (NHEAD,'(A,I4)') '\CENTER\CCORE-Plot at Frequency ', KPLOCC
      NX = '\CENTER\Depth'
      NY = '\CENTER\1-e&H-#t#&M'
      CALL PLOTANF(KANAL, NHEAD, NHEAD, NX, NY, 
     >             0., 1., FLOAT(ND), 10., 20., 0., 
     >             0., 0., 1., 0.1, 0.2, 0., 
     >             X, Y, 0, 5)

      WRITE (KANAL, '(A)') 'N=? XYTABLE'
      DO L=1, ND
        WRITE (KANAL,'(I7,E20.10)') L, WCHARM(L,KPLOCC)
      ENDDO        
      WRITE (KANAL,'(A)') 'FINISH'
      WRITE (KANAL,'(A)') 'END'


      RETURN
      END
