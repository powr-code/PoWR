      SUBROUTINE MULTSPLI(ND, N, NSUBLOW, MAXSUBL, NLOW, POPNUM, 
     >                 ELEVEL, T, WEIGHT, NNUP, NSUBNUP, LOW, NUP, 
     >                 LEVEL)

      DIMENSION NSUBLOW(MAXSUBL), NSUBNUP(MAXSUBL)
      DIMENSION POPNUM(ND, N), T(ND), ELEVEL(N), WEIGHT(N)
      CHARACTER*10 LEVEL(N) 
C***  C1 = H * C / K    ( CM * ANGSTROEM )
      DATA C1 / 1.4388 /
      NWEIGHTLOWREST = 0
      NWEIGHTNUPREST = 0      
C***  A. SPLITTING THE LOWER ENERGY LEVEL:
      IF (NLOW .GT. 0) THEN
C***     Check if sum of sublevel weights conforms with weight of original level
         WEIGHTSUM = WEIGHT(NSUBLOW(1))
          ELOWREST = WEIGHT(NSUBLOW(1)) * ELEVEL(NSUBLOW(1))
         DO J=2, NLOW
            WEIGHTSUM = WEIGHTSUM + WEIGHT(NSUBLOW(J))
            ELOWREST = ELOWREST + WEIGHT(NSUBLOW(J))*ELEVEL(NSUBLOW(J))
         ENDDO
         IF (NINT(WEIGHTSUM) .NE. NINT(WEIGHT(LOW))) THEN
            WRITE (*,'(A)') '*** WARNING: Inconsistency of stat. weights'
            WRITE (*,'(A, I3, A, I3)') 'LOWERLEVEL: ' // LEVEL(LOW)
     >           // '    WEIGHT: ', NINT(WEIGHT(LOW)),  
     >              '  =!= Sum of sublevel weights: ', NINT(WEIGHTSUM)
            WRITE (0,'(A)') 'WARNING: Inconsistency of WEIGHTS '
     >           // 'for level ' // LEVEL(LOW) // '  -- see output file'
            ! construct rest level:
            NWEIGHTLOWREST = NINT(WEIGHT(LOW) - WEIGHTSUM)
            ELOWREST=(ELEVEL(LOW) *WEIGHT(LOW)-ELOWREST)/NWEIGHTLOWREST
         ENDIF  

         DO 30 L=1,ND
         POPNUM(L,NSUBLOW(1))=1.
         DO 35 J=2,NLOW
   35    POPNUM(L,NSUBLOW(J))=EXP(C1*(ELEVEL(NSUBLOW(J-1))
     -         -ELEVEL(NSUBLOW(J)))/T(L))*WEIGHT(NSUBLOW(J))
     /         /WEIGHT(NSUBLOW(J-1))*POPNUM(L,NSUBLOW(J-1))
C***  NORMALIZATION
         SUM=0.
         IF (NWEIGHTLOWREST .GT. 0) THEN
            SUM=EXP(C1*(ELEVEL(NSUBLOW(1))-ELOWREST)/T(L))
     >           *NWEIGHTLOWREST/WEIGHT(NSUBLOW(1))*1.
            IF (L .EQ. 1) THEN
               WRITE (*,'(A)') 
     >         'Dummy SUBLEVEL inserted for normalization of popnumbers'
               WRITE (0,'(A)') 
     >         'Dummy SUBLEVEL inserted for normalization of popnumbers'
            ENDIF
         ENDIF
         DO 36 J=1,NLOW     
   36    SUM=SUM+POPNUM(L,NSUBLOW(J))
         SUMINV = POPNUM(L,LOW) / SUM
         DO 37 J=1,NLOW
   37    POPNUM(L,NSUBLOW(J))=POPNUM(L,NSUBLOW(J)) * SUMINV
   30    CONTINUE
      ENDIF

C***  B. SPLITTING THE UPPER ENERGY LEVEL: BOLTZMANN (LTE)
      IF (NNUP .GT. 0) THEN
         WEIGHTSUM = WEIGHT(NSUBNUP(1))
          ENUPREST = WEIGHT(NSUBNUP(1)) * ELEVEL(NSUBNUP(1))         
         DO J=2, NNUP
            WEIGHTSUM = WEIGHTSUM + WEIGHT(NSUBNUP(J))
            ENUPREST = ENUPREST + WEIGHT(NSUBNUP(J))*ELEVEL(NSUBNUP(J))
         ENDDO
C***     Check if sum of sublevel weights conforms with weight of original level
         IF (NINT(WEIGHTSUM) .NE. NINT(WEIGHT(NUP))) THEN
            WRITE (*,'(A)') '*** WARNING: Inconsistency of stat. weights'
            WRITE (*,'(A, I3, A, I3)') 'UPPERLEVEL: ' // LEVEL(NUP)
     >           // '    WEIGHT: ', NINT(WEIGHT(NUP)),  
     >              '  =!= Sum of sublevel weights: ', NINT(WEIGHTSUM)
            WRITE (0,'(A)') 'WARNING: Inconsistency of WEIGHTS '
     >           // 'for level ' // LEVEL(NUP) // '  -- see output file'
            NWEIGHTNUPREST = NINT(WEIGHT(NUP) - WEIGHTSUM)
            ENUPREST = (ELEVEL(NUP)*WEIGHT(NUP)-ENUPREST)/NWEIGHTNUPREST
         ENDIF  
         DO 40 L=1,ND
         POPNUM(L,NSUBNUP(1))=1.
         DO 45 J=2,NNUP
   45    POPNUM(L,NSUBNUP(J))=EXP(C1*(ELEVEL(NSUBNUP(J-1))
     -         -ELEVEL(NSUBNUP(J)))/T(L))*WEIGHT(NSUBNUP(J))
     /         /WEIGHT(NSUBNUP(J-1))*POPNUM(L,NSUBNUP(J-1))
C***     NORMALIZATION
         SUM=0.
         IF (NWEIGHTNUPREST .GT. 0) THEN
            SUM=EXP(C1*(ELEVEL(NSUBNUP(1))-ENUPREST)/T(L))
     >           *NWEIGHTNUPREST/WEIGHT(NSUBNUP(1))*1.
            IF (L .EQ. 1) THEN 
               WRITE (*,'(A)') 
     >         'Dummy SUBLEVEL inserted for normalization of popnumbers'
               WRITE (0,'(A)') 
     >         'Dummy SUBLEVEL inserted for normalization of popnumbers'
            ENDIF
         ENDIF         
         DO 46 J=1,NNUP
   46    SUM=SUM+POPNUM(L,NSUBNUP(J))
         SUMINV = POPNUM(L,NUP) / SUM
         DO 47 J=1,NNUP
   47    POPNUM(L,NSUBNUP(J))=POPNUM(L,NSUBNUP(J)) * SUMINV
   40    CONTINUE
      ENDIF

      RETURN
      
C********** Error branches ***************************************************
      
      END
