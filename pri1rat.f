      SUBROUTINE PRI1RAT (N, LEVEL, NDIM, L, CRATE, RRATE, RATCO, EN,
     $           MODHEAD, JOBNUM, NFIRST, NLAST, NATOM, 
     $           NAUTO, RDIEL, RAUTO, IONGRND, KODRLOW, LASTKDR, PRILEV, 
     $           SGAIN, SLOSS)
C***********************************************************************
C***  This Subr. gives the leading rates that populate / depopulate a 
C***  selected level
C***  The output is written into Standard-out
C***  Called from: STEAL
C***********************************************************************

      DIMENSION EN(NDIM)
      DIMENSION NFIRST(NATOM),NLAST(NATOM)
      DIMENSION CRATE(NDIM,NDIM),RRATE(NDIM,NDIM),RATCO(NDIM,NDIM)
      DIMENSION RDIEL(NDIM),RAUTO(NDIM),IONGRND(NDIM)
      DIMENSION KODRLOW(LASTKDR)
      CHARACTER LEVEL(NDIM)*10, PRILEV*10
      CHARACTER MODHEAD*100
      CHARACTER DLOSS*6, DGAIN*6, PLOSS*3, PGAIN*3
      DIMENSION SGAIN(NDIM), SLOSS(NDIM)
 

      WRITE (*,14)  MODHEAD,JOBNUM
   14 FORMAT (//,  A  ,17X,'JOB NO.',I7,/)
      WRITE (*,'(A,I3,/)') 'Leading transition rates from/to level: '
     >       // PRILEV // ' at depth index ', L

      WRITE (*,'(A, 10X, A,/)')
     >  '------------- GAINS  ----------------',
     >  '------------- LOSSES ----------------'

C***  Find level index for requested level

      DO I=1, N
         LEVIND = I
         IF (LEVEL(I) .EQ. PRILEV) GOTO 15
      ENDDO
      WRITE (*,*) '*** Requested level not found: ', PRILEV 
      WRITE (*,*) '*** Non-fatel ERROR in SUBR. PRI1RAT' 
      RETURN

   15 CONTINUE

C***  Find the appropriate element
      DO 99 NA=1, NATOM
      NFIRNA=NFIRST(NA)
      NLANA=NLAST(NA)
      IF (LEVIND .LT. NFIRNA .OR. LEVIND .GT. NLANA) CYCLE

C*******************************************************************************
C***  ADD RADIATIVE AND COLLISIONAL RATE COEFF. INTO MATRIX RATCO
      DO 11 I=NFIRNA,NLANA
      DO 11 J=NFIRNA,NLANA
   11 RATCO(I,J)=CRATE(I,J)+RRATE(I,J)
 
C***  ADD ADDITIONAL D-R TERMS INTO RATE COEFFICIENT MATRIX RATCO
      DO 12 KDR=1,LASTKDR
      LOW=KODRLOW(KDR)
      IF ((LOW .LT. NFIRNA) .OR. (LOW .GT. NLANA)) GOTO 12
      NUP=IONGRND(LOW)
      RATCO(LOW,NUP)=RATCO(LOW,NUP)+RAUTO(LOW)
      RATCO(NUP,LOW)=RATCO(NUP,LOW)+RDIEL(LOW)
   12 CONTINUE

C*******************************************************************************
C***  ADDING (ABSOLUTE) RATES INTO MATRIX RATCO
      DO 20 I=NFIRNA,NLANA
      ENI=EN(I)
      DO 20 J=NFIRNA,NLANA
   20 RATCO(I,J)=ENI*RATCO(I,J)
 
C***  Store rates
      DO J=NFIRNA, NLANA
         SGAIN(J-NFIRNA+1) = RATCO(J,LEVIND)
         SLOSS(J-NFIRNA+1) = RATCO(LEVIND,J)
      ENDDO
      SGAIN(LEVIND-NFIRNA+1) = .0
      SLOSS(LEVIND-NFIRNA+1) = .0

C***  The five leading terms
      NLEV = NLANA - NFIRNA + 1
      DO I=1, 5
C***   Gains
         JGAIN = ISMAX (NLEV, SGAIN, 1)
         IF (JGAIN+NFIRNA-1 .GT. LEVIND) THEN
            DGAIN = 'UPPER'
         ELSE
            DGAIN = 'LOWER'
         ENDIF
         IF (RRATE(JGAIN+NFIRNA-1,LEVIND) .GT. 
     >       CRATE(JGAIN+NFIRNA-1,LEVIND)) THEN
            PGAIN = 'RAD'
         ELSE
            PGAIN = 'COL'
         ENDIF
C***   Losses
         JLOSS = ISMAX (NLEV, SLOSS, 1)
         IF (JLOSS+NFIRNA-1 .GT. LEVIND) THEN
            DLOSS = 'UPPER'
         ELSE
            DLOSS = 'LOWER'
         ENDIF
         IF (RRATE(LEVIND,JLOSS+NFIRNA-1) .GT. 
     >       CRATE(LEVIND,JLOSS+NFIRNA-1)) THEN
            PLOSS = 'RAD'
         ELSE
            PLOSS = 'COL'
         ENDIF
         WRITE (*,'(2(A,2X,A,2X,1PG12.3,2X,A,10X))') 
     >     LEVEL(JGAIN+NFIRNA-1), DGAIN, SGAIN(JGAIN), PGAIN,
     >     LEVEL(JLOSS+NFIRNA-1), DLOSS, SLOSS(JLOSS), PLOSS
         SGAIN(JGAIN) = .0         
         SLOSS(JLOSS) = .0         
      ENDDO


 
   99 CONTINUE
C***  ENDLOOP  =========================================================

      RETURN
      END
