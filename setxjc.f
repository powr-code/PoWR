      SUBROUTINE SETXJC (XJCAPP,XJC,OPAC,ETAC,SCNEW,SCOLD,WCHARM,
     $                   XLAMBDA,L,NF,ND,TL,TLOLD,NOTEMP)
C***********************************************************************
C***  CALCULATE NEW CONT. SOURCE FUNCTION AND SCHARMER'S RADIATION FIELD
C***********************************************************************

      DIMENSION ETAC(NF),XJC(ND,NF),XJCAPP(NF)
      DIMENSION WCHARM(ND,NF),SCOLD(NF,ND)
      DIMENSION OPAC(NF),SCNEW(NF),XLAMBDA(NF)
      LOGICAL NOTEMP

C***  LOOP OVER ALL CONT. FREQUENCIES  ----------------------------------
      DO 1 K=1,NF
C***  PREVENT DIVIDE CHECK ERRORS
      IF (OPAC(K) .GT. .0) THEN
         SCNEW(K)=ETAC(K)/OPAC(K)
         ELSE
         SCNEW(K)=.0
         ENDIF
      XJCAPP(K)=XJC(L,K)+WCHARM(L,K)*(SCNEW(K)-SCOLD(K,L))
      IF (XJCAPP(K) .LT. 0.) THEN
         WCHARM(L,K) = 0.
c         WRITE(0,'(A,I6)') 'SETXJC: XJCAPP < 0, at K=', K
         XJCAPP(K) = XJC(L,K)
      ENDIF
   1  CONTINUE

C***  ONLY IF TEMPERATURE CORRECTIONS ARE APPLIED:
C***  ADD SPECIAL TERM AT INNER BOUNDARY, WHICH ACCOUNTS FOR THE
C***    DIRECT TEMPERATURE-DEPENDENCE OF THE RADIATION FIELD
C***    VIA THE BOUNDARY CONDITION 

      IF (.NOT. NOTEMP  .AND.  L .EQ. ND) THEN
      DO 2 K=1,NF
      XLAM=XLAMBDA(K)
      XJCAPP(K)=XJCAPP(K) +
     +     0.5*(BNUE(XLAM,TL)-BNUE(XLAM,TLOLD))
    2 CONTINUE
      ENDIF

      RETURN
      END