      SUBROUTINE SEQUIN (NDIM,N,X,K,XNEW)
C***********************************************************************
C***  X(K) MUST BE A STRICKTLY MONOTONIC SEQUENCE
C***   THIS SUBROUTINE INSERTS THE POINT XNEW AT SUITABLE POSITION K
C***  NO INSERTION IN CASE OF EXACT COINCIDENCE (K:= NEGATIVE INDEX )
C***********************************************************************

      DIMENSION X(NDIM)
 
C***  INDEX ESTIMATE BY BISECTION
      NA=1
      A=X(1)
      NB=N
      B=X(N)
      IF (XNEW.EQ.A) GOTO 9
      IF (XNEW.EQ.B) GOTO 10
      IF((XNEW-A)*(XNEW-B).GE..0) GOTO 4
    1 IF(NB-NA .EQ. 1) GOTO 2
      NH=(NA+NB)/2
      H=X(NH)
      IF (XNEW.EQ.H) GOTO 11
      IF((XNEW-A)*(XNEW-H).GT..0) GOTO 3
      NB=NH
      B=H
      GOTO 1
    3 NA=NH
      A=H
      GOTO 1
 
C***  CASES OF EXACT COINCIDENCE
    9 K=-1
      RETURN
   10 K=-N
      RETURN
   11 K=-NH
      RETURN
 
      ENTRY SEQUINE  (NDIM,N,X,K,XNEW)
C***  ENTRY POINT IF THE INSERTION INDEX K IS GIVEN *********************
      IF (K.EQ.N+1) GOTO 8
      IF (K.LE.0 .OR. K.GT.N) THEN
            CALL REMARK ('ENTRY SEQUINE: WRONG INDEX GIVEN')
            STOP 'ERROR'
            ENDIF
      NB=K
 
C***  INSERTION OF XNEW AT INDEX NB
    2 IF(N+1 .GT. NDIM) THEN
            WRITE (0,*) 'ARRAY DIMENSION TOO SMALL'
            WRITE (0,*) 'NDIM=', NDIM
            STOP 'ERROR'
            ENDIF
      K=NB
    7 DO 5 J=K,N
      I=N+K-J
    5 X(I+1)=X(I)
    8 X(K)=XNEW
      N=N+1
      RETURN
 
C***  XNEW LIES OUTSIDE X(1) ... X(N)
    4 IF(N+1 .GT. NDIM) THEN
            WRITE (0,*) 'ARRAY DIMENSION TOO SMALL'
            WRITE (0,*) 'NDIM=', NDIM
            STOP 'FATAL ERROR IN SUBR. SEQUIN'
            ENDIF
      IF((XNEW-A)*(B-A).GT..0) GOTO 6
       K = 1
       GOTO 7
    6 N=N+1
      K=N
      X(N)=XNEW
      RETURN
      END
