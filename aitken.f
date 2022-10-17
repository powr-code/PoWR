      SUBROUTINE AITKEN (ND,N,RNE,NCHARG,POPNEW,POPNUM,POP1,POP2,NATOM,
     $                   ABXYZ,NFIRST,NLAST,TNEW,T,TOLD1,TOLD2,NOTEMP)
C***********************************************************************
C***  LOGARITHMIC AITKEN EXTRAPOLATION OF POPULATION NUMBERS
C***  THE ELECTRON DENSITY IS UPDATED ACCORDING TO THE NEW POPNUMBERS
C***********************************************************************
 
      DIMENSION RNE(ND),NCHARG(N)
      DIMENSION POPNEW(ND,N),POPNUM(ND,N),POP1(ND,N),POP2(ND,N)
      DIMENSION ABXYZ(NATOM),NFIRST(NATOM),NLAST(NATOM)
      DIMENSION TNEW(ND),T(ND),TOLD1(ND),TOLD2(ND)
      LOGICAL NOTEMP
      COMMON / COMNEGT / NEGT

C**********************************************************************
C***  SET MINIMUM TEMPERATURE FOR EXTRAPOLATION
      TMIN = 3000.
C**********************************************************************
 
C***  LOOP OVER ALL DEPTH POINTS  --------------------------------------
      DO 1 L=1,ND
      RNE(L)=.0
 
C***  LOOP FOR EACH ELEMENT  -------------------------------
      DO 1 NA=1,NATOM
      NFIRNA=NFIRST(NA)
      NLANA=NLAST(NA)
      POPSUM=.0
      DO 2 J=NFIRNA,NLANA
C*** ZERO POPULATION NUMBERS ARE KEPT AT ZERO
      IF ( POPNUM(L,J).EQ.0. .OR. POP1(L,J).EQ.0. .OR. POP2(L,J)
     > .EQ.0. ) THEN
         POPNEW(L,J) = POPNUM(L,J)
         GOTO 8
         ENDIF
      A1   =ALOG10(POPNUM(L,J)/POP1(L,J))
C     THE AMPLIFICATION OF THE CORRECTIONS IS LIMITED TO A FACTOR OF 20.
      POPQ=POP1(L,J)/POP2(L,J)
      IF (POPQ .EQ. 1.) THEN 
         Q = 0.95
         ELSE
         A0   = ALOG10 (POPQ)
         Q = AMIN1(0.95,A1/A0 )
         ENDIF
      COR = A0/(1.-Q)
C***  THE CORRECTION ITSELF IS LIMITED TO A FACTOR OF 10.
      COR = AMIN1 (COR, 1. )
      COR = AMAX1 (COR,-1. )
      POPNEW(L,J)=POP2(L,J) * 10.**COR
    8 CONTINUE
      POPSUM=POPSUM+POPNEW(L,J)
    2 CONTINUE
 
C***  POPNUMBERS ARE SCALED TO ENSURE NUMBER CONSERVATION
C***  ELECTRON DENSITY IS UPDATED
      POPSUM=POPSUM/ABXYZ(NA)
      DO 3 J=NFIRNA,NLANA
      POPNEW(L,J)=POPNEW(L,J)/POPSUM
      RNE(L)=RNE(L)+NCHARG(J)*POPNEW(L,J)
    3 CONTINUE
 
C***  TEMPERATURE EXTRAPOLATION  -------------------------------
      IF ( .NOT. NOTEMP ) THEN
      IF ( T(L).EQ.0. .OR. TOLD1(L).EQ.0. .OR. TOLD2(L)
     > .EQ.0. ) THEN
         TNEWL = T(L)
         GOTO 18
         ENDIF
      A1   =ALOG10(T(L)/TOLD1(L))
C     THE AMPLIFICATION OF THE CORRECTIONS IS LIMITED TO A FACTOR OF 20.
      POPQ=TOLD1(L)/TOLD2(L)
      IF (POPQ .EQ. 1.) THEN 
         Q = 0.95
         ELSE
         A0   = ALOG10(POPQ)
         Q = AMIN1(0.95,A1/A0 )
         ENDIF
      COR = A0/(1.-Q)
C***  THE CORRECTION ITSELF IS LIMITED TO A FACTOR OF 10.
      COR = AMIN1 (COR, 1. )
      COR = AMAX1 (COR,-1. )
      TNEWL =TOLD2(L) * 10.**COR
C***  SET MINIMUM TEMPERATURE
      IF (TNEWL .LT. TMIN) THEN
         TNEW(L) = TMIN
         NEGT = NEGT + 1
         ELSE
         TNEW(L) = TNEWL
         ENDIF

   18 CONTINUE
      ENDIF
    1 CONTINUE
C***  ENDLOOPS  --------------------------------------------------------
 
      RETURN
      END
