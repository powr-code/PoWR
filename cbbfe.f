      SUBROUTINE CBBFE (OMEGA, NUP, LOW, TL, TROOT, WAVENUM, WN3, 
     >                  EINST, NDIM)
C***********************************************************************
C***  COLLISIONAL BOUND-BOUND TRANSITIONS OF IRON
C***  CALCULATION OF COLLISION STRENGTH OMEGA(UP-LOW)
C***  Nota added 25-Aug-2010:
C***  Since i line index is only assigned to supelines with 
C***  existing sublines, CBBFE is never called for "forbidden" superlines.
C***  We now insert cross-sections for such transitions in COLLI
C***  Calling tree: STEAL - COMA - COLLI - CBBFE
C***********************************************************************
 
      DIMENSION EINST(NDIM,NDIM)

      REAL, EXTERNAL :: EXPINT1EXP
      
C***  C1 = H * C / K    ( CM * KELVIN )
      DATA C1 / 1.4388 /

C***  FORMULA  OF VAN REGEMORTER (1962, APJ 136, 906)
C***           FORMULA WITH  G = 0.2
C***           GAMMA = MAX(G,.276*EXP(U0)*E1(U0))
C***            (E1: FIRST EXPONENTIAL INTEGRAL)

      U0 = C1 * WAVENUM / TL
      GAMMA = AMAX1 (0.2, 0.276*EXPINT1EXP(U0))
      OMEGA = 20.56 * GAMMA * EINST(NUP,LOW) / WN3 / TROOT

      RETURN

cccccccc stillgelegte Version von Goetz
c      SUBROUTINE CBBFE (OMEGA,IND,NUP,LOW,TL,TROOT,KEYCBB,WAVENUM,
c     $                  COCOFE,NCOMAX,NCO,WEIGHT,N,NDIM)
cC***********************************************************************
cC***  COLLISIONAL BOUND-BOUND TRANSITIONS OF HELIUM
cC***  CALCULATION OF OMEGA(UP-LOW)
cC***********************************************************************
c 
c      DIMENSION COCOFE(NCOMAX,IND),WEIGHT(NDIM)
c      CHARACTER*4 KEYCBB(IND)
c 
cC***  C1 = H * C / K    ( CM * KELVIN )
c      DATA C1 / 1.4388 /
cC***  CT = K / e (K^-1)
c      DATA CT / 8.6175E-5 /
c  
cC***  X  = LN(T) (eV)
c      X     = ALOG10(CT*TL)
cC***  CALCULATE GAMMA(X)
c      GAMMA = 0.
c      DO K=1, NCO
c         NC = NCO - K + 1
c         GAMMA = GAMMA*X + COCOFE(NC, IND)
c      ENDDO
c
c      GAMMA = 10.**GAMMA
c
cC***  CALCULATE OMEGA (DISSERTATION DREIZLER S.56)
c      OMEGA = 5.4655E-11*TROOT* GAMMA * WEIGHT(LOW)/WEIGHT(NUP)
c
c
c      RETURN
c      END
c
c
c


      END
