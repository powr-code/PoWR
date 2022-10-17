      SUBROUTINE NG4 (ND, N, RNE, NCHARG, 
     >                POPNEW, POP, POP1, POP2, POP3, 
     >                NATOM, ABXYZ, NFIRST, NLAST,
     >                TNEW, T, TOLD1, TOLD2, TOLD3, NOTEMP, TRESH, NOUT)
C***********************************************************************
C***  NG'S Acceleration method (K.C.NG 1974, CHEM. PHYS. 61, 2680):
C***  Extrapolation by considering the last four iterations
C***  === Version ===: The acceleration is now the same as suggested 
C***                   by Ivan Hubeny :
C***                   a) All popnumbers at all depths are considered
C***                      to examine one acceleration factor. If Temperature 
C***                      corrections are taken into account, T is also
C***                      included
C***                   b) No logarithmic treatment
C***                   c) Popnumbers are always weigted 
C***********************************************************************
 
      DIMENSION POPNEW(ND,N), POP(ND,N), POP1(ND,N)
      DIMENSION POP2(ND,N), POP3(ND,N)
      DIMENSION RNE(ND),NCHARG(ND)
      DIMENSION TNEW(ND), T(ND), TOLD1(ND), TOLD2(ND), TOLD3(ND)
      DIMENSION ABXYZ(NATOM), NFIRST(NATOM), NLAST(NATOM)
      LOGICAL NOTEMP
      COMMON / COMNEGT / NEGT

C***********************************************************************
C***  SET MINIMUM TEMPERATURE FOR EXTRAPOLATION OF T(R)
      TMIN=3000.
C***********************************************************************

C***  INITIALIZE COUNTER FOR NEGATIVE TEMPERATURE WARNING
      NEGT=0

C***  Popnumbers smaller than TRESH are not accounted for
      TRESH = 1.E-15

C***  Do not EXTRAP the outer points
      N1 = 5
      IF (N1 .GT. ND/2) THEN
        N1 = ND/2
      ENDIF
      NOUT = N1

      DO J=1,N
        DO L=1,ND
          IF ( POP(L,J) .LT. 1.E-99 ) POP(L,J)=1.E-99
          IF ( POP1(L,J) .LT. 1.E-99 ) POP1(L,J)=1.E-99
          IF ( POP2(L,J) .LT. 1.E-99 ) POP2(L,J)=1.E-99
          IF ( POP3(L,J) .LT. 1.E-99 ) POP3(L,J)=1.E-99
        ENDDO
      ENDDO

***  LOOP FOR ALL LEVELS AND POPNUMBERS  ---------------------------------------------
      A1=0.0
      C1=0.0
      DO J=1,N
        DO L=1,ND
          IF (POP(L,J) .LT. TRESH .OR. L .LT. N1) CYCLE
          WT = 1. / POP(L,J)
          D0 = POP(L,J) - POP1(L,J)
          D1 = D0 - POP1(L,J) + POP2(L,J)
          D2 = D0 - POP2(L,J) + POP3(L,J)
          A1 = A1 + WT*D1*D1
          B1 = B1 + WT*D1*D2
          B2 = B2 + WT*D2*D2
          C1 = C1 + WT*D0*D1
          C2 = C2 + WT*D0*D2
        ENDDO
      ENDDO
C*** Temperature corrections
      DO L=N1, ND
        WT = 1. / T(L)
        D0 = T(L) - TOLD1(L)
        D1 = D0 - TOLD1(L) + TOLD2(L)
        D2 = D0 - TOLD2(L) + TOLD3(L)
        A1 = A1 + WT*D1*D1
        B1 = B1 + WT*D1*D2
        B2 = B2 + WT*D2*D2
        C1 = C1 + WT*D0*D1
        C2 = C2 + WT*D0*D2
      ENDDO

C***  ENDLOOP  ---------------------------------------------------------

      AB = B2*A1 - B1*B1
      IF (AB .NE. 0.) THEN
        A = (B2*C1 - B1*C2) / AB
        B = (A1*C2 - B1*C1) / AB
      ELSE
        WRITE (0,*) 'No EXTRAP (NG4) executed'
        WRITE (*,*) 'No EXTRAP (NG4) executed'
        A = 0.
        B = 0.
      ENDIF

      DO J=1,N
        DO L=1,ND
          IF (POP(L,J) .LT. TRESH .OR. L .LT. N1) THEN
            POPNEW(L,J) = POP(L,J)
          ELSE
            POPNEW(L,J) = (1.-A-B)*POP(L,J) + 
     >                    A*POP1(L,J) + B*POP2(L,J)
          ENDIF
        ENDDO
      ENDDO
C*** Temperature corrections
      IF (.NOT. NOTEMP) THEN
        DO L=1, ND
          IF (L .LT. N1) THEN
            TNEW(L) = T(L)
          ELSE
            TNEW(L) = (1.-A-B)*T(L) + 
     >                  A*TOLD1(L) + B*TOLD2(L)
          ENDIF
        ENDDO
      ELSE
        DO L=1, ND
          TNEW(L) = T(L)
        ENDDO
      ENDIF

C***  Popnumbers are scaled to ensure number conservation
C***  Electron density is updated
      DO L=1,ND
        RNE(L)=0.0
        DO NA=1,NATOM
          NFIRNA=NFIRST(NA)
          NLANA=NLAST(NA)
          POPSUM=0.0
          DO J=NFIRNA,NLANA
            POPSUM=POPSUM+POPNEW(L,J)
          ENDDO
          POPSUM=POPSUM/ABXYZ(NA)
          DO J=NFIRNA,NLANA
            POPNEW(L,J)=POPNEW(L,J)/POPSUM
            RNE(L)=RNE(L)+NCHARG(J)*POPNEW(L,J)
          ENDDO
        ENDDO
      ENDDO

C*** Check for very low temperatures
      IF (.NOT. NOTEMP) THEN
        DO L=1, ND
          IF (TNEW(L) .LT. TMIN) THEN
            TNEW(L) = TMIN
            NEGT = NEGT + 1
          ENDIF
        ENDDO
      ENDIF

      RETURN
      END
