      SUBROUTINE COLIWM(Z, P, R, ND, NP, 
     >                  CWM0, CWM1, CWM2, CWM3, 
     >                  CWM1O, CWM1I, CWM3O, CWM3I, 
     >                  BSHORT_CHAR)

      DIMENSION Z(ND,NP), P(NP), R(ND)
      DIMENSION CWM0(ND,NP), CWM2(ND,NP)
      DIMENSION CWM1(ND,NP), CWM3(ND,NP)
      DIMENSION CWM1O(NP), CWM1I(NP), CWM3O(NP), CWM3I(NP)
      LOGICAL BNONCORE, BSHORT_CHAR

C***  Initialisation
      DO JP=1, NP
         LMAX=MIN0(NP+1-JP,ND+1)
         CWM1O(JP) = 0.
         CWM3O(JP) = 0.
         CWM1I(JP) = 0.
         CWM3I(JP) = 0.
         DO L=1, LMAX-1
            CWM0(L,JP) = 0.
            CWM2(L,JP) = 0.
         ENDDO
      ENDDO
      DO JP=1, NP-1
         LMAX = MIN0(NP+1-JP,ND+1)
         DO L=1, LMAX-1
            CWM1(L,JP) = 0.
            CWM3(L,JP) = 0.
         ENDDO
      ENDDO

C***  0. and 2. Moment
      DO JP=1, NP-1
         LMAX=MIN0(NP+1-JP,ND+1)
         DO L=1, LMAX-1
            A = Z(L,JP)
            AA = A*A
            AAA= A*AA
            B = Z(L,JP+1)
            BB = B*B
            BBB= B*BB
C***       Factors for Radius L
            WF0 = R(L)
            WF2 = 1./R(L)
C***       Weights for Interval [JP,JP+1]
            W0A = WF0/2.*(A-B)
            W0B = W0A
            W2A = WF2/12.*(A-B)*(3.*AA + 2.*A*B + BB)
            W2B = WF2/12.*(A-B)*(AA + 2.*A*B + 3.*BB)
C***       Second Order at the inner Boundary (non-Core)
            IF (L .EQ. (LMAX-1).AND.(JP.GT.(NP-ND))) THEN
               W0A = WF0/3.*A
               W0B = 2.*WF0/3.*A
               W2A = WF2/5.*AAA
               W2B = 2.*WF2/15.*AAA
            ENDIF
C***       Add up Weights
            CWM0(L,JP)  = CWM0(L,JP)  + W0A
            CWM0(L,JP+1)= CWM0(L,JP+1)+ W0B
            CWM2(L,JP)  = CWM2(L,JP)  + W2A
            CWM2(L,JP+1)= CWM2(L,JP+1)+ W2B
         ENDDO
      ENDDO

C***  1. und 3. Moment, COLIRAY-Version
      IF (.NOT. BSHORT_CHAR) THEN
      DO JP=1, NP-1
         LMAX = MIN0(NP+1-JP,ND+1)
         DO L=1, LMAX-2
            BNONCORE = .NOT. ((L .LT. (LMAX-2)).OR.(JP .LE. NP-ND))
            RL  = 0.5*(R(L)+R(L+1))
            RL2 = RL*RL
            A   = SQRT(RL2 - P(JP)*P(JP))
            IF (.NOT. BNONCORE) THEN
               B = SQRT(RL2 - P(JP+1)*P(JP+1))
            ELSE
C***          at the inner Boundary (non-core)
               B = 0.
            ENDIF
            AA  = A*A
            BB  = B*B
            AAA = A*AA
            BBB = B*BB
C***       Factors for Radius L
            WF1 = 1.
            WF3 = 1./RL2
C***       Weights for Interval [JP, JP+1]
            W1A = WF1/6.*(A-B)*(2.*A + B)
            W1B = WF1/6.*(A-B)*(A + 2.*B)
            W3A = WF3/20.*(A-B)*(4.*AAA + 3.*AA*B + 2.*A*BB + BBB)
            W3B = WF3/20.*(A-B)*(AAA + 2.*AA*B + 3.*A*BB + 4.*BBB)
C***       Second Order at the inner Boundary (non-Core)
            IF (BNONCORE) THEN
               W1A = WF1/3.*AA
               W1B = 0.
               W3A = WF3/5.*A*AAA
               W3B = 0.
            ENDIF
C***       Add up Weights
            CWM1(L,JP)  = CWM1(L,JP)  + W1A
            CWM3(L,JP)  = CWM3(L,JP)  + W3A
C***       no JP+1 at the inner Boundary
            IF (.NOT. BNONCORE) THEN
C***       Inner Boundary(non-core) with V=0.
            CWM1(L,JP+1)= CWM1(L,JP+1)+ W1B
            CWM3(L,JP+1)= CWM3(L,JP+1)+ W3B
            ENDIF
         ENDDO
      ENDDO

C***  1. und 3. moment at the outer boundary
      DO JP=1, NP-1
        A  = Z(1,JP)
        B  = Z(1,JP+1)
        AA = A*A
        BB = B*B
        AAA= A*AA
        BBB= B*BB
C***   Factors for Radius L=1
        WF1 = 1.
        WF3 = 1./R(1)/R(1)
C***   Weights for Interval [JP, JP+1]
        W1A = WF1/6.*(A-B)*(2.*A + B)
        W1B = WF1/6.*(A-B)*(A + 2.*B)
        W3A = WF3/20.*(A-B)*(4.*AAA + 3.*AA*B + 2.*A*BB + BBB)
        W3B = WF3/20.*(A-B)*(AAA + 2.*AA*B + 3.*A*BB + 4.*BBB)
C***   Second Order at the inner Boundary
        IF (JP .EQ. NP-1) THEN
           W1A = WF1/3.*AA
           W1B = 0.
           W3A = WF3/5.*A*AAA
           W3B = 0.
        ENDIF
C***   Add up Weights
        CWM1O(JP)  = CWM1O(JP)  + W1A
        CWM1O(JP+1)= CWM1O(JP+1)+ W1B
        CWM3O(JP)  = CWM3O(JP)  + W3A
        CWM3O(JP+1)= CWM3O(JP+1)+ W3B
      ENDDO

C***  1. und 3. moment at the inner boundary
      DO JP=1, NP-ND
        A  = Z(ND,JP)
        B  = Z(ND,JP+1)
        AA = A*A
        BB = B*B
        AAA= A*AA
        BBB= B*BB
C***   Factors for Radius L=1 are 1.
C***   Weights for Interval [JP, JP+1]
        W1A = 1./6.*(A-B)*(2.*A + B)
        W1B = 1./6.*(A-B)*(A + 2.*B)
        W3A = 1./20.*(A-B)*(4.*AAA + 3.*AA*B + 2.*A*BB + BBB)
        W3B = 1./20.*(A-B)*(AAA + 2.*AA*B + 3.*A*BB + 4.*BBB)
C***   Second Order at the inner Boundary
        IF (JP .EQ. NP-ND) THEN
           W1A = 1./3.*AA
           W1B = 1./6.*AA
           W3A = 1./5.*A*AAA
           W3B = 1./20.*A*AAA
        ENDIF
C***   Add up Weights
        CWM1I(JP)  = CWM1I(JP)  + W1A
        CWM1I(JP+1)= CWM1I(JP+1)+ W1B
        CWM3I(JP)  = CWM3I(JP)  + W3A
        CWM3I(JP+1)= CWM3I(JP+1)+ W3B
      ENDDO
      ENDIF
C***  End of COLIRAY-Version


C***  Short-Characteristics 1. and 3. Moment
      IF (BSHORT_CHAR) THEN
      DO L=1, ND
        DO JP=1, NP-L
          A  = Z(L,JP)
          B  = Z(L,JP+1)
          AA = A*A
          BB = B*B
          AAA= A*AA
          BBB= B*BB
C***   Factors for Radius L
          WF1 = 1.
          WF3 = 1./(R(L)*R(L))
C***   Weights for Interval [JP, JP+1]
          W1A = WF1/6.*(A-B)*(2.*A + B)
          W1B = WF1/6.*(A-B)*(A + 2.*B)
          W3A = WF3/20.*(A-B)*(4.*AAA + 3.*AA*B + 2.*A*BB + BBB)
          W3B = WF3/20.*(A-B)*(AAA + 2.*AA*B + 3.*A*BB + 4.*BBB)
C***   Second Order at the inner Boundary
          IF (JP .EQ. NP-L) THEN
             W1A = WF1/3.*AA
             W1B = WF1/6.*AA          
             W3A = WF3/5.*A*AAA
             W3B = WF3/20.*A*AAA      
          ENDIF
C***   Add up Weights
          CWM1(L,JP)  = CWM1(L,JP)  + W1A
          CWM1(L,JP+1)= CWM1(L,JP+1)+ W1B
          CWM3(L,JP)  = CWM3(L,JP)  + W3A
          CWM3(L,JP+1)= CWM3(L,JP+1)+ W3B
        ENDDO
      ENDDO
      ENDIF

      RETURN
      END
