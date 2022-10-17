      SUBROUTINE FREQUNORM (ND, OPASMEAN, OPASMEANTC, SMEAN, QFJMEAN, 
     >                  XJTOTL, OPAJMEAN, OPAJMEANTC, OPAPMEAN, 
     >                  QOPAHMEAN, HMEAN, EDDIHOUTJMEAN, 
     >                  RADIUS, RSTAR, DENSCON, 
     >                  FTCOLI, WJC, WJC_MIN, 
     >                  FWTEST, NF, OPC, FTFE, LASTFE, 
     >                  LPLOT_WCHARM, XLAMBDA, OPAROSS, OPALAMBDAMEAN, 
     >                  ARAD, ACONT, ATHOM, ENTOT, ABXYZ, ATMASS, NATOM,
     >                  T, GAMMAT, UNLU_TAUMAX, UNLU_TAUMAX2, TAUROSS)
C*************************************************************
C***  Some quantities which were intergrated in FREQUINT are
C***    now normalized etc. 
C*************************************************************

      IMPLICIT NONE

      INTEGER, INTENT(IN) :: ND, NF, NATOM, LASTFE

      REAL, DIMENSION(NATOM), INTENT(IN) :: ABXYZ, ATMASS
      REAL, DIMENSION(ND), INTENT(IN) :: ENTOT
      REAL, DIMENSION(ND), INTENT(OUT) :: TAUROSS
      REAL, DIMENSION(ND), INTENT(INOUT) :: ARAD, ACONT, ATHOM

      REAL, DIMENSION(ND) :: OPASMEAN, SMEAN, QFJMEAN,
     >                       XJTOTL, OPAJMEAN, QOPAHMEAN,
     >                       OPASMEANTC, OPAJMEANTC, OPAPMEAN,
     >                       OPAROSS, OPALAMBDAMEAN, T,
     >                       FTCOLI, RADIUS, HMEAN
      CHARACTER(8) :: OPC
      REAL, DIMENSION(ND,NF) :: WJC, WJC_MIN
      REAL, DIMENSION(ND,LASTFE) :: FTFE
      REAL, DIMENSION(NF) :: FWTEST, XLAMBDA
      
      INTEGER :: K, KK, L, NA, INDFE, LPLOT_WCHARM
      REAL :: ATMEAN, UNLU_TAUMAX, UNLU_TAUMAX2, TACCELERATE,
     >        TACCELERATE2, SDURCHB, FDAMP, ETAUROSS, GAMMAT,
     >        RSTAR, RINT, RL2, EWTOT, CC, EDDIHOUTJMEAN, FCOR,
     >        TACCELERATE3, BTOTL

C***  tiefenabh. clumping nach goetz
      REAL, DIMENSION(ND) :: DENSCON

C***  Maximum Value for Scharmer Weight: 1 - 1.E-10
C***  corresponding to a maximum optical depth of ca. tau = 1E10 
      REAL, PARAMETER :: WJCMAX = 0.9999999999 
C***  STEBOL = STEFAN-BOLTZMANN CONSTANT / PI (ERG/CM**2/SEC/STERAD/KELVIN**4
      REAL, PARAMETER :: STEBOL = 1.8046E-5
      REAL, PARAMETER :: AMU = 1.6605E-24     !ATOMIC MASS UNIT (in g)      
      REAL, PARAMETER :: PI4 = 12.5663706     !4 * Pi
      REAL, PARAMETER :: CLIGHT = 2.99792458E10   !C IN CGS UNITS


      ATMEAN = 0.
      DO NA=1, NATOM
          ATMEAN = ATMEAN + ABXYZ(NA) * ATMASS(NA)
      ENDDO

      DO L=1, ND
        RL2      = RADIUS(L) * RADIUS(L)
        BTOTL    = STEBOL * T(L) * T(L) * T(L) * T(L)
        OPASMEAN  (L) = DENSCON(L) * OPASMEAN  (L) / (SMEAN(L) * RSTAR)
        OPASMEANTC(L) = DENSCON(L) * OPASMEANTC(L) / (SMEAN(L) * RSTAR)
        OPALAMBDAMEAN(L) = DENSCON(L)*OPALAMBDAMEAN(L)/(SMEAN(L)* RSTAR)
        QFJMEAN (L) =           QFJMEAN(L)  /  XJTOTL(L)
        OPAJMEAN  (L) = DENSCON(L) * OPAJMEAN  (L) / (XJTOTL(L) * RSTAR)
        OPAJMEANTC(L) = DENSCON(L) * OPAJMEANTC(L) / (XJTOTL(L) * RSTAR)
        OPAPMEAN  (L) = DENSCON(L) * OPAPMEAN  (L) / (BTOTL * RSTAR)
        FTCOLI    (L) = DENSCON(L) * FTCOLI(L)
        DO INDFE=1, LASTFE
           FTFE(L,INDFE) = DENSCON(L) * FTFE(L,INDFE)
        ENDDO

        OPAROSS(L) = 4.*STEBOL*T(L)*T(L)*T(L) / OPAROSS(L) 

        IF (L == ND) CYCLE

C***    RADIATIVE ACCELERATION IN CGS UNITS
        EWTOT = 0.5*(ENTOT(L)+ENTOT(L+1)) * AMU *ATMEAN
        RINT  = 0.5*(RADIUS(L)+RADIUS(L+1))
        ARAD(L)  = PI4 * ARAD(L)  / (RSTAR*CLIGHT*RINT*RINT*EWTOT)
        ACONT(L) = PI4 * ACONT(L) / (RSTAR*CLIGHT*RINT*RINT*EWTOT)
        ATHOM(L) = PI4 * ATHOM(L) / (RSTAR*CLIGHT*RINT*RINT*EWTOT)

C***    Note: QOPAHMEAN is with the average density (not: clump), because it 
C***          enters a transfer eq. 
C***          safety check added, because floating exceptions have been
C***          encountered occasionally -  19-Oct-2015
        IF (HMEAN(L) .NE. 0.) THEN
           QOPAHMEAN(L) = QOPAHMEAN(L) / HMEAN(L)
        ELSE
           QOPAHMEAN(L) = 0.
        ENDIF

        IF (OPC .EQ. 'DIAGMIN') THEN
          DO K=1, NF
            WJC(L,K) = WJC_MIN(L,K)
          ENDDO
        ELSE
C***  Division by FWTEST now before transformation into TAU
          DO K=1, NF
            IF (FWTEST(K) .NE. 0.) THEN
              WJC(L,K) = WJC(L,K) / FWTEST(K)
            ENDIF
          ENDDO

C***    In the DIAGTAU version, the Scharmer factor was converted into a TAU 
C***       before averaging over the fine grid in FREQUINT. This is now 
C***       transformed back. 
          IF (OPC .EQ. 'DIAGTAU') THEN
            DO K=1, NF
              IF (WJC(L,K) .GE. 0.) THEN
                WJC(L,K) = 1. - EXP(-WJC(L,K))
              ELSE
                WJC(L,K) = 0.
              ENDIF
            ENDDO
          ENDIF
        ENDIF

        IF (L .EQ. LPLOT_WCHARM) THEN
          WRITE (105,'(A)') 'FINISH'
          WRITE (105,'(A)') '* Coarse Information'
          WRITE (105,'(A)') 'N=?'
          WRITE (105,'(2(E20.10,1X))')
     >      (XLAMBDA(KK), WJC(LPLOT_WCHARM,KK), KK=1, NF)
          WRITE (105,'(A)') 'FINISH'
        ENDIF

C***    In any case, the Scharmer weights are checked for valid range
        DO K=1, NF
           WJC(L,K) = AMAX1 (0.    , WJC(L,K))
           WJC(L,K) = AMIN1 (WJCMAX, WJC(L,K))
        ENDDO
 
      ENDDO

      EDDIHOUTJMEAN = EDDIHOUTJMEAN / XJTOTL(1)

      !Calculate TAUROSS now (even if not used and only stored in MODEL)
      TAUROSS(1) = 0.      
      DO L=2, ND
        TAUROSS(L) = TAUROSS(L-1) + 
     >                 0.5 * (OPAROSS(L-1) + OPAROSS(L)) *
     >                 (RADIUS(L-1) - RADIUS(L))
      ENDDO
            
      IF (GAMMAT > .0) THEN
C***  Damping of the Temperature-ALO with high Rosseland optical depth
C***    and test output > taccelerate  

C***  CPR-File
        WRITE (0, '(3(A,F6.1))') '* UNLU-PARAMETER: GAMMAT=', GAMMAT, 
     >   '  TAUMAX=', UNLU_TAUMAX, '  TAUMAX2=', UNLU_TAUMAX2

        OPEN (26, FILE='taccelerate')
        WRITE (26, '(3(A,F6.1))') '* UNLU-PARAMETER: GAMMAT=', GAMMAT, 
     >   '  TAUMAX=', UNLU_TAUMAX, '  TAUMAX2=', UNLU_TAUMAX2
        WRITE (26,'(A)') '* L    Acceler.   Accel.(damped)  ' //
     >         ' Accel.(streng.)  Tau-Ross     S/B'
        WRITE (26,'(A)') '*--------------------------------------------'
        DO L=1, ND
          TACCELERATE = OPASMEANTC(L)/(OPASMEANTC(L)-OPALAMBDAMEAN(L))
C***      Damp the amplification by exp(-tau_rosseland)
          IF (L .GT. 1) THEN
            IF (TAUROSS(L) >= 0.) THEN
              ETAUROSS = EXP(-TAUROSS(L))
            ELSE
              ETAUROSS = 1.
            ENDIF
            FDAMP = OPASMEANTC(L) * ETAUROSS /
     >           (OPASMEANTC(L) - (1.-ETAUROSS)*OPALAMBDAMEAN(L))
            OPALAMBDAMEAN(L) = FDAMP * OPALAMBDAMEAN(L)
          ENDIF
          TACCELERATE2 = OPASMEANTC(L)/(OPASMEANTC(L)-OPALAMBDAMEAN(L))

C***      Strengthen od damp corrections with LTE departure??
          SDURCHB = SMEAN(L) / (STEBOL * T(L)*T(L)*T(L)*T(L)) 
ccc       Version 1: STRENGTHEN
ccc         fstreng = 1. / amin1(1., sdurchb)
ccc ... replaced by another test version: not strengthen, but DAMP with 
cc      departure or it's inverse  ... wrh  5-Sep-2002 10:26:59
C***      NOTE THAT FCOR PLAYS THE SAME ROLE AS ETAUROSS DOES ABOVE
          IF (L .GT. 1) THEN
            FCOR = AMIN1 (SDURCHB, 1./SDURCHB)
            FDAMP = OPASMEANTC(L) * FCOR /
     >           (OPASMEANTC(L) - (1.-FCOR)*OPALAMBDAMEAN(L))
            OPALAMBDAMEAN(L) = FDAMP * OPALAMBDAMEAN(L)
          ENDIF 

          TACCELERATE3 = OPASMEANTC(L)/(OPASMEANTC(L)-OPALAMBDAMEAN(L))

          WRITE (26, '(I2,2X,5G14.4)')
     >         L, TACCELERATE, TACCELERATE2, TACCELERATE3, 
     >          TAUROSS(L), SDURCHB 
        ENDDO
      ENDIF

      RETURN
      END
