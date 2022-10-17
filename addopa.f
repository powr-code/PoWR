      SUBROUTINE ADDOPA (ND, NDDIM, MAXLIN, MAXIND, LIND, LINDS, 
     >             XK, XKMID, XKRED, DELTAX, 
     >             PARALAS, LASER, 
     >             WS, ETAL, OPAL, ETA, ETANOTH, OPA, ETAK, ETAKNOTH, 
     >             OPAK, OPAKNOTH, THOMSON, 
     >             PWEIGHT, OPAFE, ETAFE, BFECHECK, BLASERL, NUP, LOW, 
     >             N, LEVEL)
C**********************************************************************
C***  TREATMENT OF LINE BLENDS:
C***  ADD OVERLAPPING OPACITIES AND EMISSIVITIES
C***  - ALSO GENERATES THE LASER WARNING (ONLY LASER-VERSION 1)
C***  - CALLED FROM: COLI, SETXJFINE ( <- with ND=1 )
C**********************************************************************

      DIMENSION LIND(MAXLIN),LINDS(MAXLIN),XKMID(MAXIND),XKRED(MAXIND)
      DIMENSION NUP(MAXLIN), LOW(MAXLIN)
      DIMENSION OPAL(NDDIM,MAXLIN), ETAL(NDDIM,MAXLIN)
      DIMENSION ETA(ND), ETANOTH(ND), OPA(ND), WS(ND)
      DIMENSION ETAK(NDDIM), ETAKNOTH(NDDIM)
      DIMENSION OPAK(NDDIM), PWEIGHT(MAXLIN)
      DIMENSION OPAKNOTH(NDDIM), THOMSON(NDDIM)
      DIMENSION OPAFE(ND), ETAFE(ND)

      CHARACTER(10), DIMENSION(N) :: LEVEL
      LOGICAL LASER, BFECHECK, BLASERL(MAXIND)

C***  WPI = SQRT(PI)
      DATA WPI /1.772454/

      LASER = .FALSE.

C***  REAL LINE OPACITIES FROM ALL INVENTED LINES
      DO L=1,ND
        ETAK(L)     = 0.0
        ETAKNOTH(L) = 0.0
        OPAK(L)     = 0.0
        OPAKNOTH(L) = 0.0
      ENDDO

      DO NL=1, MAXLIN
        IF (LIND(NL) .EQ. 0) CYCLE
        IF (XK .GT. XKRED(LINDS(NL))) CYCLE
        DK = (XK - XKMID(LINDS(NL))) * DELTAX
        PWEIGHT(NL) = EXP(-DK*DK)
        WS(NL) = WS(NL) + PWEIGHT(NL)
        PHI = PWEIGHT(NL) / WPI

        DO L=1, ND
           ETAK(L) = ETAK(L) + ETAL(L,NL)*PHI
           OPAK(L) = OPAK(L) + OPAL(L,NL)*PHI
           OPAG = (PARALAS-1.0) * OPA(L)
           BLASERL(LIND(NL)) = 
     >          BLASERL(LIND(NL)) .OR.  OPAL(L,NL) .LT. OPAG
        ENDDO
      ENDDO


C***  IRON: ADD IRON OPACITY AND EMISSIVITY IF NECESSARY
      IF (BFECHECK) THEN
         DO L=1, ND
            ETAK(L) = ETAK(L) + ETAFE(L)
            OPAK(L) = OPAK(L) + OPAFE(L)
         ENDDO
      ENDIF

C***  Laser treatment (all versions): 
C***       restrict Total Opacity to > PARALAS * CONTINUUM
        DO L=1, ND
          OPAG = (PARALAS-1.0) * OPA(L)
          LASER = LASER .OR. (OPAK(L).LT.OPAG)
          OPAK(L)=AMAX1(OPAG,OPAK(L))
        ENDDO

C***  Add the Continuum Values to ETAK, ETAKNOTH and OPAK
      DO L=1, ND
        ETAKNOTH(L) = ETANOTH(L) + ETAK(L)
        ETAK(L)     = ETA(L) + ETAK(L)
        OPAG        = PARALAS * OPA(L)
        OPAKNOTH(L) = OPA(L) * (1.-THOMSON(L)) + OPAK(L)
        OPAK(L)     = OPA(L) + OPAK(L)
      ENDDO      
      

      RETURN
      END
