      SUBROUTINE MODELCO(ND, RADIUS, NP, P, Z, 
     >                   VELO, GRADI, NF, XLAMBDA, 
     >                   RSTAR, VDOP, NMOD, NDDIM, NPDIM, NFDIM, 
     >                   DENSCON)
C*******************************************************************
C***  Comparision of Geometric Mesh of the two Models
C***  The following quantities have to be equal for all models :
C***    ND, NP, NF, XLAMBDA, RSTAR, VDOP, BLFERR, DENSCON
C***  The relative deviation (RELDIFF) is calculated for :
C***    RADIUS, P, Z, VELO, GRADI. 
C*******************************************************************

      DIMENSION ND(NMOD), RADIUS(NDDIM, NMOD), NP(NMOD)
      DIMENSION P(NPDIM, NMOD)
      DIMENSION Z(NDDIM, NPDIM, NMOD)
      DIMENSION VELO(NDDIM,NMOD), GRADI(NDDIM,NMOD)
      DIMENSION NF(NMOD), XLAMBDA(NFDIM,NMOD), RSTAR(NMOD)
      DIMENSION VDOP(NMOD), DENSCON(ND(NMOD),NMOD)

      LOGICAL BP

      WRITE (0,*)
      WRITE (0,*) '------ COMPARISION OF THE TWO MODELS ------'

C***  ---- ND ----
      BP = .FALSE.
      DO I=1, NMOD-1
        DO J=1, NMOD
          IF (ND(I) .NE. ND(J)) BP = .TRUE.
        ENDDO
      ENDDO
      IF (BP) THEN
        WRITE (0,'(A,4(I2,1X))') 'ND : ',(ND(I),I=1, NMOD)
        STOP 'FATAL ERROR IN MODELCO'
      ENDIF

C***  ---- NP ----
      BP = .FALSE.
      DO I=1, NMOD-1
        DO J=1, NMOD
          IF (NP(I) .NE. NP(J)) BP = .TRUE.
        ENDDO
      ENDDO
      IF (BP) THEN
        WRITE (0,'(A,4(I2,1X))') 'NP : ',(NP(I),I=1, NMOD)
        STOP 'FATAL ERROR IN MODELCO'
      ENDIF

C***  ---- NF ----
      BP = .FALSE.
      DO I=1, NMOD-1
        DO J=1, NMOD
          IF (NF(I) .NE. NF(J)) BP = .TRUE.
        ENDDO
      ENDDO
      IF (BP) THEN
        WRITE (0,'(A,4(I2,1X))') 'NF : ',(NF(I),I=1, NMOD)
        STOP 'FATAL ERROR IN MODELCO'
      ENDIF

C***  ---- RADIUS ----
      RELDIFF = 0.
      DO L=1, ND(1)
        BP = .FALSE.
        DO I=1, NMOD-1
          DO J=1, NMOD
            IF (RADIUS(L,I) .NE. RADIUS(L,J)) BP = .TRUE.
            T = ABS(RADIUS(L,I) - RADIUS(L,J)) / RADIUS(L,J)
            IF (T .GT. RELDIFF) RELDIFF = T
          ENDDO
        ENDDO
        IF (.FALSE.) THEN
          WRITE (0,'(A,I2,1X,5(F20.10,1X))') 'DEPTH, RADIUS : ',
     >                             L,(RADIUS(L,I),I=1, NMOD)
        ENDIF
      ENDDO
      WRITE (0,'(A,F20.10)') 'RADIUS : RELDIFF = ',RELDIFF

C***  ---- P ----
      RELDIFF = 0.
      DO IP=1, NP(1)
        BP = .FALSE.
        DO I=1, NMOD-1
          DO J=1, NMOD
            IF (P(IP,I) .NE. P(IP,J)) BP = .TRUE.
            IF (P(IP,J) .NE. 0.) T = ABS(P(IP,I) - P(IP,J)) / P(IP,J)
            IF (T .GT. RELDIFF) RELDIFF = T
          ENDDO
        ENDDO
        IF (.FALSE.) THEN
          WRITE (0,'(A,I2,1X,5(F20.10,1X))') 'IMPACT, P : ',
     >                             IP,(P(IP,I),I=1, NMOD)
        ENDIF
      ENDDO
      WRITE (0,'(A,F20.10)') 'P      : RELDIFF = ',RELDIFF

C***  ---- Z ----
      RELDIFF = 0.
      DO L=1, ND(1)
        DO IP=1, NP(1)
          BP = .FALSE.
          DO I=1, NMOD-1
            DO J=1, NMOD
              IF (Z(L,IP,I) .NE. Z(L,IP,J)) BP = .TRUE.
              IF (Z(L,IP,J) .NE. 0.) 
     >          T = ABS(Z(L,IP,I) - Z(L,IP,J)) / Z(L,IP,J)
              IF (T .GT. RELDIFF) RELDIFF = T
            ENDDO
          ENDDO
          IF (.FALSE.) THEN
            WRITE (0,'(A,2(I2,1X),5(F20.10,1X))') 
     >            'DEPTH, IMPACT, Z : ',
     >             L, IP,(Z(L,IP,I),I=1, NMOD)
          ENDIF
        ENDDO
      ENDDO
      WRITE (0,'(A,F20.10)') 'Z      : RELDIFF = ',RELDIFF

C***  ---- VELO ----
      RELDIFF = 0.
      DO L=1, ND(1)
        BP = .FALSE.
        DO I=1, NMOD-1
          DO J=1, NMOD
            IF (VELO(L,I) .NE. VELO(L,J)) BP = .TRUE.
            IF (VELO(L,I) .NE. 0.) 
     >        T = ABS(VELO(L,I) - VELO(L,J)) / VELO(L,J)
            IF (T .GT. RELDIFF) RELDIFF = T
          ENDDO
        ENDDO
        IF (.FALSE.) THEN
          WRITE (0,'(A,I2,1X,5(F20.10,1X))') 'DEPTH, VELO : ',
     >                             L,(VELO(L,I),I=1, NMOD)
        ENDIF
      ENDDO
      WRITE (0,'(A,F20.10)') 'VELO   : RELDIFF = ',RELDIFF

C***  ---- GRADI ----
      RELDIFF = 0.
      DO L=1, ND(1)
        BP = .FALSE.
        DO I=1, NMOD-1
          DO J=1, NMOD
            IF (GRADI(L,I) .NE. GRADI(L,J)) BP = .TRUE.
            IF (GRADI(L,I) .NE. 0.) 
     >        T = ABS(GRADI(L,I) - GRADI(L,J)) / GRADI(L,J)
            IF (T .GT. RELDIFF) RELDIFF = T
          ENDDO
        ENDDO
        IF (.FALSE.) THEN
          WRITE (0,'(A,I2,1X,5(F20.10,1X))') 'DEPTH, GRADI : ',
     >                             L,(GRADI(L,I),I=1, NMOD)
        ENDIF
      ENDDO
      WRITE (0,'(A,F20.10)') 'GRADI  : RELDIFF = ',RELDIFF

C***  ---- DENSCON ----
C***  Aufpassen beim Modellvergleich, ist das so richtig???
      BP = .FALSE.
      DO L=1,ND(1)
         DO I=1, NMOD-1
            DO J=1, NMOD
               IF (DENSCON(L,I) .NE. DENSCON(L,J)) BP = .TRUE.
            ENDDO
         ENDDO
      ENDDO
      IF (BP) THEN
         WRITE (0,'(A,4(I2,1X))') 'DENSCON : ',
     >        (DENSCON(1,I),I=1, NMOD)
         STOP 'FATAL ERROR IN MODELCO'
      ENDIF

C***  ---- XLAMBDA ----
      RELDIFF = 0.
      DO K=1, ND(1)
        BP = .FALSE.
        DO I=1, NMOD-1
          DO J=1, NMOD
            IF (XLAMBDA(K,I) .NE. XLAMBDA(K,J)) BP = .TRUE.
            IF (XLAMBDA(K,I) .NE. 0.) 
     >        T = ABS(XLAMBDA(K,I) - XLAMBDA(K,J)) / XLAMBDA(K,J)
            IF (T .GT. RELDIFF) RELDIFF = T
          ENDDO
        ENDDO
        IF (.FALSE.) THEN
          WRITE (0,'(A,I2,1X,5(F20.10,1X))') 'DEPTH, XLAMBDA : ',
     >                             K,(XLAMBDA(K,I),I=1, NMOD)
        ENDIF
      ENDDO

C***  ---- RSTAR ----
      RELDIFF = 0.
      BP = .FALSE.
      DO I=1, NMOD-1
        DO J=1, NMOD
          IF (RSTAR(I) .NE. RSTAR(J)) BP = .TRUE.
          T = (RSTAR(I) - RSTAR(J)) / RSTAR(J)
          IF (T .GT. RELDIFF) RELDIFF = T
        ENDDO
      ENDDO
      IF (BP) THEN
        WRITE (0,'(A,I2,1X,5(F20.10,1X))') 
     >        'RSTAR : ',(RSTAR(I),I=1, NMOD)
        WRITE (0,'(A,F20.10)') 'RSTAR : RELDIFF = ',RELDIFF
      ENDIF

C***  ---- VDOP ----
      RELDIFF = 0.
      BP = .FALSE.
      DO I=1, NMOD-1
        DO J=1, NMOD
          IF (VDOP(I) .NE. VDOP(J)) BP = .TRUE.
          T = (VDOP(I) - VDOP(J)) / VDOP(J)
          IF (T .GT. RELDIFF) RELDIFF = T
        ENDDO
      ENDDO
      IF (BP) THEN
        WRITE (0,'(A,I2,1X,5(F20.10,1X))') 'VDOP : ',(VDOP(I),I=1, NMOD)
        STOP 'FATAL ERROR IN MODELCO'
      ENDIF

C***  ---- BLFERR ----
C *** Terminated because of the existance of the new Option 
C *** NO-BLF taken from FORMAL-CARD. Old MODELS will use new Option
C *** BNOCONT instead of BLFERR  

C      BP = .FALSE.
C      DO I=1, NMOD-1
C        DO J=1, NMOD
C          IF (BLFERR(I) .NE. BLFERR(J)) BP = .TRUE.
C        ENDDO
C      ENDDO
C      IF (BP) THEN
C        WRITE (0,'(A,5(L5,1X))') 'BLFERR : ',(BLFERR(I),I=1, NMOD)
C      ENDIF

      WRITE (0,*)

      RETURN
      END
 
 
