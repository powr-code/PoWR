      SUBROUTINE ADAPOP (POPNUM, ND, N, POPOLD, NDOLD, NOLD, NCHARG,
     >          NATOM, ABXYZ, NFIRST, NLAST, RNE, NTRANS, POPLTE, 
     >          BDEPART, ADPWEIGHT, RADIUS, ROLD, POPHELP, TAURCONT, 
     >          TAURCONTOLD, POPLTE_OLD, BTAUR, POPMIN) 
C***********************************************************************
C***  TRANSFORMATION OF POPULATION NUMBERS FROM OLD TO NEW MODEL ATOM
C***  Radically simplified version: wrh 10-Aug-2007
C***  The assignment of levels is directed by the vector NTRANS(J)
C***    IF NTRANS(J) = -1 : no action, POPNUM stays from WRSTART 
C***    IF NTRANS(J) =  0 : POPNUM set to ZERO 
C***    IF NTRANS(J) >  0 : assign POPOLD with index NTRANS
C***********************************************************************
 
      DIMENSION NCHARG(N), NTRANS(N), ADPWEIGHT(N)
      DIMENSION ABXYZ(NATOM),NFIRST(NATOM),NLAST(NATOM)
      DIMENSION RNE(ND), RADIUS(ND), TAURCONT(ND)
      DIMENSION POPNUM(ND,N), POPLTE(ND,N)
      DIMENSION POPOLD(ND,NOLD)
      DIMENSION ROLD(NDOLD), TAURCONTOLD(NDOLD)
      DIMENSION POPHELP(NDOLD,NOLD), POPLTE_OLD(NDOLD, NOLD)
      LOGICAL BDEPART, BTAUR

C***  The old popnumbers in POPHELP are interpolated with respect
C***  to the depth coordinate and stored in POPOLD (which then has
C***  still the old atomic levels, but the new radius grid)

C***  Using departure coefficients: replace old POPs by DEPARTs
      IF (BDEPART) THEN
         WRITE (0,*) 
     >      'Old DEPARTure coeficients used instead of POPNUMbers'
         DO L=1, NDOLD
            DO J=1, NOLD
               POPHELP(L,J) = POPHELP(L,J) / POPLTE_OLD(L,J)
            ENDDO
         ENDDO
         DO L=1, ND
            DO J=1, N
               POPNUM(L,J) = POPNUM(L,J) / POPLTE(L,J)
            ENDDO
         ENDDO
      ENDIF


C***  Interpolation on Tau-Grid 
      IF (BTAUR) THEN
         WRITE (0,*) 'Interpolation of Popnumbers on Tau-Grid'
         DO L=1, ND
           DO J=1, NOLD
              IF (TAURCONT(L) .GT. TAURCONTOLD(NDOLD)) THEN
                 POPOLD(L,J) = POPHELP(NDOLD,J)
              ELSE
                 CALL LIPO (POPOLD(L,J), TAURCONT(L), 
     >                    POPHELP(1,J), TAURCONTOLD, NDOLD)
              ENDIF 
           ENDDO
         ENDDO
      
      ELSE

C***  INTERPOLATION OF OLD POPNUMBERS TO THE NEW RADIUS GRID
        WRITE (0,*) 'Interpolation of Popnumbers on Radius-Grid'
        DO L=1, ND
           DO J=1, NOLD
              IF (RADIUS(L) .GT. ROLD(1)) THEN
                 POPOLD(L,J) = POPHELP(1,J)
              ELSE
                 CALL LIPO (POPOLD(L,J), RADIUS(L), 
     >              POPHELP(1,J), ROLD, NDOLD)
              ENDIF
           ENDDO
         ENDDO

      ENDIF

C*****************************************************************
C***  Now the replacement of levels according to NTRANS
C*****************************************************************

C***  Loop over all depth points --------------------------------
      DO L=1, ND

C***  Copy old POPNUMs as assigned 
         DO J=1, N
            IF (NTRANS(J) .EQ. 0) THEN
               POPNUM(L,J) = POPMIN
            ELSE IF (NTRANS(J) .GT. 0 .AND. NTRANS(J) .LE. NOLD) THEN
               POPNUM(L,J) = POPOLD(L,NTRANS(J)) * ADPWEIGHT(J)
            ENDIF
         ENDDO

C***  If POPNUM actually contains DEPARTURE coeficients,
C***   convert them now back to POPNUMs
         IF (BDEPART) THEN
            DO J=1, N
               POPNUM(L,J) = POPNUM(L,J) * POPLTE(L,J)
            ENDDO
         ENDIF

C***  Renormalization to the abundance of each element
 
C***  LOOP FOR EACH ELEMENT  -------------------------------------------
         DO NA=1, NATOM
            SUM=0.0
            NFIRNA = NFIRST(NA)
            NLANA = NLAST(NA)
            DO J = NFIRNA, NLANA
               SUM = SUM + POPNUM(L,J)
            ENDDO
            SUM = SUM / ABXYZ(NA)
            IF (SUM .NE. 0.) THEN
               DO J=NFIRNA,NLANA
                  POPNUM(L,J) = POPNUM(L,J) / SUM
               ENDDO
            ENDIF
         ENDDO
 
C***  Consistent electron density
         RNE(L) = .0
         DO J=1, N
           RNE(L) = RNE(L) + NCHARG(J) * POPNUM(L,J)
         ENDDO
 
      ENDDO   ! Deph points -------------------------------------! 

      RETURN
      END
