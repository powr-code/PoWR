      SUBROUTINE STARKBROAD (KODATIND, NOM, PHITAB,
     >    NFDIMPHITAB, NLDIMPHITAB, NLPHITAB, ND, RADIUS, 
     >    OPAC, OPAL, LINPRO, AVOIGT, NBL, MAXLAP, IPOINTERPHITAB,
     >    XMAX, XMAXBROAD, XMAXLIN, NDDIM, DXMAX, PATH_VCSSB, LOW, NUP,
     >    PHISCRATCH, XLAMNBL, ALN, T, ENTOT, RNE, LEVEL, MAINQN, 
     >    NCHARG, POPNUM, N, VDOP, ELEVEL, EION, PATH_LEMKE_DAT,
     >    DD_VDOP, DD_VDOPDU, IND_ORIGLEV, BDD_VDOP,
     >    DD_VMICDU, GRIEMPAR, TAUMAX, TAUMINBROAD, IMOD, MAXMOD, 
     >    EINST, NDIM)


C***********************************************************************      
C***  This subroutine prepares line broadening. 
C***  For H I and He II lines, the profiles 
C***  for each depth point are stored in array PHITAB.
C***  In ZONEINT, the profile fuction is interpolated from that table.
C***  
C***  For all other lines, the VOIGT function is used, 
C***  and the depth-dependent parameter AVOIGT is prepared here.
C***
C***  called from FORMAL, separately for each line (NBL) in the blend
C***********************************************************************      

      DIMENSION PHITAB(-NFDIMPHITAB:NFDIMPHITAB, NDDIM, NLDIMPHITAB)
      CHARACTER(8) :: LINPRO
      CHARACTER*256 PATH_VCSSB, PATH_LEMKE_DAT
      
      INTEGER, INTENT(IN) :: NDDIM
C* take only vector!
      REAL, DIMENSION(ND), INTENT(IN) :: DD_VDOP, DD_VDOPDU, DD_VMICDU
      DIMENSION KODATIND(2), NOM(2), MAINQN(2), NCHARG(2) 
      DIMENSION ELEVEL(2), EION(2)
      DIMENSION T(ND), ENTOT(ND), RNE(ND)
      DIMENSION IND_ORIGLEV(MAXLAP)
      DIMENSION POPNUM(NDDIM,N)
      DIMENSION AVOIGT(MAXLAP,NDDIM,MAXMOD), GRIEMPAR(MAXLAP,NDDIM,MAXMOD)
      DIMENSION EINST(NDIM,NDIM)

      CHARACTER*10 LEVEL(NDIM)
C***  Arrays for the Lemke tables 
      PARAMETER (NWS_DIM = 54)
      PARAMETER (NTS_DIM = 7)
      PARAMETER (NES_DIM = 17)
      PARAMETER (NPS_DIM = NWS_DIM * NTS_DIM * NES_DIM)
      DIMENSION STKTA_DWS(NWS_DIM), STKTA_TS (NTS_DIM)
      DIMENSION STKTA_ES (NES_DIM), STKTA_PS (NPS_DIM)
      
      LOGICAL BRESONANCE, BDD_VDOP, bHYDROGEN
C***  1 / SQRT(PI)
      DATA WPIINV / 0.564189583549 /
C***  4 * PI
      DATA PI4 / 12.56637062 /
C***  Clight in km/s
      DATA CLIGHT / 299792.458 /
C***  Prepare the profile functions for pressure broadening 
        
C***  Line is resonance line if: 
C**     1) mother lower level is the first to appear in DATOM, 
C***  or   
C***    2) level previous to mother lower level is of a different ion,
C***  or
C***    3) level previous to mother lower level is of different element

      IF (IND_ORIGLEV(LOW) .EQ. 1) THEN
        BRESONANCE = .TRUE.
      ELSE
        BRESONANCE = 
     >          (NCHARG(IND_ORIGLEV(LOW)) .NE. NCHARG(IND_ORIGLEV(LOW)-1)) 
     >          .OR.
     >          (NOM(IND_ORIGLEV(LOW)) .NE. NOM(IND_ORIGLEV(LOW)-1))
      ENDIF 
      
      bHYDROGEN = .FALSE.

C***  Pre-set AVOIGT with radiation-damping 
C***    (will not be used in all cases)

C***  IF NO EXPLICITE VALUE OF GAMMA: SUM OVER EINSTEIN COEFFICIENTS
      IF (AVOIGT(NBL,1,1) .LT. 0.) THEN
         AV=.0
C***     UPPER LEVEL:
         DO I=1, NUP-1
            IF (EINST(NUP,I) .GT. 0.) AV=AV+EINST(NUP,I)
         ENDDO
C***     LOWER LEVEL:
         DO I=1, LOW-1
            IF (EINST(LOW,I) .GT. 0.) AV=AV+EINST(LOW,I)
         ENDDO

         AVOIGT(NBL,1,1) = AV
      ENDIF     

C***  TRANSFORMATION:  GAMMA ---> A (VOIGT FUNCTION)
C***  Note: AVOIGT becomes depth dependent if VDOP is depth dependent.
      AV = AVOIGT(NBL,1,1)
      DO L=1, ND
         DNUED = 1.E13 * DD_VDOP(L) / XLAMNBL
            AVOIGT(NBL,L,IMOD) = AV / PI4 / DNUED
      ENDDO

C***  Choose the appropriate broadening mode LINPRO

C***  If LINPRO has been set already to VOIGT -> keep that setting
      IF (LINPRO .EQ. 'VOIGT   ') THEN

C***  If LINPRO='DRTRANS' was set as default by subr. DRTRANS
      ELSE IF (LINPRO .EQ. 'DRTRANS ') THEN
         LINPRO = ''
      ELSE IF (KODATIND(NOM(LOW)) .EQ. 1) THEN
C***     HYDROGEN BROADENING        
         LINPRO = 'BRD-H   '
         bHYDROGEN = .TRUE.
      ELSE IF (KODATIND(NOM(LOW)) .EQ. 2) THEN
C***     HELIUM BROADENING
         IF (NCHARG(LOW) .EQ. 0) THEN
            LINPRO = 'BRD-HeI '
         ELSE IF (NCHARG(LOW) .EQ. 1) THEN
            LINPRO = 'BRD-HeII'
         ELSE
            STOP '*** INTERNAL ERROR IN STARKPROF: NCHARG'
         ENDIF
      ELSE
C***     ALL OTHER ELEMENTS
         NELECTRON = KODATIND(NOM(LOW)) - NCHARG(LOW)
C***     Linear stark broadening for hydrogenic ions
c           (only possible if MAINQN is different!)
         IF (((NELECTRON == 1) .OR. (NELECTRON == 3)
     >        .OR. (NELECTRON == 11) .OR. (NELECTRON == 19))
     >        .AND. MAINQN(LOW) /= MAINQN(NUP)) THEN      
             LINPRO = 'L-STARK '

C***  Quadratic stark broadening, except for resonance lines  
         ELSEIF (.NOT. BRESONANCE) THEN   
            LINPRO = 'Q-STARK '
         ELSEIF (BRESONANCE) THEN   
            LINPRO = 'VOIGT   '
         ENDIF
      ENDIF
      
C******************************************************************
C********* H I ***************************************************
C******************************************************************
      IF (LINPRO .EQ. 'BRD-H   ') THEN
         NLPHITAB = NLPHITAB + 1
         IPOINTERPHITAB = NLPHITAB
         IF (NLPHITAB .GT. NLDIMPHITAB) THEN
            WRITE (0,*) '*** DIMENSION FOR LINE-BROADENING ',
     >                  ' TABLES (H and He) INSUFFICIENT'
            WRITE (0,'(A,I4)') 
     >                  'Present value: NLDIMPHITAB=', NLDIMPHITAB
            STOP ' *** FATAL ERROR'
         ENDIF

C***     Read Lemke table for the current line
         CALL READ_H_STARKDATA 
     >       (PATH_LEMKE_DAT, MAINQN(LOW), MAINQN(NUP), LEMKEFOUND,
     >       STKTA_DWS, NWS, STKTA_TS, NTS,
     >       STKTA_ES,  NES, STKTA_PS, NPS,
     >       NWS_DIM, NTS_DIM, NES_DIM, NPS_DIM)

         DO L=1, ND
            XNE = ENTOT(L) * RNE(L)

C***        Setting the profile function
            DO K=-NFDIMPHITAB, NFDIMPHITAB
               X = K * DXMAX
               XLAMK = XLAMNBL * EXP(ALN*X)
               IF (LEMKEFOUND .EQ. 1) THEN
                   PHITAB(K, L, NLPHITAB) = STARK_HI_LEMKE
     >              (XLAMK, XLAMNBL, T(L), XNE,
     >               STKTA_DWS, NWS, STKTA_TS, NTS,
     >               STKTA_ES,  NES, STKTA_PS)
               ELSE
                  PHITAB(K, L, NLPHITAB) = STARKHI
     >             (MAINQN(LOW), MAINQN(NUP), XLAMK, XLAMNBL, T(L), XNE,
     >              LEVEL(LOW), LEVEL(NUP)) 
               ENDIF 
            ENDDO

C***        Normalization
            SUM = .0
            DO K=-NFDIMPHITAB, NFDIMPHITAB
               SUM = SUM + PHITAB(K,L,NLPHITAB)
            ENDDO
            FAKNORM = 1. / (DXMAX * SUM)
            DO K=-NFDIMPHITAB, NFDIMPHITAB
               PHITAB(K,L,NLPHITAB) = PHITAB(K,L,NLPHITAB)*FAKNORM
            ENDDO

C***   We realized that the Lemke tables already include thermal broadening. 
C***   Therefore, the routine VDOP_STRUCT now provides the array 
C***   DD_VMICDU (vmic) in addition to DD_VDOPDU. ~~Tomer, 7.5.2015
C***   CONVOLGAUUS_FLEX convolves with gaussian corresponding to an arbitrary VDOP
            NDAT = 2 * NFDIMPHITAB + 1
            DD_VMICDU_SQRD = DD_VMICDU(L)**2
            CALL CONVOLGAUSS_FLEX (PHITAB(-NFDIMPHITAB,L,NLPHITAB), 
     >          PHISCRATCH, NDAT, DXMAX, DD_VDOPDU(L), DD_VMICDU_SQRD)        
         ENDDO

C******************************************************************
C********* He I ***************************************************
C******************************************************************
       ELSE IF (LINPRO .EQ. 'BRD-HeI ') THEN
          
C***   Add Stark broadening by electron impact
          DO L=1, ND
             DLAMD = XLAMNBL * DD_VDOP(L) / CLIGHT
             XNE = ENTOT(L) * RNE(L)
             CALL STARKDAMP_HEI (GAMMAHE1, SHIFT, T(L), XNE, 
     >                 XLAMNBL, NUP, LOW, LINPRO, LEVEL)
             IF (LINPRO .EQ. 'Q-STARK ') EXIT
C***            Add this GAMMA to the "a" (=Voigt function parameter)
                AVOIGT(NBL,L,IMOD) = AVOIGT(NBL,L,IMOD) + GAMMAHE1/DLAMD
             ENDDO

ccc          !!! the follwing part gives nonsense results; 
ccc          !!! Is the returned broadening parameter HEWID in diffent units???  
ccc          !!! So far, I must omit this effect         wrh 11-Jul-2008 
ccc          Attention: in case of a SECOND MODEL, POPNUM would be needed to
ccc            be scaled to te second RADIUS grid
c
cccC***      Add Stark broadening by neutral-helium impact
ccc             DO L=1, ND
cccC***            Number density of neutral helium
ccc                SUM = .0
ccc                DO J = 1, N
ccc                   IF (KODATIND(NOM(J)) .EQ. 2) THEN
ccc                      IF (NCHARG(J) .EQ. 0) SUM = SUM + POPNUM(L,J)
ccc                   ENDIF
ccc                ENDDO
ccc                HE1FRC = ENTOT(L) * SUM
ccc                CALL STARKDAMP_HEI_NEUTRAL (NUP, LOW, T(L),
ccc     >                              HE1FRC, LINPRO, HEWID,
ccc     >                              MAINQN, LEVEL)
ccc                IF (LINPRO .EQ. 'VOIGT   ') EXIT
ccc                AVOIGT(NBL,L,IMOD) = AVOIGT(NBL,L,IMOD) + HEWID/DLAMD
ccc             ENDDO

             
C******************************************************************
C********* He II **************************************************
C******************************************************************
          ELSE IF (LINPRO .EQ. 'BRD-HeII') THEN
             CALL STARKHEIIPREP 
     >          (MAINQN, NUP, LOW, LINPRO, LEVEL, PATH_VCSSB)
C***         if line not in broadening table, LINPRO was changed to L-STARK
             IF (LINPRO .EQ. 'BRD-HeII') THEN
                NLPHITAB = NLPHITAB + 1
                IPOINTERPHITAB = NLPHITAB
                IF (NLPHITAB .GT. NLDIMPHITAB) THEN
                   WRITE (0,*) '*** DIMENSION FOR LINE-BROADENING ',
     >                      ' TABLES (H and He) INSUFFICIENT'
                   WRITE (0,'(A,I4)') 
     >                      'Present value: NLDIMPHITAB=', NLDIMPHITAB
                   STOP ' *** FATAL ERROR'
                ENDIF

                DO L=1, ND
                   XNE = ENTOT(L) * RNE(L)

C***            Setting the profile function
                   DO K=-NFDIMPHITAB, NFDIMPHITAB
                      X = K * DXMAX
                      XLAMK = XLAMNBL * EXP(ALN*X)
                      PHITAB(K, L, NLPHITAB) = STARKHEII
     >                        (XLAMK, XLAMNBL, T(L), XNE) 
                   ENDDO


C***            Normalization
                   SUM = .0
                   DO K=-NFDIMPHITAB, NFDIMPHITAB
                      SUM = SUM + PHITAB(K,L,NLPHITAB)
                   ENDDO
                   FAKNORM = 1. / (DXMAX * SUM)
                   DO K=-NFDIMPHITAB, NFDIMPHITAB
                       PHITAB(K,L,NLPHITAB) = 
     >                      PHITAB(K,L,NLPHITAB)*FAKNORM
                   ENDDO

C***               CONVOLGAUUS_FLEX convolves with gaussian corresponding 
C***               to an arbitrary VDOP
                   NDAT = 2 * NFDIMPHITAB + 1
                   DD_VDOPDU_SQRD = DD_VDOPDU(L)**2
                   CALL CONVOLGAUSS_FLEX
     >                  (PHITAB(-NFDIMPHITAB,L,NLPHITAB), PHISCRATCH, 
     >                   NDAT, DXMAX, DD_VDOPDU(L), DD_VDOPDU_SQRD)  
                ENDDO
             ENDIF         
          ENDIF

C******************************************************************
C********* all other elements beside H and He: quadratic Stark  ***
C***       also fallback branch for He I if no table found
C******************************************************************
          IF (LINPRO == 'L-STARK ') THEN
C***         Linear Stark broadening for hydrogenic ions          
C            GRIEMPAR is prepared here
             DO L=1, ND
                XNE = ENTOT(L) * RNE(L)
                CALL LINSTARK (GRIEMPAR(NBL,L,IMOD), XNE, 
     >             KODATIND(NOM(LOW)),
     >             MAINQN(LOW), MAINQN(NUP), XLAMNBL, 
     >             LINPRO, LEVEL(LOW), LEVEL(NUP))
             ENDDO
          ENDIF
          
C***      Q-STARK is also fallback if no QNs are available for L-STARK
          IF (LINPRO .EQ. 'Q-STARK ') THEN

C***     Add Stark broadening by electron impact (Quadratic Stark Effect)
             DO L=1, ND
                DNUED = 1.E13 * DD_VDOP(L) / XLAMNBL
                XNE = ENTOT(L) * RNE(L)
                CALL QUADSTARK (GAMMAQUAD, T(L), XNE, 
     >             NCHARG(NUP), ELEVEL(NUP), EION(NUP), LINPRO, 
     >             LEVEL(LOW), LEVEL(NUP))
                IF (LINPRO .EQ. 'VOIGT   ') EXIT
C***            Add this GAMMA to the "a" (=Voigt function parameter)
                AVOIGT(NBL,L,IMOD) = AVOIGT(NBL,L,IMOD) + GAMMAQUAD/DNUED
             ENDDO
C***         No preset of VOIGT for DRTRANSIT lines
             IF (NCHARG(LOW) .NE. NCHARG(NUP)) LINPRO = ''
          ENDIF

C***  Bandwidth XMAXBROAD is only calculated for the first MODEL:
      IF (LINPRO .NE. '' .AND. IMOD .EQ. 1) 
     >    CALL BANDWIDTH (ND, RADIUS, OPAC, OPAL, LINPRO, AVOIGT, NDDIM, 
     >                NBL, MAXLAP, XMAX, XMAXBROAD, XMAXLIN,
     >                PHITAB(-NFDIMPHITAB, 1, NLPHITAB), 
     >                NFDIMPHITAB, DXMAX, LEVEL(NUP), LEVEL(LOW),
     >                BDD_VDOP, GRIEMPAR, VDOP, ALN, 
     >                bHYDROGEN, TAUMAX, TAUMINBROAD)

      RETURN
      END

