      SUBROUTINE TEMPCORR (TOLD, TNEW, ND, RADIUS, TEFF, 
     >           HTOTL, XJTOTL, HTOTCMF0, 
     >           FTCOLI, DTLOCAL, DTINT, DTRMAX, FLUXERR,
     >           OPASMEANTC, QFJMEAN, OPAJMEANTC, OPAPMEAN,
     >           QOPAHMEAN, EDDIHOUTJMEAN,
     >           HTOTOUTMINUS, HTOTND, OPARND, DTDRIN, UNLUTECLINE, 
     >           DUNLU_LOC, DUNLU_INT, DUNLU_RMAX, DUNLU_TB,
     >           TBTAU, TAUINT, OPALAMBDAMEAN, TAUROSS, ARAD,
     >           DTKUBAT, bTDIFFUS, HNDCORFAC, SMEAN, CORRS, FLUXEPS)          
C***********************************************************************
C***  THIS SUBROUTINE APPLIES DIRECTLY A TEMPERATURE CORRECTION, 
C***  FOLLOWING THE UNSOELD-LUCY METHOD
C***  This formalism is described in full detail in:
C***    Hamann, W.-R., Graefener, G.: 2003, Astron. Astrophys., 410, 993
C***  Called from: STEAL
C***********************************************************************

      IMPLICIT NONE

      INTEGER, INTENT(IN) :: ND
      REAL, DIMENSION(ND) :: RADIUS, TOLD, TNEW, XJTOTL,
     >                       DTLOCAL, DTINT, DTRMAX, DTKUBAT, 
     >                       TAUROSS, FTCOLI, HTOTCMF0, CORRS, FLUXERR,
     >                       QOPAHMEAN, OPASMEANTC, 
     >                       QFJMEAN, OPAJMEANTC, OPALAMBDAMEAN, 
     >                       OPAPMEAN, SMEAN
      REAL, DIMENSION(ND-1) :: ARAD, HTOTL
      
      REAL, INTENT(IN) :: TEFF, FLUXEPS, HTOTND, OPARND
      REAL, INTENT(INOUT) :: DTDRIN
      
      LOGICAL :: BSMOOTHTC, BMONOTONIC, BWARNTMIN, bKUBAT,  
     >           bOldFormat, bLOCAL, bOPAPMEAN, bINT2POINT, 
     >           bSKIPCORR, bOBEXTRAP
      LOGICAL, INTENT(IN) :: bTDIFFUS
      
      INTEGER, PARAMETER :: NPARMAX = 40
      
      CHARACTER(*) :: UNLUTECLINE
      CHARACTER(4) :: EDMETHOD
      CHARACTER(40), DIMENSION(NPARMAX) :: ACTPAR
      
      INTEGER :: I, NPAR, IPAR, L, LMINTMIN, LMAXTMIN, NTMIN, Linner,
     >           LOCMIN, LOCMAX, LCUTCORRUSE, iAUTOEXPTAU
      REAL :: TMIN, CORRMAX, CUTCORR, DUNLU_LOC, DUNLU_RMAX, DUNLU_INT, 
     >        DHINTSUM, TREF, OPAMEAN, RL2, RL2M, DHRMAX, DHRMAX1, DH,
     >        DTB1, DTB, DTBTEFF, WRH, TNEW4, TOLD4, WRADIUS, QC,
     >        QS, TLEFT, TNEWL, F, RELCORR, HTOTOUTMINUS, EDDIHOUTJMEAN,
     >        DHINT, DHCORRSUM, DT, DT1, DT2, DT3, DT4, DT5, CMAX, 
     >        DUNLU, EXPTAU, TM, DTDRTRUE, TAUINT,
     >        DUNLU_TB, rfKUBAT, taurfmin, TBTAU, FINTMIN,
     >        tempREAL, DE_LOC, DE_INT, DE_RMAX, DE_TB, HNDCORFAC, 
     >        TAVERAGE, DAMP, DAMPKUBAT, FLUXERRFAC, CLIMIT, FTEPS,
     >        AUTOFLUXTAU
       

      REAL, PARAMETER :: PI4 = 12.5663706144    !PI4 = 4*PI
      REAL, PARAMETER :: STEBOLDPI = 1.8046E-5  !STEFAN-BOLTZMANN CONSTANT (CGS-UNITS) / PI
      REAL, PARAMETER :: AMU = 1.6605E-24       !Atomic mass unit (gramm) = m_H

      !File and channel handles (=KANAL)
      INTEGER, PARAMETER :: hOUT = 6        !write to wruniqX.out (stdout)
      INTEGER, PARAMETER :: hCPR = 0        !write to wruniqX.cpr (stderr)
      
      
C***  Decode parameters of UNLUTEC option card

C***  Default values
      BSMOOTHTC  = .FALSE.
      BMONOTONIC = .FALSE.
      BWARNTMIN  = .FALSE. 
      CORRMAX = 0.
      CUTCORR = .0
      FLUXERRFAC  = .0
      TMIN = 6000.
      bSKIPCORR = .FALSE.
      EXPTAU = 0.
      iAUTOEXPTAU = -1
      FTEPS = 0.01        !Default flux epsilon for automatic flux-tau
      bOBEXTRAP = .FALSE. !Default: Do not extrapolate T for outer boundary

      bOldFormat = .TRUE.
C***  range where TBALANCE correction is gradually switched off
      TBTAU = 0.1
      taurfmin = 0.01       
C***  bOPAPMEAN = use OPAPMEAN instead of OPASMEAN in flux conservation terms 
C***     this version can be requested with the UNLUTEC option 'OPAPMEAN'
      bOPAPMEAN = .FALSE.

C***  test version: T-correction from INT term is distributed over
C***                the two depth points L-1 and L
      bINT2POINT = .FALSE.

C***  Damping of the flux (INT) term      
      TAUINT = -99.
      FINTMIN = -99.

C***  Decode parameters of UNLUTEC option cards 
C***  (already concatenated in one line
      CALL SARGC(UNLUTECLINE, NPAR)
      DO 20, I=2, MIN(NPAR,NPARMAX)
       CALL SARGV(UNLUTECLINE,I,ACTPAR(I))
   20 CONTINUE
   
C***  Check if UNLUTEC line has old or new format
      DO IPAR=2, NPAR
        SELECTCASE (ACTPAR(IPAR))
          CASE ('LOCAL', 'LOC', 'INT', 'OUT', 'RMAX')
            bOldFormat = .FALSE.
          CASE ('BALANCE', 'TBALANCE', 'TB')
            bOldFormat = .FALSE.
        ENDSELECT
      ENDDO

      IF (bOldFormat) THEN
C***    old half-fixed format *** * * * * * * * * * * * * * * * * * * * 

C***    Defaults
        DUNLU = 1.0
        DUNLU_INT = 0.5
        DUNLU_RMAX = 0.5
        DUNLU_TB = 0.

        IF (NPAR .GT. 1) THEN
          READ (ACTPAR(2),'(F10.0)', ERR=1) DUNLU
        ENDIF
        IF (NPAR .GT. 2) THEN
          READ (ACTPAR(3),'(F10.0)', ERR=1) DUNLU_INT
        ENDIF
        IF (NPAR .GT. 3) THEN
          READ (ACTPAR(4),'(F10.0)', ERR=1) DUNLU_RMAX
        ENDIF
        IPAR = 5
        DO WHILE (IPAR .LE. NPAR)
          IF (ACTPAR(IPAR) .EQ. 'SMOOTH') THEN
            BSMOOTHTC = .TRUE.
            IPAR = IPAR + 1
          ELSE IF (ACTPAR(IPAR) .EQ. 'MONOTONIC') THEN
            BMONOTONIC = .TRUE.
            IPAR = IPAR + 1

          ELSE IF (ACTPAR(IPAR) .EQ. 'CORRMAX') THEN
            IF (IPAR+1 .GT. NPAR) THEN
              WRITE (0,*) 'Error when decoding CORRMAX'
              STOP 'ERROR in Subr. TEMPCORR'
            ENDIF
            READ (ACTPAR(IPAR+1),'(F10.0)', ERR=1) CORRMAX
            IPAR = IPAR + 2

          ELSE IF (ACTPAR(IPAR) .EQ. 'CUTCORR') THEN
            IF (IPAR+1 .GT. NPAR) THEN
              WRITE (0,*) 'Error when decoding CUTCORR'
              STOP 'ERROR in Subr. TEMPCORR'
            ENDIF
            READ (ACTPAR(IPAR+1),'(F10.0)', ERR=1) CUTCORR
            IPAR = IPAR + 2

          ELSE IF (ACTPAR(IPAR) == 'KUBAT') THEN
            WRITE (0,*) 'ERROR: obsolete keyword KUBAT'
            WRITE (0,*) 'ERROR: -- use new UNLUTEC syntax!'
            STOP 'FATAL ERROR detected in subr. TEMPCORR'

          ELSE IF (ACTPAR(IPAR) .EQ. 'TMIN') THEN
            IF (IPAR+1 .GT. NPAR) THEN
              WRITE (0,*) 'Error when decoding TMIN'
              STOP 'ERROR in Subr. TEMPCORR'
            ENDIF
            READ (ACTPAR(IPAR+1),'(F10.0)', ERR=1) TMIN
            IPAR = IPAR + 2
          ELSE
            IPAR = IPAR + 1
          ENDIF
        ENDDO
        !Convert old multiple factors to new direct factors
        DUNLU_LOC = DUNLU                !local (radiative balance)
        IF (DUNLU > 0.) THEN
          DUNLU_INT  = DUNLU * DUNLU_INT     !integral term        
          DUNLU_RMAX = DUNLU * DUNLU_RMAX    !Rmax term
        ENDIF

      ELSE
C***    new, flexible format * * * * * * * * * * * * * * * * * * * * * * 
        IPAR = 2
C***    Defaults (new format)
        DUNLU_LOC  = 0.
        DUNLU_INT  = 1.
        DUNLU_RMAX = 1.        
        DUNLU_TB   = 1.

        DO WHILE (IPAR <= NPAR)
          SELECTCASE (ACTPAR(IPAR))
            CASE ('LOC', 'LOCAL')
              IF (IPAR+1 > NPAR) THEN
                WRITE (hCPR,*) 'Error when decoding LOCAL'
                STOP 'ERROR in Subr. TEMPCORR'
              ENDIF
              READ (ACTPAR(IPAR+1),'(F10.0)', ERR=1) DUNLU_LOC
              IPAR = IPAR + 2
            CASE ('INT')
              IF (IPAR+1 > NPAR) THEN
                WRITE (hCPR,*) 'Error when decoding INT'
                STOP 'ERROR in Subr. TEMPCORR'
              ENDIF
              READ (ACTPAR(IPAR+1),'(F10.0)', ERR=1) DUNLU_INT
              IPAR = IPAR + 2              
            CASE ('OUT', 'RMAX')
              READ (ACTPAR(IPAR+1),'(F10.0)', ERR=1) DUNLU_RMAX
              IPAR = IPAR + 2              

            CASE ('KUBAT', 'TBALANCE', 'BALANCE', 'TB')
              IF (IPAR+1 > NPAR) THEN
                WRITE (hCPR,*) 'Error when decoding TBALANCE'
                STOP 'ERROR in Subr. TEMPCORR'
              ENDIF
              READ (ACTPAR(IPAR+1),'(F10.0)', ERR=1) DUNLU_TB
              IPAR = IPAR + 2 
              
            CASE ('TBTAU')
C***          sets custom tau-value for switching off TB method 
              IF (IPAR+1 > NPAR) THEN
                WRITE (hCPR,*) 'Error when decoding TBTAU'
                STOP 'ERROR in Subr. TEMPCORR'
              ENDIF
              READ (ACTPAR(IPAR+1),'(F10.0)', ERR=1) tempREAL
              IF (tempREAL > 0.) THEN
                TBTAU = tempREAL
                taurfmin = 0.1 * TBTAU
              ELSE 
                WRITE (hCPR,*) 'Error when decoding TBTAU'
                STOP 'ERROR in Subr. TEMPCORR'
              ENDIF
              IPAR = IPAR + 2

            CASE ('OPAPMEAN')
              bOPAPMEAN = .TRUE.
              IPAR = IPAR + 1
            CASE ('OPASMEAN')
              bOPAPMEAN = .FALSE.
              IPAR = IPAR + 1

            CASE ('SMOOTH')
              BSMOOTHTC = .TRUE.
              IPAR = IPAR + 1

            CASE ('MONO', 'MONOTONIC')
              BMONOTONIC = .TRUE.
              IPAR = IPAR + 1
            CASE ('EXPTAU', 'COSTAU')
              IF (IPAR+1 > NPAR) THEN
                WRITE (hCPR,*) 'Error when decoding EXPTAU'
                STOP 'ERROR in Subr. TEMPCORR'
              ENDIF
              IF (ACTPAR(IPAR) == 'COSTAU') THEN
                EDMETHOD = 'COS'
              ELSE 
                EDMETHOD = 'EXP'
              ENDIF              
              IF (ACTPAR(IPAR+1) /= 'AUTO') THEN
                IF (ACTPAR(IPAR+1)(1:1) /= 'A') THEN
C***              Explicit value                
                  READ (ACTPAR(IPAR+1),'(F10.0)', ERR=1) EXPTAU
                ELSEIF (ACTPAR(IPAR+1)(1:4) == 'AMIN') THEN
C***              Explicit minimum value                
                  iAUTOEXPTAU = 1
                  READ (ACTPAR(IPAR+1)(5:),'(F10.0)', ERR=1) EXPTAU
                ELSE 
C***              Explicit maximum value                
                  iAUTOEXPTAU = 0
                  IF (ACTPAR(IPAR+1)(1:4) == 'AMAX') THEN
                    WRITE (0,*) 'AMAX'
                    READ (ACTPAR(IPAR+1)(5:),'(F10.0)', ERR=1) EXPTAU
                  ELSE
                    WRITE (0,*) 'A'
                    READ (ACTPAR(IPAR+1)(2:),'(F10.0)', ERR=1) EXPTAU
                  ENDIF
                ENDIF
              ELSE 
                iAUTOEXPTAU = 0
                EXPTAU = -1.
              ENDIF
              IPAR = IPAR + 2
            CASE ('FTEPS')
C***          sets custom epsilon for automatic flux-tau
              IF (IPAR+1 > NPAR) THEN
                WRITE (hCPR,*) 'Error when decoding FTEPS'
                STOP 'ERROR in Subr. TEMPCORR'
              ENDIF
              READ (ACTPAR(IPAR+1),'(F10.0)', ERR=1) tempREAL
              IF (tempREAL > 0.) THEN
                FTEPS = -1. * FTEPS
              ELSE 
                WRITE (hCPR,*) 'Error when decoding FTEPS'
                STOP 'ERROR in Subr. TEMPCORR'
              ENDIF
              IPAR = IPAR + 2              
            CASE ('OBEX', 'OBEXTRAP')
              bOBEXTRAP = .TRUE.
              IPAR = IPAR + 1
                                                        
            CASE ('SKIPNEGFLUX')
              IF (MINVAL(ARAD) <= 0. .OR. MINVAL(HTOTL) <= 0.) THEN
                bSKIPCORR = .TRUE.
              ENDIF
              IPAR = IPAR + 1
                            
            CASE ('CORRMAX', 'CM')
              IF (IPAR+1 > NPAR) THEN
                WRITE (hCPR,*) 'Error when decoding CORRMAX'
                STOP 'ERROR in Subr. TEMPCORR'
              ENDIF
              READ (ACTPAR(IPAR+1),'(F10.0)', ERR=1) CORRMAX
              IPAR = IPAR + 2

            CASE ('LIMIT-TO-FLUXERR')
              IF (IPAR+1 > NPAR) THEN
                WRITE (hCPR,*) 'Error: LIMIT-TO-FLUXERR needs value'
                STOP 'ERROR in Subr. TEMPCORR'
              ENDIF
              READ (ACTPAR(IPAR+1),'(F10.0)', ERR=1) FLUXERRFAC
              IPAR = IPAR + 2

            CASE ('CUTCORR', 'CC')
              IF (IPAR+1 > NPAR) THEN
                WRITE (hCPR,*) 'Error when decoding CUTCORR'
                STOP 'ERROR in Subr. TEMPCORR'
              ENDIF
              READ (ACTPAR(IPAR+1),'(F10.0)', ERR=1) CUTCORR
              IPAR = IPAR + 2

            CASE ('TMIN')
              IF (IPAR+1 > NPAR) THEN
                WRITE (hCPR,*) 'Error when decoding TMIN'
                STOP 'ERROR in Subr. TEMPCORR'
              ENDIF
              READ (ACTPAR(IPAR+1),'(F10.0)', ERR=1) TMIN
              IPAR = IPAR + 2
            CASE ('FINTMIN')
              IF (IPAR+1 > NPAR) THEN
                WRITE (hCPR,*) 'Error when decoding FINTMIN'
                STOP 'ERROR in Subr. TEMPCORR'
              ENDIF
              READ (ACTPAR(IPAR+1),'(F10.0)', ERR=1) FINTMIN
              IPAR = IPAR + 2
            CASE ('TAUINT')
              IF (IPAR+1 > NPAR) THEN
                WRITE (hCPR,*) 'Error when decoding TAUINT'
                STOP 'ERROR in Subr. TEMPCORR'
              ENDIF
              READ (ACTPAR(IPAR+1),'(F10.0)', ERR=1) TAUINT
              IPAR = IPAR + 2    
            CASE DEFAULT
              !unknown parameter => skip
              IPAR = IPAR + 1
          ENDSELECT
        ENDDO
      ENDIF
      
C***  Default values for the tau-dependent damping of the flux terms
C***    (INT and RMAX terms, i.e. 2nd and 3rd):
C***  Note that at least one of the optional parameters: TAUINT or FINTMIN
C***  must be given in order to activate this damping factor
      IF (FINTMIN < 0.)  THEN
        FINTMIN = 0.1
      ELSEIF (TAUINT < 0.) THEN
        TAUINT = 10.
      ENDIF

C***  end of decoding UNLUTEC options * * * * * * * * * * * * * * * * * 


C***********************************************************************
C***    Prepare third Term of Unsoeld-Lucy-Procedure: outer Boundary
C***********************************************************************
      DHRMAX1 = (HTOTCMF0(1) - HTOTOUTMINUS) / EDDIHOUTJMEAN -
     >           XJTOTL(1)
      DHINTSUM = .0

C**********************************************************************
C***  Loop over depth index L 
C***
C***   The temperature correction contributions from the different terms
C***   are stored in vectors DTLOCAL, DTINT, DTRMAX, DTKUBAT
C***   without their respective scaling factors, but with their
C***   depth-dependent damping, 
C***   for being plotted via the option PLOT UNLU
C**********************************************************************

C***  Note: at L=ND no correction is applied in case of TDIFFUS artistics
      Linner = ND-1
      IF (.NOT. bTDIFFUS) THEN
        Linner = ND
      ELSE
        TNEW(ND) = TOLD(ND)
        CORRS(ND) = STEBOLDPI*TOLD(ND)**4 / SMEAN(ND)
      ENDIF      
      
C***  Set epsilon for automatic flux tau-determination       
C***  This value is determined from one of the following criteratia
C***  in descending order:
C***  1. If specified, direct input parameter FTEPS on UNLUTEC card
C***  2. If specified, adapted from FLUXEPS card (i.e. desired convergence epsilon)
C***  3. Fallback: Default value given at the start of this routine
      IF (FTEPS < 0.) THEN
C***    Negative values denote a direct input via UNLUTEC FTEPS
C***    while a positive is the fallback given at the start of this routine
        FTEPS = -1. * FTEPS
      ELSEIF (FLUXEPS > 0.) THEN
C***    No direct input via UNLUTEC, but FLUXEPS card is specified:
        FTEPS = FLUXEPS
      ENDIF
      CALL TEMPCORR_FLUXERR(FLUXERR, ND, HTOTL, HTOTCMF0,
     >                      FTEPS, TAUROSS, AUTOFLUXTAU)
C***  Apply AUTOFLUXTAU is EXPTAU=AUTO in UNLUTEC card:
      IF (iAUTOEXPTAU >= 0) THEN
        IF (EXPTAU < 0.) THEN
C***      No EXPTAU input value given => just use automatic value     
          EXPTAU = AUTOFLUXTAU
        ELSEIF (iAUTOEXPTAU >= 1) THEN
C***      EXPTAU defines minimum value        
          EXPTAU = MAX(AUTOFLUXTAU, EXPTAU)
        ELSE
C***      EXPTAU defines maximum value        
          EXPTAU = MIN(AUTOFLUXTAU, EXPTAU)
        ENDIF
      ENDIF
      
C***  from numerical tests we think that the local term is not good
C***  at the inner boundary; andreas and wrh 26-Feb-2016
      FTCOLI(ND) = .0

      DO L=1, Linner
C***      (Inverse) departure of the source function from Planck
C***      (This is used in COMA for Fe rate correction damping)
          CORRS(L) = STEBOLDPI*TOLD(L)**4 / SMEAN(L)
      
C***      Under the assumption S = Bnue(Tref), the error term will be 
C***      converted into a temperature correction (in kK)

C***      TREF as mean with neighbouring points 
          IF (L == 1) THEN
            TREF = (TOLD(L) + TOLD(L+1))/2. 
          ELSEIF (L == ND) THEN 
            TREF = (TOLD(L-1) + TOLD(L))/2.
          ELSE
            TREF = (TOLD(L-1) + TOLD(L) + TOLD(L+1))/3. 
          ENDIF
                    
C***      Fill factors for optional additional exponential damping
          CALL TEMPCORR_EXPDAMP(EXPTAU, EDMETHOD, TAUROSS, 
     >                          bSKIPCORR, 
     >                          L, ND, DE_LOC, DE_INT, DE_RMAX, DE_TB)          
          
C******************************************************************
C***     First (LOCAL) term (not: thermal balance method!) 
C***     Accelerated ALO for first ("local") Term:
C***     DTB1 = correction from this term
C******************************************************************

          OPAMEAN = OPASMEANTC(L) - OPALAMBDAMEAN(L)
          DTB1    = 1. / 
     >       (TREF * TREF * TREF * 4. * STEBOLDPI * OPAMEAN * 1000.)

          DTLOCAL(L) = -FTCOLI(L) * DTB1    ! Stored for PLOT UNLU
          DT1 = DUNLU_LOC * DE_LOC * DTLOCAL(L)

C******************************************************************
C***     Second (INTEGRAL) term ) 
C***     DTB = correction from this term
C******************************************************************

C***      Note: in contrast to Hamann & Graefener, it is preabably more
C***      stable to use the Planck-weighted mean opacity (OPAPMEAN)
C***         and not the source-fuction weighted opacity OPASMEANTC
          IF (bOPAPMEAN) THEN
            OPAMEAN = OPAPMEAN(L)
          ELSE
            OPAMEAN = OPASMEANTC(L)
          ENDIF

          DTB     = 1. / 
     >       (TREF * TREF * TREF * 4. * STEBOLDPI * OPAMEAN * 1000.)


C***    The corrections of the previously treated depth points 
C***      are integrated in dhcorrsum and also taken into account
ccc this is not described in the UNLU paper, and might be questioned!
ccc  wrh  4-Jul-2016

          RL2 = RADIUS(L) * RADIUS(L)
          IF (L .EQ. 1) THEN
            DHINT = .0
            DHCORRSUM = .0
          ELSE

            DH = HTOTL(L-1) - HTOTCMF0(L-1)
            WRADIUS = RADIUS(L) - RADIUS(L-1)

C***        Already performed corrections 
            IF (L .GT. 2) THEN 
             TNEW4 = TNEW(L-1) * TNEW(L-1) * TNEW(L-1) * TNEW(L-1) 
             TOLD4 = TOLD(L-1) * TOLD(L-1) * TOLD(L-1) * TOLD(L-1)  
             RL2M = RADIUS(L-1) * RADIUS(L-1)
             WRH = 0.5 * (RADIUS(L) - RADIUS(L-2))
             IF (bOPAPMEAN) THEN
               DHCORRSUM = DHCORRSUM + 
     >           OPAPMEAN(L-1) * STEBOLDPI * (TNEW4-TOLD4) *RL2M*WRH
             ELSE
               DHCORRSUM = DHCORRSUM + 
     >           OPASMEANTC(L-1) * STEBOLDPI * (TNEW4-TOLD4) *RL2M*WRH
             ENDIF

             DH = DH + DHCORRSUM
            ENDIF 

            DHINTSUM = DHINTSUM + DH * WRADIUS * QOPAHMEAN(L-1)          
C***        Factors outside integral
            DHINT = OPAJMEANTC(L) * DHINTSUM / (RL2 * QFJMEAN(L))

          ENDIF

C***      Damping of the INT (i.e. 2nd) term (only if requested)
          IF (TAUINT > 0.) THEN 
            DAMP = (1.-FINTMIN) * (1.-EXP(-TAUROSS(L)/TAUINT)) + FINTMIN            
          ELSE
            DAMP = 1.
          ENDIF

C***      Debug Output of depth-dependent INT damping:
ccc          write (0,'(I3, 2F12.3)') L, DAMP, TAUROSS(L)

          DTINT(L) = DAMP * DHINT * DTB  ! stored for PLOT UNLU
          DT2 = DUNLU_INT * DE_INT * DTINT(L)

ccc test version wrh 22-Jul-2016
ccc The 2nd term correction is applied half to the current, and half to
ccc  the last (L-1) point  -- see also below!!
          IF (bINT2POINT) THEN
            DTINT(L) = 0.5 * DAMP * DHINT * DTB  
            IF (L > 1) DTINT(L-1) = DTINT(L-1) + 0.5 * DAMP * DHINT*DTB
          ENDIF

C***********************************************************************
C***  Third Term of Unsoeld-Lucy-Procedure: Flux at outer Boundary
C***********************************************************************
C***     Special DTB conversion for 3rd term using TEFF as reference
          DTBTEFF = 1. / 
     >       (TEFF *TEFF* TEFF * 4. * STEBOLDPI * OPAMEAN * 1000.)
          DHRMAX = OPAJMEANTC(L) * DHRMAX1 * QFJMEAN(1) 
     >           / (QFJMEAN(L) * RL2) 

C***      Damping of the RMAX (i.e. 3rd) term with same factor as 2nd
          DTRMAX(L)  = DAMP * DHRMAX  * DTBTEFF ! stored for PLOT UNLU
          DT3 = DUNLU_RMAX * DE_RMAX * DTRMAX(L)

C********************************************************************     
C***     4. term: Thermal Balance method (alternative for the LOCAL term)
C***     This correction has been already calculated by subr. PREPKUBAT
C***      in the vector DTKUBAT and is here only scaled and damped.
C***     TBALANCE term is switched off for TAUROSS > TBTAU 
C***      with the factor rfKUBAT = 1 ... 0, starting from 0.1*TBTAU
C********************************************************************
          IF (TAUROSS(L) .GE. TBTAU) THEN
             rfKUBAT = 0.
             DTKUBAT(L) = 0.
          ELSEIF (TAUROSS(L) .LE. TAURFMIN) THEN
             rfKUBAT = 1.
          ELSE
             rfKUBAT = (TBTAU-TAUROSS(L))/(TBTAU-TAURFMIN)
             DTKUBAT(L) = DTKUBAT(L) * rfKUBAT
          ENDIF

          DT4 = DUNLU_TB * DE_TB * DTKUBAT(L) 


C**********************************************************************
C***  SUM of all corrections **************************
C**********************************************************************

ccc test version wrh 22-Jul-2016
ccc The 2nd term correction is applied half to the current, and half to
ccc  the last (L-1) point  
         IF (bINT2POINT) THEN
           DT = (DT1 + 0.5*DT2 + DT3 + DT4) * 1000. 
           IF (L > 1) TNEW(L-1) = TNEW(L-1) + 0.5*DT2 * 1000.
         ELSE
           DT = (DT1 + DT2 + DT3 + DT4) * 1000. 
         ENDIF

         TNEW(L) = TOLD(L) + DT      


C***     Restrict here the T-correction to factor of two increase
C***      this is necessary to avoid a crash at the next depth point 
C***      in the local (1st) term when applying (TNEW4-TOLD4)
         TNEW(L) = MIN(TNEW(L), 2.*TOLD(L))

C***     Restrict TNEW to TMIN 
         TNEW(L) = MAX(TNEW(L), TMIN)
        
      ENDDO
C***  END OF THE LOOP OVER DEPTH INDEX L ******************************


C**********************************************************************
C***  Inner boundary
C**********************************************************************

C***  Adjustment of the temperature gradient ABS(DTDR) at the inner boundary
C***  in order to ensure the correct flux that corresponds to TEFF
C***  Correction factor HNDCORFAC is applied to the numerical gradient
C***  between the innermost two depth points
C***  DTDRIN will be written to MODEL and used by COLI - CLDIFFUS
C***  Damping: 10 times more damping than for DHINT

      HNDCORFAC = 0.25 * STEBOLDPI * TEFF**4 / HTOTND
      HNDCORFAC = 1. + (HNDCORFAC - 1.) * DUNLU_INT * 0.1
      DTDRTRUE = (TOLD(ND)-TOLD(ND-1))/(RADIUS(ND-1)-1.)
      DTDRIN = DTDRTRUE * HNDCORFAC 

C***  Option: Determine temperature of outer boundary via extrapolation
      IF (bOBEXTRAP) THEN
        TNEW(1) = TNEW(2) - (TNEW(3) - TNEW(2)) /
     >            (QOPAHMEAN(2)*(RADIUS(2) - RADIUS(3)))
     >          * (RADIUS(1) - RADIUS(2)) * QOPAHMEAN(1)
      ENDIF
      
C***  Smoothing of temperature corrections
      IF (BSMOOTHTC) THEN
        QC = 0.6
        WRITE (hCPR,'(A)') 'TEMPCORR: Temperature corrections smoothed'
      ELSE
        QC = 1.0
      ENDIF
      QS = 0.5 * (1. - QC)
      TLEFT = TNEW(1)
      TNEW(1) = TOLD(1) + (QC+QS) * (TNEW(1) - TOLD(1))  
     >               +   QS    * (TNEW(2) - TOLD(2))
      DO L=2, Linner-1
        TNEWL = TOLD(L) + QC * (TNEW(L  ) - TOLD(L  ))
     >                  + QS * (TLEFT     - TOLD(L-1))
     >                  + QS * (TNEW(L+1) - TOLD(L+1))
        TLEFT = TNEW(L)
        TNEW(L) = TNEWL
      ENDDO
      TNEW(Linner) = TOLD(Linner) 
     >               + (QC+QS) * (TNEW(Linner) - TOLD(Linner)) 
     >               +   QS    * (TLEFT - TOLD(Linner-1))

C**************************************************************
C*** Enforcing monotonic temperature structure
C**************************************************************
      IF (BMONOTONIC) THEN
      WRITE (0,*) 'TEMPCORR: Monotonic Temperature enforced'
   10    CONTINUE

C***     Find local maximum
         LOCMAX=0
         DO L=1, Linner-1
            IF (TNEW(L+1) < TNEW(L)) THEN
              LOCMAX=L
              EXIT
            ENDIF
         ENDDO
         IF (LOCMAX .EQ. 0) GOTO 11
 
C***     Find following local minimum
         IF (LOCMAX == Linner-1) THEN
C***       If the second-innermost point is a local maximum,
C***       the innermost point must have a lower temperature, thus:
           LOCMIN = Linner
         ELSE
           DO L= LOCMAX+1, Linner-1
              IF (TNEW(L+1) > TNEW(L) .OR. L .EQ. Linner-1) THEN
                  LOCMIN=L
                  EXIT
              ENDIF
           ENDDO
         ENDIF
C***     Average between local maximum and minimum
         TAVERAGE = .0
         DO L= LOCMAX, LOCMIN
             TAVERAGE = TAVERAGE + TNEW(L)
         ENDDO
         TAVERAGE = TAVERAGE / (LOCMIN - LOCMAX +1)

C***     Replace all T between with the average
         DO L= LOCMAX, LOCMIN
            TNEW(L) = TAVERAGE
         ENDDO

C***     Replace all T further out with average if yet hotter
         DO L= LOCMAX-1, 1, -1
            IF (TNEW(L) > TAVERAGE) THEN
               TNEW(L) = TAVERAGE
            ELSE
               EXIT
            ENDIF
         ENDDO

C***     Replace all T further in with average if yet cooler
         DO L= LOCMIN+1, Linner
            IF (TNEW(L) < TAVERAGE) THEN
               TNEW(L) = TAVERAGE
            ELSE
               EXIT
            ENDIF
         ENDDO
         GOTO 10  ! look for further non-monotonic parts    

   11    CONTINUE !*** no further non-monotonic parts found ****
      ENDIF

C********************************************************************
C***  Minimum temperature
C********************************************************************
      DO L=1, Linner
         IF (TNEW(L) .LT. TMIN) THEN
            TNEW(L) = TMIN
            IF (BWARNTMIN) THEN
               LMINTMIN = MIN0 (LMINTMIN, L)
               LMAXTMIN = MAX0 (LMAXTMIN, L)
               NTMIN = NTMIN + 1
            ELSE
               LMINTMIN = L
               LMAXTMIN = L
               BWARNTMIN = .TRUE.
               NTMIN = 1
            ENDIF         
         ENDIF
      ENDDO

C**************************************************************
C***  Restrict T-correction to  CORRMAX (if specified)
C**************************************************************
      IF (CORRMAX .GT. 0.) THEN
        CMAX = 0.
        DO L=1, Linner
          CMAX = AMAX1(CMAX, ABS(TNEW(L)-TOLD(L)) / TOLD(L))
        ENDDO
        IF (CMAX .GT. CORRMAX) THEN
          F = CORRMAX / CMAX
          WRITE (hCPR,'(A, F7.3)')
     >      'TEMPCORR: Delta-T was reduced by F=', 1./F
          DO L=1, Linner
            TNEW(L) = TOLD(L) + (TNEW(L) - TOLD(L)) * F
          ENDDO
        ENDIF
      ENDIF

C**************************************************************
C***  Scale rel. T-correction to a maximum of 
C***   FLUXERR*FLUXERRFAC (if specified)
C**************************************************************
      IF (MAXVAL(FLUXERR) > 0. .AND. FLUXERRFAC .GT. 0.) THEN
        CMAX = 0.
        DO L=1, Linner
          CMAX = AMAX1(CMAX, ABS(TNEW(L)-TOLD(L)) / TOLD(L))
        ENDDO
          CLIMIT = MAXVAL(FLUXERR) * FLUXERRFAC
        IF (CMAX .GT. CLIMIT) THEN
          F = CLIMIT / CMAX
          WRITE (hCPR,'(A, F7.3)') 'TEMPCORR: LIMIT-TO-FLUXERR: ' 
     >      // 'Delta-Temp was reduced by factor ', 1./F
          DO L=1, Linner
            TNEW(L) = TOLD(L) + (TNEW(L) - TOLD(L)) * F
          ENDDO
        ENDIF
      ENDIF

C**************************************************************
C***  Restrict relative T-correction to CUTCORR (if specified)
C**************************************************************
      IF (CUTCORR .GT. 0.) THEN
        LCUTCORRUSE = 0
        DO L=1, Linner
          RELCORR = (TNEW(L)-TOLD(L)) / TOLD(L)
          IF (RELCORR .GT. CUTCORR) THEN
             RELCORR = CUTCORR
             LCUTCORRUSE = L
          ELSE IF (RELCORR .LT. -CUTCORR) THEN
             RELCORR = -CUTCORR
             LCUTCORRUSE = L
          ENDIF
          TNEW(L) = TOLD(L) * (1. + RELCORR)
        ENDDO
        IF (LCUTCORRUSE > 0) THEN
          WRITE (hCPR,'(A, F7.3, A, I3)')
     >      'TEMPCORR: Corrections were cut to a relative level '
     >       // 'of ', CUTCORR, '. Innermost cut at L =', L
        ENDIF
      ENDIF
      
      IF (bINT2POINT) 
     >   WRITE (0, *) 'TEMPCORR: INT2POINT version activated'

      IF (BWARNTMIN) WRITE (*,100) TMIN, NTMIN, LMINTMIN, LMAXTMIN
  100 FORMAT ('STEAL> WARNING FROM TEMPCORR: TEMPERATURE REACHED TMIN=', 
     >        F7.0, ' K AT', I3, ' DEPTH POINTS BETWEEN L=', I3, 
     >        ' AND', I3 )

C***  Set switches to transfer information which term was used to PLOT UNLU routine     
      IF (DUNLU_TB > 0.) THEN
        bKUBAT = .TRUE.
      ELSE
        bKUBAT = .FALSE.
      ENDIF
      IF (DUNLU_LOC > 0.) THEN
        bLOCAL = .TRUE.
      ELSE
        bLOCAL = .FALSE.
      ENDIF
          
      RETURN

C***  ERROR EXIT **************************************************
1     CONTINUE
      WRITE (0,*)'TEMPCORR: ERROR WHILE DECODING THE FOLLOWING LINE:'
      WRITE (0,*) UNLUTECLINE
      STOP 'ERROR'

      END
