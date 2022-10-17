      SUBROUTINE WMODCOLI(XJCINT, FWTEST, ARAD, ACONT, ATHOM,
     >                    ND, NF, RADIUS, ENTOT, RSTAR,
     >                    XJFEMEAN, SIGMAINT, LASTFE, 
     >                    HTOTL, HTOTND, HTOTNDCOR, 
     >                    WFELOW, WFENUP, FTCOLI, NCOLIP, 
     >                    XJTOTL, XKTOTL, XNTOTL, WJC,
     >                    DBDTINT, DBDTOPAINT, DBDTINT_M, DBDTOPAINT_M,
     >                    OPASMEAN, OPASMEANTC, OPAPMEAN, QFJMEAN, 
     >                    SMEAN, OPAJMEAN, OPAJMEANTC, 
     >                    QOPAHMEAN, EDDIHOUTJMEAN, 
     >                    HTOTOUTMINUS, LASTIND, 
     >                    FERATUL, FERATLU, BCOLIP, EPSGMAX, 
     >                    FTFE, EMCOLI, FF_INFO, IFF_DK, IFF_MAX_MS,
     >                    IFF_WCHARM, OPALAMBDAMEAN, TOTOUT, bKALPHA,
     >                    ALPHAF, RHO, XMU, TAUROSS, OPAROSS)
C****************************************************************
C*** Handles all MODEL write (MSWRIT) statements from COLI
C***  only called by COLI main program
C****************************************************************

      IMPLICIT NONE

      INTEGER, PARAMETER :: TINYINT = SELECTED_INT_KIND(2)

      INTEGER, INTENT(IN) :: ND, NF, LASTFE, IFF_MAX_MS

      REAL, DIMENSION(ND) :: ARAD, ACONT, ATHOM,
     >                       RADIUS, ENTOT,
     >                       HTOTL, FTCOLI, 
     >                       XJTOTL, XKTOTL, XNTOTL,
     >                       EPSGMAX, OPALAMBDAMEAN, SMEAN,
     >                       ALPHAF, RHO, XMU,          !used for hydro
     >                       QOPAHMEAN, 
     >                       OPAJMEANTC, OPAJMEAN, OPASMEAN,
     >                       OPASMEANTC, OPAPMEAN, QFJMEAN

      REAL, DIMENSION(LASTFE) :: SIGMAINT
      REAL, DIMENSION(ND,NF) :: XJCINT, WJC
      REAL, DIMENSION(LASTFE,ND) :: FERATLU, FERATUL
      REAL, DIMENSION(ND,LASTFE) :: XJFEMEAN, FTFE, WFELOW, WFENUP
      REAL, DIMENSION(NF) :: EMCOLI, FWTEST
      REAL, DIMENSION(10) :: FF_INFO
C      INTEGER, DIMENSION(IFF_MAX_MS) ::  IFF_DK
      INTEGER(KIND=TINYINT), DIMENSION(IFF_MAX_MS) :: IFF_DK   !erst, wenn MS-Storage aktualisiert
      INTEGER(KIND=TINYINT), DIMENSION(IFF_MAX_MS,ND) :: IFF_WCHARM
      
      REAL, DIMENSION(ND), INTENT(IN) :: TAUROSS, OPAROSS

      INTEGER :: L, INDF, NCOLIP, LASTIND, IND, INDFE, IDUMMY, IERR,
     >           IFF_N_MS
      REAL :: EDDIHOUTJMEAN, DBDTINT_M, DBDTOPAINT_M, OPARNDM,
     >        DBDTOPAINT, DBDTINT, OPARND, RSTAR, HTOTOUTMINUS,
     >        TOTOUT, HTOTND, HTOTNDCOR

      LOGICAL :: BCOLIP, bKALPHA
      
      CHARACTER(8) :: NAME, CKIND

C***  NORMALIZE AND SAVE XJC'S
      DO INDF=1, NF
         DO L=1, ND
            IF (FWTEST(INDF) .NE. 0.) THEN
               XJCINT(L, INDF) = XJCINT(L, INDF) / FWTEST(INDF)
     >                                 /RADIUS(L)/RADIUS(L)
C!!!  The Division by FWTEST is now done in FREQUNORM
C!!!               WJC(L,INDF) = WJC(L,INDF)/FWTEST(INDF)
            ELSE
               XJCINT(L, INDF) = 0.
               WJC(L,INDF) = 0.
            ENDIF
         ENDDO
         WRITE (NAME,'(A3,I4)') 'XJC', INDF
         CALL WRITMS 
     >        (3, XJCINT(1, INDF), ND, NAME, -1, IDUMMY, IERR)
         WRITE (NAME,'(A3,I4)') 'WJC', INDF
         CALL WRITMS 
     >        (3,    WJC(1, INDF), ND, NAME, -1, IDUMMY, IERR)
      ENDDO

      DO INDF=1, NF
         EMCOLI(INDF)=EMCOLI(INDF)/FWTEST(INDF)
      ENDDO
      CALL WRITMS (3,EMCOLI,NF,'EMCOLI  ',-1, IDUMMY, IERR)


      CALL WRITMS(3, TAUROSS,ND, 'TAUROSS ', -1, IDUMMY, IERR)
      CALL WRITMS(3, OPAROSS,ND, 'OPAROSS ', -1, IDUMMY, IERR)

      CALL WRITMS(3, ARAD, ND-1, 'ARAD    ', -1, IDUMMY, IERR)
      CALL WRITMS(3, ACONT, ND-1, 'ACONT   ', -1, IDUMMY, IERR)
      CALL WRITMS(3, ATHOM, ND-1, 'ATHOM   ', -1, IDUMMY, IERR)
C      CALL WRITMS(3, RHO,    ND, 'RHO     ', -1, IDUMMY, IERR)
C      CALL WRITMS(3, XMU,    ND, 'XMU     ', -1, IDUMMY, IERR)
C***  save radiative acceleration and fractions in internal units
C      (this will be used in STEAL->HYDROSOLVE)
      IF (bKALPHA) THEN
        CALL WRITMS(3, ALPHAF, ND, 'ALPHAF  ', -1, IDUMMY, IERR)
      ENDIF



C***  CALCULATE AND SAVE LOCAL RADIATIVE ENERGY LOSS IN CGS
      DO L=1, ND
         FTCOLI(L) = FTCOLI(L)/RSTAR
      ENDDO
      CALL WRITMS(3, FTCOLI, ND, 'FTCOLI  ', -1, IDUMMY, IERR)

C***  IRON: NORMALIZE AND SAVE MEAN INTENSITIES 'XJFEMEAN'
C           AND DIAGONAL WEIGHTS 'WFELOW' AND 'WFENUP'
      DO INDFE = 1, LASTFE
         DO L=1, ND
            XJFEMEAN(L, INDFE) = XJFEMEAN(L, INDFE) / SIGMAINT(INDFE)
     >                           /RADIUS(L)/RADIUS(L)
            WFELOW(L, INDFE)   = WFELOW(L, INDFE) / SIGMAINT(INDFE)
            WFENUP(L, INDFE)   = WFENUP(L, INDFE) / SIGMAINT(INDFE)
            FTFE(L,INDFE)      = FTFE(L,INDFE)/RSTAR
         ENDDO
         IND = LASTIND + INDFE
         IF (IND <= 9999) THEN
           WRITE (NAME,'(A3,I4,A1)') 'XJL', IND, ' '
         ELSE
           WRITE (NAME,'(A3,I5)') 'XJL', IND
         ENDIF
         CALL WRITMS 
     >        (3, XJFEMEAN(1, INDFE), ND, NAME, -1, IDUMMY, IERR)
         WRITE (NAME,'(A3,I4)') 'WFL', INDFE
         CALL WRITMS 
     >        (3, WFELOW(1, INDFE), ND, NAME, -1, IDUMMY, IERR)
         WRITE (NAME,'(A3,I4)') 'WFU', INDFE
         CALL WRITMS 
     >        (3, WFENUP(1, INDFE), ND, NAME, -1, IDUMMY, IERR)
         WRITE (NAME,'(A3,I4)') 'FTF', INDFE
         CALL WRITMS
     >        (3, FTFE(1, INDFE), ND, NAME, -1, IDUMMY, IERR)
      ENDDO

C***  Radiative Rates for Iron superlevels
      IF (LASTFE .GT. 0) THEN
         DO L=1, ND
            WRITE (NAME, '(A5, I3)') 'FERLU', L
            CALL WRITMS 
     >           (3, FERATLU(1,L), LASTFE, NAME, -1, IDUMMY, IERR)
            WRITE (NAME, '(A5, I3)') 'FERUL', L
            CALL WRITMS 
     >           (3, FERATUL(1,L), LASTFE, NAME, -1, IDUMMY, IERR)
         ENDDO
      ENDIF


C***  Save the integrated Flux of all Lines
C***  Note that HTOTL and NTOTL are on interstices (i.e. only valid ND-1 entries)
C!!!  Unsinnige Laenge beibehalten wegen Inkompatibilitaet mit aelteren Modellen
      CALL WRITMS(3, HTOTL, ND, 'HTOTL   ', -1, IDUMMY, IERR)
      CALL WRITMS(3,XJTOTL, ND, 'JTOTL   ', -1, IDUMMY, IERR)
      CALL WRITMS(3,XKTOTL, ND, 'KTOTL   ', -1, IDUMMY, IERR)
      CALL WRITMS(3,XNTOTL, ND, 'NTOTL   ', -1, IDUMMY, IERR)

C***  Calculate and Save Rosseland Opacity at the inner Boundary
      OPARND = DBDTINT / DBDTOPAINT
C      WRITE (0,*) 'OPAROSS(ND) ', OPARND, OPAROSS(ND)
      CALL WRITMS(3, OPARND, 1, 'OPARND  ', -1, IDUMMY, IERR)
      
C***  Save Eddington Flux at inner boundary 
C***  calculated in FREQUINT
      CALL WRITMS(3, HTOTND, 1, 'HTOTND  ', -1, IDUMMY, IERR)

C***  Diffusion/LTE correction for Eddington Flux at inner boundary
C***  HTOTNDCOR = int( B_nu(T_ND)/2 - h_nu * J_nu,ND , d-nu )
C***  calculated in FREQUINT
      CALL WRITMS(3, HTOTNDCOR, 1, 'HTNDCOR ', -1, IDUMMY, IERR)

C***  Write NCOLIP
      CALL WRITMS(3, NCOLIP,    1, 'NCOLIP  ', -1, IDUMMY, IERR)

C***  Unsoeld-Lucy Terms for TEMPEQ
      CALL WRITMS (3,OPASMEAN,     ND, 'OPASMEAN', -1, IDUMMY, IERR)
      CALL WRITMS (3,OPAJMEAN,     ND, 'OPAJMEAN', -1, IDUMMY, IERR)
      CALL WRITMS (3,OPAPMEAN,     ND, 'OPAPMEAN', -1, IDUMMY, IERR)      
      CALL WRITMS (3,QFJMEAN,      ND, 'QFJMEAN ', -1, IDUMMY, IERR)
      CALL WRITMS (3,OPASMEANTC,   ND, 'OPASMTC ', -1, IDUMMY, IERR)
      CALL WRITMS (3,OPAJMEANTC,   ND, 'OPAJMTC ', -1, IDUMMY, IERR)
      CALL WRITMS (3,SMEAN,        ND, 'SMEAN   ', -1, IDUMMY, IERR)
      CALL WRITMS (3,QOPAHMEAN,  ND-1, 'QOPAHMEA', -1, IDUMMY, IERR)
      CALL WRITMS (3,EDDIHOUTJMEAN, 1, 'EDDIHOJM', -1, IDUMMY, IERR)
      CALL WRITMS (3,OPALAMBDAMEAN,ND, 'OPALMEAN', -1, IDUMMY, IERR)
      
      IF (BCOLIP) THEN
        CALL WRITMS(3, EPSGMAX, ND-1, 'EPSGMAX ', -1, IDUMMY, IERR)
      ENDIF

      CALL WRITMS (3,TOTOUT,1, 'TOTOUT  ',-1, IDUMMY, IERR)

      RETURN
      END
