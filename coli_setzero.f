      SUBROUTINE COLI_SETZERO(ND, NDDIM, NIT, NPDIM, NFDIM, MAXFEIND, 
     >             MAXLIN, DBDTINT, DBDTOPAINT, EDDIHOUTJMEAN, 
     >             HTOTOUTMINUS, HTOTND, HTOTNDCOR,
     >             DBDTINT_M, DBDTOPAINT_M,
C***  with ND
     >             OPA, ETA,
     >             XJTOTL, HTOTL, XKTOTL, XNTOTL, ARAD, ACONT, ATHOM, 
     >             FTCOLI, OPAKOLD, ETAKOLD, OPAKNOTHO, ETAKNOTHO, 
     >             OPAO, THOMSONO, ETANOTHO, DJDSMO_OLD, DJDOMO_OLD, 
     >             OPASMEAN, QFJMEAN, OPAJMEAN, OPASMEANTC, OPAJMEANTC, 
     >             OPAPMEAN, SMEAN, QLFOLD, EPSGMAX, OPAROSS, OPALAMBDAMEAN,
C***  with ND-1
     >             QOPAHMEAN, HMEAN, QLHOLD, OPAKHOLD, 
C***  with NDDIM,NIT
     >             XJLOLD, XJLMO_OLD, EDDIFO, S_OLD, OPAK_OLD, EPSG,
C***  with NDDIM,NPDIM
     >             CWM0, CWM2, CWM1, CWM3,
C***  with NDDIM,NFDIM,
     >             WJC,
C***  with NDDIM,NPDIM,NIT
     >             XIPLUS_OLD, XIMINUS_OLD,
C***  with NDDIM-1, NIT; In COLI it is also NDDIM,NIT
     >             XHLOLD, XHLMO_OLD, EDDIGO, 
C***  with MAXFEIND, NDDIM
     >             XJFEMEAN, FERATLU, FERATUL, FTFE, WFELOW, WFENUP,
C***  with MAXLIN
     >             LIND, LINDS, WS,
C***  with MAXIND
     >             MAXIND, BLASERL,
C***  with NFDIM
     >             EMCOLI,
C***  no Arrays
     >             OPAMAX1, OPAMAX1_LAMBDA, IOPAMAX1_K)

C****************************************************************
C***  Presets all given variables to zero
C***    Called by COLI
C****************************************************************

      IMPLICIT NONE

      INTEGER, INTENT(IN) :: ND, NDDIM, NIT, NPDIM, NFDIM,
     >                       MAXFEIND, MAXLIN, MAXIND
      REAL, DIMENSION(ND) :: XJTOTL, HTOTL, XKTOTL, XNTOTL, 
     >                       ARAD, ACONT, ATHOM, FTCOLI, 
     >                       OPAKOLD, ETAKOLD, OPAKNOTHO, ETAKNOTHO, 
     >                       OPAO, THOMSONO, ETANOTHO, 
     >                       DJDSMO_OLD, DJDOMO_OLD, OPASMEAN, 
     >                       QFJMEAN, OPAJMEAN, OPASMEANTC, 
     >                       OPAJMEANTC, OPAPMEAN, SMEAN, QLFOLD,
     >                       EPSGMAX, OPAROSS, OPALAMBDAMEAN
      REAL, DIMENSION(NDDIM) :: OPA, ETA
      REAL, DIMENSION(NDDIM,NIT) :: S_OLD, OPAK_OLD, EPSG
      REAL, DIMENSION(NDDIM,NPDIM) :: CWM0, CWM2, CWM1, CWM3
      REAL, DIMENSION(NDDIM, MAXFEIND) :: WFELOW, WFENUP, FTFE, 
     >                                    XJFEMEAN
      REAL, DIMENSION(NDDIM,NFDIM) :: WJC
      
      REAL, DIMENSION(ND-1) :: QOPAHMEAN, HMEAN, QLHOLD, OPAKHOLD

      REAL, DIMENSION(NDDIM,NIT) :: XJLOLD, XJLMO_OLD, EDDIFO, 
     >                              XHLOLD, XHLMO_OLD, EDDIGO

      REAL, DIMENSION(NDDIM,NPDIM,NIT) :: XIPLUS_OLD, XIMINUS_OLD

      REAL, DIMENSION(NFDIM) :: EMCOLI
      
      REAL, DIMENSION(MAXFEIND,NDDIM) :: FERATLU, FERATUL

      INTEGER, DIMENSION(MAXLIN) :: LIND, LINDS
      REAL, DIMENSION(MAXLIN) :: WS
      
      LOGICAL, DIMENSION(MAXIND) :: BLASERL
      
      INTEGER :: IOPAMAX1_K
      REAL :: OPAMAX1_LAMBDA, OPAMAX1, HTOTOUTMINUS, EDDIHOUTJMEAN,
     >        DBDTOPAINT_M, DBDTINT_M, DBDTOPAINT, DBDTINT,
     >        HTOTND, HTOTNDCOR

C***  Set Zero
      DBDTINT        = 0.
      DBDTOPAINT     = 0.
      DBDTINT_M      = 0.
      DBDTOPAINT_M   = 0.
      EDDIHOUTJMEAN  = 0.
      HTOTOUTMINUS   = 0.
      HTOTND         = 0.
      HTOTNDCOR      = 0.
 
      OPA            = 0.
      ETA            = 0.
      CWM0           = 0.
      CWM1           = 0.
      CWM2           = 0.
      CWM3           = 0.

      XJTOTL         = 0.
      HTOTL          = 0.
      XKTOTL         = 0.
      XNTOTL         = 0.
      ARAD           = 0.
      ACONT          = 0.
      ATHOM          = 0.
      FTCOLI         = 0.
      OPAKOLD        = 0.
      ETAKOLD        = 0.
      OPAKNOTHO      = 0.
      ETAKNOTHO      = 0.
      OPAO           = 0.
      THOMSONO       = 0.
      ETANOTHO       = 0.
      DJDSMO_OLD     = 0.
      DJDOMO_OLD     = 0.
      OPASMEAN       = 0.
      QFJMEAN        = 0.
      OPAJMEAN       = 0.
      OPASMEANTC     = 0.
      OPAJMEANTC     = 0.
      OPAPMEAN       = 0.
      OPAROSS        = 0.
      OPALAMBDAMEAN  = 0.
      SMEAN          = 0.
      QLFOLD         = 0.
      EPSGMAX        = 0.
 
      QOPAHMEAN      = 0.
      HMEAN          = 0.
      QLHOLD         = 0.
      OPAKHOLD       = 0.
 
      XJLOLD         = 0.
      XJLMO_OLD      = 0.
      EDDIFO         = 0.
      S_OLD          = 0.
      OPAK_OLD       = 0.
      EPSG           = 0.
C***  BLUE BOUNDARY CONDITION FOR SHORTRAY 
      XIPLUS_OLD     = 0.
      XIMINUS_OLD    = 0.

      XHLOLD         = 0.
      XHLMO_OLD      = 0.
      EDDIGO         = 0.
 
      XJFEMEAN       = 0.
      FERATLU        = 0.
      FERATUL        = 0.
      FTFE           = 0.
      WFELOW         = 0.
      WFENUP         = 0.

      LIND           = 0.
      LINDS          = 0.
      WS             = 0.
      BLASERL        = .FALSE.

      OPAMAX1        = 0.
      OPAMAX1_LAMBDA = 0.
      IOPAMAX1_K     = 0

      EMCOLI         = 0.
      

      RETURN
      END
