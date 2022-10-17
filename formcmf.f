       SUBROUTINE FORMCMF (TEFF, U, UK, NDDIM, NFLDIM,
     $        Z, OPA, ETA, ETANOTH, ETACK, ETACCK, OPACK, ETANCK, 
     $        XCMFRED, XCMFBLUE, DXCMF, XMAX, RMAX, VMAXDU, 
     $        THOMSON, THOMCK, OPAL, ETAL, 
     $        RADIUS, ND, NP, P,
     $        TA, TB, TC, UB, GA, H, QQ, S, V, VA, VB, PP,
     $        BCORE, DBDR, VDU, GRADI, VDOP,
     $        W0, NBLINE, DELXLAP, XRANGERED, XRANGEBLUE,
     $        XJNUE, OPAK, ETAK,
     $        XPLOT, YPLOT, XJC, NFDIM, XJCIND,
     $        LINE, MODHEAD, JOBNUM, NREDMAX, WREDI, 
     $        WREDI0, VDUEL, BWESEX,
     $        XLAMBDA, XLAM, T, RNE, POPNUM, ENTOT, RSTAR, MAINPRO, 
     $        MAINLEV, NOM, KODAT, NDIM, N, MAXATOM, LEVEL, NCHARG, 
     $        WEIGHT, ELEVEL, EION, EINST, ALPHA, SEXPO, ADDCON1, 
     $        ADDCON2, ADDCON3, IGAUNT, SIGMATHK, SEXPOK, EDGEK, NF, 
     $        KONTNUP, KONTLOW, LASTKON, XDATA, REDIS, LBREF, INDREF, 
     $        NPDIM, A, B, C, W, BX, WX, EDDI, IWARNJ0, 
     $        IWARN, FNUEC, LSDWL, ST, BELIFI, 
     >        DENSCON, FILLFAC, ENTOTDENS, bBIGBANDLIMIT,
C *** IRON: Additional Parameter used in SUBROUTINE FECHECK called from FORMCMF
     >        INDRB, IFRBSTA, IFRBEND, LASTFE,
     >        CLIGHT, VDOPFE, DXFE, XLAM0FE,
     >        INDFEACT, BFECHECK,
C *** IRON: Additional Parameter used in SUBROUTINE CMFFEOP called from FORMCMF
     >        SIGMAFE, OPAFE, ETAFE, 
     >        IFENUP, IFELOW, BIRONLINES,
     >        SIGMAACT, OPAFEI, ETAFEI,NFL,BAIRWAVELENGTHSET, 
     >        BAIRWAVELENGTH, VSIDU,
     >        DD_VDOPDU, NATOM, YSCRATCH, MAXLAP, IND_ELLINE, 
     >        bDDOPAFORMCMF, bDDFECONVOL)

C***********************************************************************
C***  CALLED FROM: FORMAL
C***  LINE RADIATION TRANSFER IN THE COMOVING FRAME, ACCOUNTING ITERATIVELY
C***     FOR PHOTON REDISTRIBUTION BY ELECTRON SCATTERING.
C***     THE RESULTING FREQUENCY-DEPENDING CONTINUUM EMISSIVITY IS STORED
C***     IN THE ARRAY ETACK(L,K). 
C***     THE ORIGINAL CONTINUUM EMISSIVITY IS STORED IN THE ARRAY
C***     ETACCK(L,K). THIS ARRAY IS USED IN PREPRAY TO CALCULATE THE RAY
C***     EMISSIVITY OF THE CONTINUUM
C***     THE CMF FREQUENCY OF INDEX K IS: XK = XCMFBLUE - (K-1) * DXCMF
C***     K RANGES FROM 1 TO NFL. L IS THE DEPTH INDEX.
C***
C***  The main part runs twice: 
C***     1) without lines (BWITHLINES = .FALSE.)
C***     2) with    lines (BWITHLINES = .TRUE.)
C***
C***  Output arrays:
C***     ETACCK(L,K) = continuum emissivity incl. redistribution, no lines
C***     ETACK (L,K) = continuum emissivity incl. redistributed line photons 
C***     OPACK (L,K) = continuum opacity    incl. Thomson term 
C***
C***  Only onternally used:
C***     ETAK  (L,K) = true emissivity: continuum, in 2) + sum of lines 
C***     ETANCK(L,K) = internal storage of emissivity without Thomson
C***********************************************************************

      INTEGER, INTENT(IN) :: ND, NATOM, NDDIM, NFLDIM, MAXLAP

      INTEGER, DIMENSION(MAXLAP) :: IND_ELLINE
      REAL, DIMENSION(ND) :: RADIUS, OPA, ETA, ETANOTH, THOMSON
      REAL, DIMENSION (NDDIM, NATOM), INTENT(IN) :: DD_VDOPDU
      DIMENSION ETACK(NDDIM,NFLDIM), OPACK(NDDIM,NFLDIM)
      DIMENSION ETACCK(NDDIM,NFLDIM)
      DIMENSION ETANCK(NDDIM,NFLDIM),THOMCK(NDDIM,NFLDIM)
      DIMENSION BCORE(NFLDIM), DBDR(NFLDIM)
      DIMENSION XLAMBDA(NFDIM), XJC(ND,NF)
      DIMENSION T(NDDIM), ENTOT(NDDIM), ENTOTDENS(NDDIM)
      DIMENSION W0(ND)
      DIMENSION U(ND,NP), UK(ND,NP), V(ND,NP), Z(ND,NP)
      DIMENSION PP(ND), VDU(ND), GRADI(ND)
      DIMENSION XJNUE(NDDIM,NFLDIM)
      DIMENSION OPAL(NDDIM,NBLINE), ETAL(NDDIM,NBLINE)
      DIMENSION DELXLAP(NBLINE)
      DIMENSION OPAK(NDDIM,NFLDIM), ETAK(NDDIM,NFLDIM)
      DIMENSION XPLOT(NFLDIM), YPLOT(NFLDIM), XJCIND(ND)
      DIMENSION INDFEACT(LASTFE) 
      REAL, DIMENSION(NFLDIM) :: YSCRATCH
C***  ARRAYS TO HANDLE DEPTH-DEPENDENT REDISTRIBUTION INTEGRAL
      DIMENSION VDUEL(NDDIM), WREDI0(NDDIM), WREDI(NREDMAX,NDDIM)
      REAL, DIMENSION(NDDIM, NFLDIM) :: OPAFE, ETAFE
      REAL, DIMENSION(NDDIM, NFLDIM) :: OPACFEK, ETACFEK
      INTEGER, DIMENSION(MAXATOM) :: KODAT

      CHARACTER*80  HEADER, XTEXT, YTEXT
      CHARACTER MODHEAD*100
      LOGICAL REDIS, bBIGBANDLIMIT, 
     >        bDDOPAFORMCMF, bDDFECONVOL
      LOGICAL BFECHECK, BFEWING, BIRONLINES, BWITHLINES
      LOGICAL BAIRWAVELENGTHSET, BAIRWAVELENGTH
      
      INTEGER :: NAFE

      DIMENSION QQ(ND), S(ND)

      DIMENSION DENSCON(ND),FILLFAC(ND)

C***  1 / SQRT(PI)
      DATA WPIINV / 0.564189583549 /

C***  TEST PLOT FACILITY: J-NUE AT DEPTH POINT LJPLOT (0 = DISABLED)
      LJPLOT = 0

C-----------------------------------------------------------------
C***  MAXIMUM NUMBER OF ITERATION
      MAXITER = 100
C-----------------------------------------------------------------
C***  EPSILON CONVERGENCE-CRITERION
      EPS = 0.001
C-----------------------------------------------------------------
C***  MAXIMUM HALF-BANDWIDTH FOR CMF LINE TRANSFER IN DOPPLER UNITS
      CMFBAND = 4.5
C-----------------------------------------------------------------
C***  SPACING OF CMF FREQUENCY POINTS
      DXCMF = 0.3
C-----------------------------------------------------------------
C***  HALF BANDWIDTH OF ELECTRON REDISTRIBUTION INTEGRAL
C***     IN UNITS OF THE ELECTRON DOPPLER VELOCITY
      BWES = 2. * BWESEX
C-----------------------------------------------------------------
C***  THE FREQUENCY RANGES (BOTH, CMF AND OBS. FRAME) ARE EXTENDED
C***      BY ONE HALF BANDWIDTH, MULTIPLIED WITH THE FOLLOWING FACTORS:
      BWEXRED = 2.
      BWEXBLU = 1.5
C-----------------------------------------------------------------
 
C***  Switch for adding lines (incl. iron bands) - else only continuum 
      BWITHLINES = .FALSE.
      WRITE (0,'(A)') "FORMCMF: Continuum Thomson redistribution"

C***  CALCULATE THE RANDOM VELOCITY OF THE ELECTRONS IN DOPPLER UNITS
C***       ATOMIC MASS OF HELIUM
           ATMASS = 4.
      DO 27 L=1, ND
        VTHHE2 = 0.01651 * T(L) / ATMASS
        VMICRO2 = VDOP * VDOP - VTHHE2
        VTHEL2 = 30.3165 * T(L)
        VDUEL(L) = SQRT (VMICRO2 + VTHEL2) / VDOP
   27 CONTINUE
C----------------------------------------------------------------------
C***  Prepare Clump density
      DO L=1, ND
         ENTOTDENS(L) = ENTOT(L) * DENSCON(L)
      ENDDO

      ALN = ALOG (1. - VDOP / CLIGHT)
      XLAMLN = ALOG (XLAM)

C***  DEFINE CMF FREQUENCY RANGE
C***  Since for the continuum, electron redistribution is now calculated
C***  in any case,, the range is always extended independent from the
C***  REDIS option - wrh  3-Jul-2020  
C***  The extension for electron reditribution is also returned to
C***  the calling FORMAL
cc      IF (REDIS) THEN
        FREMINEX = BWEXRED * BWES * VDUEL(ND)
        FREMAXEX = BWEXBLU * BWES * VDUEL(ND)
cc      ELSE
cc        FREMINEX = .0
cc        FREMAXEX = .0
cc      ENDIF

C***  The Band (XCMFRED,XBLUE) must cover the whole range needed in 
C***      OBSFRAM --> ZONEINT
C***  It is not clear to me why *twice* VDU(1) is needed (wrh 27-May-2008)
      BIGGERBAND =  AMAX1(XMAX,CMFBAND)
C***  Range should not grow by more than 0.2c due to XMAX (Tomer, 8.1.2015)
C***  This range cutting can be suppressed by the CARDS option NO-BIGBANDCUT     
      IF (bBIGBANDLIMIT) THEN
        BIGGERBAND = AMIN1(BIGGERBAND, 0.2 * CLIGHT / VDOP)
      ENDIF  
      XCMFBLUE = DELXLAP(NBLINE) + BIGGERBAND + FREMAXEX + 2.1*VMAXDU
     >                                                   + VSIDU
      XCMFRED  = DELXLAP(1)      - BIGGERBAND - FREMINEX - 2.1*VMAXDU
     >                                                   - VSIDU

C***  The range may be extended by the RANGE option (but not shrinked)
      IF (XRANGERED .NE. XRANGEBLUE) THEN
         XCMFBLUE = AMAX1 (XCMFBLUE, XRANGEBLUE+ BIGGERBAND + 
     >                              FREMAXEX + 1.1*VMAXDU) 
         XCMFRED  = AMIN1 (XCMFRED , XRANGERED - BIGGERBAND - 
     >                              FREMAXEX - 1.1*VMAXDU) 
         XLAMBDAC = XLAM * ( EXP( XRANGEBLUE * ALN ) 
     >               +        EXP( XRANGERED  * ALN ) ) / 2. 
      ELSE
         XLAMBDAC = XLAM * ( EXP( XCMFBLUE * ALN ) 
     >               +        EXP( XCMFRED  * ALN ) ) / 2. 
      ENDIF

C***  Calculate Iron OPAs/ETAs for Air, lower limit = 3200 Ang
      BAIRWAVELENGTH = ( (.NOT. BAIRWAVELENGTHSET) .AND. 
     >      ( XLAMBDAC .GE. 3200. .AND. XLAMBDAC .LE. 10000. ) )
     >  .OR. ( BAIRWAVELENGTHSET .AND. BAIRWAVELENGTH ) 

      NFL = 1 + (XCMFBLUE - XCMFRED) / DXCMF

      IF (NFL .GT. NFLDIM) THEN
         CALL REMARK ('NFLDIM INSUFFICIENT')
         WRITE (0,*) 'CMFBLUE, RED: ', XCMFBLUE, XCMFRED
         WRITE (0,*) 'XMAX, TOMERTHRESH: ', XMAX, 0.2 * CLIGHT / VDOP
   99    FORMAT (' NFLDIM GIVEN:', I8, ', BUT REQIERED:', I8)
         WRITE (0,99) NFLDIM, NFL
         PRINT 99, NFLDIM, NFL
         STOP 'ERROR'
         ENDIF

C*********************************************************
C***  DEFINE WEIGHTS FOR REDISTRIBUTION INTEGRAL
C***     THE REDISTRIBUTION FUNCTION IS THE ANGLE-AVERAGED
C***      REDISTRIBUTION FUNCTION FOR ELECTRON SCATTERING, R(e,A)
C***      (cf. MIHALAS: STELLAR ATMOSPHERES, P. 432)
C***      R(e,a) IS CODED AS FUNCTION FIERFC
C*********************************************************
C***      NUMBER OF INTEGRATION STEPS TO ONE SIDE : NREDI
      NREDI = NINT (BWES * 2. * VDUEL(ND) / DXCMF)
      IF (NREDI .GT. NREDMAX) THEN
        CALL REMARK ('DIMENSION NREDMAX INSUFFICIENT')
        WRITE (0,*) 'NREDI=',NREDI, '  NREDMAX=',NREDMAX
        STOP 'ERROR'
      ENDIF
      DO 28 L=1, ND
        WREDI0(L) = FIERFC(0.)
        SUM = WREDI0(L)
        DO 40 I=1,NREDI
          X = I * DXCMF/ (2. * VDUEL(L))
          WREDI(I,L) = FIERFC(X)
          SUM = SUM + 2. * WREDI(I,L)
   40   CONTINUE
C***  NORMALISATION
        WREDI0(L) = WREDI0(L) / SUM
        DO 41 I=1,NREDI
   41     WREDI(I,L) = WREDI(I,L) / SUM
   28 CONTINUE
 
C***  VERSION OF OUTER BOUNDARY CONDITION : NO INCIDENT RADIATION
      XIMINUS=0.
C***  FREQUENCY DERIVATIVE OF IMINUS
      DXI=0.

C***  IRON: WIDTH OF RED LINE WING IN IRON-DELTAX-UNITS
      IF (BIRONLINES) THEN
         DFEDUMMY= 1.
      ENDIF
  
C***  IRON: SET LOWEST LINE TO '1', NO ACTIVE FE-LINES 
      INDFEACT(1) = 1
      MAXFEACT    = 0
      BFECHECK    = .FALSE.
      BFEWING     = .FALSE.
      VDOPMAX = MAXVAL(DD_VDOPDU)

C***  CALCULATION OF OPAK, ETAK = ALL BLENDING LINES
      WRITE (0,*) 'Calculating total opacities...'
      DO 12 K=1,NFL
        XK = XCMFBLUE - (K-1) * DXCMF
        XLAMREF = EXP (XLAMLN+ XK * ALN)
        DO 15 L=1,ND
          OPAK(L,K) = 0.
          ETAK(L,K) = 0.
   15   CONTINUE
C***  LOOP EINSCHRAENKEN?
        DO 14 IBLEND=1,NBLINE
          XKBL = XK - DELXLAP(IBLEND)
C***      Worst Doppler width appoximation          
          XKLBMAXDU = XKBL*XKBL / VDOPMAX/VDOPMAX
C***      Do not evaluate PHI for more than 10 Dopplerwidths          
          IF (XKLBMAXDU > 4.5) CYCLE
          IF (bDDOPAFORMCMF) THEN
            NA = IND_ELLINE(IBLEND)          
          ELSE
            PHIX = EXP(-XKBL*XKBL) * WPIINV
          ENDIF
          DO 13 L=1,ND
            IF (bDDOPAFORMCMF) THEN
              XKBLL = XKBL / DD_VDOPDU(L,NA)
              PHIX = EXP(-XKBLL*XKBLL) * WPIINV / DD_VDOPDU(L,NA)
            ENDIF
            OPAK(L,K) = OPAK(L,K) + OPAL(L,IBLEND) * PHIX
            ETAK(L,K) = ETAK(L,K) + ETAL(L,IBLEND) * PHIX
              
   13     CONTINUE
   14   CONTINUE

       
        IF (BIRONLINES) THEN
C***       IRON: CHECK FOR ACTIVE BOUND-BOUND TRANSITIONS OF GENERIC ION
           CALL FECHECK (XLAMREF, INDRB, IFRBSTA, IFRBEND, LASTFE,
     >                    CLIGHT, VDOPFE, DXFE, XLAM0FE,
     >                    INDFEACT, MAXFEACT, BFECHECK, BFEWING,
     >                    DFEDUMMY)

C***       IRON: NON-LTE IRON OPACITY AT GIVEN FREQUENCY FOR ALL DEPTH POINTS
           CALL CMFFEOP (XLAMREF, ND, N, INDFEACT, MAXFEACT, LASTFE,
     >                    SIGMAFE, OPAFE(1,K), ETAFE(1,K), INDRB,
     >                    IFRBSTA, IFRBEND, IFENUP, IFELOW,
     >                    CLIGHT, VDOPFE, DXFE, XLAM0FE,
     >                    ELEVEL, WEIGHT, RSTAR, POPNUM, ENTOT,
     >                    SIGMAACT, OPAFEI, ETAFEI, T, 1, TEFF) 

        ENDIF

C***    Inner Boundary values at each CMF frequency (diffusion approx.):
C***    BCORE = B_nue    DBDR = dB/dr  
C***    BCOREL and DBDRL re stored in the vectors BCORE and DBDR. 
C***    The local variables are used later in CALL BACKJCU
        CALL DIFFUS (XLAMREF, T, RADIUS, ND,
     >               BCOREL, DBDRL, DUMMY, DUMMY, .TRUE.)
        BCORE(K) = BCOREL
        DBDR(K) = DBDRL
 
        CALL COOP (XLAMREF,ND,T,RNE,POPNUM,ENTOTDENS,RSTAR,
     $             OPA,ETANOTH,THOMSON,IWARN,MAINPRO,MAINLEV,NOM,KODAT,
     $             NDIM,N,MAXATOM,LEVEL,NCHARG,WEIGHT,ELEVEL,EION,EINST,
     $             ALPHA,SEXPO,
     $             ADDCON1, ADDCON2, ADDCON3,
     $             IGAUNT,SIGMATHK,SEXPOK,EDGEK,0,NF,DUMMY,
     $             RADIUS,KONTNUP,KONTLOW,LASTKON,XDATA)
C***    Opacity is scaled down from clump to average value
        DO L=1, ND
          OPA(L) = OPA(L) * FILLFAC(L)
          ETANOTH(L) = ETANOTH(L) * FILLFAC(L)
        ENDDO

C***    U IS NEEDED AS BLUE WING BOUNDARY FOR THE RADIATION TRANSFER;
C***    IN CASE OF DWL PLOT (LSDWL > 0) FNUEC IS ALSO NEEDED
        IF (K .EQ. 1 .OR. LSDWL .GT. 0) THEN
C+++       abarnisk 06/2004: "ELIMIN" ersetzt hier "BACKJCU"-Aufruf 
C***       Note: EDDI is output parameter, only needs dimensioned space
           CALL ELIMIN (XLAMREF,FNUEC,DUMMY2,U,Z,A,B,C,W,BX,WX,XJCIND,
     $          RADIUS,P,BCORE,DBDR,OPA,ETANOTH,THOMSON,EDDI,ND,NP,NPDIM,
     $          ENTOT,INDREF, IDUMMY, ST, BELIFI, 0)
           IF (LSDWL .GT. 0) GOTO 31
        ELSE
            CALL BACKJC (XJC, ND, NF, XJCIND, XLAMBDA, XLAMREF, RADIUS)
        ENDIF

C***    FILL FREQUENCY-INDEXED ARRAYS 
        DO 17 L=1,ND
          ETANCK(L,K) = ETANOTH(L)
          OPACK (L,K) = OPA(L)
          THOMCK(L,K) = THOMSON(L)
C***      Initialize ETACK with coherent scattering
          ETACK(L,K) = ETANOTH(L) + OPA(L) * THOMSON(L) * XJCIND(L)
   17   CONTINUE

   12 CONTINUE

      IF (BIRONLINES .AND. bDDFECONVOL) THEN
C***      Convolve iron line opacities and emissivities with
C***       the correct DD_VDOPDU (if larger than VDOPFE)
        WRITE (0,*) 'Convolving iron line opacities with specified'
     >              // ' VDOP stratification...'
        NAFE = KODAT(26)
        CALL CONVOLOPAFE(OPAFE, YSCRATCH, NDDIM, ND, NFL, DXCMF,
     >                   VDOP, VDOPFE, DD_VDOPDU(1,NAFE))
        CALL CONVOLOPAFE(ETAFE, YSCRATCH, NDDIM, ND, NFL, DXCMF,
     >                   VDOP, VDOPFE, DD_VDOPDU(1,NAFE))
      ENDIF   
      WRITE(0,'(A)') 'STARTING FORMCMF WITH CONTINUUM ONLY'
   
C*** Originally, NOREDIS would cause to exit the routine here, adding 
C*** only coherent Thomson term to ETA
C*** Now ETA and J are calculated iteratively, with the continuum Thomson 
C*** term redistributed
C*** Tomer, 27/5/14
     
 301  CONTINUE 
C*** ITERATION LOOP =================================================
      DO 30 ITER=1, MAXITER

C***  SET XJNUE = 0.
      DO 16 L=1,ND
      DO 16 K=1,NFL
   16   XJNUE(L,K) = 0.


C***  LOOP OVER IMPACT PARAMETERS
      DO 11 JP=1,NP-1

      LMAX=MIN0(NP+1-JP,ND)
      LZ=LMAX-1

C***  *****************************************************************
C***  The following formalism is a difference scheme for the ray-by-ray 
C***  radiative transfer, as descibed in my Thesis (Hamann 1979). 
C***  U and V vectors are the intensity-like and flux-like Feautrier
C***  variables
C***  - In program COLI we have meanwhile replaced this method by a
C***    short-characteristic integration, which guarentees positive
C***    intensities. One might consider to apply shortchar here, too
C***    (wrh, 16-Apr-2013)
C***  *****************************************************************

C***  INITIALIZE UK = FEAUTRIER U AT BLUE WING
      DO 6 L=1,LMAX
    6   UK(L,JP) = U(L,JP)

C***  INITIALIZE FEAUTRIER V AS DERIVATIVE OF U
       DO 5 L=1,LZ
        X=0.5*(OPACK(L,1)+OPACK(L+1,1))
        V(L,JP)=(UK(L+1,JP)-UK(L,JP)) / (X*(Z(L,JP)-Z(L+1,JP)))
    5  CONTINUE

C***  GENERATE WEIGHTS FOR INTEGRATION OF THE 0.-MOMENT Jnue
      CALL GENW0 (JP, LMAX, ND, NP, Z, RADIUS, P, W0)

C***  PP(L) = VELOCITY GRADIENT, PROJECTED ON THE PRESENT RAY J, /DELTAX
C***  GDU = GRADI, CONVERTED IN DOPPLER UNITS
C***  VDU = VELO IN DOPPLER UNITS, WAS CONVERTED ALREADY IN FORMAL 
      DO 4 L=1,LMAX
      RL = RADIUS(L)
      Y  = Z(L,JP) / RL
      YY = Y*Y
      GDU = GRADI(L) / VDOP
    4 PP(L) = (YY * GDU + (1.-YY) * VDU(L)/RL) / DXCMF


C***  LOOP FOR ALL (ORIGINAL) FREQUENCY POINTS
C***  THE CURRENT INDEX K REFERS TO THE RESULTING U,V , WHERE ALL QUANTITIES
C***  ARE TAKEN. THE FREQUENCY DERIVATIVE STRETCHES BACK TO K-1

      DO 1 K=1,NFL
         
C+++ abarnisk 06/2004: Verbesserung des Startwertes fuer U (bzw. UK) durch
C+++ mehrfaches "auf der Stelle iterieren" von "CMFSET_FORMAL" mit OPACK(1,1)
C+++ ETACK(1,1). Laenge der Schleife willkuerlich gewaehlt!
         
         IF (K.EQ.1) THEN 
            DO 311 I=1,100      ! oder 50,20,10???
               CALL CMFSET_FORMAL(Z(1,JP),
     $           ND,LMAX,TA,TB,TC,UB,VA,VB,GA,H,S,
     $           OPACK(1,K),ETACK(1,K),PP,BCORE(K),DBDR(K),
     $           RADIUS,XIMINUS,DXI,OPAK(1,K),ETAK(1,K),
     >           OPAFE(1,K),ETAFE(1,K),BWITHLINES)
               CALL VMALV (VA,VB,V(1,JP),QQ,LMAX)
               DO 501 L=1, LMAX
 501              QQ(L) = QQ(L) + S(L)
                  CALL MDV (UB,UK(1,JP),LMAX)
                  CALL VADD (UK(1,JP),QQ,LMAX)
                  CALL INVTRI (TA,TB,TC,UK(1,JP),LMAX)
                  CALL MDV (H,V(1,JP),LZ)
                  CALL GMALU (GA,UK(1,JP),S,LMAX)
                  CALL VADD (V(1,JP),S,LZ)
 311           CONTINUE
               GOTO 2
         ENDIF
            
            CALL CMFSET_FORMAL(Z(1,JP),
     $           ND,LMAX,TA,TB,TC,UB,VA,VB,GA,H,S,
     $           OPACK(1,K),ETACK(1,K),PP,BCORE(K),DBDR(K),
     $           RADIUS,XIMINUS,DXI,OPAK(1,K),ETAK(1,K),
     >           OPAFE(1,K),ETAFE(1,K),BWITHLINES)
            
      CALL VMALV (VA,VB,V(1,JP),QQ,LMAX)

C***  THE FOLLOWING INLINING SAVES CPU TIME:
C***   (NOTE: IT SEEMS THAT INLINING OF THE OTHER CALLS WITH MATRIX COLUMNS
C***          AS FORMAL PARAMETERS DOES NOT SAVE TIME.)
C!!      CALL VADD (QQ,S,LMAX)
      DO 50 L=1, LMAX
   50 QQ(L) = QQ(L) + S(L)

      CALL MDV (UB,UK(1,JP),LMAX)

      CALL VADD (UK(1,JP),QQ,LMAX)

      CALL INVTRI (TA,TB,TC,UK(1,JP),LMAX)
C***  NOW U IS THE FIELD AT THE NEW INDEX K
      CALL MDV (H,V(1,JP),LZ)
      CALL GMALU (GA,UK(1,JP),S,LMAX)

      CALL VADD (V(1,JP),S,LZ)
 
2     CONTINUE
C***  end of the Feautrier difference scheme ************************

C***  ADDING THE NEW UK(L) TO THE ANGEL-AVERAGED INTENSITY XJNUE
      DO 20 L=1,LMAX
20    XJNUE(L,K) = XJNUE(L,K) + UK(L,JP) * W0(L)

    1 CONTINUE
   11 CONTINUE

C***  UPDATING THE FREQUENZ DEPENDENT CONTINUUM EMISSIVITY
C***  ANGLE AVERAGED FREQUENCY-REDISTRIBUTION
C***  CORMAX = MAX. RELATIVE CORRECTION OF ETA
      CORMAX = 0.
      DO 18 L=1,ND
        DO 18 K=1,NFL
C***  INITIALIZING OF REDISTRIBUTION INTEGRAL
          OPATH = OPACK(L,K) * THOMCK(L,K)
          ETACKNEW = ETANCK(L,K) + OPATH * XJNUE(L,K) * WREDI0(L)
          DO 19 I=1,NREDI
            KPI = K+I
            IF (KPI .GT. NFL) KPI = NFL
            KMI = K-I
            IF (KMI .LT.1 ) KMI = 1
            ETACKNEW = ETACKNEW + OPATH * (XJNUE(L,KPI) +
     >                 XJNUE(L,KMI)) * WREDI(I,L)
   19     CONTINUE

            RELCOR = ABS (ETACKNEW / ETACK(L,K) - 1.) 
            IF (RELCOR .GT. CORMAX) CORMAX = RELCOR
          ETACK(L,K) = ETACKNEW
   18 CONTINUE

C***  TEST PLOT FACILITY ********************************************
      IF (LJPLOT .GT. 0) THEN

C***  DEFINE THE Y-VECTOR
      DO K=1, NFL
         YPLOT(K) = XJNUE(LJPLOT,K) 
      ENDDO

      IF (ITER .EQ. 1) THEN
C***  DEFINE X-VECTOR: Wavelength
      DO  K=1, NFL
        XK = XCMFBLUE - (K-1) * DXCMF
          XPLOT(K) = EXP (XLAMLN+ XK * ALN)
      ENDDO
      XTEXT   = '\CENTER\#l# [\A]'

C***  Label Y-AXIS
      YTEXT = 'INTENSITY J&T#n#&M'
      WRITE (HEADER, 10) BWITHLINES, LJPLOT, MODHEAD(13:33)
   10 FORMAT ('J&T#n#&M:  WITHLINES=', L1, '  L=', I2,
     >        '   MODEL ', A )

        CALL PLOTANFS (1, HEADER, HEADER, XTEXT, YTEXT,
     $             0., 0., 0., 0., 0., 0.,
     $             0., 0., 0., 0., 0., 0.,
     $             XPLOT, YPLOT, NFL, 'COLOR=2')
      ELSE
        ISYMBOL = 5
        IF (ITER .EQ. 2) ISYMBOL = 9
        IF (ITER .EQ. 3) ISYMBOL = 10
        CALL PLOTCON (1,XPLOT,YPLOT,NFL,ISYMBOL)
      ENDIF
      ENDIF
C********************************************************************

      WRITE (0,'(A,I3,2X,A,F11.5)') 
     >      'FORMCMF : Iteration-Nr.:', ITER, 'Cormax=', CORMAX

C***  Check for convergence 
      IF (CORMAX .LT. EPS .AND. ITER .GT. 1) THEN
         WRITE(0,'(A,F6.4)')'Converged: Cormax <',  EPS
         IF (.NOT.BWITHLINES) THEN 
C***        This was the first run without lines
C***        -> save the continuum emissivity (whole array)
         DO L=1, ND
            DO K=1, NFL
               ETACCK(L,K) = ETACK(L,K)
            ENDDO
         ENDDO
C***     If NOREDIS=TRUE, then skip the Thomson redistribution of lines 
C***     and exit routine
            IF (.NOT. REDIS) THEN 
                WRITE(0,'(A)') 'Line redistribution due to electron '
     >                         // ' scattering skipped! (NOREDIS)'
                GOTO 31
            ENDIF
            BWITHLINES = .TRUE.
            WRITE(0,'(A)')'NOW FORMCMF WITH LINES'
            GOTO 301
         ELSE
C***     This was the second run **with** lines
            GOTO 31
         ENDIF
      ENDIF

C***  Check for divvergence 
      IF (CORMAX .GT. 1.E6 .AND. ITER .GT. 1 .AND. BWITHLINES) THEN
         WRITE (0,'(A)')  
     >   '*** WARNING: Iteration of electron-scattering redistribution',
     >   ' terminated early because of apparent divergence'
         WRITE (*,'(2A)')  
     >   '*** WARNING: Iteration of electron-scattering redistribution',
     >   ' terminated early because of apparent divergence'
          GOTO 31

         STOP '*** FATAL ERROR in FORMCMF'
      ENDIF

C***  END OF ITERRATION LOOP ===================================
30    CONTINUE
      WRITE (*,'(A,I3,2A)')  'MAX. NUMBER (',MAXITER,') OF ITERATIONS ',
     >                       'IN FORMCMF INSUFFICIENT'
      WRITE (0,'(A,I3,2A)')  'MAX. NUMBER (',MAXITER,') OF ITERATIONS ',
     >                       'IN FORMCMF INSUFFICIENT'

31    CONTINUE

      RETURN
      END
 
 





