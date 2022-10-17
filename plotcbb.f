C***  MAIN PROGRAM   ******************************************************
      SUBROUTINE PLOTCBB
C*******************************************************************************
C***  STEAL-LIKE PROGRAM for generating Plots of the 
C***   Collision strength OMEGA for bound-bound transitions
C*******************************************************************************

C***  SET ARRAY DIMENSION PARAMETERS
C***  IRON: ADD GENERIC ION TO MAXATOM
      PARAMETER ( MAXATOM =          26 )
      PARAMETER ( NDIM    =        1560 )
      PARAMETER ( MAXAUTO =        2400 )
      PARAMETER ( MAXIND  =       20000 )
      PARAMETER ( MAXINDE = MAXAUTO+MAXIND )
      PARAMETER ( MAXNDR  =         200 )
      PARAMETER ( NDDIM   =          89 )
      PARAMETER ( NFLDIM  =          40 )
      PARAMETER ( MAXXDAT =          11 )

C***  MAXIMUM ION CHARGE WHICH MAY OCCUR (SEE ALSO SUBR. GAUNTFF)
      PARAMETER ( MAXION = 17 )

      PARAMETER ( NDIMP2  = NDIM + 2 )
      PARAMETER ( NFDIM   = 2*NDIM + 120 )
      PARAMETER ( MAXKONT = NFDIM/2 )
      PARAMETER ( MAXKODR = NDIM )

      PARAMETER (MAXLAP = 40)

C***  NUMBER OF POINTS PER PLOTTED DATASET
      PARAMETER ( MAXPLOTN   = 1000)

C***  HANDLING OF DIELECTRONIC RECOMBINATION / AUTOIONIZATION (SUBR. DATOM)
      COMMON / COMAUTO / LOWAUTO(MAXAUTO),WAUTO(MAXAUTO)
     $                  ,EAUTO(MAXAUTO),AAUTO(MAXAUTO),IONAUTO(MAXAUTO)
     $                  ,KRUDAUT(MAXAUTO),DRRATEN(NDIM)
     $                  ,RDIEL(NDIM),RAUTO(NDIM)
     $                  ,DRJLW(NDIM),DRJLWE(NDIM),DRLJW(NDIM)
     $                  ,KODRNUP(MAXKODR),KODRLOW(MAXKODR)

      COMMON // NCHARG(NDIM),WEIGHT(NDIM),ELEVEL(NDIM),EION(NDIM)
     $ ,MAINQN(NDIM),EINST(NDIM,NDIM),ENLTE(NDIM),NOM(NDIM)
     $ ,IONGRND(NDIM), ENOLD(NDIM)
     $ ,ABXYZ(MAXATOM),KODAT(MAXATOM),ATMASS(MAXATOM),STAGE(MAXATOM)
     $ ,NFIRST(MAXATOM),NLAST(MAXATOM)
     $ ,ALTESUM(4,NDIM)
     $ ,NCHARGDR(MAXNDR),ELEVDR(MAXNDR),NOMDR(MAXNDR)
     $ ,WEIGHTDR(MAXNDR),INDNUPDR(MAXAUTO)
      DIMENSION XDATA(MAXXDAT)
      DIMENSION SIGMATHK(MAXATOM),SEXPOK(MAXATOM),EDGEK(MAXATOM)

      DIMENSION XPLOT(MAXPLOTN), YPLOT(MAXPLOTN)

      COMMON /COMKONT/ SIGMAKI(NFDIM,MAXKONT)
     $ ,ALPHA(MAXKONT),SEXPO(MAXKONT)
     $ ,ADDCON1(MAXKONT), ADDCON2(MAXKONT), ADDCON3(MAXKONT)
     $ ,KONTNUP(MAXKONT),KONTLOW(MAXKONT)

      CHARACTER*8 IGAUNT(MAXKONT), KEYCBF(MAXKONT)

      COMMON /COMIND/  INDNUP(MAXIND),INDLOW(MAXIND),RUDLINE(MAXIND)
     $ ,COCO(4,MAXIND)

      CHARACTER LEVEL(NDIM)*10
      CHARACTER*10 MAINPRO(NDDIM),MAINLEV(NDDIM)
      CHARACTER*10 ELEMENT(MAXATOM)
      CHARACTER*10 IONDR(MAXNDR),LEVELDR(MAXNDR)
      CHARACTER*4 KEYCBB(MAXIND)
      CHARACTER*2 SYMBOL(MAXATOM)
      CHARACTER   KARTE*80, ACTPAR*10, CSTRING*8

C***  IRON: COMMON BLOCK FOR IRON-SPECIFIC DATA
C      INTEGER, PARAMETER :: INDEXMAX = 1E7, NFEREADMAX = 3E5    !std
      INTEGER, PARAMETER :: INDEXMAX = 4E7, NFEREADMAX = 5E5     !vd20
C      INTEGER, PARAMETER :: INDEXMAX = 1E8, NFEREADMAX = 6E5     !xxl

      PARAMETER ( MAXFEIND  =       1500 ) 
      COMMON /IRON/ INDRB(MAXFEIND), INDRF(MAXFEIND),
     >              IFRBSTA(MAXFEIND), IFRBEND(MAXFEIND),
     >              IFENUP(MAXFEIND), IFELOW(MAXFEIND),
     >              INDFEACT(MAXFEIND), SIGMAACT(MAXFEIND),
     >              SIGMAINT(MAXFEIND)

      DIMENSION FEDUMMY(NFEREADMAX)
      DIMENSION SIGMAFE(INDEXMAX)

      LOGICAL BFEMODEL

      CHARACTER TIM1*10

C***  C1 = H * C / K    ( CM * ANGSTROEM )
      DATA C1 / 1.4388 /

C***  Operating system:
      COMMON / COMOS / OPSYS
      CHARACTER*8 OPSYS

      CALL INSTALL
                                 
      IF (OPSYS .EQ. 'CRAY' .OR. OPSYS .EQ. 'SGI') THEN
        CALL CLOCK(TIM1)
      ELSE
        CALL TIME(TIM1)
      ENDIF

C***  DEFINE CHANNEL 42 FOR PLOTS
      KANAL = 42
      OPEN (KANAL, FILE='PLOT', STATUS='UNKNOWN')

C***  READING THE ATOMIC DATA FROM FILE DATOM
      CALL       DATOM (NDIM,N,LEVEL,NCHARG,WEIGHT,ELEVEL,EION,MAINQN,
     $                  EINST,ALPHA,SEXPO,
     $                  ADDCON1, ADDCON2, ADDCON3, 
     $                  IGAUNT,COCO,KEYCBB,ALTESUM,
     $                  INDNUP,INDLOW,LASTIND,MAXIND,MAXATOM,NATOM,
     $                  ELEMENT,SYMBOL,NOM,KODAT,ATMASS,STAGE,
     $                  SIGMATHK,SEXPOK,EDGEK,NFIRST,
     $                  NLAST,NAUTO,MAXAUTO,LOWAUTO,WAUTO,EAUTO,AAUTO,
     $                  IONAUTO,KRUDAUT,KONTNUP,KONTLOW,LASTKON,MAXKONT,
     $                  IONGRND,KODRNUP,KODRLOW,LASTKDR,MAXKODR,KEYCBF,
C***  IRON: ADDITIONAL PARAMETERS FOR IRON-GROUP LINE BLANKETING
     >                  'STEAL   ', INDEXMAX, NFEREADMAX, MAXFEIND,
     >                  LASTFE, SIGMAFE, INDRB, INDRF,
     >                  IFENUP, IFELOW, IFRBSTA, IFRBEND, FEDUMMY, 
     >                  VDOPFE, DXFE, XLAM0FE, SIGMAINT, BFEMODEL)
 

      CALL PLOTANF (42, '', 'Collision strength #W#&Tup-low&M',
     >              '\CENTER\log T / K', 
     >              '\CENTER\#W# [10&H-6&M cgs]', 
     >              0,0,0,0,0,0,
     >              0,0,0,0,0,0,
     >              XPLOT, YPLOT, 0, 5) 

      POSLABEL = 14
      BACKSPACE (42)
      WRITE (42,'(A)') '\FONT HELVET'
      WRITE (42,'(A)') 'END'


C*** Read range of line indices that shall be plotted from input 

      INDSTART = 0
      INDEND   = 0
     
      OPEN (1, FILE='CARDS', STATUS='UNKNOWN')
    8 READ (1, '(A)', END=10) KARTE
      IF (KARTE(:4) .NE. 'LINE') GOTO 8

      CALL SARGC (KARTE, NPAR)
      IF (NPAR .LT. 2) GOTO 990
      CALL SARGV (KARTE, 2, ACTPAR)
      READ (ACTPAR,'(I10)',ERR=990) INDSTART

      IF (NPAR .LT. 4) GOTO 10
      CALL SARGV (KARTE, 4, ACTPAR)
      READ (ACTPAR,'(I10)',ERR=990) INDEND

      GOTO 8

C**** Reading of input lines completed
   10 CONTINUE

      INDSTART = AMAX0 (      1 , INDSTART)
      INDEND   = AMAX0 (INDSTART, INDEND  )
      INDEND   = AMIN0 (LASTIND , INDEND  )

      write (0,*) '**** CARDS reading completed' 
      write (0,*) '**** INDSTART=', INDSTART 
      write (0,*) '**** INDEND  =', INDEND 

      DO 99 IND=INDSTART, INDEND

      NUP=INDNUP(IND)
      LOW=INDLOW(IND)
      WAVENUM=ELEVEL(NUP)-ELEVEL(LOW)
      WN2=WAVENUM*WAVENUM
      WN3=WN2*WAVENUM

      TMIN = 5000.
      TMAX = 600 000.
      
      TMINLOG = ALOG10(TMIN)
      TMAXLOG = ALOG10(TMAX)
      NPOINTS = 300
      DTLOG = (TMAXLOG - TMINLOG) / FLOAT(NPOINTS-1)

      CSTRING = 'COLOR=4'

      DO II=1,NPOINTS
         XLOGTL = TMINLOG + FLOAT(II) * DTLOG
         TL = 10.**XLOGTL
         TROOT = SQRT(TL)
         T32 = TL * TROOT

C***  'NONE': TRANSITIONS WITH UNKNOWN COLLISIONAL COEFFICIENTS
C***          (COLLISIONAL CROSS SECTION SIGMA(LOW,UP) IS SET TO  PI*A0**2)
      IF (KEYCBB(IND) .EQ. 'NONE') THEN
C                           ====
              OMEGA=5.465E-11*TROOT*(1.+C1*WAVENUM/TL)*WEIGHT(LOW)/
     /               WEIGHT(NUP)

C***  'ZERO': NO COLLISIONAL TRANSITION
      ELSE IF (KEYCBB(IND) .EQ. 'ZERO') THEN
C                                ====
               OMEGA=0.

C***  HELIUM  ==========================================================
      ELSE IF (NOM(LOW) .EQ. KODAT(2)) THEN
            CALL CBBHE (OMEGA,IND,NUP,LOW,TL,TROOT,T32,NDIM,N,NCHARG,
     $                  EINST,COCO,KEYCBB,WEIGHT,WAVENUM,WN2,WN3)
 
C***  HYDROGEN  ========================================================
      ELSE IF (NOM(LOW) .EQ. KODAT(1)) THEN
            CALL CBBH  (OMEGA,IND,NUP,LOW,TL,TROOT,T32,NDIM,N,NCHARG,
     $                  EINST,COCO,KEYCBB,WEIGHT,WAVENUM,WN2,WN3)
 
C***  NITROGEN  ========================================================
      ELSE IF (NOM(LOW) .EQ. KODAT(7)) THEN
            CALL CBBN  (OMEGA,IND,NUP,LOW,TL,TROOT,T32,NDIM,N,NCHARG,
     $                  EINST,COCO,KEYCBB,WEIGHT,WAVENUM,WN2,WN3)
 
C***  GENERIC ION ======================================================
      ELSE IF (NOM(LOW) .EQ. KODAT(26)) THEN
            CALL CBBFE (OMEGA, NUP, LOW, TL, TROOT, WAVENUM, WN3, 
     >                  EINST, NDIM)
 
      ELSE
         CALL CBBMORE  (OMEGA,IND,NUP,LOW,TL,TROOT,T32,NDIM,N,NCHARG,
     $                  EINST,COCO,KEYCBB,WEIGHT,WAVENUM,WN2,WN3)
      ENDIF

      XPLOT(II) = XLOGTL

cc      IF (OMEGA .GT. OMEGAMIN) THEN
cc         YPLOT(II) = ALOG10(OMEGA)
         YPLOT(II) = OMEGA * 1.E6
      ENDDO

      BACKSPACE (42)      
      WRITE (42,'(A,F6.2,A,I5,2X,A,2X,A)') 
     >   '\LUN XMAX YMIN -7 ', POSLABEL, ' .2 LINE', IND,  
     >   LEVEL(NUP) // ' - ' // LEVEL(LOW), KEYCBB(IND)
      WRITE (42,'(A,G12.5,1X,G12.5,A,I5)') 
     >      '\LUN ', XPLOT(1), YPLOT(1), ' M0 .2 .2 ', IND  
      WRITE (42,'(A)') 'END'
      POSLABEL = POSLABEL - 0.4




      CALL PLOTCON (42, XPLOT, YPLOT, NPOINTS, 5)

   99 CONTINUE

      CLOSE (42)

      STOP 'O.K.'

C********* Error stop
  990 WRITE (0,*) 'INVALID INPUT LINE: ', KARTE(:IDX(KARTE))
      STOP '*** ERROR in PLOTCBB'

      END
