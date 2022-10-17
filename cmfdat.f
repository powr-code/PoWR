      SUBROUTINE CMFDAT      
C**********************************************************************
C     READS IN CMFGEN ATOMIC DATA AND CONVERTS IT INTO DATOM FILE
C     PER DEFAULT THE SUPERLEVELS ARE USED TO SUM LEVELS
C**********************************************************************
      
      IMPLICIT NONE

      INTEGER   NLMAX     ! MAX. NUMBER OF CMFGEN LEVELS
      INTEGER   MAXELEM   ! MAX. NUMBER OF ELEMENTS
      INTEGER   NLINMAX   ! MAX. NUMBER OF CMFGEN TRANSITIONS
      PARAMETER (NLMAX=999, MAXELEM=118, NLINMAX=999999)
      INTEGER NCMF      ! ACTUAL NUMBER OF CMFGEN LEVELS 
      INTEGER NCMFSUPER ! NUMBER OF CMF SUPERLEVELS
      INTEGER N, I, NPAR, NION, NSUPERLEVEL, NTEMP, NL, K, NZ
      INTEGER POS1, POS2
      INTEGER ILOW, IUP, NTRANSFOUND, QNTEST, IZ
      
      REAL ENERGY, WEIGHT, EION, FMEAN, FSUM, GSUM !,XKOEFF
      REAL  XNU, CC, SIGMA 
      CHARACTER LINE*132, ACTPAR*80, ACTPAR2*80, CION*10,test*5
      CHARACTER LNAME*10, FLEVNAME*80, FLINNAME*80, ELEMSYM*2
      CHARACTER MEDIUM*3, FORMUP*10, FORMLOW*10, CXKOEFF*4
      CHARACTER(LEN=10) LEVNAME(NLMAX)


!***  LEVELS FROM CMFGEN 
      CHARACTER(LEN=30) CMFNAME(NLMAX)  ! NAME OF CMFGEN SUBLEVELS
      REAL,DIMENSION(NLMAX)::CMFWEIGHT, ! WEIGHT OF CMFGEN SUBLEVEL
     >                       CMFENERGY  ! ENERGY OF CMFGEN SUBLEVEL 
     >                      ,ELEVEL     ! ENERGY OF CMFGEN SUPERLEVEL
 
      REAL ELEVELM  ! FOR TESTING CORRECT ENERGY ORDER
      INTEGER,DIMENSION(NLMAX)::CMFSUPER ! ID OF CMFGEN SUPERLEVEL
     >                         , ID      ! ID OF CMFGEN SUBLEVEL 
     >                         , QN      ! MAIN QUANTUM NUMBER OF SUPERLEVEL
     >                         , QNLEV   ! MAIN QUANTUM NUMBER OF LEVEL

      LOGICAL,DIMENSION(NLMAX)::MIXQN   ! Number of the SUPERLEVELS with mixed MAIN QUANTUM NUMBER
      LOGICAL TRUE

      INTEGER NSPLITSUPER(NLMAX)         ! NUMBER OF SUBLEVELS IN SUPERLEVEL

!***  LINES FROM CMFGEN
      CHARACTER(LEN=30) LOWERLEVEL(NLINMAX) ! NAME OF LOWER CMFGEN SUBLEVELS
      CHARACTER(LEN=30) UPPERLEVEL(NLINMAX) ! NAME OF UPPER CMFGEN SUBLEVELS
      REAL,DIMENSION(NLINMAX)::  FVALUE     ! F-VALUE OF TRANSITION
     >                         ,LAM         ! WAVELENGTH OF TRANSITION
      INTEGER,DIMENSION(NLINMAX):: IDLOW,IDUP ! ID OF LOWER AND UPPER SUBLEVEL
      INTEGER NTRANS                        ! NUMBER OF TRANSITIONS

      CHARACTER CX*1, CKB*4, CFVAL*10, CHION*1
      CHARACTER CHLANG*1,CHPAR*1,CHSPIN*1,CHTOTLANG*1

      REAL,DIMENSION(NLMAX) :: DATOMWEIGHT

      LOGICAL BBEGINREAD, BFIRSTENTRY, BDATOM, BFORMAL, BQNTEST
      LOGICAL BCONT, BLINES, BLEVELS
      LOGICAL BTRANS(NLMAX,NLMAX) 

      CHARACTER, EXTERNAL :: ARABIC2ROMAN*5
      INTEGER, EXTERNAL :: ROMAN2ARABIC, GETQN, IDX

!***  FILE WITH OLD LEVEL NAMES
      CHARACTER OLDLEVNAMFILE*80 ! FILE CONTAINING NAMES FOR LEVELS
      LOGICAL   BOLDNAMES
      CHARACTER (LEN=10) OLDNAMES(NLMAX)
      INTEGER NOLDLEVNUM, NOLDNAMES


C***  Link data to identify program version
      CHARACTER LINK_DATE*30, LINK_USER*10, LINK_HOST*60
      COMMON / COM_LINKINFO / LINK_DATE, LINK_USER, LINK_HOST
      CHARACTER dat*8,tim*10,zone*5,stim*8,sdat*9,username*20
      INTEGER value(8)


      CHARACTER(LEN=10) ELEMENT(MAXELEM) ! ELEMENT name
      CHARACTER(LEN=2)  SYMB(MAXELEM)    ! ELEMENT SYMBOL
      REAL              ATMASS(MAXELEM)  ! ATOM MASS

      CHARACTER(LEN=132) KSHELLFILE
      CHARACTER(LEN=60) KSHELLPATH, KSHELLFNAME
 
      REAL        FMINALLOWED                                                                                                   
      DATA        FMINALLOWED / 0.01 /   

      DATA        CC /2.997929E10/


      DATA (ELEMENT(I),I=1,20) / 
     >     'HYDROGEN', 'HELIUM', 'LITHIUM', 'BERYLLIUM',
     >     'BORON', 'CARBON', 'NITROGEN', 'OXYGEN', 
     >     'FLUORINE', 'NEON', 'SODIUM', 'MAGNESIUM', 'ALUMINIUM', 
     >     'SILICON', 'PHOSPHORUS', 'SULFUR', 'CHLORINE', 'ARGON', 
     >     'POTASSIUM', 'CALCIUM' /

      DATA (SYMB(I),I=1,20) / 
     >     'H ', 'HE', 'LI', 'BE', 'B ', 'C ', 'N ', 'O ', 'F ', 'NE', 
     >     'NA', 'MG', 'AL', 'SI', 'P ', 'S ', 'CL', 'AR', 'K ', 'CA' /

      DATA (ATMASS(I),I=1,20) / 
     >      1.00,  4.00,  6.93,  9.01, 10.81, 12.01, 14.01, 16.00, 
     >     19.00, 20.18, 22.99, 24.31, 26.98, 28.09, 30.97, 32.07,
     >     35.45, 39.95, 39.10, 40.08 /


      KSHELLPATH='/home/corona/wrh/work/wrdata-archive/'
      KSHELLFNAME='DATOM_K-SHELL.allIons+averages_Verner-sexpo'
      KSHELLFILE=TRIM(KSHELLPATH) // KSHELLFNAME


!*** ---------------------------------------------------------------------

      IF (IARGC() .GE. 1) THEN
         CALL GETARG(1,LINE)
         IF ( LINE (1:2) .EQ. '-h' .OR. LINE(1:6) .EQ. '--help') THEN
          WRITE (*,'(/,A)') "-------------------------------------------------------"
          WRITE (*,'(A)')   "Converts CMF atomic data to PoWR DATOM and FORMAL_CARDS"
          WRITE (*,'(A)') "-------------------------------------------------------"
          WRITE (*,'(A)')   "Required file: INPUT"
          WRITE (*,'(A)') "Required data in INPUT: "
          WRITE (*,'(A)') "ELEMENT=    ! element symbol, e.g., ELEMENT=SI"
          WRITE (*,'(A)') "ION=        ! ionization stage, e.g., ION=IV"  
          WRITE (*,'(A)') "LEVELS=     ! file with CMF superlevels, e.g., LEVELS=f_to_s_split.dat" 
          WRITE (*,'(A)') "LINES=      ! file with CMF lines, e.g., LINES=osc_op_split.dat" 
          WRITE (*,'(A)')   "-------------------------------------------------------"
          WRITE (*,'(A)') "Options for file INPUT:"
          WRITE (*,'(A)') "FORMAL_CARDS"
          WRITE (*,'(A)') "DATOM"
          WRITE (*,'(A)') "DATOM_LEVELS ! print out only DATOM levels"
          WRITE (*,'(A)') "DATOM_LINES  ! print out only DATOM lines"
          WRITE (*,'(A)') "DATOM_CONT   ! print out only DATOM continua"
          WRITE (*,'(A)') "NLEVEL=      ! restrict max. CMF-superlevel, e.g., NLEVEL=50"
          WRITE (*,'(A)') "OLDNAMES=    ! file with own levelnames, instead of automatic names"
          WRITE (*,'(A)') "             ! format of file is 1st col: number   2nd col: name, e.g, "
          WRITE (*,'(A)') '             ! 1 = "SI IV 8P21" '
          WRITE (*,'(A)') '             ! 3 = "SI IV 3D.3" '
          WRITE (*,'(A)')   "-------------------------------------------------------"
          STOP          
          ELSE
          WRITE (0,'(A)') "Get info with option -help"
          STOP
         ENDIF
         
      ENDIF


c     Identify user of opdat and date and time of usage
      CALL GETENV ('USER', USERNAME)
      CALL DATE_AND_TIME (dat,tim,zone,value)      
          WRITE (sdat,'( I2.2,"-",I2.2,"-",A2)')
     >    value(3),value(2),dat(3:4)
          WRITE (stim,'(I2.2,":",I2.2,":",I2.2)')
     >    value(5),value(6),value(7)

C     Version output
      WRITE(*,'(2A)') '*---- cmfdat Version from ', LINK_DATE
      WRITE(*,'(4A)') '*---- linked by ', LINK_USER(:IDX(LINK_USER)),
     >      ' at host ', LINK_HOST(:IDX(LINK_HOST))
c      WRITE (0,'(A)') '*---- get info:  opdat -help '
c     Producer output
      WRITE(*,'(6A)') "*---- DATOM created on: ", sdat , ' ', stim
     > ,' by ', USERNAME

      BTRANS        = .FALSE.
      NCMFSUPER     = 0
      FLEVNAME      = ''
      FLINNAME      = ''
      BBEGINREAD    = .FALSE.
      BFIRSTENTRY   = .TRUE.
      MEDIUM        = 'VAC'
      BDATOM        = .FALSE.
      BCONT         = .FALSE.
      BLINES        = .FALSE.
      BLEVELS       = .FALSE.
      BLINES        = .FALSE. 
      BFORMAL       = .FALSE. 
      NSUPERLEVEL   = NLMAX
      OLDLEVNAMFILE = ''
      BOLDNAMES     = .FALSE.
      OLDNAMES      = ''
      TRUE          = .TRUE.

!***  READ INPUT STEERING FILE -----------------------------------------------

      OPEN(2,FILE='INPUT',STATUS='OLD',ERR=992)
      DO 
       READ(2,400,END=12) LINE
       CALL SARGV(LINE,1,ACTPAR)
       IF (ACTPAR .EQ. 'ELEMENT') THEN
        CALL SARGV (LINE,2,ACTPAR2)
        ELEMSYM=ACTPAR2
       ELSEIF (ACTPAR .EQ. 'ION') THEN
        CALL SARGV (LINE,2,ACTPAR2)
        CION=ACTPAR2
        NION=ROMAN2ARABIC(CION)
       ELSEIF (ACTPAR .EQ. 'LEVELS') THEN
        CALL SARGV (LINE,2,ACTPAR2)
        READ(ACTPAR2,'(A)') FLEVNAME
       ELSEIF (ACTPAR .EQ. 'LINES') THEN
        CALL SARGV (LINE,2,ACTPAR2)
        READ(ACTPAR2,'(A)') FLINNAME
       ELSEIF (ACTPAR .EQ. 'OLDNAMES') THEN
        CALL SARGV (LINE,2,ACTPAR2)
        READ(ACTPAR2,'(A)') OLDLEVNAMFILE
        BOLDNAMES = .TRUE.
       ELSEIF (ACTPAR .EQ. 'NLEVEL') THEN
        CALL SARGV (LINE,2,ACTPAR2)
        READ(ACTPAR2,'(I4)') NSUPERLEVEL
       ELSEIF (ACTPAR .EQ. 'DATOM') THEN
         BDATOM = .TRUE.
       ELSEIF (ACTPAR .EQ. 'DATOM_LEVELS') THEN
         BLEVELS = .TRUE.
       ELSEIF (ACTPAR .EQ. 'DATOM_LINES') THEN
         BLINES = .TRUE.
       ELSEIF (ACTPAR .EQ. 'DATOM_CONT') THEN
         BCONT = .TRUE.
       ELSEIF (ACTPAR .EQ. 'FORMAL_CARDS') THEN
         BFORMAL = .TRUE.
       ENDIF       
      ENDDO

 12   CLOSE(2)

       IF (NION .LT. 10 .AND. NION .GT. 1) THEN
        WRITE (CHION,'(I1)') NION
       ELSEIF (NION .EQ. 1) THEN
        CHION = "I"
       ELSEIF (NION .EQ. 10) THEN
        CHION = "X"
       ELSE
        WRITE (*,*) "ION = ", NION
        STOP 'ION .GT. X NOT SUPPORTED'
       ENDIF
      
      
      WRITE (*,'(2A)')  "*---- LEVELS   = ", TRIM(FLEVNAME)
      WRITE (*,'(2A)')  "*---- LINES    = ", TRIM(FLINNAME)
      CALL GETENV ('PWD', LINE)
      WRITE (*,'(2A)')  "*---- DIRECTORY= ", TRIM(LINE)


!***  READ FILE WITH OLD LEVEL NAMES -----------------------------------------

      IF (BOLDNAMES) THEN
      OPEN (2,FILE=OLDLEVNAMFILE,FORM='FORMATTED',STATUS='OLD',ERR=990)
      
      ! EACH LINE IS FOR A LEVEL NAME
      NOLDNAMES = 0
      LOOP_OLDNAMES : DO N = 1, NLMAX
         READ(2,400,END=22) LINE
         CALL SARGC(LINE, NPAR)
         IF (NPAR .NE. 0) THEN
            CALL SARGV (LINE,1,ACTPAR)
            IF (ACTPAR(1:1) .EQ. '#' .OR. ACTPAR(1:1) .EQ. '*')
     >           CYCLE LOOP_OLDNAMES
            READ (ACTPAR,'(I4)',ERR=993) NOLDLEVNUM
            CALL SARGV (LINE,2,ACTPAR)
            OLDNAMES(NOLDLEVNUM) = ACTPAR
            NOLDNAMES = NOLDNAMES + 1
         ENDIF   
      ENDDO LOOP_OLDNAMES
      
 22   CLOSE(2)
      ! DO N = 1, NOLDNAMES+3
      !   IF (OLDNAMES(N) .NE. '') THEN
      !      WRITE (*,'(A,A,A)') "=",OLDNAMES(N),"="
      !   ELSE
      !      WRITE (*,*) n , " IS EMPTY"
      !   ENDIF
      !ENDDO
 
      ENDIF


!***  READ CMFGEN LEVELS WITH SUPERLEVEL DESIGNATIONS ------------------------

      OPEN(1,FILE=FLEVNAME,
     > FORM='FORMATTED',STATUS='OLD',ERR=991)

      NL = 0
      LOOP_LEVELREAD : DO N=1, NLMAX
       READ(1,400,END=11) LINE
       CALL SARGC(LINE, NPAR) 
       IF (NPAR .NE. 0 .AND. .NOT. BBEGINREAD) THEN
        DO I=1, NPAR
         CALL SARGV (LINE,I,ACTPAR)
         IF (TRIM(ACTPAR) .EQ. '!Entry') THEN
          BBEGINREAD = .TRUE.
          CYCLE LOOP_LEVELREAD
         ENDIF
        ENDDO
       ENDIF

       IF (BBEGINREAD) THEN
        IF (NPAR .EQ. 0) CYCLE
        CALL SARGV(LINE,6,ACTPAR)
        NL = NL + 1
        READ (ACTPAR,'(I4)') NTEMP 
        IF (NTEMP .GT. NSUPERLEVEL) CYCLE
        CMFSUPER(NL) = NTEMP
        CALL SARGV(LINE,8,ACTPAR)
        READ (ACTPAR,'(I4)') ID(NL)
        CALL SARGV(LINE,1,ACTPAR)
        READ (ACTPAR,'(A)') CMFNAME(NL)
        CALL SARGV(LINE,2,ACTPAR)
        READ (ACTPAR,'(F10.0)',ERR=993) CMFWEIGHT(NL)
        CALL SARGV(LINE,3,ACTPAR)
        READ (ACTPAR,'(F18.0)') CMFENERGY(NL)
        IF(BFIRSTENTRY) THEN
         CALL SARGV(LINE,5,ACTPAR)
          READ (ACTPAR,'(F12.0)') EION
          EION = 1.E8 / EION 
        ENDIF
        NCMFSUPER = MAX0(CMFSUPER(NL),NCMFSUPER) ! NUMBER OF SUPERLEVELS
        NCMF = NL ! NUMBER OF CMF-SUBLEVELS
        BFIRSTENTRY=.FALSE.
       ENDIF
      ENDDO LOOP_LEVELREAD
      

 11   CLOSE(1)

      WRITE (0,'(A,I5,A,I5,A)') 
     >   "* FOUND " , NCMF, " LEVELS IN ", NCMFSUPER, " SUPERLEVELS"

      ! CREATE VECTOR WITH NUMBER OF SUBLEVES IN SUPERLEVELS
      DO I = 1, NCMF
         IF ( CMFSUPER(I) .LE. NSUPERLEVEL ) 
     >     NSPLITSUPER (CMFSUPER(I)) = NSPLITSUPER(CMFSUPER(I)) + 1
      ENDDO
      
!***  READ IN LINES -------------------------------------------------

      OPEN(3,FILE=FLINNAME,
     > FORM='FORMATTED',STATUS='OLD',ERR=994)

      BBEGINREAD = .FALSE.
      NL = 0 

      LOOP_LINEREAD : DO N=1, NLINMAX

       READ(3,400,END=13) LINE
       CALL MSARGC(LINE, NPAR) 
       IF (NPAR .NE. 0 .AND. .NOT. BBEGINREAD) THEN
        DO I=1, NPAR
         CALL MSARGV (LINE,I,ACTPAR)
         IF (ACTPAR .EQ. 'Transition') THEN
          BBEGINREAD = .TRUE.
          CYCLE LOOP_LINEREAD
         ENDIF
        ENDDO
       ENDIF

       IF (BBEGINREAD) THEN
        IF (NPAR .EQ. 0) CYCLE
        IF (NL .GT. 2 .AND. NPAR .LE. 1) EXIT LOOP_LINEREAD
        NL = NL + 1
        CALL MSARGV(LINE,1,ACTPAR)
        LOWERLEVEL(NL) = ACTPAR
        CALL MSARGV(LINE,2,ACTPAR)
        UPPERLEVEL(NL) = ACTPAR
        CALL MSARGV(LINE,3,ACTPAR)
        READ (ACTPAR,'(F20.0)',ERR=201) FVALUE(NL) 
        CALL MSARGV(LINE,5,ACTPAR)
        READ (ACTPAR,'(F20.0)') LAM(NL)
        CALL MSARGV(LINE,6,ACTPAR)
        READ (ACTPAR,'(I4)') IDLOW(NL)
        CALL MSARGV(LINE,7,ACTPAR)
        READ (ACTPAR,'(I4)') IDUP(NL)
        NTRANS = NL ! NUMBER OF FOUND TRANSITIONS
       ENDIF  
      ENDDO LOOP_LINEREAD

 13   CLOSE(3)


!***  WRITE OUT DATOM LEVELS ----------------------------------------

      DO I = 1 , 20
        IF ( ELEMSYM .EQ. SYMB(I)) NZ = I
      ENDDO
      IZ = ROMAN2ARABIC(CION) - 1

      IF (BDATOM .OR. BLEVELS) THEN

C     write first element information
      WRITE (*,'(A)') '*KEYWORD--  ---NAME--- SYMB   ATMASS   STAGE'
      WRITE (*,'(A)') '********************************************'
      WRITE (*,'(A,A,1X,A,A,A,3X,F6.2,3X,I3)') 
     >     'ELEMENT     ', ELEMENT(NZ) , '(', SYMB(NZ) , ')', ATMASS(NZ)
     >     , MAX(IZ,1)
      ! WRITE (*,'(12X,A)') '=========='
      WRITE (*,'(11X,A,$)') ' '
      DO I=1, LEN_TRIM(ELEMENT(NZ))
         WRITE (*,'(A,$)') '='
      ENDDO
      WRITE (*,'(1X)')
      WRITE (*,'(A)') '********************************************'
!      WRITE (*,'(A)') 


      WRITE (*,'(A)') 
     > '*KEYWORD--  ---NAME--- CH WEIG--ENERGY-- EION----- QN'

      ENDIF ! BDATOM 

      ELEVELM = -1.

      DO I=1, NCMFSUPER
         
       ENERGY = 0.
       
       DATOMWEIGHT(I) = 0. 
       ! EXTRACTING AND CHECKING MAIN QUANTUM NUMBER FOR SUPERLEVEL
       QN(I) = 0
       QNTEST = QN(I) 
       BQNTEST = .FALSE. 

       DO N=1, NCMF
        IF (CMFSUPER(N) .EQ. I) THEN
         QN(I) = GETQN(TRIM(CMFNAME(N)),CHLANG,CHPAR,CHSPIN,CHTOTLANG)
         QNLEV(N) = QN(I)
         IF (QNTEST .NE. QN(I) .AND. BQNTEST) THEN
          WRITE (0,'(A,I3,A,I3,1X,I3)') 
     >  "WARNING : MIXING OF QN IN SUPERLEVEL ", I, 
     >   " MAIN QUANTUM NUMBERS FOUND: ", QNTEST, QN(I)
         MIXQN(I) = .TRUE.
         ELSE
         MIXQN(I) = .FALSE.
         ENDIF
         QNTEST = QN(I)
         BQNTEST = .TRUE.
         ENERGY =  ENERGY + CMFWEIGHT(N) * CMFENERGY(N)
         DATOMWEIGHT(I) = DATOMWEIGHT(I) + CMFWEIGHT(N) 
        ENDIF
       ENDDO
       ELEVEL(I) = ENERGY / DATOMWEIGHT(I)

       LNAME(1:10) = ".........."
       LNAME(1:2) = ELEMSYM
       LNAME(3:3) = CHION

       IF (QN(I) .LT. 10) THEN
        WRITE (LNAME(4:4),'(I1)')  QN(I)
        LNAME(5:5) = CHLANG
       ELSE
        WRITE (LNAME(4:5),'(I2)')  QN(I)
       ENDIF

       LNAME(6:6) = CHPAR
       LNAME(7:7) = CHSPIN
       LNAME(8:8) = CHTOTLANG

       IF ( I .LT. 100) THEN
        WRITE (LNAME(9:10),'(I2.2)') I
       ELSE
        WRITE (LNAME(8:10),'(I3.3)') I
       ENDIF

       LEVNAME(I) = LNAME
       IF (OLDNAMES(I) .NE. '') LEVNAME(I) = OLDNAMES(I)
       
       IF (ELEVEL(I) .LE. ELEVELM) THEN
        WRITE (0,'(A,I4,A,I4)') 
     >    "==> WRONG ENERGY ORDER FOR CMF-SUPERLEVELS "
     >    ,I-1, " AND ", I
          WRITE (0,'(A,1X,F10.1)')  LNAME, ELEVEL(I)
          STOP 
       ENDIF
       ELEVELM = ELEVEL(I)

       IF (BDATOM .OR. BLEVELS) THEN
        IF (I .EQ. 1) THEN
         WRITE (*,101) 
     >     LEVNAME(I), IZ,INT(DATOMWEIGHT(I)),ELEVEL(I),EION,QN(I)
        ELSE 
         WRITE (*,102) 
     >     LEVNAME(I) ,IZ,INT(DATOMWEIGHT(I)),ELEVEL(I),QN(I)
        ENDIF
       ENDIF ! BDATOM
      ENDDO

      IF (BDATOM .OR. BLEVELS) THEN
      WRITE (*,'(A)') 
     > '*KEYWORD--  ---NAME--- CH WEIG--ENERGY-- EION----- QN'
      ENDIF  ! BDATOM

!***  SUPERLINES ---------------------------------------------------

      IF (BDATOM .OR. BLINES) THEN
       WRITE (*,'(A,A)')
     >'*KEYWORD--UPPERLEVEL  LOWERLEVEL--EINSTEIN  RUD-CEY -',
     >'-COLLISIONAL COEFFICIENTS--'
      ENDIF

      DO IUP = 2, NCMFSUPER     ! LOOP OVER UPPER SUPERLEVELS
       DO ILOW = 1, IUP - 1 ! LOOP OVER LOWER SUPERLEVELS
        FSUM = 0.
        GSUM = 0. 
        FMEAN = 0. 
        NTRANSFOUND = 0
        DO K = 1, NTRANS        ! LOOP OVER TRANSITIONS
         IF ( CMFSUPER(IDUP(K)) .EQ. IUP 
     >    .AND. CMFSUPER(IDLOW(K)) .EQ. ILOW) THEN
          NTRANSFOUND = NTRANSFOUND + 1
          FSUM = FSUM + CMFWEIGHT(IDLOW(K)) * FVALUE(K)
         ENDIF
        ENDDO   ! LOOP OVER TRANSITIONS

        IF (NTRANSFOUND .GT. 0) THEN 
           BTRANS(IUP,ILOW) = .TRUE.
           FMEAN = FSUM / DATOMWEIGHT(ILOW)
        ENDIF

C   ...   ALLOWED TRANSITIONS                                                                                                   
          ! NOTE: THE KEYWORDS KB22, KB24 ARE NOT ALLOWED FOR HE I

          IF (FMEAN  .GE. FMINALLOWED) THEN                                                                                
           WRITE (CFVAL,'(1X,1PG9.2)') -FMEAN
           CX=' '
           IF (ELEMSYM .EQ. 'HE' .AND. CION .EQ. 'I') THEN
            CKB='JEFF'    
            CXKOEFF=""
           ELSE                                                                                                             
            CKB='KB22'
            IF (ABS(QN(IUP) - QN(ILOW)) .EQ. 0) THEN
              CXKOEFF = "0.70"  
            ELSE                                                                                                               
              CXKOEFF = "0.20"                                                                                                   
            ENDIF         

           ENDIF


C   ...   FORBIDDEN TRANSITIONS                                                                                                 

          ELSEIF (FMEAN .LT. FMINALLOWED) THEN                                                                            
             WRITE (CFVAL,'(1X,1PG9.2)') -FMEAN                                                                           
             CX='X'             
C            RESONANCE AND PSEUDO-RESONANCE LINES                                                                               
             IF (ILOW .EQ. 1 .AND. FMEAN .GT. 1.E-98) THEN                                                                  
                CX=' '                                                                                                          
             ELSEIF (FMEAN .LT. 1.E-98) THEN                                                                              
                CFVAL='          '                                                                                              
             ENDIF   

C            COLLISIONAL CROSS SECTIONS DEPEND ON DELTA NQ                                                                      
             IF (ABS(QN(IUP) - QN(ILOW)) .EQ. 0) THEN                                                                    
                IF (ELEMSYM .EQ. 'HE' .AND. CION .EQ. 'I') THEN
                   CKB='UPS2'
                   CXKOEFF=""
                ELSE
                   CKB='KB24'                                                                                                      
                   CXKOEFF="1.00"
                ENDIF
             ELSEIF (ABS(QN(IUP) - QN(ILOW) ) .EQ. 1) THEN                                                                
                IF (ELEMSYM .EQ. 'HE' .AND. CION .EQ. 'I') THEN
                   CKB='UPS1'
                   CXKOEFF=""
                ELSE
                   CKB='KB24'                                                                                                      
                   CXKOEFF="0.05"
                ENDIF   
             ELSE                                                                                                               
                CKB='NONE'                                                                                                      
                CXKOEFF="0.00"                                                                                                     
C                removed 17.09.2013:                                                                                            
c                CFVAL='          '                                                                                             
             ENDIF                                                                                                              
                                                                                                                                
          ENDIF             

         
         IF (BDATOM .OR. BLINES) THEN
          
          WRITE(*,103) LEVNAME(IUP),LEVNAME(ILOW),CFVAL,CX,CKB,CXKOEFF
 103      FORMAT('LINE      ',A10,'  ',A10,A10,'   ',A1,' ',A4,'    ',
     >        A)
          ENDIF

       ENDDO ! LOOP OVER LOWER SUPERLVELS
      ENDDO ! LOOP OVER UPPER SUPERLEVELS

      IF (BDATOM) THEN
       WRITE (*,'(A,A)')
     >'*KEYWORD--UPPERLEVEL  LOWERLEVEL--EINSTEIN  RUD-CEY -',
     >'-COLLISIONAL COEFFICIENTS--'
      ENDIF

!*** WRITE OUT CONTINUA -----------------------------------------------
      
      IF (BDATOM .OR. BCONT) THEN
       WRITE(*,'(A,A)')
     >'*KEYWORD  LOWERLEVEL ----SIGMA ----ALPHA ----SEXPO -IGAUNT- -K'
     >, 'EYCBF- -IONLEV---' 

       DO I=1, NCMFSUPER
         XNU   = (EION - ELEVEL(I)) * CC
         SIGMA = 4.541E-10 / FLOAT(IZ+1) / SQRT(XNU)
         SIGMA = SIGMA*1.E18

            WRITE(*,104) LEVNAME(I), SIGMA
 104        FORMAT ('CONTINUUM ',A10,' ',1PE9.2E2)
       ENDDO

       WRITE(*,'(A,A)')
     >'*KEYWORD  LOWERLEVEL ----SIGMA ----ALPHA ----SEXPO -IGAUNT- -K'
     >, 'EYCBF- -IONLEV---'


C***  WRITE out K-SHELL data ---------------------------------------
       
      OPEN (42,FILE=KSHELLFILE,STATUS='UNKNOWN',ERR=843)
!     !'/home/corona/wrh/work/wrdata-archive/DATOM_K-SHELL.
!     !>allIons+averages_Verner-sexpo'
!     >      ,STATUS='UNKNOWN',ERR=843)
      DO 
       READ (42,'(A132)',END=842) LINE
       IF (LINE(16:17) .EQ. SYMB(NZ)) THEN
          READ (LINE(18:20),'(I3)') I
          IF (I .EQ. IZ+1) THEN
           WRITE(*,'(A)')
     >     '*KEYWORD--*****SY*I*<-K-SIGMA><-K-SEXPO>***<-K-EION->'
           WRITE (*,'(A132)') LINE
           WRITE(*,'(A)')
     >     '*KEYWORD--*****SY*I*<-K-SIGMA><-K-SEXPO>***<-K-EION->'
          ENDIF
       ENDIF   
      ENDDO
 842  CLOSE (42)

      ENDIF ! END OF BDATOM OR BCONT

!***  FORMAL CARDS ---------------------------------------------------

      IF (BFORMAL) THEN

      WRITE (*,'(A)') "* FORMAL "
      WRITE (*,*) ""

      DO IUP = 2, NCMFSUPER     ! LOOP OVER UPPER SUPERLEVELS

       DO ILOW = 1, IUP - 1     ! LOOP OVER LOWER SUPERLEVELS

        IF (BTRANS(IUP,ILOW)) THEN
         IF ( NSPLITSUPER(IUP) .GT. 1 .OR.
     >        NSPLITSUPER(ILOW) .GT. 1 ) THEN
           
         WRITE(*,'(A)') '+MULTIPLET ?'
         WRITE (*,'(A,A,1X,A,A)') "UPPERLEVEL=", LEVNAME(IUP), 
     >               "LOWERLEVEL=", LEVNAME(ILOW)  
      
         DO I = 1, NCMF   ! split lower level
            IF (CMFSUPER(I) .EQ. ILOW) THEN
               IF (NSPLITSUPER(ILOW) .GT. 1) THEN
                  FORMLOW = '..........'
                  FORMLOW(1:2) = ELEMSYM
                  FORMLOW(3:3) = CHION
                  FORMLOW(5:6) = LEVNAME(ILOW)(9:10)
                  IF (MIXQN(ILOW) .EQV. TRUE) THEN   ! check if cmf superlevels are for the same J values 
                     FORMLOW(8:8) =   "n"            ! or angular momentums of different main quantum numbers
                     IF (QNLEV(I) .LT. 10) THEN           
                        WRITE (FORMLOW(10:10), '(I1)') QNLEV(I)
                        ELSEIF ( QNLEV(I) .LT. 100) THEN
                        WRITE (FORMLOW(9:10), '(I2)')  QNLEV(I)
                     ENDIF    
                  ELSE
                     POS1 = index(CMFNAME(I),"[")
                     POS2 = POS1 + 3
                     IF ( POS1 .NE. 0 ) THEN
                        FORMLOW(7:7) =   "J"
                        FORMLOW(8:10) =   CMFNAME(I)(POS1+1:POS2) 
                     ENDIF
                  ENDIF                  
                  WRITE(*,107) FORMLOW, INT(CMFWEIGHT(I)), 
     >                 CMFENERGY(I)
 107              FORMAT ('/LOWERLEVEL ',A10, ' ', I4, ' ', F18.3)
               ENDIF
            ENDIF
         ENDDO

         DO I = 1, NCMF  ! split upper level
            IF (CMFSUPER(I) .EQ. IUP) THEN  
             IF (NSPLITSUPER(IUP) .GT. 1) THEN
                FORMUP = '..........'
                FORMUP(1:2) = ELEMSYM
                FORMUP(3:3) = CHION
                FORMUP(5:6) = LEVNAME(IUP)(9:10)
                IF (MIXQN(IUP) .EQV. TRUE) THEN   ! check if cmf superlevels are for the same J values 
                   FORMUP(8:8) =   "n"            ! or angular momentums of different main quantum numbers
                   IF (QNLEV(I) .LT. 10) THEN           
                      WRITE (FORMUP(10:10), '(I1)') QNLEV(I)
                      ELSEIF ( QNLEV(I) .LT. 100) THEN
                      WRITE (FORMUP(9:10), '(I2)')  QNLEV(I)
                   ENDIF    
                ELSE
                   POS1 = index(CMFNAME(I),"[")
                   POS2 = POS1 + 3
                   IF ( POS1 .NE. 0 ) THEN
                      FORMUP(7:7) =   "J"
                      FORMUP(8:10) =   CMFNAME(I)(POS1+1:POS2) 
                   ENDIF
                ENDIF   
                WRITE(*,108) FORMUP, INT(CMFWEIGHT(I)), 
     >                    CMFENERGY(I)
 108            FORMAT ('/UPPERLEVEL ',A10, ' ', I4, ' ', F18.3)
             ENDIF
          ENDIF
         ENDDO

         DO K = 1, NTRANS        ! LOOP OVER TRANSITIONS
          IF ( CMFSUPER(IDUP(K)) .EQ. IUP 
     >     .AND. CMFSUPER(IDLOW(K)) .EQ. ILOW) THEN

             IF (NSPLITSUPER(ILOW) .GT. 1) THEN
                FORMLOW = '..........'
                FORMLOW(1:2) = ELEMSYM
                FORMLOW(3:3) = CHION
                FORMLOW(5:6) = LEVNAME(ILOW)(9:10)
                IF (MIXQN(ILOW) .EQV. TRUE) THEN   ! check if cmf superlevels are for the same J values 
                   FORMLOW(8:8) =   "n"            ! or angular momentums of different main quantum numbers
                   IF (QNLEV(IDLOW(K)) .LT. 10) THEN           
                      WRITE (FORMLOW(10:10), '(I1)') QNLEV(IDLOW(K))
                      ELSEIF ( QNLEV(IDLOW(K)) .LT. 100) THEN
                      WRITE (FORMLOW(9:10), '(I2)')  QNLEV(IDLOW(K))
                   ENDIF    
                ELSE
                   POS1 = index(CMFNAME(IDLOW(K)),"[")
                   POS2 = POS1 + 3
                   IF ( POS1 .NE. 0 ) THEN
                      FORMLOW(7:7) =   "J"
                      FORMLOW(8:10) =   CMFNAME(IDLOW(K))(POS1+1:POS2) 
                   ENDIF
                ENDIF   
             ELSE
                FORMLOW =  LEVNAME(ILOW)
             ENDIF

             IF (NSPLITSUPER(IUP) .GT. 1) THEN
                FORMUP = '..........'
                FORMUP(1:2) = ELEMSYM
                FORMUP(3:3) = CHION
                FORMUP(5:6) = LEVNAME(IUP)(9:10)
                  IF (MIXQN(IUP) .EQV. TRUE) THEN   ! check if cmf superlevels are for the same J values 
                     FORMUP(8:8) =   "n"            ! or angular momentums of different main quantum numbers
                    IF (QNLEV(IDUP(K)) .LT. 10) THEN           
                     WRITE (FORMUP(10:10), '(I1)') QNLEV(IDUP(K))
                    ELSEIF ( QNLEV(IDUP(K)) .LT. 100) THEN
                     WRITE (FORMUP(9:10), '(I2)')  QNLEV(IDUP(K))
                    ENDIF    
                  ELSE

                POS1 = index(CMFNAME(IDUP(K)),"[")
                POS2 = POS1 + 3
                IF ( POS1 .NE. 0 ) THEN
                   FORMUP(7:7) =   "J"
                   FORMUP(8:10) =   CMFNAME(IDUP(K))(POS1+1:POS2) 
                ENDIF
                  ENDIF   
             ELSE   
                FORMUP =  LEVNAME(IUP)
             ENDIF   
           MEDIUM = 'VAC'
           IF (LAM(K) .GE. 2000. ) MEDIUM = 'AIR'
  
           IF (LAM(K) .LT. 1.E+07) THEN
            WRITE(*,109) FORMUP, 
     >                  FORMLOW,
     >                                  -FVALUE(K), LAM(K), MEDIUM
            ! different output format for different ranges of wavelength
            ! multiple.f reads 10 digits for f-value, but we use one for a blank
 109        FORMAT ('/SUBLINE ',A10,'  ',A10,1X,1PE9.2, '  ',0PF10.2,5X,A3)
           ELSEIF (LAM(K) .GE. 1.E+10) THEN
            WRITE (0,'(A,A,A,A)') 
     >       "WARNING: LAMBDA >= 1.E10 ANG FOR SUBLINE", FORMUP,' ' 
     >        , FORMLOW  
            WRITE(*,111) FORMUP, 
     >                  FORMLOW,
     >                                  -FVALUE(K), LAM(K), MEDIUM
 111        FORMAT ('/SUBLINE ',A10,'  ',A10,1X,1PE9.2, '  ',0PG10.2E2,5X,A3)

           ELSE
            WRITE(*,110) FORMUP, 
     >                  FORMLOW,
     >                                  -FVALUE(K), LAM(K), MEDIUM
 110        FORMAT ('/SUBLINE ',A10,'  ',A10,1X,1PE9.2, '  ',0PG10.4E1,5X,A3)
           ENDIF

          ENDIF
         
         ENDDO   ! LOOP OVER TRANSITIONS

        WRITE(*,'(A)') '-MULTIPLET'
        WRITE (*,*) '' 

        ELSE ! -------- LINE WITHOUT SPLITTING --------------------------
         DO K = 1, NTRANS
          IF ( CMFSUPER(IDUP(K)) .EQ. IUP 
     >     .AND. CMFSUPER(IDLOW(K)) .EQ. ILOW) THEN
           MEDIUM = 'VAC'
           IF (LAM(K) .GE. 2000. ) MEDIUM = 'AIR'
           WRITE (*,'(A10,F14.2,5X,3A)') 
     >      "+LINE ??? ", LAM(K), MEDIUM
           WRITE(*,112)
     >      LEVNAME(IUP), LEVNAME(ILOW)      
 112       FORMAT ('UPPERLEVEL=',A10,' LOWERLEVEL=',A10)      
           WRITE (*,*) ""
          ENDIF
         ENDDO ! LOOP OVER TRANSITIONS

        ENDIF ! MULITPLET / LINE
        ENDIF ! BTRANS

       ENDDO ! LOOP OVER LOWER SUPERLVELS
      ENDDO ! LOOP OVER UPPER SUPERLEVELS

      ENDIF  ! END OF BFORMAL -------------------------------------------
      RETURN
C *** error codes ******************************************************
 201  WRITE (*,*) "PROBLEMS WHILE READING FVALUE ", ACTPAR, NL, N, LINE
      STOP
 843  STOP 'ERROR: FILE FOR K-SHELL DATA MISSING'
 990  STOP 'ERROR WHILE TRYING TO OPEN FILE OLDNAMES'
 991  STOP 'ERROR while trying to open file f_to_s.dat ... '
 992  STOP 'Steering file  INPUT missing'
 994  STOP 'ERROR while trying to open file osc.dat ... '
 993  WRITE (*,*) "ERROR IN LINE", N
      WRITE (*,*) LINE
      STOP

C *** Formats
  400 FORMAT (A132)
  101 FORMAT('LEVEL     ','  ',A10,' ',I2,' ',I4,2F10.1,1X,I2)    
  102 FORMAT('LEVEL     ','  ',A10,' ',I2,' ',I4,F10.1,11X,I2)    

      END
