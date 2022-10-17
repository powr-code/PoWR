      SUBROUTINE PLOTOPA (ND, MODHEAD, NATOM, ATMASS, ABXYZ, 
     >           FILLFAC, OPTIONPLOTOPA, NPLOTOPA, 
     >           T,RNE,POPNUM,ENTOT,RSTAR,
     $           OPA,ETA,THOMSON,IWARN,MAINPRO,MAINLEV,
     $           NOM,KODAT,NDIM,N,MAXATOM,
     $           LEVEL,NCHARG,WEIGHT,ELEVEL,EION,EINST,
     $           SIGMATHK,SEXPOK,EDGEK,
     $           K,NF,SIGMAKI,RADIUS,
     $           KONTNUP,KONTLOW,LASTKON,XDATA)
C***********************************************************************
C***  PLOT of opacity over radius at given wavelength
C***  CALLED FROM: COMO
C***  Option: PLOT OPA LAMBDA = x.x [PERMASS]
C***  wrh + adriane, 10-Nov-2005 18:26:39
C***********************************************************************

      DIMENSION OPA(ND), RADIUS(ND), FILLFAC(ND), ENTOT(ND)
      DIMENSION ABXYZ(MAXATOM), ATMASS(MAXATOM)
      CHARACTER MODHEAD*100, YTEXT*60
      CHARACTER OPTIONPLOTOPA*(80)(NPLOTOPA), XLAMSTR*20, ACTPAR*20
      LOGICAL BPERMASS
C***  Atomic Mass Unit
      DATA AMU /1.66E-24/

C**  OPACITY PLOTS
        OPEN (2, FILE='PLOT', STATUS='UNKNOWN', POSITION='APPEND')
        CALL JSYMSET ('G2','TRANSFER')
      DO IPLOTOPA=1, NPLOTOPA 
        XLAMSTR = 'UNDEF'
        BPERMASS = .FALSE.
        CALL SARGC (OPTIONPLOTOPA(IPLOTOPA), NPAR)
        DO I=3, NPAR
          CALL SARGV (OPTIONPLOTOPA(IPLOTOPA), I, ACTPAR)
          IF (ACTPAR .EQ. 'LAMBDA') THEN
             IF (NPAR .GT. I) THEN
                CALL SARGV (OPTIONPLOTOPA(IPLOTOPA), I+1, XLAMSTR)
             ELSE
                GOTO 99
             ENDIF
          ELSEIF  (ACTPAR .EQ. 'PERMASS') THEN
             BPERMASS = .TRUE.
          ENDIF
          READ (XLAMSTR,'(F10.0)',ERR=99) XLAM
        ENDDO

        GOTO 100
   99   WRITE (*,*) 'ERROR in option PLOT OPA: ', 
     >     'wavelength could not be decoded, plot option skipped' 
        WRITE (*,*) 'The option was: ' // 
     >     OPTIONPLOTOPA(IPLOTOPA)(:IDX(OPTIONPLOTOPA))
       WRITE (0,*) 'ERROR in option PLOT OPA: ', 
     >     'wavelength could not be decoded, plot option skipped' 
        WRITE (0,*) 'The option was: ' // 
     >     OPTIONPLOTOPA(IPLOTOPA)(:IDX(OPTIONPLOTOPA))
        CYCLE
  100   CONTINUE  

      CALL COOP (XLAM, ND,T,RNE,POPNUM,ENTOT,RSTAR,
     $           OPA,ETA,THOMSON,IWARN,MAINPRO,MAINLEV,
     $           NOM,KODAT,NDIM,N,MAXATOM,
     $           LEVEL,NCHARG,WEIGHT,ELEVEL,EION,EINST,
     $           DUMMY,DUMMY,DUMMY,
     $           DUMMY,DUMMY,DUMMY,
     $           SIGMATHK,SEXPOK,EDGEK,
     $           K,NF,SIGMAKI,RADIUS,
     $           KONTNUP,KONTLOW,LASTKON,XDATA)

C*** Mass absorption coefficient
        IF (BPERMASS) THEN
C***       ATMEAN = MEAN ATOMIC WEIGHT ( IN AMU )
           ATMEAN=0.
           DO NA=1,NATOM
              ATMEAN=ATMEAN+ABXYZ(NA)*ATMASS(NA)
           ENDDO
C***       Opacity was calculated with ENTOT = density in the clumps
C***       --> must be divided by mass density in the clumps
C***       Units of OPA must be converted from [per Rstar] into [per cm]
           DO L=1, ND
              RHO = AMU * ATMEAN * ENTOT(L) 
              OPA(L) = OPA(L) / RSTAR / RHO
           ENDDO
        ELSE
C***       Opacity was calculated with ENTOT = density in the clumps
C***       --> average opacity must account for the voids
           OPA(L) = OPA(L) * FILLFAC(L) 
        ENDIF

 
C***     AUTO-SCALING by WRplot
         XMIN = .0
         XMAX = .0
 
         XTICK = 0.1
         XABST = 1.
         YTICK = 10.
         YABST = 50.
         KANAL = 2
         IF (BPERMASS) THEN
          YTEXT = '\CENTER\Mass absorption coeff. [cm&H2&M/g]' 
         ELSE
          YTEXT = '\CENTER\Opacity [per R\*]' 
         ENDIF
         CALL PLOTANFS (KANAL, 
     >        OPTIONPLOTOPA(IPLOTOPA)(:IDX(OPTIONPLOTOPA(IPLOTOPA))), 
     >        '&E'//MODHEAD,
     >        '\CENTER\Radius [R\*]',
     >        YTEXT(:IDX(YTEXT)) // ' at ' 
     >           // XLAMSTR(:IDX(XLAMSTR)) // ' Ang',
     >        0., XMIN, XMAX, XTICK, XABST,.0,
     >        0., 1.,   RADIUS(1), YTICK, YABST, .0,
     >        RADIUS, OPA, ND, 'PEN=3, COLOR=2')

      ENDDO




      RETURN
      END
