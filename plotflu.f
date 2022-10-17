      SUBROUTINE PLOTFLU (NF, XLAMBDA, EMFLUX, EMCOLI, MODHEAD, JOBNUM, 
     >                    KANAL, RSTAR, TOTOUT, KARTE)
C***********************************************************************
C***  DIRECT TRANSFER OF EMERGENT FLUX PLOT
C***    - DIFFERENT VERSIONS ACCORDING TO THE PARAMETER "SWITCH"
C***********************************************************************
      DIMENSION XLAMBDA(NF), EMFLUX(NF), EMCOLI(NF)
      INTEGER, PARAMETER ::  NFMAX = 5000
      DIMENSION X(NFMAX), Y(NFMAX)
      CHARACTER MODHEAD*100,HEADER*60, YTEXT*100, CENTER*8, ANG*2
      CHARACTER ACTPAR*10, CFAC*9, KARTE*(*)

C***  SWITCH BETWEEN PLOT OF LOG F-NUE / LOG F-LAMBDA / T-RAD
      CHARACTER(4) :: SWITCH
C***  SWITCH BETWEEN FLUX AT RSTAR OR AT 10pc distance 
      LOGICAL B10PC, BXMIN, BXMINERR, BXMAX, BXMAXERR

      REAL, PARAMETER :: CLIGHT = 2.9979E18     !CLIGHT = SPEED OF LIGHT (IN ANGSTROEM PER SECOND) 
      REAL, PARAMETER :: PI = 3.14159

      REAL, PARAMETER :: PARSEC = 3.08561E18    !1 PARSEC IN CM

      IF (TOTOUT == 6HUNDEF. ) THEN
         PRINT 7
    7    FORMAT (//'WARNING: OPTION "PLOT FLUX" IGNORED - ',
     >           'EMERGENT CONT. FLUX NOT YET CALCULATED ',//)
         RETURN
      ENDIF

C***  Default values
      SWITCH = 'TRAD'
      B10PC = .FALSE.
      BXMIN = .FALSE.
      BXMAX = .FALSE.
      XMIN = 2.
      XMAX = 4.

C***  Read further parameters 
      CALL SARGC(KARTE, NVAR)
      DO IVAR = 3, NVAR
         CALL SARGV(KARTE, IVAR, ACTPAR)
         IF (BXMIN) THEN
            READ (ACTPAR,'(F10.0)',ERR=90) VALUE
            XMIN = VALUE
            BXMIN    = .FALSE.
            BXMINERR = .FALSE.
   90       CONTINUE
         ELSE IF (BXMAX) THEN
            READ (ACTPAR,'(F10.0)',ERR=91) VALUE
            XMAX = VALUE
            BXMAX    = .FALSE.
            BXMAXERR = .FALSE.
   91       CONTINUE
         ELSE IF (ACTPAR .EQ. 'FLAM') THEN
            SWITCH = 'FLAM'
         ELSE IF (ACTPAR .EQ. 'LOGF') THEN
            SWITCH = 'LOGF'
         ELSE IF (ACTPAR .EQ. 'TRAD') THEN
            SWITCH = 'TRAD'
         ELSE IF (ACTPAR .EQ. '10PC') THEN
            B10PC = .TRUE.
         ELSE IF (ACTPAR .EQ. 'XMIN') THEN
            BXMIN = .TRUE.
            BXMINERR = .TRUE.
         ELSE IF (ACTPAR .EQ. 'XMAX') THEN
            BXMAX = .TRUE.
            BXMAXERR = .TRUE.
         ELSE
            WRITE (*,95) ACTPAR(:IDX(ACTPAR)), KARTE(:IDX(KARTE))
            WRITE (0,95) ACTPAR(:IDX(ACTPAR)), KARTE(:IDX(KARTE))
   95       FORMAT ('WARNING: UNRECOGNIZED PARAMETER: ', A, 
     >         ' ON INPUT LINE: ', A) 
         ENDIF
      ENDDO

      IF (BXMINERR) THEN
            WRITE (*,96) 'XMIN', KARTE(:IDX(KARTE))
            WRITE (0,96) 'XMIN', KARTE(:IDX(KARTE))
      ENDIF
      IF (BXMAXERR) THEN
            WRITE (*,96) 'XMAX', KARTE(:IDX(KARTE))
            WRITE (0,96) 'XMAX', KARTE(:IDX(KARTE))
      ENDIF
   96 FORMAT ('WARNING: ERROR WHEN DECODING VALUE OF ', A, 
     >        ' ON INPUT LINE: ', A) 
C************************

      CALL JSYMSET ('G2','TRANSFER')
      CALL REMARK ('FLUX PLOT DATA TO BE ROUTED')
      CENTER = CHAR(92) // 'CENTER' // CHAR(92)
      ANG    = CHAR(92) // 'A' 

      IF (NF .GT. NFMAX) THEN
        WRITE (KANAL, *)
     >   '* WARNING: DIMENSION NFMAX INSUFFICIENT IN PLOTFLU'
        WRITE (0, *)
     >   '* WARNING: DIMENSION NFMAX INSUFFICIENT IN PLOTFLU'
      ENDIF

      NFF = MIN0(NF, NFMAX)

C***  Log(Flux) value repacing for negative fluxes
      YMINLOG = -100. 


C***  CONSTRUCT HEADER LINE
      WRITE(UNIT=HEADER,FMT=3) MODHEAD(13:36), JOBNUM
    3 FORMAT ('EMERGENT FLUX OF MODEL ',A20,3X,'JOB ',I7)

C***  WRITE "INBOX" OPTION
      WRITE (KANAL,*) 'PLOT: ', HEADER
      WRITE (KANAL,*) 'KASDEF INBOX'

      IF (B10PC) THEN
C***  FACTOR OF CONVERSION FROM ASTROPHYSICAL FLUX AT RSTAR (EMFLUX) TO 
C***  PHYSICAL FLUX IN 10 PC
        FAC = PI * RSTAR*RSTAR / (PARSEC*PARSEC) / 100.
        FACLOG = ALOG10(FAC)
        CFAC = '  in 10pc'
      ELSE
        FAC = 1.
        FACLOG = .0
        CFAC = '    '
      ENDIF

C***  PLOT OF LOG F-NUE VERSUS LOG LAMBDA  -----------------------------
      IF (SWITCH .EQ. 'LOGF') THEN
      XSCALE=20./(XMAX-XMIN)
      XTICK=0.1
      XABST=0.5
 
      YMIN=-10.
      YMAX=0.
      IF (B10PC) THEN
        YMIN = YMIN - 16.
        YMAX = YMAX - 16.
      ENDIF
      YSCALE=15./(YMAX-YMIN)
      YTICK=1.
      YABST=5.

      DO 4 K=1, NFF
      X(K) = ALOG10(XLAMBDA(K))
      IF (EMFLUX(K) .LE. .0) THEN
         Y(K) = YMINLOG
         ELSE
         Y(K) = ALOG10(EMFLUX(K)) + FACLOG
         ENDIF
    4 CONTINUE

       YTEXT =  CENTER // 
     >     'log  F&T#n#&M / (erg cm&H-2&M s&H-1&M Hz&H-1&M)' // CFAC

      CALL PLOTANF (KANAL,HEADER,HEADER
     $ ,CENTER//'log (#l#/'//ANG//')'
     > ,YTEXT
     $ ,XSCALE,XMIN,XMAX,XTICK,XABST,.0
     $ ,YSCALE,YMIN,YMAX,YTICK,YABST,.0
     $ ,X, Y, NFF, 5)

      DO K=1, NFF
         IF (EMCOLI(K) .LE. .0) THEN
            Y(K) = YMINLOG
         ELSE
            Y(K) = ALOG10(EMCOLI(K)) + FACLOG
         ENDIF
      ENDDO
      CALL PLOTCONS (KANAL, X, Y, NFF, 'COLOR=2')       


C***  PLOT OF LOG F-LAMBDA VERSUS LOG LAMBDA  -----------------------------
      ELSE IF (SWITCH .EQ. 'FLAM') THEN
      XSCALE=20./(XMAX-XMIN)
      XTICK= 0.1
      XABST= 0.5
 
      YMIN = 3.
      YMAX = 13.
      IF (B10PC) THEN
        YMIN = YMIN - 16.
        YMAX = YMAX - 16.
      ENDIF
      YSCALE=15./(YMAX-YMIN)
      YTICK=1.
      YABST=5.

      DO 5 K=1, NFF
      X(K) = ALOG10 (XLAMBDA(K))
      IF (EMFLUX(K) .LE. .0) THEN
         Y(K) = YMINLOG
         ELSE
         Y(K) = ALOG10 (EMFLUX(K)) 
     >          + ALOG10(CLIGHT / XLAMBDA(K) / XLAMBDA(K)) + FACLOG
         ENDIF
    5 CONTINUE 
 
      YTEXT = CENTER // 
     > 'log  F&T#l#&M / (erg cm&H-2&M s&H-1&M '//ANG//'&H-1&M)' // CFAC

      CALL PLOTANF (KANAL,HEADER,HEADER
     $ ,CENTER//'log (#l# / '//ANG//')'
     > ,YTEXT
     $ ,XSCALE,XMIN,XMAX,XTICK,XABST,.0
     $ ,YSCALE,YMIN,YMAX,YTICK,YABST,.0
     $ , X, Y, NFF, 5)

      DO K=1, NFF
         IF (EMCOLI(K) .LE. .0) THEN
            Y(K) = YMINLOG
         ELSE
            Y(K) = ALOG10(EMCOLI(K)) 
     >             + FACLOG + ALOG10(CLIGHT / XLAMBDA(K) / XLAMBDA(K)) 
         ENDIF
      ENDDO
      CALL PLOTCONS (KANAL, X, Y, NFF, 'COLOR=2')       


      ELSE

C***  PLOT OF T-RAD VERSUS LOG LAMBDA  ---------------------------------
      XSCALE=20./(XMAX-XMIN)
      XTICK=0.5
      XABST=1.

      DO 1 K=1, NFF
      X(K)=ALOG10(XLAMBDA(K))
      Y(K)=TRADFUN(XLAMBDA(K),EMFLUX(K))/1000.
    1 CONTINUE
 
      YMIN=Y(1)
      DO 2 K=2,NFF
      IF (YMIN .GT. Y(K)) YMIN=Y(K)
    2 CONTINUE

      YMIN=INT(YMIN/5.)*5.
      YMAX=YMIN+75.
      YSCALE=.2
      YTICK=5.
      YABST=10.
 
      CALL PLOTANF (KANAL,HEADER,HEADER
     $ ,CENTER//'log (#l# / '//ANG//')'
     $ ,CENTER//'T&Trad&M / kK'
     $ ,XSCALE,XMIN,XMAX,XTICK,XABST,.0
     $ ,YSCALE,YMIN,YMAX,YTICK,YABST,.0
     $ ,X, Y, NFF, 5)

      DO K=1, NFF
         IF (EMCOLI(K) .LE. .0) THEN
            Y(K) = YMINLOG
         ELSE
            Y(K)=TRADFUN(XLAMBDA(K),EMCOLI(K))/1000.
         ENDIF
      ENDDO
      CALL PLOTCONS (KANAL, X, Y, NFF, 'COLOR=2')       

      ENDIF

      RETURN
      END
