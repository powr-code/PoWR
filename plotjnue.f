      SUBROUTINE PLOTJNUE (KANAL, MODHEAD, JOBNUM, ND, NF, XLAMBDA,
     >     PLOTOPT, MAXPLOTN, XPLOT, YPLOT, XJC )
 
C*****************************************************************
C***  PLOT JNUE 
C*****************************************************************
C***  STEBOL = STEFAN-BOLTZMANN CONSTANT / PI (ERG/CM**2/SEC/STERAD/KELVIN**4
      DATA STEBOL /1.8046E-5/
      CHARACTER MODHEAD*100, HEADLINE*110, CENTER*8 

      DIMENSION XPLOT(MAXPLOTN), YPLOT(MAXPLOTN) 
      DIMENSION XLAMBDA(NF), XJC(ND,NF)
      CHARACTER*(*) PLOTOPT
      CHARACTER ACTPAR*20, XTEXT*100, YTEXT*100, MODE*20
      LOGICAL BTRAD

      CENTER = CHAR(92) // 'CENTER' // CHAR(92)
      HEADLINE = 'M'//MODHEAD(13:)
      WRITE (HEADLINE(90:), FMT='(A8,I7)') ' JOB No.', JOBNUM

      BTRAD = .FALSE.

C***  JNUE PLOT (Continuum)
       CALL SARGC (PLOTOPT,NPAR)
       IF (NPAR .LT. 4) GOTO 91
       CALL SARGV (PLOTOPT,3,MODE)

       IF (MODE .EQ. 'L') THEN
          IF (NF .GT. MAXPLOTN) GOTO 92
          CALL SARGV (PLOTOPT,4,ACTPAR)
          READ (ACTPAR, '(I10)') L
          IF (L .LT. 1 .OR. L .GT. ND) GOTO 91           
       ELSEIF (MODE .EQ. 'K') THEN
          IF (ND .GT. MAXPLOTN) GOTO 92
          CALL SARGV (PLOTOPT,4,ACTPAR)
          READ (ACTPAR, '(I10)') K
          IF (K .LT. 1 .OR. K .GT. NF) GOTO 91           
       ELSE
          GOTO 91
       ENDIF

       IF (NPAR .GE. 5) THEN
          CALL SARGV (PLOTOPT,5,ACTPAR)
          BTRAD = ACTPAR .EQ. 'TRAD'
       ENDIF

       IF (BTRAD) THEN 
          YTEXT = CENTER // 'J&T#n#&M expressed as T&Trad&M / kK' 
       ELSE
          YTEXT = CENTER // 'J&T#n#&M (erg cm&H-2&M s&H-1&M Hz&H-1&M)' 
       ENDIF

       IF (MODE .EQ. 'L') THEN 
          XTEXT = CENTER // 'log #l# / \A'
          WRITE (YTEXT, '(2A,I5)') YTEXT(:IDX(YTEXT)), 
     >                             ' at depth index L =', L
          DO K=1, NF
             LAST = K-1
             IF (XLAMBDA(K) .GT. 1.E5) EXIT
             XPLOT(K) = ALOG10(XLAMBDA(K))
             IF (BTRAD) THEN 
                YPLOT(K) = TRADFUN (XLAMBDA(K), XJC(L,K)) / 1000.
             ELSE
                YPLOT(K) = XJC(L,K)
             ENDIF
          ENDDO
       ELSE
          XTEXT = CENTER // 'Depth Index L'
          WRITE (YTEXT, '(2A,I5)') YTEXT(:IDX(YTEXT)), 
     >                             ' at freq. index K =', K
          DO L=1, ND
             XPLOT(L) = FLOAT(L)
             IF (BTRAD) THEN 
                YPLOT(L) = TRADFUN (XLAMBDA(K), XJC(L,K)) / 1000.
             ELSE
                YPLOT(L) = XJC(L,K)
             ENDIF
          ENDDO
          LAST = ND
       ENDIF
       CALL PLOTANFS (KANAL,HEADLINE, '&E'//HEADLINE,
     >        XTEXT, YTEXT, 
     >        .0, .0, .0, .0, .0, .0,
     >        .0, .0, .0, .0, .0, .0,
     $        XPLOT, YPLOT, LAST, 'PEN=1 COLOR=2')
       CALL JSYMSET ('G2','TRANSFER')



      GOTO 10

C***  ERROR BRANCHES
   91 WRITE (0, '(A,/,A,/,A)') 'PLOTOPTIONS: WARNING: ' // 
     >      'L=... OR K=... (DEPTH OR FREQ. INDEX) MUST BE SPECIFIED',
     >      'THE FOLLOWING PLOT OPTION IS IGNORED:',
     >      PLOTOPT(:IDX(PLOTOPT))
      GOTO 10

   92 WRITE (0, '(A,/,A,/,A)') 'PLOTOPTIONS: WARNING: ' //
     >      'DIMENSION MAXPLOTN INSUFFICIENT',
     >      'THE FOLLOWING PLOT OPTION IS IGNORED:',
     >      PLOTOPT(:IDX(PLOTOPT))
      GOTO 10

   10 CONTINUE

      RETURN
      END
