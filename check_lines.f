      SUBROUTINE CHECK_LINES(XK, XKMIN, XKMID, XKMAX, 
     >             LINECHECK, ILINECHECK, NLINE, MAXLIN, LEVEL,
     >             LIND, LINDS, WS, BPLOT, RADIUS, NLACT, LINE, 
C***  for PRELINECL
     >             NUP, LOW, N, XLAM, NDIM, ND, XJLMEAN, ELEVEL, 
     >             INDNUP, INDLOW, NDDIM, 
C***  and also for LIOP
     >             EINST, WEIGHT, XLAMSOR, ENTOT, POPNUM, RSTAR, 
     >             OPAL, ETAL, VDOP)
C****************************************************************
C***  Adds new Lines to the List of Active Lines
C***    The Line Opacities are prepared
C***  Deletes Lines which are finished from the List
C***    The J-bars are writted to the Model File
C***
C***    Called by COLI
C****************************************************************

      DIMENSION XKMIN(NLINE), XKMID(NLINE), XKMAX(NLINE), LINE(NLINE)
      DIMENSION LIND(MAXLIN), LINDS(MAXLIN), WS(MAXLIN)
      DIMENSION NUP(MAXLIN), LOW(MAXLIN), XLAM(MAXLIN)
      DIMENSION XLAMSOR(NLINE)
      DIMENSION RADIUS(ND)
      DIMENSION EINST(NDIM,NDIM), WEIGHT(NDIM)
      DIMENSION XJLMEAN(NDDIM,MAXLIN)
      DIMENSION OPAL(NDDIM,MAXLIN), ETAL(NDDIM,MAXLIN)

      LOGICAL BPLOT

      CHARACTER NAME*8
      CHARACTER LEVEL(NDDIM)*10

C***  Which lines are active?
C***    Add new lines
        DO
          IF (XK .LT. XKMIN(LINECHECK) .OR. LINECHECK .GT. NLINE) EXIT
          IF (XK .GE. XKMIN(LINECHECK) .AND.
     >        XK .LE. XKMAX(LINECHECK)) THEN
C***  Search first free index
            DO NLNEW=1, MAXLIN
              IF (LIND(NLNEW) .EQ. 0) THEN
                LIND(NLNEW) = ILINECHECK
                LINDS(NLNEW) = LINECHECK
                WS(NLNEW) = 0.
C***  Prepare line quantities
                CALL PRELINE (NUP(NLNEW), LOW(NLNEW), ILINECHECK, N, 
     >                 XLAM(NLNEW), ND, XJLMEAN(1,NLNEW), 
     >                 ELEVEL, NDIM, INDNUP, INDLOW)
                CALL LIOP (EINST(NUP(NLNEW),LOW(NLNEW)), 
     >                 WEIGHT(LOW(NLNEW)), 
     >                 WEIGHT(NUP(NLNEW)), LOW(NLNEW), NUP(NLNEW), ND, 
     >                 XLAMSOR(LINECHECK), ENTOT, POPNUM, RSTAR, 
     >                 OPAL(1,NLNEW), ETAL(1,NLNEW), VDOP)
                IF (BPLOT) THEN
                  WRITE (38,'(A,F9.2,1X,F3.0,1X,A,I4,1X,I4,F11.3)') 
     >              'KASDEF LUN ', XK, FLOAT(NLNEW), '0. 0.1 0.2 &E', 
     >              LINECHECK, LINE(LINECHECK), 
     >              XLAMSOR(LINECHECK)
                  WRITE (38,'(A,F9.2,A)')
     >              'KASDEF SYM ', XKMID(LINECHECK), 
     >              ' 0.0 0. 0.2 -0.2 8'
                  WRITE (38,'(A,F9.2,1X,F5.1,A)')
     >              'KASDEF SYM ', XKMID(LINECHECK), 
     >              FLOAT(NLNEW), ' 0. -0.2 -0.2 8'
                ENDIF
                NLACT = NLACT + 1
                IF (NLACT .GT. MAXLIN) THEN
                  WRITE (0,*) 'Capacity of MAXLIN exceeded'
                  STOP 'ERROR in Subr. COLI'
                ENDIF
                EXIT
              ENDIF
            ENDDO
          ENDIF
          LINECHECK = LINECHECK + 1
          ILINECHECK = LINE(LINECHECK)
        ENDDO
C***    Check the old lines
        DO NL=1, MAXLIN
          LACT = LIND(NL)
          LACTS = LINDS(NL)
          IF (LACT .EQ. 0) CYCLE
          IF (XKMAX(LACTS) .LT. XK) THEN
C***  Line has been finished -> Store XJL
            DO L=1, ND
              XJLMEAN(L,NL) = XJLMEAN(L,NL) / WS(NL)
     >                         /RADIUS(L)/RADIUS(L)
            ENDDO
            IF (LACT <= 9999) THEN
              WRITE (NAME,'(A3,I4,A1)') 'XJL',LACT,' '
            ELSE
              WRITE (NAME,'(A3,I5)') 'XJL',LACT
            ENDIF
            CALL WRITMS 
     >        (3, XJLMEAN(1,NL), ND, NAME, -1, IDUMMY, IERR)
            LIND(NL) = 0
            NLACT = NLACT - 1
          ELSE
            IF (BPLOT) THEN
              WRITE (40+NL,'(I8,1X,F3.0)') NINT(XK), FLOAT(NL)
            ENDIF
          ENDIF
        ENDDO

      RETURN
      END
