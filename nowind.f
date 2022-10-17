      SUBROUTINE NOWIND (NOWIND_LINE, RCON, NATOM, ATMASS, ABXYZ, 
     >                   ND, RNE, T, RADIUS, VELO, VDOP, TAURCONT, 
     >                   ND_MERGED, RADIUS_MERGED, PGRID_MERGED,
     >                   NP_MERGED, ZGRID_MERGED, IERR)

C***********************************************************************
C*** "NOWIND" CARD allows to calculate the mergent spectrum 
C*** as if the (outer part of) the wind would not exist
C***     The wind is either cut due to VELO, RADIUS or TAU 
C***  In case of SECONDMODEL, the criteria refer to the the first model 
C***  Note that the _MERGED input parameters are over-written
C***********************************************************************

      IMPLICIT NONE
      CHARACTER NOWIND_LINE*(*), ACTPAR*20
      INTEGER NPAR, NATOM, ND, L, NA, IDX, IERR
      INTEGER ND_MERGED, NP_MERGED, NCORE, JP, JMAX, I
      INTEGER INDCUT, ISRCHFLT

      REAL RCUT, RCON, XMUL, VSONIC, VDOP, VCUT, ATMEAN, TAUCUT
      REAL RR, PJ, PJPJ
      REAL ATMASS(NATOM), ABXYZ(NATOM)
      REAL RNE(ND), T(ND), RADIUS(ND), VELO(ND), TAURCONT(ND)
      REAL RADIUS_MERGED(ND_MERGED), PGRID_MERGED(NP_MERGED)
      REAL ZGRID_MERGED(ND_MERGED*NP_MERGED)

      REAL, PARAMETER :: RGAS = 8.3145E7        !Gas Constant (CGS)

C****************************************************************
C***  First part: defining RCUT based on structure of first MODEL 
C****************************************************************

      CALL SARGC (NOWIND_LINE, NPAR)

C***  No parameters (default): cutting at connection radius 
      IF (NPAR == 1) THEN                  
         RCUT = RCON
         GOTO 20
      ENDIF

      CALL SARGV(NOWIND_LINE, 2, ACTPAR)
      IF (ACTPAR == 'OFF') THEN
         NOWIND_LINE = 'NONE'
         RETURN

      ELSEIF (ACTPAR == 'RADIUS') THEN
      IF (NPAR .LT.  3) GOTO 90                  
         CALL SARGV(NOWIND_LINE, 3, ACTPAR)
         READ (ACTPAR, '(F10.0)', ERR=91) RCUT
         GOTO 20
      ELSEIF (ACTPAR == 'VDOP') THEN
         VCUT = VDOP
         GOTO 10

      ELSEIF (ACTPAR == 'TAU') THEN
      IF (NPAR .LT.  3) GOTO 90                  
         CALL SARGV(NOWIND_LINE, 3, ACTPAR)
         READ (ACTPAR, '(F10.0)', ERR=91) TAUCUT
         GOTO 25

C***  For compatibility with older syntax, allow keywords
C***      VDOP or SONIC without preceeding VEL0=
      ELSEIF (ACTPAR == 'VDOP') THEN
         VCUT = VDOP
         GOTO 10

      ELSEIF (ACTPAR == 'SONIC') THEN
C***     Need to calculate sonic velocity  
         ATMEAN = 0.
         DO NA=1, NATOM
            ATMEAN = ATMEAN + ABXYZ(NA) * ATMASS(NA)
         ENDDO
         DO L=1, ND
            XMUL = ATMEAN / (1. + RNE(L))
            VSONIC = SQRT(RGAS * T(L)/ XMUL) * 1.E-5
            IF (VELO(L) < VSONIC) EXIT
         ENDDO
         VCUT = VSONIC
         GOTO 10

      ELSEIF (ACTPAR == 'VELO') THEN
         IF (NPAR .LT.  3) GOTO 90                  
         CALL SARGV(NOWIND_LINE, 3, ACTPAR)
         IF (ACTPAR == 'VDOP') THEN
            VCUT = VDOP
         ELSEIF (ACTPAR == 'SONIC') THEN
C***        Need to calculate sonic velocity  
            ATMEAN = 0.
            DO NA=1, NATOM
               ATMEAN = ATMEAN + ABXYZ(NA) * ATMASS(NA)
            ENDDO
            DO L=1, ND
               XMUL = ATMEAN / (1. + RNE(L))
               VSONIC = SQRT(RGAS * T(L)/ XMUL) * 1.E-5
               IF (VELO(L) < VSONIC) EXIT
            ENDDO
            VCUT = VSONIC

         ELSE
            READ (ACTPAR, '(F10.0)', ERR=91) VCUT
         ENDIF
         GOTO 10

C***  Assume that the second parameter is the cut velocity value
      ELSE
         READ (ACTPAR, '(F10.0)', ERR=91) VCUT
         GOTO 10
      ENDIF

C***  velocity VCUT was specified
   10 CONTINUE
C***  Test if VCUT is within valid range
      IF (VCUT .GE. VELO(1)) THEN
         WRITE (0,*) '**** ERROR: VCUT .GE. max. velocity'
         GOTO 99
      ENDIF 
      IF (VCUT .LE. VELO(ND)) THEN
         WRITE (0,*) '**** ERROR: VCUT .LE. min. velocity'
         GOTO 99
      ENDIF 
      CALL SPLINPO (RCUT,   VCUT, RADIUS,   VELO, ND) 
      CALL SPLINPO (TAUCUT, VCUT, TAURCONT, VELO, ND) 
      GOTO 30

C***  If RADIUS was specified:
   20 CONTINUE
C***  Test if RCUT is within valid range
      IF (RCUT .GE. RADIUS(1)) THEN
         WRITE (0,*) '**** ERROR: RCUT .GE. RMAX'
         GOTO 99
      ENDIF 
      IF (RCUT .LE. 1.) THEN
         WRITE (0,*) '**** ERROR: RCUT .LE. 1.'
         GOTO 99
      ENDIF 
      CALL SPLINPO (VCUT,   RCUT, VELO,     RADIUS, ND) 
      CALL SPLINPO (TAUCUT, RCUT, TAURCONT, RADIUS, ND) 
      GOTO 30

C***  If TAU was specified:
   25 CONTINUE
C***  Test if TAUCUT is within valid range
      IF (TAUCUT .LE. TAURCONT(1)) THEN
         WRITE (0,*) '**** ERROR: TAUCUT .LE. TAURCONT(1)'
         GOTO 99
      ENDIF 
      IF (TAUCUT .GE. TAURCONT(ND)) THEN
         WRITE (0,*) '**** ERROR: TAUCUT .GE. TAUMAX'
         GOTO 99
      ENDIF 
      CALL SPLINPO (VCUT, TAUCUT, VELO,     TAURCONT, ND) 
      CALL SPLINPO (RCUT, TAUCUT, RADIUS,   TAURCONT, ND) 
      GOTO 30

    
   30 Continue
      WRITE(0,80) NOWIND_LINE(:IDX(NOWIND_LINE)), VCUT, RCUT, TAUCUT
      WRITE(*,80) NOWIND_LINE(:IDX(NOWIND_LINE)), VCUT, RCUT, TAUCUT
   80 FORMAT (A,/,'Thus removing layers above:',/,
     >        '     VCUT   =', F9.3, ' km/s',/,
     >        '     RCUT   =', F9.3, ' Rstar',/,
     >        '     TAUCUT =', F9.3)

C****************************************************************
C***  Second part: usining RCUT to shorten RADIUS_MERGED 
C****************************************************************

      IF (RCUT .GE. RADIUS_MERGED(1)) THEN
         WRITE (0,*) '**** ERROR: RCUT .GE. RMAX of second model'
         GOTO 99
      ENDIF

C***  Now define the new RADIUS_MERGED vector
      NCORE = NP_MERGED - ND_MERGED
      RADIUS_MERGED(1) = RCUT
      INDCUT =  ISRCHFLT(ND_MERGED, RADIUS_MERGED, 1, RCUT) - 1
      DO L= INDCUT+1, ND_MERGED
         RADIUS_MERGED(L-INDCUT+1) = RADIUS_MERGED(L)
      ENDDO  
      WRITE (0,'(A,I4)') 
     >     'Number of radial points (total wind) ND=', ND_MERGED
      ND_MERGED = ND_MERGED - INDCUT + 1
      NP_MERGED = NP_MERGED - INDCUT + 1
      WRITE (0,'(A,I4)') 
     >     'Reduced by NOWIND to                 ND=', ND_MERGED

C***  Re-establish PGRID and ZGRID according to the new RADIUS vector
      DO L=1, ND_MERGED
         JP = NP_MERGED + 1  - L
         PGRID_MERGED(JP) = RADIUS_MERGED(L)
      ENDDO

      DO L = 1, ND_MERGED
        RR = RADIUS_MERGED(L) * RADIUS_MERGED(L)
        JMAX = NP_MERGED + 1 - L
        DO JP = 1, JMAX
          PJ = PGRID_MERGED(JP)
          PJPJ = PJ * PJ
          I = (JP-1)*ND_MERGED + L
          IF ( (RR-PJPJ) .GT. .0) THEN
             ZGRID_MERGED(I) = SQRT(RR-PJPJ)
          ELSE
             ZGRID_MERGED(I) = .0
          ENDIF
        ENDDO
c        write (0,*) L, RADIUS_MERGED(L)
      ENDDO

c      do jp=1, np_merged 
c        write (0,*) jp, pgrid_MERGED(jp)
c      enddo

      RETURN

C***  Error branches
   90 WRITE (0,*) '*** ERROR: This option requieres a value'
      GOTO 99

   91 WRITE (0,*) '*** ERROR when decoding floating-point number'
      GOTO 99

   99 WRITE (0,*) 'THE ERROR OCCURED IN THE FOLLOWING LINE:'
      WRITE (0,*) NOWIND_LINE(:IDX(NOWIND_LINE))
      WRITE (0,*) 'Therefore, this spectral range must be skipped!'
      IERR = 1
      RETURN

      END
