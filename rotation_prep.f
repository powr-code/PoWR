      SUBROUTINE ROTATION_PREP (KARTE, RADIUS, VELO, TAURCONT, 
     >                    ND, RCOROT, NPHI, VSIDU, XMAX,
     >                    ND_MERGED, RADIUS_MERGED, 
     >                    NP_MERGED, PGRID_MERGED, ZGRID_MERGED, 
     >                    NDDIM, NPDIM, NPHIMAX, PHIWEIGHT, PHIARR,
     >                    BVSINI_AT_RCOROT, DX_WINDROT)

C*************************************************************************     
C***  Prepares wind-rotation formalism 
C***  This subroutine can be skipped if wind-rotation is not requested, 
C***    i.e. if no VSINI > .0 has been specified
C***  In the first part, the input line is analyzed on basis of model1
C*************************************************************************   

      IMPLICIT NONE     
      CHARACTER KARTE*(*), ACTPAR*20, LINE*200, ACTPAR2*40
      REAL, DIMENSION(NDDIM) ::  RADIUS, VELO, TAURCONT
      REAL RADIUS_MERGED(NDDIM)
      INTEGER, DIMENSION(NPDIM) :: NPHI, LPHIEND
      REAL, DIMENSION(NPDIM) :: PGRID_MERGED
      REAL, DIMENSION(NDDIM*NPDIM) :: ZGRID_MERGED
      REAL, DIMENSION(NPHIMAX,NPDIM) :: PHIWEIGHT, PHIARR
      LOGICAL BVSINI_AT_RCOROT      

C***  because of IMPLICIT NONE:
      REAL PI, RCOROT, RCOROT_V, RCOROT_TAU, XMAX, VSIDU, DPMAX
      REAL DX_WINDROT, RR, COSPHI, PJ, PJPJ, VSIDU_AT_P, DELTACOSPHI
      INTEGER IDX, NCORE, NCORE_ADD, NCORE_NEW, NP_MERGED, NP_NEW
      INTEGER NDDIM, NPDIM, NPHIMAX, ND, ND_MERGED, JP, L, I, J, JMAX 
      INTEGER LPHI, MAXNPHI, NPHISUM, MAXJP, NPHI_MEAN

      DATA PI / 3.14159265358979 /
      
C***********************************************************************
C***  Part 1:interpret the RCOROT line; 
C***    this is done on basis of MODEL1 
C***  alternative parameters are:
C***    RSTAR = x.x
C***    TAUROSS = x.x
C***    VELOKMS = x.x 
C***********************************************************************
      CALL SARGV(KARTE, 1, ACTPAR)
C***  In every case, RCOROT in all units (velocity, optical depth, 
C***    star radius) is interpolated from the given one.      
      IF (ACTPAR .NE. 'RCOROT') THEN
        WRITE (0,'(2A)') 'RCOROT not specified in FORMAL_CARDS, ', 
     >            'adopting the default: RCOROT RSTAR=1.0' 
        RCOROT = 1.        
        CALL SPLINPO (RCOROT_V, RCOROT, VELO, RADIUS, ND)
        CALL SPLINPO (RCOROT_TAU, RCOROT, TAURCONT, RADIUS, ND)   
      ELSE
        CALL SARGV(KARTE, 2, ACTPAR)
        IF (ACTPAR .EQ. 'VELOKMS') THEN            
                CALL SARGV(KARTE, 3, ACTPAR)
                READ (ACTPAR, '(F20.0)', ERR = 99) RCOROT_V
                CALL SPLINPO (RCOROT, RCOROT_V, RADIUS, VELO, ND)
                CALL SPLINPO (RCOROT_TAU, RCOROT_V, TAURCONT, VELO, ND)
        ELSE IF (ACTPAR .EQ. 'TAUROSS') THEN
                CALL SARGV(KARTE, 3, ACTPAR)
                READ (ACTPAR, '(F20.0)', ERR = 99) RCOROT_TAU
                CALL SPLINPO (RCOROT, RCOROT_TAU, RADIUS, TAURCONT, ND)
                CALL SPLINPO (RCOROT_V, RCOROT_TAU, VELO, TAURCONT, ND)
        ELSE IF (ACTPAR .EQ. 'RSTAR') THEN
                CALL SARGV(KARTE, 3, ACTPAR)
                READ (ACTPAR, '(F20.0)', ERR = 99) RCOROT
                CALL SPLINPO (RCOROT_V, RCOROT, VELO, RADIUS, ND)
                CALL SPLINPO (RCOROT_TAU, RCOROT, TAURCONT, RADIUS, ND)
        ELSE
                WRITE(0,*) '*** ERROR: RCOROT specification invalid'
                WRITE(0,*) '*** Allowed keywords:  VELOKMS | RSTAR |', 
     >                     'TAUROSS'
                GOTO 99
        ENDIF
  
       WRITE (0, '(A)') 'Corotation radius specified by input line:'
       WRITE (0, '(A)') KARTE(:IDX(KARTE))
       WRITE (0, '(A,G15.5)') 
     >      'Corresponding radius [in Rstar]: ', RCOROT
       WRITE (0, '(A,G15.5)') 
     >      'Corresponding velocity: [in km/s]: ', RCOROT_V
       WRITE (0, '(A,G15.5,/)') 
     >      'Corresponding optical depth (continuum Rosseland mean):  ',  
     >      RCOROT_TAU  

      ENDIF



C***********************************************************************
C***  Part 2:
C***  Add more impact parameters if necessary to resolve rotation 
C***    profile with DX_WINDROT Doppler units
C***  This id done on basis of the _MERGED arrays, which may differ
C***    from model 1 in case of SECOND_MODEL
C***********************************************************************

C***  XMAX should be big enough to contain potential rotational broadening 
      XMAX = AMAX1(XMAX, VSIDU)

C***  First, express Delta-P in Doppler units
      DPMAX = VSIDU / DX_WINDROT 
      IF (BVSINI_AT_RCOROT) DPMAX = DPMAX / RCOROT 

C***  Requested number of core rays
      NCORE_NEW = INT(DPMAX) + 1  

C***  Was the original number of core-rays sufficient?
      NCORE = NP_MERGED - ND_MERGED
      IF (NCORE_NEW .LE. NCORE) THEN
        WRITE (0,'(A,I4)') 
     >   'Number of core impact parameters is sufficient: NCORE=', NCORE
        GOTO 1
      ENDIF

C***  Not sufficient: determine the necessary number NCORE_NEW
      NCORE_ADD =  NCORE_NEW - NCORE
      NP_NEW = NP_MERGED + NCORE_ADD
      IF (NP_NEW .GT. NPDIM) THEN
          WRITE(0,'(A,I4)') '*** ERROR: insufficient NPDIM=', NPDIM
          WRITE(0,'(A,I4)') 'Needed for wind rotation: NPDIM=', NP_NEW
          STOP '*** FATAL ERROR in subroutine ROTATION_PREP'
      ENDIF

      WRITE (0,'(A,I4,A,I4)') 
     >       'Number of core impact parameters enhanced from ', NCORE,
     >       ' to', NCORE_NEW

C***  Shift the non-core impact parameters 
      DO JP=NP_MERGED, NCORE, -1
         PGRID_MERGED(JP+NCORE_ADD) = PGRID_MERGED(JP)
      ENDDO

C***  Calculate equally spaced core impact parameters 
      DO JP=2, NCORE_NEW
         PGRID_MERGED(JP) = (JP-1)/FLOAT(NCORE_NEW)
      ENDDO
         
      NP_MERGED = NP_NEW

C***  Z-array must be updated with new P-vector
      DO L=1,ND_MERGED
        RR = RADIUS_MERGED(L) * RADIUS_MERGED(L)
        JMAX = NP_MERGED +1 -L
        DO J=1,JMAX
          PJ=PGRID_MERGED(J)
          PJPJ=PJ*PJ
          I=(J-1)*ND_MERGED+L
          IF ( (RR-PJPJ) .GT. .0) THEN
             ZGRID_MERGED(I)=SQRT(RR-PJPJ)
          ELSE
             ZGRID_MERGED(I)=.0
          ENDIF
        ENDDO
      ENDDO
C***  End of branch for insertion of more core impact parameters 

    1 CONTINUE

C***********************************************************************
C***  Determine the number of phi angles NPHI
C***    needed to resolve line profiles with DX_WINDROT Doppler units
C***  Note: NPHI is individual for each impact parameter, depending on
C***    the maximum rotation of the shell with radius P(JP)  
C***********************************************************************
      MAXNPHI = 0
      NPHISUM = 0
      DO JP=1, NP_MERGED 
        IF (PGRID_MERGED(JP) .LE. RCOROT) THEN
          VSIDU_AT_P = VSIDU * PGRID_MERGED(JP)
        ELSE 
          VSIDU_AT_P = VSIDU * RCOROT * RCOROT / PGRID_MERGED(JP)
        ENDIF
        IF (BVSINI_AT_RCOROT) VSIDU_AT_P = VSIDU_AT_P / RCOROT 

        NPHI(JP) = 2*INT(VSIDU_AT_P/DX_WINDROT) + 1
        NPHI(JP) = MAX(2, NPHI(JP))
        IF (NPHI(JP) .GT. MAXNPHI) THEN
           MAXNPHI = NPHI(JP)
           MAXJP   = JP
        ENDIF
        NPHISUM = NPHISUM + NPHI(JP)
ccc        WRITE(0,*) JP, P(JP), NPHI(JP)
      ENDDO  
      
      WRITE(0,'(A,I4,A,I4,A,F10.3)') 
     > 'Maximum number of phi angles:', MAXNPHI, 
     >     ' at P(', MAXJP, ') =', PGRID_MERGED(MAXJP)

      NPHI_MEAN = NINT(FLOAT(NPHISUM)/FLOAT(NP_MERGED))
      WRITE(0,'(A,I4,/)') 
     > 'Average number of phi angles per impact parameter:', NPHI_MEAN 

      DO JP=1, NP_MERGED  

C***     Define angles: equally spaced in cosine
         DELTACOSPHI = 2. / (NPHI(JP)-1)
         DO LPHI = 1, NPHI(JP)
           COSPHI = 1. - DELTACOSPHI*(LPHI-1) 
           COSPHI = MAX(COSPHI, -1.)
           PHIARR(LPHI,JP) = ACOS(COSPHI)
         ENDDO

C***    Calculate the angle integration weights by trapezoidal rule  
C***      but omitting the factors 0.5 at all weights
         PHIWEIGHT(1,      JP) = PHIARR(2,JP) - PHIARR(1,JP)
         DO LPHI=2, NPHI(JP)-1
           PHIWEIGHT(LPHI,JP) = PHIARR(LPHI+1,JP) - PHIARR(LPHI-1,JP)
         ENDDO
         PHIWEIGHT(NPHI(JP),JP) = PHIARR(NPHI(JP),JP)
     >                           - PHIARR(NPHI(JP)-1,JP)

C***     Normalization: since phi = 0 ... pi, the sum of 
C***       PHIWEIGHT is normalized to 2*pi
         DO LPHI=1, NPHI(JP)
            PHIWEIGHT(LPHI,JP) =  PHIWEIGHT(LPHI,JP) / (2.*PI) 
         ENDDO

      ENDDO

      RETURN 
      
   99 WRITE (0,'(A)') '*** ERROR when decoding parameter ' // ACTPAR
      WRITE (0,'(A)') '*** The error occured in the following line:'
      WRITE (0,'(A)') KARTE(:IDX(KARTE))
      STOP '*** Fatal Error in rotation_prep'
   
      END
      
