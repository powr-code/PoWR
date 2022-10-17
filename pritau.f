      SUBROUTINE PRITAU (LSTAU,MODHEAD,JOBNUM,ND,RADIUS,RNE,
     >                   ENTOT,T,TAUTHOM,TAUROSS,TAUROSSCONT,
     >                   DENSCON,FILLFAC)
C***********************************************************************
C***  PRINTOUT OF THE NLTE OPTICAL DEPTH SCALES (ROSSELAND, THOMSON)
C***
C***  LSTAU = n if CARDS option PRINT TAU n is set (default = 1)
C***   only every nth depth point is printed (1st and last are always written)
C***   if no PRINT TAU is found in the cards file, steal sets  LSTAU = 1
C***   when a model is converged and the final printout is forced
C***********************************************************************
 
      IMPLICIT NONE

      INTEGER, INTENT(IN) :: ND, LSTAU, JOBNUM
      REAL, DIMENSION(ND), INTENT(IN) :: RADIUS, RNE, ENTOT ,T ,
     >                                   TAUTHOM, TAUROSS,
     >                                   TAUROSSCONT,
C***  tiefenabh. clumping nach goetz
     >                                   DENSCON, FILLFAC

      CHARACTER(100), INTENT(IN) :: MODHEAD
 
      REAL :: RL1, DENS, ENELOG, R13, R23, T13, T23, R1, T1, 
     >        TAU13, TAU23, TAU1, RLOG
      INTEGER :: L

      PRINT 1,MODHEAD,JOBNUM
    1 FORMAT (1X,  A  ,20X,'JOB NO.',I7,
     $       ///,20X,'O P T I C A L   D E P T H   S C A L E S',
     $       /,20X,39('='),//,
     $       1X,'DEPTH',7X,'R-1',6X,'LOG(R-1)',3X,'EL. TEMPERATURE',2X,
     $       'LOG(PARTICLE DENS.)',2X,'LOG(EL. DENSITY)',2X,
     $       'TAU-THOMSON',2X,'TAU-ROSSELAND',4X,'TAU-ROSSCONT',/,
     $       1X,'INDEX',31X,'(KELVIN)',7X,'(ATOMS PER CM+3)',4X,
     $       '(EL. PER CM+3)',4X,'(NLTE)',7X,'(NLTE)',/)
 
C***  LOOP OVER ALL DEPTH POINTS  --------------------------------------
      DO L=1, ND
        IF ((((L-1)/LSTAU)*LSTAU /= (L-1)) .AND. (L /= ND)) CYCLE
        RL1=RADIUS(L)-1.
        IF (L < ND) THEN
          RLOG=ALOG10(RL1)
        ELSE
          RLOG=-9.99E+99
        ENDIF 
        DENS=ALOG10(ENTOT(L))             !particle density
        ENELOG=ALOG10(RNE(L)*ENTOT(L))    !electron density
        PRINT 3, L, RL1, RLOG, T(L), DENS, ENELOG, 
     >           TAUTHOM(L), TAUROSS(L), TAUROSSCONT(L)
    3   FORMAT (2X,I3,4X,G10.3,2X,G10.3,5X,F9.0,10X,F7.3,13X,F8.3,8X,
     >          F7.3,2X,F15.10,1X,F15.10)
      ENDDO
C***  ENDLOOP  ---------------------------------------------------------
 
c      IF (DENSCON_FIX .NE. 1.) THEN
c         FILLFAC_FIX=1./DENSCON_FIX
c         PRINT 6, DENSCON_FIX, FILLFAC_FIX
c 6       FORMAT (/, ' MODEL WITH CLUMPING: DENSCON_FIX =', F6.2, 
c     >        '    FILLFAC_FIX =', F6.4)  
c      ENDIF

      PRINT 17
   17 FORMAT (/////,1X,'RADII AND TEMPERATURES FOR DIFFERENT OPTICAL',
     $       ' DEPTHS:      RADIUS',5X,'TEMPERATURE',/,58X,24('-'))
C***  CALCULATE RADII AND CORRESPONDING TEMPERATURES FOR TAU-THOMSON = 1/3,
C***                                                                 = 2/3,
C***                                                                 = 1.0
      TAU13=0.333333333333
      IF (TAUTHOM(ND) < TAU13) THEN
         R13=1.
         T13=T(ND)
      ELSE
         CALL LIPO (R13,TAU13,RADIUS,TAUTHOM,ND)
         CALL LIPO (T13,R13,T,RADIUS,ND)
      ENDIF
      TAU23=0.666666666666
      IF (TAUTHOM(ND) < TAU23) THEN
         R23=1.
         T23=T(ND)
      ELSE
         CALL LIPO (R23,TAU23,RADIUS,TAUTHOM,ND)
         CALL LIPO (T23,R23,T,RADIUS,ND)
      ENDIF
      TAU1=1.
      IF (TAUTHOM(ND) < TAU1 ) THEN
         R1 =1.
         T1 =T(ND)
      ELSE
         CALL LIPO (R1 ,TAU1 ,RADIUS,TAUTHOM,ND)
         CALL LIPO (T1 ,R1 ,T,RADIUS,ND)
      ENDIF
      PRINT 18, R13,T13,R23,T23,R1,T1
   18 FORMAT (35X,'TAU-THOMSON = 1/3:',5X,F7.3,6X,F9.0,/,
     >        35X,'TAU-THOMSON = 2/3:',5X,F7.3,6X,F9.0,/,
     >        35X,'TAU-THOMSON = 1.0:',5X,F7.3,6X,F9.0,/)
C***  CALCULATE RADII AND CORRESPONDING TEMPERATURES FOR TAU-ROSSELAND = 1/3,
C***                                                                   = 2/3,
C***                                                                   = 1.0
      IF (TAUROSS(ND) < TAU13) THEN
         R13=1.
         T13=T(ND)
      ELSE
         CALL LIPO (R13,TAU13,RADIUS,TAUROSS,ND)
         CALL LIPO (T13,R13,T,RADIUS,ND)
      ENDIF
      IF (TAUROSS(ND) < TAU23) THEN
         R23=1.
         T23=T(ND)
      ELSE
         CALL LIPO (R23,TAU23,RADIUS,TAUROSS,ND)
         CALL LIPO (T23,R23,T,RADIUS,ND)
      ENDIF
      IF (TAUROSS(ND) < TAU1 ) THEN
         R1 =1.
         T1 =T(ND)
      ELSE
         CALL LIPO (R1 ,TAU1 ,RADIUS,TAUROSS,ND)
         CALL LIPO (T1 ,R1 ,T,RADIUS,ND)
      ENDIF

      PRINT 19, R13,T13,R23,T23,R1,T1
   19 FORMAT (33X,'TAU-ROSSELAND = 1/3:',5X,F7.3,6X,F9.0,/,
     >        33X,'TAU-ROSSELAND = 2/3:',5X,F7.3,6X,F9.0,/,
     >        33X,'TAU-ROSSELAND = 1.0:',5X,F7.3,6X,F9.0)
 

      RETURN
      END
