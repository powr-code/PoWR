      SUBROUTINE CONNECT_VELOTABLE (T, RADIUS, VELO, ND, RSTAR, RMAX, 
     >       GEFFL, XMU, VTURB, IMIN, TEFF, XMASS, CALLED_FROM)
C******************************************************************************
C***  CONNECTIOM OF THE INNER VELOCITY-FIELD WITH TABULATED WIND PART 
C***  CALLED FROM: INITVEL, VELTHIN if bVELOTABLE = .TRUE.
C***  written by wrh 11-Dec-2024
C***  See velthin.f for meaning of the variables
C***
C***  As in the beta-law case, this subroutine has the task to determine
C***  the radius RCON where the photosphere and wind are connected. 
C***
C***  This subroutine can be called from INITVEL or from VELTHIN, and
C***    internally adapts to these two cases:
C***  If called from INITVEL, 
C***    the radius grid is not yet established. 
C***    Hence, all depth-indexed variables are DUMMY in the 
C***    CALL CONNECT_VELOTABLE, and cannot be used inside. 
C***  If called from VELTHIN, in contrast, the radius grid has been 
C***     established already, and the hydrostatic equation has been 
C***     integrated from the inner boundary (L=ND) outward till L=IMIN.
C***     IMIN = smallest depth index for which the hydrostatic velocity is
C***         defined (avoiding overflow due to exponential growth of VELO)
C***  The condition for the connection point is that both, velocity and
C***  velocity-gradient, are continuous. 
C***
C***  The connection point RCON is a virtual point that must not coincide 
C***  with a depth point of the radius grid. 
C***
C***  The connection point is placed such that the velocity there is 
C***  VSOUND * FSONICTA, 
C***  where the latter factor is defined by the user on the VELOTABLE  
C***  line in the CARDS input file. However, FSONICTA is automatically 
C***  reduced if necessary such that the gradient at the connection 
C***  point is not steeper than half of the maximum gradient that is
C***  present in the tabulated wind velocity field. 
C***
C***  The second part of this subroutine is to re-scale the tabulated
C***  velocity field such that the connection-point conditions are met. 
C***  This requires an internal iteration. 
C******************************************************************************
      USE COMVELOTABLE
  
      IMPLICIT NONE
      INCLUDE 'interfacebib.inc'

      INTEGER, INTENT(IN) :: ND
      REAL, INTENT(IN) :: RSTAR, RMAX, VTURB, TEFF, XMASS
      REAL, DIMENSION(ND), INTENT(IN) :: T, RADIUS, XMU, GEFFL
      REAL, DIMENSION(ND), INTENT(INOUT) :: VELO
      CHARACTER*(*) CALLED_FROM

      INTEGER ISRCHFLT, NTABLE_WIND

      INTEGER :: IMIN, ICON, L,
     >           LGRADMAX_TABLE,LGRAD_HALF, ICON_SAFE, I1, I2,ITER
      REAL :: VSOUND, GRADICON, RTABLE_AT_VCON, SCALE_R, SCALE_V,
     >        GRADITABLE_MAX, GRADIL, RTABLE_CON, VTABLE_CON,
     >        VTABLE_RMAX, HSCALE_ICON, P, Q, GRADICON_MINUS,
     >        HSCALE_ICON_MINUS,QVELO, PVELO, SCALE, GRADICON_ICON,
     >        VCON
      REAL, EXTERNAL :: WRVEL

      !Common block VELPAR needed to transfer velocity parameters
      ! This routine updates RCON, VCON
      REAL :: VFINAL, VMIN, BETA, VPAR1, VPAR2, RCON, HSCALE,
     >                 BETA2, BETA2FRACTION, VPAR1_2, VPAR2_2
      COMMON /VELPAR/ VFINAL,VMIN,BETA,VPAR1,VPAR2,RCON,HSCALE,
     >                 BETA2, BETA2FRACTION, VPAR1_2, VPAR2_2

      !Physical Constants
      REAL, PARAMETER :: BOLTZ = 1.38E-16   !BOLTZMANN CONSTANT (ERG/DEG)
      REAL, PARAMETER :: AMU = 1.66E-24     !ATOMIC MASS UNIT (GRAMM)
      REAL, PARAMETER :: RGAS = BOLTZ / AMU !Gas Constant (CGS) = BOLTZ / AMU

C***  Always start from CARDS value:
      FSONICTA = FSONICTA_CARDS

C***  Iteration loop: 
C***    re-scaling of the velotable in RTABLE and VTABLE will also
C***    modify the velocity gradient; therefore, the whole allgorithm 
C***    is iteratively repeated
      DO 100 ITER = 1, 20

C***  Check for valid parameter CALLED_FROM
      IF (CALLED_FROM .NE. 'INITVEL' .AND. 
     >    CALLED_FROM .NE. 'VELTHIN') GOTO 91

C**********************************************************
C***  Find maximum gradient in the tabulated velocity field
C**********************************************************
      GRADITABLE_MAX = .0
      DO L=NTABLE, 1, -1
       IF (L .EQ. 1) THEN
         GRADIL = (VTABLE(2)-VTABLE(1))/(RTABLE(2)-RTABLE(1))
       ELSEIF (L .EQ. NTABLE) THEN
         GRADIL = (VTABLE(L)-VTABLE(L-1))/(RTABLE(L)-RTABLE(L-1))
       ELSE
         GRADIL = (VTABLE(L+1)-VTABLE(L-1))/(RTABLE(L+1)-RTABLE(L-1))
       ENDIF
C*     Find maximum and its index
       IF (GRADIL .GT. GRADITABLE_MAX) THEN
          GRADITABLE_MAX = GRADIL
       ENDIF 
      ENDDO

C**********************************************************
C***  Establish RCON and VCON
C**********************************************************

C***  Branching
      IF (CALLED_FROM .EQ. 'VELTHIN') GOTO 50

C************************************************************
C***  INITVEL branch
C***    in this stage the radius grid is not yet established
C***********************************************************
      VSOUND = SQRT(RGAS * TEFF / XMASS) * 1.E-5 
      VCON = FSONICTA*VSOUND
      GRADICON = VCON / HSCALE
      IF (VCON .LT. VMIN) THEN
         WRITE (0,*) CALLED_FROM // 
     >               '/CONNECT_VELOTABLE: NO HYDROSTATIC DOMAIN'
         VCON = VMIN
         RCON = 1.0
         GRADICON = .0
         GOTO 60
      ENDIF

C***  Problem case: Maximum gradient in TABLE lower then requested       
C***  --> reduce VCON such that the gradient is small enough 
      IF (GRADICON .GE. GRADITABLE_MAX/2) THEN
         VCON = GRADITABLE_MAX/2 * HSCALE
         GRADICON = VCON / HSCALE
         FSONICTA = MIN (FSONICTA, VCON/VSOUND)
      ENDIF

      RCON = HSCALE * LOG(VCON/VMIN) + 1.
      RCON = MAX(RCON, 1.0)

C***  End of INITVEL branch
      GOTO 60

C***************************************************************
C*** VELTHIN branch
C***  if called from VELTHIN, radius grid is already established
C***************************************************************
   50 CONTINUE

C***  find last subsonic index, ICON
      DO L=ND, IMIN, -1
C***  We neglect here that VSOUND should be interpolated 
         VSOUND = SQRT(RGAS * T(L) / XMU(L)) * 1.E-5   
         ICON = MIN(L+1,ND)
         IF (VELO(L) .GT. FSONICTA*VSOUND) EXIT
      ENDDO

      VCON = FSONICTA*VSOUND

C***  No hydrostatic domain
      IF (ICON .EQ. ND) THEN
         WRITE (0,*) CALLED_FROM // 
     >               '/CONNECT_VELOTABLE: NO HYDROSTATIC DOMAIN'
         VCON = VMIN
         RCON = 1.0
         GOTO 60
      ENDIF

C***  Define V-gradient at connection point
      HSCALE_ICON = (BOLTZ*T(ICON)/(XMU(ICON)*AMU) + (VTURB*1.E5)**2)
     >               /GEFFL(ICON) /RSTAR
      GRADICON = VELO(ICON) / HSCALE_ICON 
C***  Interpolation of the velocity gradient over V if possible
      IF (ICON .GT. IMIN) THEN
         HSCALE_ICON_MINUS = (BOLTZ*T(ICON-1)/(XMU(ICON-1)*AMU) 
     >        + (VTURB*1.E5)**2) /GEFFL(ICON-1) /RSTAR
         GRADICON_MINUS = VELO(ICON-1) / HSCALE_ICON_MINUS
         P = (FSONICTA*VSOUND-VELO(ICON))/(VELO(ICON-1)-VELO(ICON))
         Q = 1-P
         GRADICON = Q*GRADICON + P*GRADICON_MINUS
      ENDIF

C***  Problem case: Maximum gradient in TABLE lower then requested       
C***  --> reduce VCON such that the gradient is small enough 
      IF (GRADICON .GE. GRADITABLE_MAX/2) THEN

C***     Go inward with ICON till gradient is below GRADITABLE_MAX/2
         ICON_SAFE = ND
         DO L=ICON, ND
            HSCALE_ICON = (BOLTZ*T(L)/(XMU(L)*AMU) + (VTURB*1.E5)**2)
     >                   /GEFFL(L) /RSTAR
            GRADICON_ICON = VELO(L) / HSCALE_ICON 
            ICON_SAFE = L
            IF (GRADICON_ICON .LT. GRADITABLE_MAX/2) EXIT
         ENDDO
         ICON = ICON_SAFE

C***     Interpolation of the velocity gradient over V if possible
         IF (ICON .GT. IMIN) THEN
            HSCALE_ICON_MINUS = (BOLTZ*T(ICON-1)/(XMU(ICON-1)*AMU) 
     >                   + (VTURB*1.E5)**2) /GEFFL(ICON-1) /RSTAR
            GRADICON_MINUS = VELO(ICON-1) / HSCALE_ICON_MINUS
            P = (FSONICTA*VSOUND-VELO(ICON))/(VELO(ICON-1)-VELO(ICON))
            Q = 1. - P
            GRADICON = Q*GRADICON_ICON + P*GRADICON_MINUS
            PVELO = (GRADITABLE_MAX/2-GRADICON_ICON) / 
     >              (GRADICON_MINUS-GRADICON_ICON)
            QVELO= 1 - PVELO
            VCON = QVELO*VELO(ICON)+PVELO*VELO(ICON-1)
C***        Make sure that no extrapolation happened 
            VCON = MAX (VCON, VELO(ICON))
            VCON = MIN (VCON, VELO(ICON-1))
         ELSE
            GRADICON = GRADICON_ICON
            VCON = VELO(ICON)
         ENDIF

         FSONICTA = MIN (FSONICTA, VCON/VSOUND)
      ENDIF
C***  End of problem case

      VCON = FSONICTA*VSOUND

      CALL SPLINPOX (RCON, VCON, RADIUS(IMIN), VELO(IMIN),ND-IMIN+1)


C***  end of VELTHIN branch

C***  from here, both branches are united
   60 CONTINUE
      IF (FSONICTA .NE. FSONICTA_CARDS)
     >   WRITE (0,'(A,F7.4)') CALLED_FROM // '/CONNECT_VELOTABLE: '
     >         // 'VCON reduced to VSOUND *', FSONICTA  

C**************************************************************
C***  The tabulated velocity field is how re-scaled in order to 
C***  meet the connection point with velocity and gradient 
C**************************************************************

C***  Branch if NO HYDROSTATIC DOMAIN 
C***  Do not care about the gradient when no hydrostatic domain

      IF (RCON .LE. 1.0) THEN
         VCON = VMIN

C***     shorten table to first radius where VTABLE(L) .LT. VCON
         NTABLE_WIND = ISRCHFLT(NTABLE,VTABLE, 1, VCON) - 1
         RTABLE_CON = RTABLE(NTABLE_WIND)

C***     Rescale RTABLE such that RTABLE_CON --> RCON = 1.0) 
         SCALE = (RCON-RTABLE(1)) / (RTABLE_CON - RTABLE(1))
         DO L=1, NTABLE
            RTABLE(L) = RTABLE(1) + SCALE*(RTABLE(L)-RTABLE(1))
         ENDDO
         GOTO 61
      ENDIF

C***  Branch if HYDROSTATIC DOMAIN exists:

C***  Find radius where the tabulated velo has same gradient
C***   note that RTABLE, VTABLE are in descending order!
      DO L=NTABLE, 1, -1
       IF (L .EQ. 1) THEN
         GRADIL = (VTABLE(2)-VTABLE(1))/(RTABLE(2)-RTABLE(1))
       ELSEIF (L .EQ. NTABLE) THEN
         GRADIL = (VTABLE(L)-VTABLE(L-1))/(RTABLE(L)-RTABLE(L-1))
       ELSE
         GRADIL = (VTABLE(L+1)-VTABLE(L-1))/(RTABLE(L+1)-RTABLE(L-1))
       ENDIF
C*     Find first table entry where gradient matches
       IF (GRADIL .GT. GRADICON) THEN
          RTABLE_CON = RTABLE(L)
          GOTO 21
       ENDIF
      ENDDO
   21 CONTINUE

C***  re-scale the radius points such that RTABLE_CON becomes RCON
C***  but RTABLE(1) stays unchanged
   61 CONTINUE
      SCALE_R = (RCON-RTABLE(1)) / (RTABLE_CON - RTABLE(1))
      DO L=1, NTABLE
         RTABLE(L) = RTABLE(1) + SCALE_R*(RTABLE(L)-RTABLE(1))
      ENDDO
      IF (RCON .LT. RTABLE(NTABLE)) RTABLE(NTABLE) = RCON

C*** Now re-scale the velocities such that VTABLE_CON becomes VCON
C*** but VTABLE(RMAX) stays VFINAL
ccc   61 CONTINUE
      IF (RCON .GT. RTABLE(NTABLE)) THEN
         CALL SPLINPO (VTABLE_CON, RCON, VTABLE, RTABLE, NTABLE)
      ELSE
        VTABLE_CON = VTABLE(NTABLE)
      ENDIF
      CALL SPLINPO (VTABLE_RMAX, RMAX, VTABLE, RTABLE, NTABLE)
      SCALE_V = (VCON-VTABLE_RMAX) / (VTABLE_CON - VTABLE_RMAX)
      DO L=1, NTABLE
         VTABLE(L) = VTABLE_RMAX + SCALE_V*(VTABLE(L)-VTABLE_RMAX)
      ENDDO

      IF ((ABS(SCALE_R-1) .LT. 1e-6) .AND. 
     >    (ABS(SCALE_V-1) .LT. 1e-6)) GOTO 200
  100 CONTINUE
C***  End of iteration loop ****************************************

      WRITE (0,*) 
     >   '*** WARNING from CONNECT_VELOTABLE: LOOP NOT CONVERGED'

  200 CONTINUE

C***  Re-establish wind velocity (only VELTHIN branch)
      IF (CALLED_FROM .EQ. 'VELTHIN') THEN
         DO L=1, ICON-1
            VELO(L) = WRVEL(RADIUS(L))
         ENDDO  
      ENDIF

ccc   test: check if monotonic
      DO L=2, ND
       if (velo(l-1) .le. velo(l)) then
         write (0,'(A,I4,X,F8.3,X,1PG14.6)') 
     >         'L, VELO(L), R(L) =', l-1, velo(l-1), radius(l-1)
         write (0,'(A, I4,X,F8.3,X,1PG14.6)') 
     >         'L, VELO(L), R(L) =', l, velo(l), radius(l)
         write (0,*) 'ICON=', ICON
         write (0,*) 'RCON=', RCON
         write (0,*) 'VCON=', VCON
         stop 'ERROR: velo not monotonic'
       endif
      enddo

      RETURN

C******************* ERROR BRANCH ***********************
   91 WRITE (0,*) '*** INVALID CALLED_FROM parameter:' // CALLED_FROM
      STOP '*** FATAL INTERNAL ERROR detected in CONNECT_VELOTABLE'

      END


