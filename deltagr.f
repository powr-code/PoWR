      FUNCTION DELTAGR(R)
C***  DIFFERENCE OF VELOCITY GRADIENTS FROM BOTH VELOCITY LAWS: 
C***  THE BETA LAW (OUTER REGION) AND THE EXPONENTIAL LAW (INNER REGION)
C***  THIS ROUTINE SERVES FOR ESTIMATING THE CONNECTION POINT BETWEEN BOTH
C***  REGIONS AND IS CALLED FROM SUBROUTINE INITVEL, MAIN PROGRAM WRSTART
C***  BRANCH FOR BETA .LE. 0 (SQRT-LOG-LAW) REMOVED (wrh 14-Aug-2009) 

      COMMON/VELPAR/ VFINAL,VMIN,BETA,VPAR1,VPAR2,RCON,HSCALE, 
     >       BETA2, BETA2FRACTION, VPAR1_2, VPAR2_2

C***  New parametrization of beta law, wrh  5-Mar-2001 
C***  Two-beta-law: no change, since no contribution from second term
      RP = VPAR2 + R

      GRADOUT = VPAR1*BETA*(1.-1./RP)**(BETA-1.) /(RP*RP)      

C***  Prevent overflow of exp
      ARG = AMIN1 (100., (R-1.)/HSCALE)       
      GRADIN = VMIN * EXP(ARG) / HSCALE

      DELTAGR= GRADOUT - GRADIN 

      RETURN
      END
