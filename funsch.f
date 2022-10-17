      FUNCTION FUNSCH (X)
C***********************************************************************
C***  THIS FUNCTION IS USED TO ESTIMATE THE SCHARMER LINE CORES IN THE
C***  BOUNDARY ZONES
C***  --  USED AS FORMAL PARAMETER FOR "REGULA"
C***  REGULA IS CALLED BY SUBROUTINE "COFREQ"
C***   COMMON / COMFUN / IS ALSO DEFINED IN CORFREQ
C***********************************************************************

      COMMON / COMFUN / DELTAV,XMIN
      EXTERNAL ERF

      XA=AMAX1(X-DELTAV,XMIN)
      FUNSCH=ERF(X)-ERF(XA)

      RETURN
      END
