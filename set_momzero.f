      SUBROUTINE SET_MOMZERO(ND, XJL, XHL, XKL, XNL, 
     >                           XHO, XHI, XNO, XNI, 
     >                           XHOM, XHOP, XNOM, XNOP)
C****************************************************************
C***  The Moments, integrated by SHORTRAY are set to Zero
C***    Called by COLI
C****************************************************************

      DIMENSION XJL(ND), XHL(ND), XKL(ND), XNL(ND)

      XJL = 0.
      XHL = 0.
      XKL = 0.
      XNL = 0.
      XHO = 0.
      XHI = 0.
      XNO = 0.
      XNI = 0.
C***  For special treatment of outer boundary
      XHOM = 0.
      XHOP = 0.
      XNOM = 0.
      XNOP = 0.

      RETURN
      END
