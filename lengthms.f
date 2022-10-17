      SUBROUTINE LENGTHMS(ICHANNEL, NDIM, NAME, IERR)
C************************************************************
C***  ROUTINE VON wrh
C************************************************************

      CALL CMSSTORE (ICHANNEL, IDUMMY, IDUMMY, NAME, NDUMMY, XDUMMY, 
     >               NDIM, 'LENGTH', IERR)

      RETURN
      END
