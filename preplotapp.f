      SUBROUTINE PREPLOTAPP(ND, NF, L, 
     >    XJL, XJLAPP, XJL_PLOTDATA, XRED, XBLUE, 
     >    XJC, XJCAPP, WCHARM, 
     >    XJC_PLOTDATA_I, XJC_PLOTDATA_L, 
     >    IPLOT_XJLAPP, IPLOT_XJCAPP, LPLOT_XJCAPP, 
     >    MODE)
C********************************************************************
C***  Prepares Plot of FS- and Approximate Radiation Fields
C***    also the cores of lines and the WCHARM of the Continua are 
C***    plotted. 
C***  
C***  The Line Plot shows a specific Line over depth
C***  For the Continua two plots are prepared
C***    Over depth for one specific frequency and
C***    over frequency at one specific depth
C***
C***  For Mode=0: XJLAPP and XJCAPP is stored before they are
C***              recalculated in SETXJFINE
C***  For Mode=1: The recalculated quantities are stored
C********************************************************************

      DIMENSION XJL(ND,2), XJLAPP(2)
      DIMENSION XJL_PLOTDATA(ND,5), XRED(2), XBLUE(2)
      DIMENSION XJC(ND,NF), XJCAPP(NF), WCHARM(ND,NF)
      DIMENSION XJC_PLOTDATA_I(ND,4), XJC_PLOTDATA_L(NF,4) 

      IF (IPLOT_XJLAPP .GT. 0) THEN
        IF (MODE .EQ. 0) THEN
          XJL_PLOTDATA(L,1) = XJL(L,IPLOT_XJLAPP)
          XJL_PLOTDATA(L,2) = XJLAPP(IPLOT_XJLAPP)
          XJL_PLOTDATA(L,4) = XRED(IPLOT_XJLAPP)
          XJL_PLOTDATA(L,5) = XBLUE(IPLOT_XJLAPP)
        ELSE
          XJL_PLOTDATA(L,3) = XJLAPP(IPLOT_XJLAPP)
        ENDIF
      ENDIF

      IF (IPLOT_XJCAPP .GT. 0) THEN
        IF (MODE .EQ. 0) THEN
          XJC_PLOTDATA_I(L,1) = XJC(L,IPLOT_XJCAPP)
          XJC_PLOTDATA_I(L,2) = XJCAPP(IPLOT_XJCAPP)
          XJC_PLOTDATA_I(L,3) = XJCAPP(IPLOT_XJCAPP)
          XJC_PLOTDATA_I(L,4) = WCHARM(L,IPLOT_XJCAPP)
        ELSE
          XJC_PLOTDATA_I(L,3) = XJCAPP(IPLOT_XJCAPP)
        ENDIF
      ENDIF

      IF (L .EQ. LPLOT_XJCAPP) THEN
        IF (MODE .EQ. 0) THEN
          DO K=1, NF
            XJC_PLOTDATA_L(K,1) = XJC(L,K)
            XJC_PLOTDATA_L(K,2) = XJCAPP(K)
            XJC_PLOTDATA_L(K,4) = WCHARM(L,K)
          ENDDO
        ELSE
          DO K=1, NF
            XJC_PLOTDATA_L(K,3) = XJCAPP(K)
          ENDDO
        ENDIF
      ENDIF

      RETURN
      END
