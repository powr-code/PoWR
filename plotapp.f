      SUBROUTINE PLOTAPP(KANAL, ND, NF, 
     >    XJL_PLOTDATA, XJC_PLOTDATA_I, XJC_PLOTDATA_L, 
     >    IPLOT_XJLAPP, IPLOT_XJCAPP, LPLOT_XJCAPP, NITER_PLOT_JAPP)
C********************************************************************
C***  Plots the data prepared in PREPLOTAPP
C********************************************************************

      DIMENSION XJL_PLOTDATA(ND,5)
      DIMENSION XJC_PLOTDATA_I(ND,4), XJC_PLOTDATA_L(NF,4) 

      CALL JSYMSET ('G2','TRANSFER')

C***  First Plot (XJL vs Depth), Data Table also added
      WRITE (KANAL,'(A,I5)') 
     >  ' PLOT   : XJL vs XJLAPP for Line index ', IPLOT_XJLAPP
      WRITE (KANAL,'(2A,I5)') 
     >  ' HEADER : XJL, &2XJLAPP&Told&M&1 and &4XJLAPP&Tnew&M&1 ', 
     >  'for Line index ', IPLOT_XJLAPP
      WRITE (KANAL,'(A)') 
     >  ' X-ACHSE:\CENTER\Depth L'
      WRITE (KANAL,'(A)') 
     >  ' Y-ACHSE:\CENTER\J'
      WRITE (KANAL,'(A)') '     MASSTAB'
      WRITE (KANAL,'(A)') ' X: AUTO'
      WRITE (KANAL,'(A)') ' Y: '
C***  Write Plotdata for Depth-Plot at file and Plot first Graph
      WRITE (KANAL,'(A)') '* PLOTAPP DATA DEPTH'
      WRITE (KANAL,'(A3,3A21,3X,2A8,3A21,1X,A8)') 
     >  '* L', 'XJL', 'XJLAPP', 'XJLAPPNEW', 'XRED', 'XBLUE', 
     >  'XJC', 'XJCAPP', 'XJCAPPNEW', 'WCHARM'
      WRITE (KANAL,'(A)') 'N=? XYTABLE SELECT 1 2'
      DO L=1, ND
        WRITE (KANAL,'(I3,3(1X,E20.10),3X,2(1X,F7.4), 
     >                     3(1X,E20.10),1X,F8.6)') 
     >    L, (XJL_PLOTDATA(L,II), II=1, 5), 
     >       (XJC_PLOTDATA_I(L,II), II=1,4)
      ENDDO
      WRITE (KANAL,'(A)') 'FINISH'
      WRITE (KANAL,'(A)') ' N=? COLOR=2 XYTABLE SELECT 1 3'
      WRITE (KANAL,'(A)') 
     >  'COMMAND INCLUDE SELF INCKEY="* PLOTAPP DATA DEPTH"'
      WRITE (KANAL,'(A)') ' N=? COLOR=4 XYTABLE SELECT 1 4'
      WRITE (KANAL,'(A)') 
     >  'COMMAND INCLUDE SELF INCKEY="* PLOTAPP DATA DEPTH"'
      WRITE (KANAL,'(A)') ' END'


C***  Second Plot (XRED, XBLUE vs Depth)
      WRITE (KANAL,'(A,I5)') 
     >  ' PLOT   : XBLUE, XRED for Line index ', IPLOT_XJLAPP
      WRITE (KANAL,'(A,I5)') 
     >  ' HEADER : XBLUE, &2-XRED&1 for Line index ', IPLOT_XJLAPP
      WRITE (KANAL,'(A)') 
     >  ' X-ACHSE:\CENTER\Depth L'
      WRITE (KANAL,'(A)') 
     >  ' Y-ACHSE:\CENTER\J'
      WRITE (KANAL,'(A)') '     MASSTAB'
      WRITE (KANAL,'(A)') ' X: AUTO'
      WRITE (KANAL,'(A)') ' Y: '
      WRITE (KANAL,'(A)') ' N=? COLOR=2 XYTABLE SELECT 1 5'
      WRITE (KANAL,'(A)') 'COMMAND Y* -1.'
      WRITE (KANAL,'(A)') 
     >  'COMMAND INCLUDE SELF INCKEY="* PLOTAPP DATA DEPTH"'
      WRITE (KANAL,'(A)') ' N=? COLOR=1 XYTABLE SELECT 1 6'
      WRITE (KANAL,'(A)') 
     >  'COMMAND INCLUDE SELF INCKEY="* PLOTAPP DATA DEPTH"'
      WRITE (KANAL,'(A)') ' END'

C***  Third Plot (XJC vs Depth)
      WRITE (KANAL,'(A,I5)') 
     >  ' PLOT   : XJC vs XJCAPP for Frequency index ', IPLOT_XJCAPP
      WRITE (KANAL,'(2A,I5)') 
     >  ' HEADER : XJC, &2XJCAPP&Told&M&1 and &4XJCAPP&Tnew&M&1 ', 
     >  'for Frequency index ', IPLOT_XJCAPP
      WRITE (KANAL,'(A)') 
     >  ' X-ACHSE:\CENTER\Depth L'
      WRITE (KANAL,'(A)') 
     >  ' Y-ACHSE:\CENTER\J'
      WRITE (KANAL,'(A)') '     MASSTAB'
      WRITE (KANAL,'(A)') ' X: AUTO'
      WRITE (KANAL,'(A)') ' Y: '
      WRITE (KANAL,'(A)') ' N=? COLOR=1 XYTABLE SELECT 1 7'
      WRITE (KANAL,'(A)') 
     >  'COMMAND INCLUDE SELF INCKEY="* PLOTAPP DATA DEPTH"'
      WRITE (KANAL,'(A)') ' N=? COLOR=2 XYTABLE SELECT 1 8'
      WRITE (KANAL,'(A)') 
     >  'COMMAND INCLUDE SELF INCKEY="* PLOTAPP DATA DEPTH"'
      WRITE (KANAL,'(A)') ' N=? COLOR=4 XYTABLE SELECT 1 9'
      WRITE (KANAL,'(A)') 
     >  'COMMAND INCLUDE SELF INCKEY="* PLOTAPP DATA DEPTH"'
      WRITE (KANAL,'(A)') ' END'

C***  Fourth Plot (WCHARM vs Depth)
      WRITE (KANAL,'(A,I5)') 
     >  ' PLOT   : WCHARM for Frequency index ', IPLOT_XJCAPP
      WRITE (KANAL,'(A,I5)') 
     >  ' HEADER : WCHARM for Frequency index ', IPLOT_XJCAPP
      WRITE (KANAL,'(A)') 
     >  ' X-ACHSE:\CENTER\Depth L'
      WRITE (KANAL,'(A)') 
     >  ' Y-ACHSE:\CENTER\ '
      WRITE (KANAL,'(A)') '     MASSTAB'
      WRITE (KANAL,'(A)') ' X: AUTO'
      WRITE (KANAL,'(A)') ' Y: '
      WRITE (KANAL,'(A)') ' N=? COLOR=1 XYTABLE SELECT 1 10'
      WRITE (KANAL,'(A)') 
     >  'COMMAND INCLUDE SELF INCKEY="* PLOTAPP DATA DEPTH"'
      WRITE (KANAL,'(A)') ' END'

C***  Fifth Plot (XJC vs INDEX)
      WRITE (KANAL,'(A,I5)') 
     >  ' PLOT   : XJC vs XJCAPP for Depth ', LPLOT_XJCAPP
      WRITE (KANAL,'(2A,I5)') 
     >  ' HEADER : XJC, &2XJCAPP&Told&M&1 and &4XJCAPP&Tnew&M&1 ', 
     >  'for Depth ', LPLOT_XJCAPP
      WRITE (KANAL,'(A)') 
     >  ' X-ACHSE:\CENTER\Frequency K'
      WRITE (KANAL,'(A)') 
     >  ' Y-ACHSE:\CENTER\J'
      WRITE (KANAL,'(A)') '     MASSTAB'
      WRITE (KANAL,'(A)') ' X: AUTO'
      WRITE (KANAL,'(A)') ' Y: '
C***  Write Plotdata for Index-Plot at file and Plot first graph
      WRITE (KANAL,'(A)') '* PLOTAPP DATA INDEX'
      WRITE (KANAL,'(A3,3A21,A11)') 
     >  '* K', 'XJC', 'XJCAPP', 'XJCAPPNEW', 'WCHARM'
      WRITE (KANAL,'(A)') 'N=? XYTABLE SELECT 1 2'
      DO K=1, NF
        WRITE (KANAL,'(I3, 
     >                     2(1X,E20.10),3X,F8.6)') 
     >    K, (XJC_PLOTDATA_L(K,II), II=1,3)
      ENDDO
      WRITE (KANAL,'(A)') 'FINISH'
      WRITE (KANAL,'(A)') ' N=? COLOR=2 XYTABLE SELECT 1 3'
      WRITE (KANAL,'(A)') 
     >  'COMMAND INCLUDE SELF INCKEY="* PLOTAPP DATA INDEX"'
      WRITE (KANAL,'(A)') 'FINISH'
      WRITE (KANAL,'(A)') ' N=? COLOR=4 XYTABLE SELECT 1 4'
      WRITE (KANAL,'(A)') 
     >  'COMMAND INCLUDE SELF INCKEY="* PLOTAPP DATA INDEX"'
      WRITE (KANAL,'(A)') ' END'

C***  Fourth Plot (WCHARM vs Index)
      WRITE (KANAL,'(A,I5)') 
     >  ' PLOT   : WCHARM for Depth ', LPLOT_XJCAPP
      WRITE (KANAL,'(A,I5)') 
     >  ' HEADER : WCHARM for Depth ', LPLOT_XJCAPP
      WRITE (KANAL,'(A)') 
     >  ' X-ACHSE:\CENTER\Frequency K'
      WRITE (KANAL,'(A)') 
     >  ' Y-ACHSE:\CENTER\ '
      WRITE (KANAL,'(A)') '     MASSTAB'
      WRITE (KANAL,'(A)') ' X: AUTO'
      WRITE (KANAL,'(A)') ' Y: '
      WRITE (KANAL,'(A)') ' N=? COLOR=1 XYTABLE SELECT 1 5'
      WRITE (KANAL,'(A)') 
     >  'COMMAND INCLUDE SELF INCKEY="* PLOTAPP DATA INDEX"'
      WRITE (KANAL,'(A)') ' END'


      RETURN
      END
