      SUBROUTINE XRUDI (XJ,WAVENUM,XJC,XLAMBDA,ND,NF,L)
C***********************************************************************
C***  INTERPOLATION OF THE CONTINUUM RADIATION FIELD AT WAVENUM
C***  LINEAR INTERPOLATION OF THE RADIATION TEMPERATURE
C***********************************************************************
      DIMENSION XJC(ND,NF),XLAMBDA(NF)

      WLENG=1.E8/WAVENUM
      NA=1
      A=XLAMBDA(1)
      NB=NF
      B=XLAMBDA(NF)
      IF ((WLENG-A)*(WLENG-B) .GT. .0) THEN
         WRITE (0,*) 'RUDIMENTAL TRANSITION OUTSIDE WAVELENGTH GRID',
     >               ' AT ', WLENG, ' ANGSTROEM'
         WRITE (0,*) '*** FATAL ERROR in SUBROUTINE XRUDI'
         STOP '*** FATAL ERROR in SUBROUTINE XRUDI'
         ENDIF
   10 IF ( NB-NA .EQ. 1) GOTO 12
      NH=(NA+NB)/2
      H=XLAMBDA(NH)
      IF ((WLENG-A)*(WLENG-H) .GT. .0) GOTO 13
      NB=NH
      B=H
      GOTO 10
   13 NA=NH
      A=H
      GOTO 10
   12 P=(WLENG-A)/(B-A)
C***  LINEAR INTERPOLATION OF THE RADIATION TEMPERATURE
      TRADA=TRADFUN(XLAMBDA(NA),XJC(L,NA))
      TRADB=TRADFUN(XLAMBDA(NB),XJC(L,NB))
      TRAD=P*TRADB+(1.-P)*TRADA
      XJ=BNUE(WLENG,TRAD)
      RETURN
      END
