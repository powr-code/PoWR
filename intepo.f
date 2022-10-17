      SUBROUTINE INTEPO (ND,N,RNE,NCHARG,POPNUM,ENTOT,NATOM,ABXYZ,
     $                   NFIRST,NLAST,IFRO,ITO,MODE,NOPOP)
C***********************************************************************
C***  LINEAR INTERPOLATION BETWEEN GIVEN DENSITY POINTS OR 
C***  EXTRAPOLATION FROM A GIVEN DENSITY POINT TO INNER OR OUTER BOUNDARY
C***
C***  THE ELECTRON DENSITY IS UPDATED ACCORDING TO THE NEW POPNUMBERS
C***********************************************************************
 
      DIMENSION RNE(ND),NCHARG(N),POPNUM(ND,N),ENTOT(ND)
      DIMENSION ABXYZ(NATOM),NFIRST(NATOM),NLAST(NATOM)

      CHARACTER*7 MODE
      LOGICAL NOPOP

      IF (IFRO.LT.1.OR.IFRO.GT.ND) STOP 'IFRO'
      IF (ITO .LT.1.OR.ITO .GT.ND) STOP 'ITO'
      IF ((MODE .EQ. 'INTERPO').AND.(IFRO.GT.ITO)) THEN
         IHELP=IFRO
         IFRO=ITO
         ITO=IHELP
      ENDIF

      IF (NOPOP) RETURN
  
      IF (IFRO .NE. ITO) THEN
        ADENFR=LOG(ENTOT(IFRO))
        ADENTO=LOG(ENTOT(ITO))
        ADIFF=ADENFR-ADENTO
        IF (MODE .EQ. 'EXTR IN') ADIFF = -ADIFF
      ENDIF 

C***  LOOP OVER ALL DEPTH POINTS  --------------------------------------
      DO 1 L=1,ND
      RNE(L)=.0
      IF (IFRO .NE. ITO) THEN
        IF (MODE .EQ. 'EXTR IN') THEN
          ADEN=(ADENTO-LOG(ENTOT(L)))/ADIFF
        ELSE
          ADEN=(ADENFR-LOG(ENTOT(L)))/ADIFF
        ENDIF
      ENDIF

C***  LOOP FOR EACH ELEMENT  ------------------------------
      DO 1 NA=1,NATOM
      NFIRNA=NFIRST(NA)
      NLANA=NLAST(NA)
      POPSUM=.0
      DO 2 J=NFIRNA,NLANA

      IF (MODE .EQ. 'INTERPO') THEN
        IF (L.LE.IFRO.OR.L.GE.ITO) GOTO 4
        if (POPNUM(IFRO,J) .gt. 0.) then
        APOPJ=LOG(POPNUM(IFRO,J))
        else
        apopj = -1000
        endif
        if (POPNUM(ITO,J) .gt. 0.) then
        APOPD=APOPJ-LOG(POPNUM(ITO,J))
        else
        apopd = -1000
        endif
        POPNUM(L,J)=EXP(APOPJ-APOPD*ADEN)
      ENDIF

      IF (MODE .EQ. 'EXTR IN') THEN
        IF (L .LE. IFRO) GOTO 4
        IF (IFRO .EQ. ITO) THEN
          POPNUM(L,J) = POPNUM(IFRO,J)
        ELSE
          APOPJ=LOG(POPNUM(ITO,J))
          APOPD=APOPJ-LOG(POPNUM(IFRO,J))
C***      Extrapolated popnumber might not exceed the element abundance
          POPNUM(L,J)=AMIN1(EXP(APOPJ-APOPD*ADEN), ABXYZ(NA))
        ENDIF
      ENDIF

      IF (MODE .EQ. 'EXTROUT') THEN
        IF (L .GE. IFRO) GOTO 4
        IF (IFRO .EQ. ITO) THEN 
           POPNUM(L,J) = POPNUM(IFRO,J)
        ELSE
           APOPJ=LOG(POPNUM(IFRO,J))
           APOPD=APOPJ-LOG(POPNUM(ITO,J))
C***       Extrapolated popnumber might not exceed the element abundance
           POPNUM(L,J)=AMIN1(EXP(APOPJ-APOPD*ADEN), ABXYZ(NA))
c          ENDIF
        ENDIF
      ENDIF

    4 POPSUM=POPSUM+POPNUM(L,J)

    2   CONTINUE
 
C***  POPNUMBERS ARE SCALED TO ENSURE NUMBER CONSERVATION
C***  ELECTRON DENSITY IS UPDATED
      POPSUM=POPSUM/ABXYZ(NA)
      DO 3 J=NFIRNA,NLANA
      POPNUM(L,J)=POPNUM(L,J)/POPSUM
      RNE(L)=RNE(L)+NCHARG(J)*POPNUM(L,J)
    3 CONTINUE
 
    1 CONTINUE
C***  ENDLOOPS  --------------------------------------------------------
 
      RETURN
      END
