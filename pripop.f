      SUBROUTINE PRIPOP (LSPOP,WEIGHT,NCHARG,NOM,TNEW,BUNLU,
     $             ND,N,RNE,ITNE,LEVEL,POPNUM,DEPART,JOBNUM,MODHEAD,
     $             ABXYZ,NFIRST,NLAST,NATOM,NATOUT, SMPOP)
C***********************************************************************
C***  OUTPUT OF RELATIVE POPULATION NUMBERS AND DEPARTURE COEFFICIENTS
C***********************************************************************
 
      DIMENSION POPNUM(ND,N),DEPART(ND,N)
      DIMENSION RNE(ND),ITNE(ND),TNEW(ND)
      DIMENSION WEIGHT(N),NCHARG(N),NOM(N)
      DIMENSION ABXYZ(NATOM),NFIRST(NATOM),NLAST(NATOM)
      LOGICAL BUNLU
      CHARACTER LEVEL(N)*10
      CHARACTER MODHEAD*100
      CHARACTER PRILINE*130,NUMBER*16

C***  SMACH1 IS THE SMALLEST MACHINE CONSTANT WITH 1.+(-SMACH1) .NE. 1.
      SMACH1=SMACH(1)

C**********************************************************************
C***  SET SMALLPOP = THRESHOLD FOR SMALL-POPNUMBER WARNING "SMA"
C***  NOTE: THIS SHOULD AGREE WITH THE WARNINGS ISSUED BY SUBR. PRICORR
C**********************************************************************
      SMALLPOP = SMPOP
 
      PRINT 1,MODHEAD,JOBNUM
    1 FORMAT (1X,  A  ,20X,'JOB NO.',I7,
     $           //,20X,'RELATIVE NON-LTE POPULATION NUMBERS (LOG) PER',
     $ ' ELEMENT AND DEPARTURE COEFFICIENTS',/,20X,80('-'))
 
C***  OUTPUT ONLY FOR ONE OR ALL ELEMENTS
      IF (NATOUT .NE. 0) THEN
         NASTART=NATOUT
         NASTOP=NATOUT
      ELSE
         NASTART=1
         NASTOP=NATOM
      ENDIF
         
C***  LOOP OVER SPECIFIED ELEMENTS  ******************************************
      DO 99 NA=NASTART,NASTOP
      ABNA=ABXYZ(NA)
      ABLOG=ALOG10(ABNA)
C!!!      RELMACH=SMACH1*ABNA
      RELMACH = SMALLPOP
      NFIRNA=NFIRST(NA)
      NLANA=NLAST(NA)
      PRINT 21, LEVEL(NFIRNA)(1:2),ABLOG
   21 FORMAT (///,1X,39('+'),'  LOG OF REL. ',A2,'-ABUNDANCE (BY ',
     $        'NUMBER): ',F7.3,2X,39('+'))
      J1=NFIRNA
    4 J2=MIN0(NLANA,J1+4)
      PRINT 2,(LEVEL(J),J=J1,J2)
    2 FORMAT (//,4X,'EL.   IT',4X,5(7X,A10,6X))
      PRINT 11
   11 FORMAT(1X)
 
C***  LOOP OVER ALL DEPTH POINTS  --------------------------------------
      DO 3 L=1,ND
      IF(((L-1)/LSPOP)*LSPOP.NE.(L-1) .AND. L.NE.ND) GOTO 3
      RNESUM=0.0
      DO 23 I=NFIRNA,NLANA
   23 RNESUM=RNESUM+NCHARG(I)*POPNUM(L,I)
      RELRNE=RNESUM/ABNA
      ENCODE (130,9,PRILINE) L,RELRNE,ITNE(L)
    9 FORMAT(I3,F7.3,I3)
 
      DO 5 J=J1,J2
      IF (POPNUM(L,J) .NE. .0) THEN
          ENCODE(16,8,NUMBER) ALOG10(ABS(POPNUM(L,J)))-ABLOG,DEPART(L,J)
    8     FORMAT (F6.2,G10.3)
      ELSE
          NUMBER='****  ZERO  ****'
      ENDIF
      I=18+(J-J1)*23
      PRILINE(I+4:I+19)=NUMBER
C***  WARNINGS: NEGATIVE OR SMALL ABS. POP.NUMBERS
      IF (POPNUM(L,J) .LE. RELMACH) PRILINE(I:I+3) = 'SMA'
      IF (POPNUM(L,J) .LT. .0     ) PRILINE(I:I+3) = 'NEG'
C***  LASER WARNING BETWEEN BOUND LEVELS
      IF (J .GT. NFIRNA) THEN
      IF ((NCHARG(J-1) .EQ. NCHARG(J)) .AND. (NOM(J-1) .EQ. NOM(J))
     $    .AND. (WEIGHT(J)*POPNUM(L,J-1) .LT. WEIGHT(J-1)*POPNUM(L,J)))
     $      PRILINE(I-2:I-2)='*'
      ENDIF
    5 CONTINUE
 
      PRINT 6,PRILINE
    6 FORMAT (A)
    3 CONTINUE
C***  ENDLOOP  ---------------------------------------------------------
 
      IF(J2.EQ.NLANA) GOTO 99
      J1=J1+5
      GOTO 4
   99 CONTINUE
C***  ENDLOOP  *********************************************************

C***  PRINTOUT OF THE TEMPERATURE STRATIFICATION
      IF (.NOT. BUNLU) THEN
      PRINT 10
   10 FORMAT (///,40X,'TEMPERATURE STRATIFICATION',/,40X,26('='),
     $        //,40X,'DEPTH INDEX      T (KELVIN)',/)
      DO 12 L=1,ND
      IF(((L-1)/LSPOP)*LSPOP.NE.(L-1) .AND. L.NE.ND) GOTO 12
      PRINT 13,L,TNEW(L)
   13 FORMAT (40X,I10,F15.0)
   12 CONTINUE
      ENDIF
 
      RETURN
      END
