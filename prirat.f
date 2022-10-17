      SUBROUTINE PRIRAT (ITNE,N,LEVEL,NDIM,L,CRATE,RRATE,RATCO,EN,
     $           IFRRA,MODHEAD,JOBNUM,NETTO,NFIRST,NLAST,NATOM,NATOUT,
     $           NAUTO,RDIEL,RAUTO,IONGRND,KODRLOW,LASTKDR)
C***********************************************************************
C***  OUTPUT OF RATE MATRIX ON FT11   **********************************
C***********************************************************************

      DIMENSION EN(NDIM)
      DIMENSION NFIRST(NATOM),NLAST(NATOM)
      DIMENSION CRATE(NDIM,NDIM),RRATE(NDIM,NDIM),RATCO(NDIM,NDIM)
      DIMENSION RDIEL(NDIM),RAUTO(NDIM),IONGRND(NDIM)
      DIMENSION KODRLOW(LASTKDR)
      CHARACTER LEVEL(NDIM)*10
      CHARACTER MODHEAD*100
      CHARACTER*12 NUMBERS(10)
 
 
C***  OPEN OUTPUT FILE AND WRITE HEADLINE FOR THE RATE PRINTOUT
      IF (L.EQ.1.OR.L.EQ.IFRRA) THEN
            OPEN (11,FILE='RATES', STATUS='UNKNOWN')
      ENDIF
      WRITE (11,14)  MODHEAD,JOBNUM
   14 FORMAT (1X,  A  ,20X,'JOB NO.',I7,//)
 
C***  OUTPUT ONLY FOR ONE OR ALL ELEMENTS
      IF (NATOUT .NE. 0) THEN
         NASTART=NATOUT
         NASTOP=NATOUT
      ELSE
         NASTART=1
         NASTOP=NATOM
      ENDIF
         
C***  LOOP OVER SPECIFIED ELEMENTS  ==========================================
      DO 99 NA=NASTART,NASTOP
      NFIRNA=NFIRST(NA)
      NLANA=NLAST(NA)

C*******************************************************************************
C***  ADD RADIATIVE AND COLLISIONAL RATE COEFF. INTO MATRIX RATCO
      DO 11 I=NFIRNA,NLANA
      DO 11 J=NFIRNA,NLANA
   11 RATCO(I,J)=CRATE(I,J)+RRATE(I,J)
 
C***  ADD ADDITIONAL D-R TERMS INTO RATE COEFFICIENT MATRIX RATCO
      DO 12 KDR=1,LASTKDR
      LOW=KODRLOW(KDR)
      IF ((LOW .LT. NFIRNA) .OR. (LOW .GT. NLANA)) GOTO 12
      NUP=IONGRND(LOW)
      RATCO(LOW,NUP)=RATCO(LOW,NUP)+RAUTO(LOW)
      RATCO(NUP,LOW)=RATCO(NUP,LOW)+RDIEL(LOW)
   12 CONTINUE

      IF (NETTO .NE. 1) THEN
      DO 50 J1=NFIRNA,NLANA,10
      J2=MIN0(NLANA,J1+9)
      WRITE (11,1) L, (LEVEL(J),J=J1,J2)
    1 FORMAT(//,20X,
     $       'TOTAL (COLL.+RAD.+D-R) RATE COEFF. AT DEPTH POINT NO.',
     $       I3,//,11X,10(2X,A10))
      WRITE (11,4)
    4 FORMAT (1X)
      DO 5 I=NFIRNA,NLANA
      ENCODE (120,40,NUMBERS) (RATCO(I,J),J=J1,J2)
      DO 43 M=1,10
      IF (NUMBERS(M) .EQ. '    0.00E+00') NUMBERS(M)='    0       '
   43 CONTINUE
      WRITE (11,42) LEVEL(I),(NUMBERS(M),M=1,10)
    5 CONTINUE
   50 CONTINUE
      ENDIF
 
C*******************************************************************************
C***  RATIO COLLISIONAL / RADIATIVE RATE COEFFICIENTS
      IF (NETTO .NE. 1) THEN
      DO 60 J1=NFIRNA,NLANA,10
      J2=MIN0(NLANA,J1+9)
      WRITE (11,61) L, (LEVEL(J),J=J1,J2)
   61 FORMAT (//,20X,'RATIO COLLISIONAL / RADIATIVE RATE COEFF. AT L=',
     $        I3,//,11X,10(2X,A10))
      WRITE (11,4)
      DO 65 I=NFIRNA,NLANA
      DO 66 J=J1,J2
      M=1+J-J1
      IF (RRATE(I,J) .NE. .0) THEN
         ENCODE (12,40,NUMBERS(M)) CRATE(I,J)/RRATE(I,J)
         IF (NUMBERS(M) .EQ. '    0.00E+00') NUMBERS(M)='    0       '
         ELSE
         NUMBERS(M)=' '
         ENDIF
   66 CONTINUE
      MMAX=1+J2-J1
      WRITE (11,42) LEVEL(I),(NUMBERS(M),M=1,MMAX)
   65 CONTINUE
   60 CONTINUE
      ENDIF

C*******************************************************************************
C***  RATIO D-R / (RADIATIVE + COLLISIONAL) RATE COEFFICIENTS
      IF ((NETTO .NE. 1) .AND. (NAUTO .GT. 0)) THEN
      DO 70 J1=NFIRNA,NLANA,10
      J2=MIN0(NLANA,J1+9)
      WRITE (11,71) L, (LEVEL(J),J=J1,J2)
   71 FORMAT (//,20X,
     $        'RATIO D-R / (COLLISIONAL + RADIATIVE) RATE COEFF. AT L=',
     $        I3,//,11X,10(2X,A10))
      WRITE (11,4)
      DO 75 I=NFIRNA,NLANA
      DO 76 J=J1,J2
      DRRATE=.0
      IF ((I .LT. J) .AND. (J .EQ. IONGRND(I))) DRRATE=RAUTO(I)
      IF ((I .GT. J) .AND. (I .EQ. IONGRND(J))) DRRATE=RDIEL(J)
      M=1+J-J1
      IF (RRATE(I,J)+CRATE(I,J) .NE. .0) THEN
         ENCODE (12,40,NUMBERS(M)) DRRATE/(CRATE(I,J)+RRATE(I,J))
         IF (NUMBERS(M) .EQ. '    0.00E+00') NUMBERS(M)='    0       '
         ELSE
         NUMBERS(M)=' '
         ENDIF
   76 CONTINUE
      MMAX=1+J2-J1
      WRITE (11,42) LEVEL(I),(NUMBERS(M),M=1,MMAX)
   75 CONTINUE
   70 CONTINUE
      ENDIF
 
C*******************************************************************************
C***  ADDING (ABSOLUTE) RATES INTO MATRIX RATCO
      DO 20 I=NFIRNA,NLANA
      ENI=EN(I)
      DO 20 J=NFIRNA,NLANA
   20 RATCO(I,J)=ENI*RATCO(I,J)
 
      IF (NETTO .NE. 1) THEN
      DO 51 J1=NFIRNA,NLANA,10
      J2=MIN0(NLANA,J1+9)
      WRITE (11,21) L, (LEVEL(J),J=J1,J2)
   21 FORMAT(//,20X,'  TOTAL RATES                  AT DEPTH POINT NO.',
     $       I3,//,11X,10(2X,A10))
      WRITE (11,4)
      DO 22 I=NFIRNA,NLANA
      ENCODE (120,40,NUMBERS) (RATCO(I,J),J=J1,J2)
      DO 44 M=1,10
      IF (NUMBERS(M) .EQ. '    0.00E+00') NUMBERS(M)='    0       '
   44 CONTINUE
      WRITE (11,42) LEVEL(I),(NUMBERS(M),M=1,10)
   22 CONTINUE
   51 CONTINUE
      ENDIF
 
C*******************************************************************************
      DO 52 J1=NFIRNA,NLANA,10
      J2=MIN0(NLANA,J1+9)
      WRITE (11,3) L, (LEVEL(J),J=J1,J2)
    3 FORMAT(//,20X,'  NET RATES                    AT DEPTH POINT NO.',
     $       I3,//,11X,10(2X,A10))
      WRITE (11,4)
      DO 8 I=NFIRNA,NLANA
      ENCODE (120,40,NUMBERS) (RATCO(I,J)-RATCO(J,I),J=J1,J2)
   40 FORMAT (1P,10E12.2)
      DO 41 M=1,10
      IF (NUMBERS(M) .EQ. '    0.00E+00') NUMBERS(M)='    0       '
   41 CONTINUE
      WRITE (11,42) LEVEL(I),(NUMBERS(M),M=1,10)
   42 FORMAT (1X,A10,10A12)
    8 CONTINUE
   52 CONTINUE
 
C*******************************************************************************
      DO 53 J1=NFIRNA,NLANA,10
      J2=MIN0(NLANA,J1+9)
      WRITE (11,9) L, (LEVEL(J),J=J1,J2)
    9 FORMAT(//,20X,'   RELATIVE NET RATES          AT DEPTH POINT NO.',
     $       I3,//,11X,10(2X,A10))
      WRITE (11,4)
      DO 10 I=NFIRNA,NLANA
         SCALE=0.
         DO 30 J=NFIRNA,NLANA
   30       SCALE=SCALE+ABS(RATCO(I,J)-RATCO(J,I))
         SCALE=SCALE/2.
         IF (SCALE.LE.0.) SCALE=1.
      ENCODE (120,40,NUMBERS) ((RATCO(I,J)-RATCO(J,I))/SCALE,J=J1,J2)
      DO 45 M=1,10
      IF (NUMBERS(M) .EQ. '    0.00E+00') NUMBERS(M)='    0       '
   45 CONTINUE
      WRITE (11,42) LEVEL(I),(NUMBERS(M),M=1,10)
   10 CONTINUE
   53 CONTINUE
 
   99 CONTINUE
C***  ENDLOOP  =========================================================

      RETURN
      END
