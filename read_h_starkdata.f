      SUBROUTINE READ_H_STARKDATA (PATH_LEMKE_DAT,LOW,NUP,LEMKEFOUND,       
     >      STKTA_DWS, NWS, STKTA_TS, NTS,
     >      STKTA_ES,  NES, STKTA_PS, NPS,
     >      NWS_DIM, NTS_DIM, NES_DIM, NPS_DIM)
C*****************************************************************
C     Preparatory subroutine for stark-broadening of hydrogen lines
C     This routines reads for the specified line the broadening table 
C      LEMKE_HI.DAT (default path: /home/corona/wrh/work/wrdata/) 
C      from Simbad, created by Michael Lemke
C      see -->  Lemke M.: 1997, A&A Suppl. 122, 285
C     The table comprises the first4 spectral series
C      (Lyman, Balmer, Bracket, Paschen) fur upper levels till n=22
C
C     Content of this subroutine:
C     The file is opened, the block for the line specified by the 
C      principle quantum numbers (nup,low) is searched, and the 
C      data arrays are read and stored:
C       - One vector each for the three axes 
C          ("scaled wavelength", temperature, el. density)
C       - One matrix (here as vector) with the three indices   
C         for the broadening profile
C     * Note: The profile is NOT NORMALIZED 
C     The subroutine returns LEMKEFOUND=1 when the line was found, 
C         and LEMKEFOUND=0 else. 
C******************************************************************

      CHARACTER*256 PATH_LEMKE_DAT, FILENAME 
      CHARACTER*3 LOCAL_SPECIES

      DIMENSION STKTA_DWS(NWS), STKTA_TS (NTS), STKTA_ES (NES)
      DIMENSION STKTA_PS (NPS)

      LOGICAL STKTA_QHALF      !If TRUE, on half of profile tabulated
      LOGICAL STKTA_WSCA       ! Lambda (in A) o
C***  Note: Both are always true in the used table

      IF (PATH_LEMKE_DAT .EQ. 'default') THEN
         FILENAME = '/home/corona/wrh/work/wrdata/LEMKE_HI.DAT'
      ELSE
         FILENAME = 
     >    PATH_LEMKE_DAT(:IDX(PATH_LEMKE_DAT)) // '/' // 'LEMKE_HI.DAT'
      ENDIF
      KANAL=13


C***  Fortran channel
      LUSTK = 80

      LOCAL_SPECIES=' '
      OPEN(UNIT=LUSTK,FILE=FILENAME, FORM='FORMATTED',
     &               STATUS='OLD',ACTION='READ',ERR=99)
      DO 
           READ(LUSTK,*,END=90) LOCAL_SPECIES(1:3),
     &               LOW_STRK,NUP_STRK,STKTA_QHALF,STKTA_WSCA

            READ(LUSTK,*,END=90) NWS, NTS, NES, NPS
C*** check if arrays are dimension sufficiently large
            IF (NWS .GT. NWS_DIM) GOTO 91
            IF (NTS .GT. NTS_DIM) GOTO 92
            IF (NES .GT. NES_DIM) GOTO 93
            IF (NPS .GT. NPS_DIM) GOTO 94

            READ(LUSTK,*, END=90) (STKTA_DWS(I),I=1,NWS)
            READ(LUSTK,*, END=90) (STKTA_TS(I), I=1,NTS)
            READ(LUSTK,*, END=90) (STKTA_ES(I), I=1,NES)
            READ(LUSTK,*, END=90) (STKTA_PS(I), I=1,NPS)

C***  Check if wanted profile was found, then exit reading
          IF (LOCAL_SPECIES .EQ. 'HI ' .AND.
     >        LOW_STRK .EQ. LOW .AND. NUP_STRK .EQ. NUP) EXIT

      END DO

      CLOSE(UNIT=LUSTK)
      LEMKEFOUND = 1
      RETURN


C***  Set flag if this line was not found in the data table
   90 CONTINUE
      LEMKEFOUND = 0
      RETURN

C*******************************************************
   99 WRITE (0,*) '*** ERROR: cannot open file ', 
     >        FILENAME(:IDX(FILENAME))
      WRITE (0,*) '*** ERROR: Broadening data for H I not available'
      GOTO 100

   91 WRITE (0,*) 'DIMENSION NWS_DIM TOO SMALL:',
     >  ' Needed: ', NWS 
      GOTO 100

   92 WRITE (0,*) 'DIMENSION NTS_DIM TOO SMALL:',
     >  ' Needed: ', NTS
      GOTO 100

   93 WRITE (0,*) 'DIMENSION NES_DIM TOO SMALL:',
     >  ' Needed: ', NES 
      GOTO 100

   94 WRITE (0,*) 'DIMENSION NPS_DIM TOO SMALL:',
     >  ' Needed: ', NPS 
      GOTO 100


  100 STOP '*** FATAL ERROR in READ_H_STARKDATA'

      END
