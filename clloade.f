      SUBROUTINE CLLOADE (NCHANE, NZE1, NZE2, NFRO, EDDIA, NDEDDIA, 
     >                    EDDIF, EDDIG, ND, 
     >                    EDDIHOUT, EDDIHIN, EDDINOUT, EDDININ, 
     >                    BCLEERR, BCOLIP, XHI, XHO, EPSG, 
     >                    XHOM, XNOM, EDDIHOUTP, EDDINOUTP)
C**********************************************************************
C***  Reads the old EDDIEs from file fort.<NCHANE>
C***  If this file does not exist or is inappropriate (BCLEERR = .TRUE.)
C***    then the EDDIEs are not read but the counters NZE1 and NZE2 are
C***    updated
C**********************************************************************

      DIMENSION EDDIA(NDEDDIA, NFRO)
      DIMENSION EDDIF(ND), EDDIG(ND-1)
      DIMENSION EPSG(ND)
      CHARACTER NAME*8
      LOGICAL BCLEERR, BCOLIP

C***  In the first call of this routine, the EDDIMIX-Factors are read
C***    (or set ZERO, if not present on old EDDI file)
      IF (NZE1 .EQ. 0) THEN
         IF (.NOT. BCLEERR) THEN
            CALL READMS (NCHANE, EPSG, ND-1, 'EPSG    ', IERR)
         ENDIF
         IF (BCLEERR .OR. IERR .EQ. -10) THEN
            DO L=1, ND-1
               EPSG(L) = .0
            ENDDO
         ENDIF
      ENDIF

C***
      IF (NZE2 .EQ. NFRO) THEN
        NZE2 = 0
      ENDIF
      NZE2 = NZE2 + 1
      IF (NZE2 .EQ. 1) THEN
C***  Read next array from file
        NZE1 = NZE1 + 1
        IF (.NOT. BCLEERR) THEN
          WRITE (NAME,'(A2,I5,A1)') 'ED', NZE1, ' '
          CALL READMS (NCHANE, EDDIA, NDEDDIA*NFRO, NAME, IERR)
          IF (IERR .NE. 0) THEN
            WRITE (0,'(A,I2)') 
     >        'Error when reading from file fort.', NCHANE
            STOP 'ERROR in Subr. CLLOADE'
          ENDIF
        ENDIF
      ENDIF

C***  Restore EDDIs if data is present
      IF (.NOT. BCLEERR) THEN
        NDA = 2*ND - 1
        DO L=1, ND
          EDDIF(L) = EDDIA(L, NZE2)
        ENDDO
        DO L=1, ND-1
          EDDIG(L) = EDDIA(ND+L, NZE2)
        ENDDO
        EDDIHOUT = EDDIA(NDA+1, NZE2)
        EDDIHIN  = EDDIA(NDA+2, NZE2)
        EDDINOUT = EDDIA(NDA+3, NZE2)
        EDDININ  = EDDIA(NDA+4, NZE2)
        XHI      = EDDIA(NDA+5, NZE2)
        IF (.NOT. BCOLIP) THEN
          XHO      = EDDIA(NDA+6, NZE2)
        ENDIF
        XHOM     = EDDIA(NDA+7, NZE2)
        XNOM     = EDDIA(NDA+8, NZE2)
        EDDIHOUTP= EDDIA(NDA+9, NZE2)
        EDDINOUTP= EDDIA(NDA+10, NZE2)
      ENDIF

      RETURN
      END
