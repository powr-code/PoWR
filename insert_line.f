      SUBROUTINE INSERT_LINE 
     >                 (LINE, INDBL, NBLINE, INDLAP, XLAMLAP, DELXLAP,
     $                   XLAMBL, LINPROBL, AVOIGTBL, XLAM, MAXLAP, ALN,
     >                   ND, LINPRO, AVOIGT, NMOD, NDDIM, MAXMOD )

C*******************************************************************************
C***  This subroutine inserts an additional line in the list
C***  such that it stays sorted by ascending DELTAX
C***  Various line parameters are inserted accordingly in
C***  their respective arrays:
C***     INDLAP, XLAMLAP, LINPRO, AVOIGT 
C***  CALLED FROM FORMAL - PREFORM  
C*******************************************************************************

      IMPLICIT NONE

      INTEGER I, LOCDXGT, NDDIM, MAXLAP, MAXMOD, L, INDBL, NBL, LINE
      INTEGER ISRCHFGT, IMOD, NMOD, NBLINE, NPAR
      REAL ALN, AVOIGTBL, XLAMBL, XLAM, DELTAX 

      INTEGER, DIMENSION (MAXMOD) :: ND
      INTEGER, DIMENSION (MAXLAP) :: INDLAP
      REAL,    DIMENSION (MAXLAP) :: XLAMLAP, DELXLAP
      REAL, DIMENSION (MAXLAP,NDDIM,MAXMOD) :: AVOIGT
      CHARACTER KARTE*80
      CHARACTER*8 LINPRO(MAXLAP), PROF, PROFDEFAULT, LINPROBL
      INTEGER, PARAMETER :: NPARMAX = 4 
      CHARACTER*10 ACTPAR(NPARMAX)

C***  If this is the very first line: initialize the vectors
      IF (NBLINE .EQ. 0) THEN
         DELTAX = .0
         LOCDXGT = 1
         LINE = INDBL
      ELSE
         DELTAX = ALOG(XLAMBL/XLAM) / ALN
         LOCDXGT=ISRCHFGT(NBLINE,DELXLAP,1,DELTAX)
         IF (LOCDXGT .LE. NBLINE) THEN
            CALL SHIFT (INDLAP,LOCDXGT,NBLINE)
            CALL SHIFT (XLAMLAP,LOCDXGT,NBLINE)
            CALL SHIFT (DELXLAP,LOCDXGT,NBLINE)
            CALL SHIFT (LINPRO,LOCDXGT,NBLINE)
            DO IMOD=1, NMOD
               DO L=1, ND(IMOD)
                  CALL SHIFT (AVOIGT(1,L,IMOD),LOCDXGT,NBLINE)
               ENDDO
            ENDDO
          ENDIF
      ENDIF

      NBLINE=NBLINE+1

      IF (NBLINE .GE. MAXLAP) THEN
               PRINT *,'*** ERROR: MAX. NUMBER OF BLEND'
     $                ,' COMPONENTS EXCEEDED: MAXLAP=', MAXLAP
               WRITE (0,*) '*** ERROR: MAX. NUMBER OF BLEND'
     >                     ,' COMPONENTS EXCEEDED: MAXLAP=', MAXLAP
               STOP '*** FATAL ERROR detected by subd. INSERT_LINE' 
      ENDIF

      INDLAP(LOCDXGT)=INDBL
      XLAMLAP(LOCDXGT)=XLAMBL
      DELXLAP(LOCDXGT)=DELTAX
      LINPRO(LOCDXGT)=LINPROBL
      DO IMOD=1, NMOD
         AVOIGT(LOCDXGT,1,IMOD)=AVOIGTBL
      ENDDO

      RETURN
      END
