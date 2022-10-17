      SUBROUTINE PLOTHSUM (HTOTL, HTOTM, HTOTG, HTOTOBS, HTOTCMF0,
     >                     ND, MODHEAD, JOBNUM, KANAL, TEFF, BHTOTERR,
     >                     FLUXEPS)
C******************************************************************************
C***  DIRECT TRANSFER OF HSUM PLOT
C***  TOTAL (FREQUENCY-INTEGRATED) FLUX versus DEPTH INDEX
C***  Note: there are some more curves available, but de-activated
C***        wrh 14-Apr-2003 11:58:53
C******************************************************************************
    
      IMPLICIT NONE

      INTEGER, PARAMETER :: NDMAX = 200

      INTEGER, INTENT(IN) :: ND, KANAL, JOBNUM
      REAL, INTENT(IN) :: TEFF, FLUXEPS

      CHARACTER(100) :: MODHEAD 
      CHARACTER(110) :: HEADLINE
      CHARACTER(8) :: CENTER

      REAL, DIMENSION(NDMAX) :: X, Y
      REAL, DIMENSION(ND) :: HTOTL, HTOTM, HTOTG, HTOTOBS, HTOTCMF0
      LOGICAL :: BHTOTERR

      INTEGER :: L
      REAL :: XMIN, XMAX, YMIN, YMAX, TRAD, HSUM

C***  STEBOL = STEFAN-BOLTZMANN CONSTANT / PI (ERG/CM**2/SEC/STERAD/KELVIN**4
      REAL, PARAMETER :: STEBOL = 1.8046E-5

      IF (BHTOTERR) THEN
         WRITE (*, '(A)')
     >         'PLOTHTOT: WARNING: HTOTL  not present in MODEL file'
     >         //' --  PLOT disabled!'
         RETURN
      ENDIF 


      IF (ND .GT. NDMAX) THEN
        WRITE (0, '(A)') 'DIMENSION INSUFFICIENT - HSUM PLOT SUPPRESSED'       
        WRITE (0, '(A)') 'NON-FATAL ERROR IN SUBROUTINE PLOTHSUM'
        RETURN
      ENDIF
 
      CENTER = CHAR(92) // 'CENTER' // CHAR(92)
      CALL JSYMSET ('G2','TRANSFER')
 
      HEADLINE = 'HSUM: M'//MODHEAD(13:)
      WRITE (HEADLINE(90:), '(A8,I7)') ' JOB No.', JOBNUM


C***  HTOTOBS: Conversion into radiation temperatures
      DO L=1, ND-1
         X(L) = FLOAT(L)
         HSUM = HTOTOBS(L)
         TRAD = (4. * ABS(HSUM) / STEBOL)**0.25
         IF (HSUM .LT. .0) TRAD = -TRAD
         Y(L) = TRAD / 1000.
      ENDDO

cc      CALL PLOTANF (KANAL,HEADLINE, '&E'//HEADLINE,
cc     $        CENTER//'DEPTH INDEX L',
cc     $        CENTER//'T&Trad&M / KK',
cc     >             0., XMIN, XMAX, 5., 10., 0.,
cc     >             0., YMIN, YMAX, .2,  1., 0.,
cc     $        X,Y,ND-1,5)

C***  HTOT+HMECH+HGRAV: Conversion into radiation temperatures
cc      DO L=1, ND-1
cc         HSUM = HTOTOBS(L)+HTOTM(L)+HTOTG(L)
cc         TRAD = (4. * ABS(HSUM) / STEBOL)**0.25
cc         IF (HSUM .LT. .0) TRAD = -TRAD
cc         Y(L) = TRAD / 1000.
cc      ENDDO
cc      CALL PLOTCONS (KANAL,X,Y,ND-1,'COLOR= 3') 

C***  HTOT+HGRAV: Conversion into radiation temperatures
cc      DO L=1, ND-1
cc         HSUM = HTOTOBS(L)+HTOTG(L)
cc         TRAD = (4. * ABS(HSUM) / STEBOL)**0.25
cc         IF (HSUM .LT. .0) TRAD = -TRAD
cc         Y(L) = TRAD / 1000.
cc      ENDDO
cc      CALL PLOTCONS (KANAL,X,Y,ND-1,'COLOR= 3') 

      WRITE (KANAL,'(A)') 'PLOT: ' // HEADLINE
      WRITE (KANAL,'(A)') '\DEFINECOLOR 8 = 0.8 0.8 0.8'      
C***  HTOTCMF0: Conversion into radiation temperatures
      DO L=1, ND-1
         Y(L) = HTOTCMF0(L) 
         TRAD = (4. * ABS(Y(L)) / STEBOL)**0.25
         IF (Y(L) .LE. .0) TRAD = -TRAD
         Y(L) = TRAD / 1000.
      ENDDO      

      CALL PLOTANFS (KANAL,HEADLINE, '&E'//HEADLINE,
     $        CENTER//'DEPTH INDEX L',
     $        CENTER//'T&Trad&M / kK',
     >             0., 0., 0.,0., 0., 0.,
     >             0., 0., 0.,0., 0., 0.,
     $        X,Y,ND-1, 'PEN=8')
      
      IF (FLUXEPS .GT. 0.) THEN
         DO L=1, ND-1
            X(L) = FLOAT(L)
            Y(L) = Y(L) * ((1. + FLUXEPS)**0.25)
         ENDDO
         CALL PLOTCONS (KANAL,X,Y,ND-1,'COLOR= 8 PEN = 3') 
         DO L=1, ND-1
            X(L) = FLOAT(L)
            Y(L) = Y(L) *(((1. - FLUXEPS)/(1. + FLUXEPS))**0.25)
         ENDDO
         CALL PLOTCONS (KANAL,X,Y,ND-1,'COLOR= 8 PEN = 3') 
      ENDIF

C***  HTOTL: Conversion into radiation temperatures
      DO L=1, ND-1
         X(L) = FLOAT(L)
         Y(L) = HTOTL(L) 
         TRAD = (4. * ABS(Y(L)) / STEBOL)**0.25
         IF (Y(L) .LE. .0) TRAD = -TRAD
         Y(L) = TRAD / 1000.
      ENDDO
      CALL PLOTCONS (KANAL,X,Y,ND-1,'COLOR= 2 PEN = 3') 

C***  Continuum
      X(1) = 1.
      X(2) = FLOAT(ND-1)
      Y(1) = TEFF / 1000.
      Y(2) = TEFF / 1000. 
      CALL PLOTCONS (KANAL,X,Y,2, 'COLOR=3') 

C***  Add an invisible dataset to span TEFF +- 5000
      Y(1) = TEFF/1000. + 5
      Y(2) = TEFF/1000. - 5
      CALL PLOTCON (KANAL,X,Y,2,0) 
      
      RETURN
      END
