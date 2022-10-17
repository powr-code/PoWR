      FUNCTION KHOLTSMARK(IQNLOW, IQNUP, IZ)
C***********************************************************************      
C***  returns the Holtsmark wing broadening coefficient K_low,up
C***  (required for the treatment of linear Stark broadening)
C***
C***  IZ = NUCLEAR CHARGE 
C***       (note: this is sometimes called "ion charge" in the literature)
C***
C***  table and calculation taken from H. Griem, 1960, ApJ 132, p.883
C***   (The Griem calculations are an approximative method 
C***    for the results described in Unsoeld, Page 320ff)
C***
C***  Please note that K_low,up is NOT dimensonless, but instead has
C***   the unit of [ c*h^7 / (m^3 * e^9) ]  = cm^(3/2) * s * g^(-1/2)
C***   [All constant units have to evaluated in cgs]
C***
C***  called by LINSTARK
C***********************************************************************      

      IMPLICIT NONE
      
      INTEGER, INTENT(IN) :: IQNLOW, IQNUP, IZ
      REAL :: KHOLTSMARK
      INTEGER :: L1, U1

C**** Table 1 from Griem for KHOLTSMARK 
C        ( -1  = no valid entry )
      REAL, DIMENSION(2:18,1:4) :: KTAB = 
     >   RESHAPE(  (/
     >      0.293E-3,       -1.,     -1.,    -1.,
     >      0.560E-3,  0.143E-1,     -1.,    -1.,
     >      0.940E-3,  0.188E-1,   0.163,    -1.,
     >      1.430E-3,  0.262E-1,   0.174,    0.98,
     >      2.040E-3,  0.356E-1,   0.213,    0.91,
     >      2.750E-3,  0.470E-1,   0.267,    1.02,
     >      3.580E-3,  0.600E-1,   0.332,    1.20,
     >      4.510E-3,  0.720E-1,   0.406,    1.42,
     >      5.600E-3,  0.880E-1,   0.49 ,    1.68,
     >      6.700E-3,  1.100E-1,   0.58 ,    1.96,
     >      8.000E-3,  1.290E-1,   0.68 ,    2.28,
     >      9.400E-3,  1.520E-1,   0.80 ,    2.63,
     >     10.800E-3,  1.760E-1,   0.92 ,    3.01,
     >     12.400E-3,  2.010E-1,   1.04 ,    3.41,
     >     14.100E-3,  2.280E-1,   1.18 ,    3.84,
     >     16.000E-3,  2.570E-1,   1.33 ,    4.31,
     >     17.900E-3,  2.880E-1,   1.48 ,    4.8
     >  /) , (/ 17, 4 /), ORDER=(/ 2, 1 /) )
      
      !File and channel handles (=KANAL)
      INTEGER, PARAMETER :: hCPR = 0        !write to wruniqX.cpr (stderr)
      INTEGER, PARAMETER :: hOUT = 6        !write to wruniqX.out (stdout)
      
      
C***  STOP if called with invalid quantum numbers      
      IF (IQNLOW <= 0) THEN
        WRITE (hCPR,*) 'ERROR: IQNLOW <= 0'
        WRITE (hCPR,'(2(A,I3))') 'IQNLOW = ', IQNLOW, ' IQNUP = ', IQNUP
        STOP 'FATAL ERROR IN KHOLTSMARK'
      ENDIF
      IF (IQNUP <= 0) THEN
        WRITE (hCPR,*) 'ERROR: IQNUP <= 0'
        WRITE (hCPR,'(2(A,I3))') 'IQNLOW = ', IQNLOW, ' IQNUP = ', IQNUP
        STOP 'FATAL ERROR IN KHOLTSMARK'
      ENDIF
C***  STOP if called with reserve QN roles
      IF (IQNLOW >= IQNUP) THEN
        WRITE (hCPR,*) 'ERROR: IQNLOW >= IQNUP NOT ALLOWED'
        WRITE (hCPR,'(2(A,I3))') 'IQNLOW = ', IQNLOW, ' IQNUP = ', IQNUP
        STOP 'FATAL ERROR IN KHOLTSMARK'
      ENDIF

C***  Determine if table values can be used
      IF (IQNUP <= 18 .AND. IQNLOW <= 4) THEN
c        KHOLTSMARK = KTAB(IQNLOW, IQNUP)
        KHOLTSMARK = KTAB(IQNUP, IQNLOW)
      ELSE
C***    asymptotic formula for everything else
C        (note that Hubeny incorrectly writes 5.5 * 1.E-4 in his 1994 paper, 
C         but recalculating Griem's equations leads to 5.5 * 1.E-5 as he writes in his 1960 paper)
        KHOLTSMARK = 5.5 * 1.E-5 
     >        * (IQNUP*IQNLOW)**4 / (IQNUP*IQNUP - IQNLOW*IQNLOW)
      ENDIF
      
C***  Correction for hydrogenic ions beyond hydrogen
C***   after Hubeny & Mihalas (Book, 2015), page 260, Eq. (8.151)
C***   see also Eq. (33) in Griem (1960)
C***  Note: This factor can either be put on K_low,up or Holtsmark normal field strength F0
C***        SYNSPEC (formal integral for TLUSTY models) puts it in F0, Griem puts it on K_low,up
      KHOLTSMARK = KHOLTSMARK / (FLOAT(IZ)**5)
      
      
      RETURN
      END
