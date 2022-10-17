      SUBROUTINE QUADSTARK (GAMMAQUAD, T, XNE, NCHARG, ELEVEL, EION, 
     >                      LINPRO, LEVELLOW, LEVELNUP)
C**********************************************************************
C***  Quadratic Stark effect
C***     following Cowley (1971: Observatory 91, 139)
C***  GAMMAQUAD     - damping parameter - Lorentz profil
C***  T             - temperature
C***  XNE           - electron density
C***  ELEVEL        - Excitation energy of the upper level  
C***  NCHARG        - Ion charge
C***   Note: Z in Cowley's formula is "Charge seen by the active electron"  
C***           i.e. Z = NCHARG + 1
C**********************************************************************

      CHARACTER*10 LEVELLOW, LEVELNUP
      CHARACTER*8  LINPRO

C***  Rydberg Energy (infinite mass) in Kaiser
      DATA RYD / 109 737.3 / 
C** The constant is: 8/sqrt(3) * pi  * hbar^2 * m_e^(-3/2) * k_B^(-1/2)  in cgs (See also Eqs. (11) & (23) in Freudenstein+77,  224)
      DATA COWLEY_FAC / 5.0E-5 /
C***  The effective quantum number is calculated with reference to the
C***  ionization energy. However, EION may be not known if there is no higher 
C***  ion in the data than the considered one. 
      IF (EION .LE. .0) THEN
           LINPRO = 'VOIGT   '
           GAMMAQUAD = .0
           WRITE (0,*) '*** WARNING: QUADSTARK cannot handle ', 
     >           LEVELLOW, ' - ', LEVELNUP, 
     >           '  (ionization energy not known)'
           RETURN
      ENDIF 

C***  "distance from continuum" in kaiser
      DELTAE = EION - ELEVEL


C***    WARNING !!!!!!!!
C***    eigentlich muesste man die Energiedifferenz bei doppelt 
C***    angeregten Zustaenden zu dem entsprechenden angeregten Level des 
C***    oberen Ions bilden, siehe Cowley. Diese Information zu beschaffen
C***    erfordert aber einiges Nachdenken.

C***  If energy is close or above 100 Kayser: skip QUADSTARK
C***  use VOIGT as fallback 
      IF (DELTAE .LT. 100) THEN
           LINPRO = 'VOIGT   '
           GAMMAQUAD = .0
           WRITE (0,*) '*** WARNING: QUADSTARK cannot handle ', 
     >           LEVELLOW, ' - ', LEVELNUP, 
     >           '  (upper level too high)'
           RETURN
      ENDIF 
        
C***  primary effective quantum number n* of the upper level
      EFFQN = (NCHARG+1) * SQRT (RYD / DELTAE)

C**   This Gamma is the FWHM of the Lorentz profile 
C***  (see Freudenstein+77,  224, 2w = GAMMAQUAD = FWHM)
C**   Note: It might be more accurate to use T = const = 10kK, see 
C***  discussion by Cowley + Freudenstein. -- Tomer, 18.11.2014
      GAMMAQUAD = COWLEY_FAC * XNE / SQRT(T) * (EFFQN**2/(NCHARG+2))**2

      RETURN
      END

