      SUBROUTINE STARKDAMP_HEI_NEUTRAL 
     >       (NUP, LOW, TEMP, HE1FRC, LINPRO, HEWID, MAINQN, LEVEL)
C***************************************************************************  
C
C  Ca;lculation of the Stark-effect damping parameter for HeI lines
C
C  This subroutine only accounts for the collisions with other *neutral*
C  helium particles. This is probably a minor effect, compared to the 
C  contribution by electrons. Thesefore I do not implement this routine.
C                                                     wrh 10-Jul-2008  
C 
C  Compute approximate broadening by He collisions (p-d approx)
C  
C  NUP,NLOW - pinciple quantum numbers
C  TEMP     - Temperatur
C  HE1FRC   - neutral helium number density in cm^-3
C  HEWID    - lorentz type broadening parameter simply add to
C             radiation-damping parameter 
C
C  mario parade may 2008
C
C  Output in usual rad /s cm^3
C   
C  27/4/00 PB
C

      REAL PI, C, HE1FRC, GAMMAF
      REAL SIGMAH(4), ALPHAH(4), X, GX, A0, SIGMA, VBAR, KB,M0,TEMP
      REAL HEWID, GVW
      DIMENSION MAINQN(2)
      CHARACTER*8 LINPRO
      CHARACTER*10 LEVEL(2)

      PARAMETER (C=2.997925E+18,PI=3.14159265,A0=5.29177249E-11)    
      PARAMETER (KB=1.380658E-23,M0=1.660540E-27)

      DATA SIGMAH/ 834., 1998., 3837., 5700. /
      DATA ALPHAH/ .280, .330, .260, .250 /


      NNUP = MAINQN(NUP) 
      NLOW = MAINQN(LOW)
      LINE = NNUP - NLOW

C***  This routine obviously only covers differences 
C***  in the principal quantum up to 4  - wrh  9-Jul-2008 
C***  Moreover, it is checked whether both MAINQN are known
      IF (LINE .GT. 4 .OR. NNUP .LE. 0 .OR. NLOW .LE. 0) THEN
        IF (LINE .GT. 4) THEN
           WRITE (*,'(A)') '*** WARNING: ' // 
     >      'Stark broadening cancelled for '
     >      // LEVEL(LOW) // ' - ' // LEVEL(NUP) // ': delta-n > 4'
           WRITE (0,'(A)') '*** WARNING: ' // 
     >      'Stark broadening cancelled for '
     >      // LEVEL(LOW) // ' - ' // LEVEL(NUP) // ': delta-n > 4'
        ELSE IF (NNUP .LE. 0) THEN
           WRITE (*,'(A)') '*** WARNING: ' // 
     >      'Stark broadening cancelled for '
     >      // LEVEL(LOW) // ' - ' // LEVEL(NUP) // ': no MAINQN(NUP)'
           WRITE (0,'(A)') '*** WARNING: ' // 
     >      'Stark broadening cancelled for '
     >      // LEVEL(LOW) // ' - ' // LEVEL(NUP) // ': no MAINQN(NUP)'
        ELSE IF (NLOW .LE. 0) THEN
           WRITE (*,'(A)') '*** WARNING: ' // 
     >      'Stark broadening cancelled for '
     >      // LEVEL(LOW) // ' - ' // LEVEL(NUP) // ': no MAINQN(LOW)'
           WRITE (0,'(A)') '*** WARNING: ' // 
     >      'Stark broadening cancelled for '
     >      // LEVEL(LOW) // ' - ' // LEVEL(NUP) // ': no MAINQN(LOW)'
        ENDIF
	HEWID = 0.0  
        LINPRO = 'VOIGT   '
	GOTO 2
      ENDIF	



      X = 2. - ALPHAH(LINE) * .5
      GX = X - 1.0
      GAMMAF = 1+(-.5748646+(.9512363+(-.6998588+(.4245549-.1010678*GX
     ;           )*GX)*GX)*GX)*GX
      SIGMA = SIGMAH(LINE) * A0 * A0
      GVW = (4./PI)**(ALPHAH(LINE)*0.5)*GAMMAF*1.E4*SIGMA
      VBAR = SQRT(8.*KB*TEMP/PI/M0 * (1./1.008 + 1./4.003))
      HEWID = GVW * ((VBAR/1.E4)**(1.-ALPHAH(LINE))) 
      HEWID = HEWID * 0.625            !polarisability difference term
      HEWID = HEWID * HE1FRC * 1.E6    !in cgs


 2    CONTINUE

      RETURN

      END
