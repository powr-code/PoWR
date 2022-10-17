	FUNCTION VCSE1F(X)
C
C  E1 function calculator for VCS approximation. It's rough, but 
C  arranged to be fast. X must be >=0.
C
C  From Atlas9 & Kurucz code (mario)
C  Approximation der vidal-cooper-smith daten
C

      VCSE1F=0.0
      IF(X.LE.0.0) RETURN
      IF(X.LE.0.01) THEN
        VCSE1F=-LOG(X)-0.577215+X
      ELSE IF(X.LE.1.0) THEN
        VCSE1F=-LOG(X)-0.57721566+X*(0.99999193+X*(-0.24991055+
     +                            X*(0.05519968+X*(-0.00976004+
     +                            X*0.00107857))))
      ELSE IF(X.LE.30.) THEN
        VCSE1F=(X*(X+2.334733)+0.25062)/(X*(X+3.330657)+
     +         1.681534)/X*EXP(-X)
      END IF
C

      RETURN
      END
