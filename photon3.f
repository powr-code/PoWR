      FUNCTION PHOTON3 (X,NO)
C***********************************************************************
C***  CALCULATION OF PART OF THE SUM (I= 3 TO 6) IN THE EXTENDED VERSION
C***  OF FORMULA 12 (SEE DESCRIPTION OF PROGRAM "DETAIL") FOR RBF CROSS
C***  SECTIONS
C***  ---  CALLED FROM PHOTOCS: IF (IGAUNT(LOW) .EQ. 'DETAILN3')  ---
C***  NECESSARY FOR PHOTOIONISATION CROSSSECTIONS OF 10 LEVELS OF THE
C***  QUARTET SYSTEM IN NITROGEN III
C***  THE LAST 4 COEFFICIENTS ARE STORED IN DATA STATEMENTS!!!
C***********************************************************************

C***  NUMBER OF CONSIDERED LEVELS:
      PARAMETER ( N3QUART = 10 )

      DIMENSION AI(3:6,N3QUART)
C***  DETAIL LEVEL:     A4P1  ==>  N 32P2P4.2               (1)
      DATA (AI(I,1),I=3,6) /
     $ 6.37367E-2, 0.23060, 0.24562, 6.99181E-2 /
C***  DETAIL LEVEL:     A4S1  ==>  N 32P3S4.6               (2)
      DATA (AI(I,2),I=3,6) /
     $ -0.22891, -0.13430, -3.79655E-3, 5.93594E-3 /
C***  DETAIL LEVEL:     A4P2  ==>  N 33S'P412               (3)
      DATA (AI(I,3),I=3,6) /
     $ -0.56510, -0.30266, -0.10805, 3.69164E-3 /
C***  DETAIL LEVEL:     A4D1  ==>  N 33P'D415               (4)
      DATA (AI(I,4),I=3,6) /
     $ -0.50098, -0.32843, -0.15584, -3.08269E-2 /
C***  DETAIL LEVEL:     A4S2  ==>  N 33P'S417               (5)
      DATA (AI(I,5),I=3,6) /
     $ -0.55232, -0.40950, -0.20078, -3.97182E-2 /
C***  DETAIL LEVEL:     A4P3  ==>  N 33P'P418               (6)
      DATA (AI(I,6),I=3,6) /
     $ -0.52906, -0.41945, -0.21380, -4.24112E-2 /
C***  DETAIL LEVEL:     A4F1  ==>  N 33D'F422               (7)
      DATA (AI(I,7),I=3,6) /
     $ 4.06244E-2, -5.55455E-2, 4.59197E-2, 1.38044E-2 /
C***  DETAIL LEVEL:     A4D2  ==>  N 33D'D423               (8)
      DATA (AI(I,8),I=3,6) /
     $ 7.02743E-2, 2.16001E-2, 7.56214E-2, 1.59807E-2 /
C***  DETAIL LEVEL:     A4P4  ==>  N 33D'P424               (9)
      DATA (AI(I,9),I=3,6) /
     $ 4.44634E-2, 8.28837E-2, 9.21828E-2, 1.73687E-2 /
C***  DETAIL LEVEL:     A4P5  ==>  N 34S'P431               (10)
      DATA (AI(I,10),I=3,6) /
     $ -0.13676, 0.17396, 0.12893, 2.99775E-2 /

C***  ERROR STOP:
      IF ((NO .LT. 1) .OR. (NO .GT. N3QUART)) STOP 'ERROR'

C***  CALCULATION OF PART OF THE SUM (I= 3 TO 6):
      XLN=ALOG(X)
      PHOTON3=XLN*(AI(3,NO)+XLN*(AI(4,NO)+XLN*(AI(5,NO)+AI(6,NO)*XLN)))

      RETURN
      END
