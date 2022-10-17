      FUNCTION VOIGTH(A,V)
C***********************************************************************
C***  normalized (!!!) VOIGT function H(a,v)
C***  ----- nach DETLEF KOESTER -----
C***********************************************************************
      DATA WPI /1.77245385/

      V=ABS(V)
      VV=V*V
      IF ((A.LT.0.2) .AND. (V.GT.5.0)) GO TO 70
      IF ((A.GE.0.2) .AND. (A.GT.1.4 .OR. (A+V).GT.3.2)) GO TO 80
      H0=EXP(-VV)
      H2=(1.-2.*VV)*H0
      IF (V.GT.2.4) GO TO 30
      IF (V.LE.1.3) GO TO 20
      H1=(-.220416*VV+1.989196*V-6.61487)*VV+9.39456*V-4.4848
      GO TO 40
   20 H1=(.42139*VV-2.34358*V+3.28868)*VV-.15517*V-1.1247
      GO TO 40
   30 H1=((-.0032783*VV+.0429913*V-.188326)*VV+.278712*V+.55415)/
     *   (VV-1.5)
   40 CONTINUE
      IF (A.GE..2) GO TO 60
   50 H=((H2*A+H1)*A+H0)
      VOIGTH=H/WPI
      RETURN

   60 HH1=H1+H0*1.12838
      HH2=H2+HH1*1.12838-H0
      HH3=(1.-H2)*0.37613-HH1*0.66667*VV+HH2*1.12838
      HH4=(3.*HH3-HH1)*0.37613+H0*0.66667*VV*VV
      H=((((HH4*A+HH3)*A+HH2)*A+HH1)*A+H0)*
     *  (((-.122727278*A+.532770573)*A-.96284325)*A+.979895032)
      VOIGTH=H/WPI
      RETURN

   70 H=((2.12/VV+.8463)/VV+.5642)*A/VV
      VOIGTH=H/WPI
      RETURN

   80 AA=A*A
      U=(AA+VV)*1.4142
      UU=U*U
      H=((((AA-10.*VV)*AA*3.+15.*VV*VV)/UU+3.*VV-AA)/UU+1.)*A*.79788/U
      VOIGTH=H/WPI
      RETURN

      END
