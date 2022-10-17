      FUNCTION TAUCMU (XMU,L,R,OPA,ND)
C***********************************************************************
C***  CALCULATION OF CONTINUUM OPTICAL DEPTH AT RADIUS R(L)
C***  AS FUNCTION OF MU (ANGLE COSINE)
C***  CORE RAYS: CONTINUUM OPTICAL DEPTH TO INNER BOUNDARY
C***********************************************************************
      DIMENSION R(ND),OPA(ND)

      RL=R(L)
      XSIN2=1.-XMU*XMU
      XSIN=SQRT(XSIN2)
      RMIN=RL*XSIN
      RMIN2=RMIN*RMIN
      TCMU=0.0

C***  INTEGRATION OF TAUC(MU) INSIDE RL (I.E. FOR RADII .LT. RL)
      IF (XMU .GT. 0.0) THEN
         RLMU=RL
         OPALM=OPA(L)
         SXOLD=0.0
         LMAX=ISRCHFLT(ND,R(1),1,RMIN)-1
C***     INTEGRATION FROM RL TO R(LMAX)
         DO 11 LMU=L+1,LMAX
         RLMU1=R(LMU)
         OPALM1=OPA(LMU)
         SX=RL*XMU-SQRT(RLMU1*RLMU1-RMIN2)
         DS=SX-SXOLD
         TCMU=TCMU+(OPALM+OPALM1)/2.*DS
         RLMU=RLMU1
         OPALM=OPALM1
   11    SXOLD=SX
C***     STEP FROM R(LMAX) TO RMIN
         IF (LMAX .GE. ND) THEN

C***        CORE RAY
            TAUCMU=TCMU
            RETURN

         ELSE
            OPARMIN=OPALM+(OPA(LMAX+1)-OPALM)/(RLMU-R(LMAX+1))*
     *                    (RLMU-RMIN)
            SX=RL*XMU
            DS=SX-SXOLD
            TCMU=TCMU+(OPALM+OPARMIN)/2.*DS
         ENDIF
C***     FACTOR 2. FOR SYMMETRICAL RUN OF MU-RAY INSIDE RL
         TCMU=2.*TCMU
      ENDIF

C***  INTEGRATION OF TAUC(MU) OUTSIDE RL (I.E. FOR RADII .GT. RL)
      RLMU=RL
      OPALM=OPA(L)
      SXOLD=0.0
      DO 12 LMU=L-1,1,-1
      RLMU1=R(LMU)
      OPALM1=OPA(LMU)
      SX=RL*XMU+SQRT(RLMU1*RLMU1-RMIN2)
      DS=SX-SXOLD
      TCMU=TCMU+(OPALM+OPALM1)/2.*DS
      RLMU=RLMU1
      OPALM=OPALM1
   12 SXOLD=SX

      TAUCMU=TCMU

      RETURN
      END
