      SUBROUTINE CHECK_CONT(BKONCHECK, BKONCHECK2, 
     >             KONCHECK, KONCHECK2, NF, NF2, K, XK, 
     >             XLAMKOLD, XKC, XKC2, BANDM, BANDP, 
C***  For Continuum Interpolation
     >             CMFBAND, CMFBANDR, DXMAX, 
     >             KCCHECK, KCL, KCU, KCONT, KCDMAX, KCDELTA, 
     >             XLAMBDA, KONTACT, KONTAUP, DFKONT, BPLOT)
C****************************************************************
C***  It is checked wether Continua ar near the actual 
C***    frequency Point. Additionally the indices for the Continuum
C***    Interpolation are derived. 
C***    Called by COLI
C****************************************************************

      DIMENSION XKC(NF), XKC2(NF2), XLAMBDA(NF)
      LOGICAL BKONCHECK, BKONCHECK2, BPLOT

C***  Continuum point near actual frequency?
        BKONCHECK = .FALSE.
        IF (KONCHECK .LE. NF) THEN
          DO
            IF (XK .LT. XKC(KONCHECK)-BANDM) THEN
              EXIT
            ENDIF
            IF (XK .GE. XKC(KONCHECK)-BANDM .AND.
     >          XK .LE. XKC(KONCHECK)+BANDP) THEN
              BKONCHECK = .TRUE.
              EXIT
            ENDIF
            IF (XK .GT. XKC(KONCHECK)+BANDP) THEN
              KONCHECK = KONCHECK + 1
              IF (KONCHECK .GT. NF) EXIT
            ENDIF
          ENDDO
        ENDIF
        BKONCHECK2 =  .FALSE.
        IF (KONCHECK2 .LE. NF2) THEN
          DO
            IF (XK .LT. XKC2(KONCHECK2)-CMFBAND/DXMAX) THEN
              EXIT
            ENDIF
            IF (XK .GE. XKC2(KONCHECK2)-CMFBAND/DXMAX .AND.
     >          XK .LE. XKC2(KONCHECK2)+CMFBANDR/DXMAX) THEN
              BKONCHECK2 = .TRUE.
              EXIT
            ENDIF
            IF (XK .GT. XKC2(KONCHECK2)+CMFBANDR/DXMAX) THEN
              KONCHECK2 = KONCHECK2 + 1
              IF (KONCHECK2 .GT. NF2) EXIT
            ENDIF
          ENDDO
        ENDIF

C***    Test for Continuum Point in Interpolation Interval
        IF (KCCHECK .LE. NF) THEN
           DO
C***         First check, if Continuum is just passed
              IF (XK .GE. XKC(KCCHECK)) THEN
C***            Cont. Interval Finished
C***            ->  new Continuum, newstart (KCL = KCU)
                 KCCHECK = KCCHECK+1
                 KCL = K
                 KCU = K
                 IF (KCCHECK .GT. NF) EXIT
C***         Then find the Length of the Interpolation Interval
              ELSE
                 KCONT = INT(XKC(KCCHECK))
                 IF (KCU .LE. (KCONT-KCDMAX)) THEN
C***            Normal Case -> Length = Maximum Length
                    KCDELTA = KCDMAX
                    EXIT
                 ELSE
C***            Continuum detected -> Shorter Interval is cosen
                    KCDELTA = KCONT-KCU
                    EXIT
                 ENDIF
              ENDIF
           ENDDO
        ENDIF
 
C***    New index for Array 'XJC' reached ?
        IF (XLAMKOLD .GT. XLAMBDA(KONTAUP)) THEN
           KONTACT = KONTACT + 1
           KONTAUP = KONTACT + 1
           DFKONT  = 1./XLAMBDA(KONTACT) - 1./XLAMBDA(KONTAUP)
        ENDIF

        IF (BPLOT .AND. BKONCHECK) THEN
          WRITE (39,'(I6,1X,F3.0)') K, -2.
        ENDIF
        IF (BPLOT .AND. BKONCHECK2) THEN
          WRITE (39,'(I6,1X,F3.0)') K, -1.
        ENDIF

      RETURN
      END
