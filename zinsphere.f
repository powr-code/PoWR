      SUBROUTINE ZINSPHERE(ZINTER, ALPHA, SDIST, SRAD, 
     >                     LPSTA, LPEND, NPHI, 
     >                     P, NP)
C******************************************************************
C***  CALCULATES THE POINTS OF INTERSECTION FOR THE SPHERE MODEL
C***  ALPHA IS THE ANGLE OF THE INCLINATION OF THE CONE
C***  SRAD IS THE RADIUS OF THE SPHERE AND SDIST THE DISPLACEMENT
C******************************************************************

      DIMENSION ZINTER(2, NP, NPHI), P(NP)

      DATA PI / 3.14159265358979 /

c      iphitest = 1
c      write (*,*) 'iphitest, nphi=', iphitest, nphi

c      write (*,*) 'ZINSPHERE: SDIST, SRAD=', SDIST, SRAD
c      write (*,*) 'ZINSPHERE: ALPHA=',ALPHA

      SA = SIN(ALPHA)
      CA = COS(ALPHA)

      DO LPHI=LPSTA, LPEND
c      do lphi=iphitest, iphitest
        IF (NPHI .GT. 1) THEN
          PHI=PI * (LPHI-1) / (NPHI-1)
        ENDIF
c      write (*,*) 'phi=', phi
        SPHI  = SIN(PHI)
        SPHI2 = SPHI * SPHI
        CPHI  = COS(PHI)
        CPHI2 = CPHI * CPHI
        DO JP=1, NP
          XP = P(JP)
c      write (*,*) 'xp=', xp
C***  NOTE RAD IS THE RADICANT
          RAD = SRAD*SRAD - (SDIST*SPHI*CA)*(SDIST*SPHI*CA) -
     >            (XP-SDIST*CA*CPHI)*(XP-SDIST*CA*CPHI)
c          write (*,*) 'xnenn=',xnenn
c          stop 'test in zinsphere'
c         write (*,*) 'rad=', rad
          IF (RAD .LT. 0.) THEN
c            write (0,*) 'radicant negative',rad
            ZINTER(1, JP, LPHI) = 0.
            ZINTER(2, JP, LPHI) = 0.
          ELSE
c            write (0,*) 'radicant positiv',rad
            SQRRAD = SQRT(RAD)
            Z0 = SDIST * SA
c        write (*,*) 'z0=', z0
            ZINTER(1, JP, LPHI) = Z0 + SQRRAD
c        write (*,*) 'zinter=', zinter(1,jp,lphi), zinter(2,jp,lphi)
            ZINTER(2, JP, LPHI) = Z0 - SQRRAD
          ENDIF

        ENDDO
      ENDDO

c      write (*,'(3f15.5)') 
c     >    (p(ip),zinter(1,ip,iphitest),zinter(2,ip,iphitest),ip=1, np)

c      stop 'test nach zinsphere'

      RETURN
      END
