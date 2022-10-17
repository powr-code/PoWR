      SUBROUTINE RUKU(RSTART, REND, T, ISTEP)

c      write (*,*) 'rstart, rend, t, istep:',rstart, rend, t, istep

      H = (REND - RSTART) / ISTEP
      H2 = H / 2.

c      write (*,*) 'h,h2=',h,h2

      TL = T
      RL = RSTART

c      stop 'stop in ruku'

      DO I = 1, ISTEP

        RIN = RL
        TIN = TL
        CALL DTDR(1, RIN, TIN, TPRIME)
        X1 = H * TPRIME

        RIN = RL + H2
        TIN = TL + X1/2.
        CALL DTDR(1, RIN, TIN, TPRIME)
        X2 = H * TPRIME

        RIN = RL + H2
        TIN = TL + X2/2.
        CALL DTDR(1, RIN, TIN, TPRIME)
        X3 = H * TPRIME

        RIN = RL + H
        TIN = TL + X3
        CALL DTDR(1, RIN, TIN, TPRIME)
        X4 = H * TPRIME

        RL = RL + H
        TL = TL + ((X1/2.) + X2 + X3 + (X4/2.)) / 3.

c        write (*,*) 'rl, tl=',rl,tl

      ENDDO

      T = TL

      RETURN
      END
