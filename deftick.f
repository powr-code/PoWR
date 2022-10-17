       SUBROUTINE DEFTICK(UVAL,NTMARK,NSMARK,BUV,BTM,BSM)

C***  Calculate "good" values for an axis-layout
C***  so that distances between ticks are multiples of 1,2 or 5
C***  This - incomplete - Version assumes 0.0 as lower value 
C***   and uval is greater than 0.0

        IMPLICIT NONE
C***    Parameters
C***    IN:
        REAL UVAL           ! the greatest upper value to appear
        INTEGER NTMARK      ! the number of tickmarks to approximate
        INTEGER NSMARK      ! the number of s-marks to approximate

C***    OUT:
        real buv            ! best upper value (an "even" number)
        real btm            ! best ticksmarks        
        real bsm            ! best s-marks

C**     locals

        real    amantissa   ! a..., because this is a real value
                            ! ...mantisse, because this is ...
        integer nexponent   ! n..., because ...
        real    aloga       ! see above
        real    bum         ! means: best upper mantissa
        integer n
        real    x
        real d10,d5,d2

C so nothing will go wrong
        if (uval .le. 0.0) uval=1.0     

C Describe uval in the form amantissa * 10**nexponent         
        aloga=alog10(uval)
        nexponent = AINT(aloga)
        if (aloga .lt. 0) nexponent=nexponent-1
        amantissa = 10**(aloga-nexponent)

C  Here we decide, what the best mantissa for the upper value is
        if (amantissa .le. 10.0) bum=10.0
        if (amantissa .le.  5.0) bum= 5.0
        if (amantissa .le.  2.0) bum= 2.0

C And here comes the best upper value
        buv=bum*(10.0**nexponent)

C This distance has a tickmark, if exactly ntmark's would appear
        x=bum/ntmark

C And these are the errors, if a good distance is used
        d10=abs( 1.0-x)
        d5 =abs( 0.5-x)
        d2 =abs( 0.2-x)

C This gets the distance with the smallest error
        n=10
        if (d5 .lt. d10) n=5
        if ((d2 .lt. d5) .and. (d2 .lt. d10)) n=2

C Don't forget the exponent
        btm=n*(10.0**(nexponent-1))

C Same procedure for smarks
        x=bum/nsmark
        d10=abs( 1.0-x)
        d5 =abs( 0.5-x)
        d2 =abs( 0.2-x)

        n=10
        if (d5 .lt. d10) n=5
        if ((d2 .lt. d5) .and. (d2 .lt. d10)) n=2

        bsm=n*(10.0**(nexponent-1))

C=================================
	return
        END
