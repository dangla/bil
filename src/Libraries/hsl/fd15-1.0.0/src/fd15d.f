* COPYRIGHT (c) 1988 AEA Technology
* Original date 17 Feb 2005

C 17th February 2005 Version 1.0.0. Replacement for FD05.

      DOUBLE PRECISION FUNCTION FD15AD(T)
C----------------------------------------------------------------
C  Fortran 77 implementation of the Fortran 90 intrinsic
C    functions: EPSILON, TINY, HUGE and RADIX.  Note that
C    the RADIX result is returned as DOUBLE PRECISION.
C
C  The CHARACTER argument specifies the type of result:
C       
C   'E'  smallest positive real number: 1.0 + DC(1) > 1.0, i.e.
C          EPSILON(DOUBLE PRECISION)
C   'T'  smallest full precision positive real number, i.e.
C          TINY(DOUBLE PRECISION)
C   'H'  largest finite positive real number, i.e.
C          HUGE(DOUBLE PRECISION)
C   'R'  the base of the floating point arithematic, i.e.
C          RADIX(DOUBLE PRECISION)
C
C    any other value gives a result of zero.
C----------------------------------------------------------------
      CHARACTER T

      IF ( T.EQ.'E' ) THEN
         FD15AD = EPSILON(1.0D0)
      ELSE IF ( T.EQ.'T' ) THEN
         FD15AD = TINY(1.0D0)
      ELSE IF ( T.EQ.'H' ) THEN
         FD15AD = HUGE(1.0D0)
      ELSE IF ( T.EQ.'R' ) THEN
         FD15AD = DBLE(RADIX(1.0D0))
      ELSE
         FD15AD = 0.0D0
      ENDIF
      RETURN
      END
